package workspace

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/lishi/chemvcs/internal/objectstore"
)

func TestRestoreDirectory_DeletesExtraFiles(t *testing.T) {
	// Create a temporary directory for the test repository
	tmpDir := t.TempDir()
	repoDir := filepath.Join(tmpDir, "repo")
	os.MkdirAll(repoDir, 0755)

	// Initialize object store
	store, err := objectstore.NewStore(repoDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	scanner := NewScanner(store)

	// Create working directory with initial files
	workDir := filepath.Join(tmpDir, "work")
	os.MkdirAll(workDir, 0755)

	// Create initial snapshot with file1.txt and file2.txt
	os.WriteFile(filepath.Join(workDir, "file1.txt"), []byte("content1"), 0644)
	os.WriteFile(filepath.Join(workDir, "file2.txt"), []byte("content2"), 0644)
	os.MkdirAll(filepath.Join(workDir, "subdir"), 0755)
	os.WriteFile(filepath.Join(workDir, "subdir", "file3.txt"), []byte("content3"), 0644)

	// Now modify working directory - add extra files
	os.WriteFile(filepath.Join(workDir, "extra.txt"), []byte("should be deleted"), 0644)
	os.MkdirAll(filepath.Join(workDir, "extradir"), 0755)
	os.WriteFile(filepath.Join(workDir, "extradir", "extra2.txt"), []byte("also deleted"), 0644)

	// Create a new snapshot with only file1.txt
	workDir2 := filepath.Join(tmpDir, "work2")
	os.MkdirAll(workDir2, 0755)
	os.WriteFile(filepath.Join(workDir2, "file1.txt"), []byte("content1"), 0644)

	rootObj2, err := scanner.ScanDirectory(workDir2)
	if err != nil {
		t.Fatalf("Failed to scan directory: %v", err)
	}

	rootHash2, err := store.PutObject(rootObj2)
	if err != nil {
		t.Fatalf("Failed to store root object: %v", err)
	}

	// Restore snapshot2 into workDir (should delete extra files and file2.txt, subdir)
	if err := scanner.RestoreDirectory(rootHash2, workDir); err != nil {
		t.Fatalf("Failed to restore directory: %v", err)
	}

	// Verify that only file1.txt exists
	if _, err := os.Stat(filepath.Join(workDir, "file1.txt")); os.IsNotExist(err) {
		t.Error("file1.txt should exist but doesn't")
	}

	// Verify that file2.txt was deleted
	if _, err := os.Stat(filepath.Join(workDir, "file2.txt")); !os.IsNotExist(err) {
		t.Error("file2.txt should be deleted but still exists")
	}

	// Verify that extra.txt was deleted
	if _, err := os.Stat(filepath.Join(workDir, "extra.txt")); !os.IsNotExist(err) {
		t.Error("extra.txt should be deleted but still exists")
	}

	// Verify that subdir was deleted
	if _, err := os.Stat(filepath.Join(workDir, "subdir")); !os.IsNotExist(err) {
		t.Error("subdir should be deleted but still exists")
	}

	// Verify that extradir was deleted
	if _, err := os.Stat(filepath.Join(workDir, "extradir")); !os.IsNotExist(err) {
		t.Error("extradir should be deleted but still exists")
	}
}

func TestRestoreDirectory_PreservesChemVCS(t *testing.T) {
	// Create a temporary directory for the test repository
	tmpDir := t.TempDir()
	repoDir := filepath.Join(tmpDir, "repo")
	os.MkdirAll(repoDir, 0755)

	// Initialize object store
	store, err := objectstore.NewStore(repoDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	scanner := NewScanner(store)

	// Create working directory
	workDir := filepath.Join(tmpDir, "work")
	os.MkdirAll(workDir, 0755)

	// Create .chemvcs directory
	chemvcsDir := filepath.Join(workDir, ".chemvcs")
	os.MkdirAll(chemvcsDir, 0755)
	os.WriteFile(filepath.Join(chemvcsDir, "config"), []byte("test config"), 0644)

	// Create a file
	os.WriteFile(filepath.Join(workDir, "file1.txt"), []byte("content1"), 0644)

	// Scan to create snapshot
	rootObj, err := scanner.ScanDirectory(workDir)
	if err != nil {
		t.Fatalf("Failed to scan directory: %v", err)
	}

	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("Failed to store root object: %v", err)
	}

	// Add extra file
	os.WriteFile(filepath.Join(workDir, "extra.txt"), []byte("extra"), 0644)

	// Restore directory (should delete extra.txt but preserve .chemvcs)
	if err := scanner.RestoreDirectory(rootHash, workDir); err != nil {
		t.Fatalf("Failed to restore directory: %v", err)
	}

	// Verify .chemvcs still exists
	if _, err := os.Stat(filepath.Join(chemvcsDir, "config")); os.IsNotExist(err) {
		t.Error(".chemvcs/config should be preserved but was deleted")
	}

	// Verify extra.txt was deleted
	if _, err := os.Stat(filepath.Join(workDir, "extra.txt")); !os.IsNotExist(err) {
		t.Error("extra.txt should be deleted but still exists")
	}
}

func TestRestoreDirectory_EmptyDirectory(t *testing.T) {
	// Create a temporary directory for the test repository
	tmpDir := t.TempDir()
	repoDir := filepath.Join(tmpDir, "repo")
	os.MkdirAll(repoDir, 0755)

	// Initialize object store
	store, err := objectstore.NewStore(repoDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	scanner := NewScanner(store)

	// Create working directory with files
	workDir := filepath.Join(tmpDir, "work")
	os.MkdirAll(workDir, 0755)
	os.WriteFile(filepath.Join(workDir, "file1.txt"), []byte("content1"), 0644)
	os.WriteFile(filepath.Join(workDir, "file2.txt"), []byte("content2"), 0644)

	// Create empty directory snapshot
	emptyDir := filepath.Join(tmpDir, "empty")
	os.MkdirAll(emptyDir, 0755)

	rootObj, err := scanner.ScanDirectory(emptyDir)
	if err != nil {
		t.Fatalf("Failed to scan empty directory: %v", err)
	}

	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("Failed to store root object: %v", err)
	}

	// Restore empty snapshot (should delete all files)
	if err := scanner.RestoreDirectory(rootHash, workDir); err != nil {
		t.Fatalf("Failed to restore directory: %v", err)
	}

	// Verify all files are deleted
	entries, err := os.ReadDir(workDir)
	if err != nil {
		t.Fatalf("Failed to read directory: %v", err)
	}

	// Should only contain .chemvcs if it exists
	for _, entry := range entries {
		if entry.Name() != ".chemvcs" {
			t.Errorf("Found unexpected entry after cleanup: %s", entry.Name())
		}
	}
}
