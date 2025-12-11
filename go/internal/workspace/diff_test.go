package workspace

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/lishi/chemvcs/internal/objectstore"
)

func TestComputeChanges_NoChanges(t *testing.T) {
	// Create test directory
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	// Create test file
	if err := os.WriteFile(filepath.Join(tmpDir, "test.txt"), []byte("content"), 0644); err != nil {
		t.Fatalf("failed to create test file: %v", err)
	}

	// Scan and store
	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	rootObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	// Compare with itself
	changes, err := ComputeChanges(rootHash, rootHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Should have no changes
	if len(changes.Added) != 0 {
		t.Errorf("expected 0 added, got %d", len(changes.Added))
	}
	if len(changes.Modified) != 0 {
		t.Errorf("expected 0 modified, got %d", len(changes.Modified))
	}
	if len(changes.Deleted) != 0 {
		t.Errorf("expected 0 deleted, got %d", len(changes.Deleted))
	}
}

func TestComputeChanges_FileAdded(t *testing.T) {
	// Create initial directory
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan empty directory
	oldObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	oldHash, err := store.PutObject(oldObj)
	if err != nil {
		t.Fatalf("failed to store old root: %v", err)
	}

	// Add a file
	if err := os.WriteFile(filepath.Join(tmpDir, "new.txt"), []byte("new content"), 0644); err != nil {
		t.Fatalf("failed to create new file: %v", err)
	}

	// Scan again
	newObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	newHash, err := store.PutObject(newObj)
	if err != nil {
		t.Fatalf("failed to store new root: %v", err)
	}

	// Compute changes
	changes, err := ComputeChanges(oldHash, newHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Should have 1 added file
	if len(changes.Added) != 1 {
		t.Errorf("expected 1 added, got %d", len(changes.Added))
	} else if changes.Added[0] != "new.txt" {
		t.Errorf("expected 'new.txt' added, got %q", changes.Added[0])
	}

	if len(changes.Modified) != 0 {
		t.Errorf("expected 0 modified, got %d", len(changes.Modified))
	}
	if len(changes.Deleted) != 0 {
		t.Errorf("expected 0 deleted, got %d", len(changes.Deleted))
	}
}

func TestComputeChanges_FileDeleted(t *testing.T) {
	// Create directory with file
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	testFile := filepath.Join(tmpDir, "delete-me.txt")
	if err := os.WriteFile(testFile, []byte("will be deleted"), 0644); err != nil {
		t.Fatalf("failed to create test file: %v", err)
	}

	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan with file
	oldObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	oldHash, err := store.PutObject(oldObj)
	if err != nil {
		t.Fatalf("failed to store old root: %v", err)
	}

	// Delete file
	if err := os.Remove(testFile); err != nil {
		t.Fatalf("failed to delete file: %v", err)
	}

	// Scan again
	newObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	newHash, err := store.PutObject(newObj)
	if err != nil {
		t.Fatalf("failed to store new root: %v", err)
	}

	// Compute changes
	changes, err := ComputeChanges(oldHash, newHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Should have 1 deleted file
	if len(changes.Deleted) != 1 {
		t.Errorf("expected 1 deleted, got %d", len(changes.Deleted))
	} else if changes.Deleted[0] != "delete-me.txt" {
		t.Errorf("expected 'delete-me.txt' deleted, got %q", changes.Deleted[0])
	}

	if len(changes.Added) != 0 {
		t.Errorf("expected 0 added, got %d", len(changes.Added))
	}
	if len(changes.Modified) != 0 {
		t.Errorf("expected 0 modified, got %d", len(changes.Modified))
	}
}

func TestComputeChanges_FileModified(t *testing.T) {
	// Create directory with file
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	testFile := filepath.Join(tmpDir, "modify-me.txt")
	if err := os.WriteFile(testFile, []byte("original content"), 0644); err != nil {
		t.Fatalf("failed to create test file: %v", err)
	}

	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan initial state
	oldObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	oldHash, err := store.PutObject(oldObj)
	if err != nil {
		t.Fatalf("failed to store old root: %v", err)
	}

	// Modify file
	if err := os.WriteFile(testFile, []byte("modified content"), 0644); err != nil {
		t.Fatalf("failed to modify file: %v", err)
	}

	// Scan again
	newObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	newHash, err := store.PutObject(newObj)
	if err != nil {
		t.Fatalf("failed to store new root: %v", err)
	}

	// Compute changes
	changes, err := ComputeChanges(oldHash, newHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Should have 1 modified file
	if len(changes.Modified) != 1 {
		t.Errorf("expected 1 modified, got %d", len(changes.Modified))
	} else if changes.Modified[0] != "modify-me.txt" {
		t.Errorf("expected 'modify-me.txt' modified, got %q", changes.Modified[0])
	}

	if len(changes.Added) != 0 {
		t.Errorf("expected 0 added, got %d", len(changes.Added))
	}
	if len(changes.Deleted) != 0 {
		t.Errorf("expected 0 deleted, got %d", len(changes.Deleted))
	}
}

func TestComputeChanges_ComplexScenario(t *testing.T) {
	// Create directory structure
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	subDir := filepath.Join(tmpDir, "subdir")
	if err := os.MkdirAll(subDir, 0755); err != nil {
		t.Fatalf("failed to create subdir: %v", err)
	}

	// Create initial files
	files := map[string]string{
		"keep.txt":          "unchanged",
		"modify.txt":        "original",
		"delete.txt":        "will be deleted",
		"subdir/nested.txt": "nested file",
	}

	for relPath, content := range files {
		path := filepath.Join(tmpDir, relPath)
		if err := os.WriteFile(path, []byte(content), 0644); err != nil {
			t.Fatalf("failed to create file %s: %v", relPath, err)
		}
	}

	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan initial state
	oldObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	oldHash, err := store.PutObject(oldObj)
	if err != nil {
		t.Fatalf("failed to store old root: %v", err)
	}

	// Make changes
	// - Modify modify.txt
	if err := os.WriteFile(filepath.Join(tmpDir, "modify.txt"), []byte("modified"), 0644); err != nil {
		t.Fatalf("failed to modify file: %v", err)
	}
	// - Delete delete.txt
	if err := os.Remove(filepath.Join(tmpDir, "delete.txt")); err != nil {
		t.Fatalf("failed to delete file: %v", err)
	}
	// - Add new.txt
	if err := os.WriteFile(filepath.Join(tmpDir, "new.txt"), []byte("new file"), 0644); err != nil {
		t.Fatalf("failed to create new file: %v", err)
	}
	// - Modify nested file
	if err := os.WriteFile(filepath.Join(subDir, "nested.txt"), []byte("modified nested"), 0644); err != nil {
		t.Fatalf("failed to modify nested file: %v", err)
	}

	// Scan new state
	newObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	newHash, err := store.PutObject(newObj)
	if err != nil {
		t.Fatalf("failed to store new root: %v", err)
	}

	// Compute changes
	changes, err := ComputeChanges(oldHash, newHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Verify changes
	// Added: new.txt
	if len(changes.Added) != 1 {
		t.Errorf("expected 1 added, got %d: %v", len(changes.Added), changes.Added)
	} else if changes.Added[0] != "new.txt" {
		t.Errorf("expected 'new.txt' added, got %q", changes.Added[0])
	}

	// Modified: modify.txt, subdir/nested.txt
	if len(changes.Modified) != 2 {
		t.Errorf("expected 2 modified, got %d: %v", len(changes.Modified), changes.Modified)
	} else {
		// Check both are present (order may vary)
		foundModify := false
		foundNested := false
		for _, path := range changes.Modified {
			if path == "modify.txt" {
				foundModify = true
			}
			if path == "subdir/nested.txt" {
				foundNested = true
			}
		}
		if !foundModify {
			t.Errorf("'modify.txt' not in modified list")
		}
		if !foundNested {
			t.Errorf("'subdir/nested.txt' not in modified list")
		}
	}

	// Deleted: delete.txt
	if len(changes.Deleted) != 1 {
		t.Errorf("expected 1 deleted, got %d: %v", len(changes.Deleted), changes.Deleted)
	} else if changes.Deleted[0] != "delete.txt" {
		t.Errorf("expected 'delete.txt' deleted, got %q", changes.Deleted[0])
	}
}

func TestComputeChanges_DirectoryAdded(t *testing.T) {
	// Create initial empty directory
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan empty directory
	oldObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	oldHash, err := store.PutObject(oldObj)
	if err != nil {
		t.Fatalf("failed to store old root: %v", err)
	}

	// Add new directory with files
	newDir := filepath.Join(tmpDir, "newdir")
	if err := os.MkdirAll(newDir, 0755); err != nil {
		t.Fatalf("failed to create new dir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(newDir, "file1.txt"), []byte("content1"), 0644); err != nil {
		t.Fatalf("failed to create file1: %v", err)
	}
	if err := os.WriteFile(filepath.Join(newDir, "file2.txt"), []byte("content2"), 0644); err != nil {
		t.Fatalf("failed to create file2: %v", err)
	}

	// Scan again
	newObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}
	newHash, err := store.PutObject(newObj)
	if err != nil {
		t.Fatalf("failed to store new root: %v", err)
	}

	// Compute changes
	changes, err := ComputeChanges(oldHash, newHash, store)
	if err != nil {
		t.Fatalf("ComputeChanges failed: %v", err)
	}

	// Should have 2 added files
	if len(changes.Added) != 2 {
		t.Errorf("expected 2 added, got %d: %v", len(changes.Added), changes.Added)
	}

	// Both files should be in added list
	expected := map[string]bool{
		"newdir/file1.txt": false,
		"newdir/file2.txt": false,
	}
	for _, path := range changes.Added {
		if _, ok := expected[path]; ok {
			expected[path] = true
		}
	}
	for path, found := range expected {
		if !found {
			t.Errorf("expected %q in added list", path)
		}
	}
}
