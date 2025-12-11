package workspace

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/lishi/chemvcs/internal/objectstore"
)

func TestScanDirectory_Empty(t *testing.T) {
	// Create temporary directory structure
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	// Create scanner
	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan empty directory
	rootObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	// Verify root object
	if rootObj.Type != "folder" {
		t.Errorf("expected Type='folder', got %q", rootObj.Type)
	}
	if len(rootObj.Refs) != 0 {
		t.Errorf("expected empty folder, got %d refs", len(rootObj.Refs))
	}
}

func TestScanDirectory_WithFiles(t *testing.T) {
	// Create temporary directory structure
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	// Create test files
	testFiles := map[string]string{
		"file1.txt": "content of file 1",
		"file2.txt": "content of file 2",
	}

	for name, content := range testFiles {
		path := filepath.Join(tmpDir, name)
		if err := os.WriteFile(path, []byte(content), 0644); err != nil {
			t.Fatalf("failed to create test file: %v", err)
		}
	}

	// Create scanner
	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan directory
	rootObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	// Verify root object
	if rootObj.Type != "folder" {
		t.Errorf("expected Type='folder', got %q", rootObj.Type)
	}
	if len(rootObj.Refs) != 2 {
		t.Errorf("expected 2 file refs, got %d", len(rootObj.Refs))
	}

	// Verify all refs are objects
	for i, ref := range rootObj.Refs {
		if ref.Kind != "object" {
			t.Errorf("ref[%d]: expected Kind='object', got %q", i, ref.Kind)
		}
		if ref.ID == "" {
			t.Errorf("ref[%d]: ID is empty", i)
		}

		// Verify file object exists
		fileObj, err := store.GetObject(ref.ID)
		if err != nil {
			t.Errorf("ref[%d]: failed to get object: %v", i, err)
			continue
		}

		if fileObj.Type != "file" {
			t.Errorf("ref[%d]: expected Type='file', got %q", i, fileObj.Type)
		}

		// Verify file has blob reference
		if len(fileObj.Refs) != 1 {
			t.Errorf("ref[%d]: expected 1 blob ref, got %d", i, len(fileObj.Refs))
			continue
		}

		blobRef := fileObj.Refs[0]
		if blobRef.Kind != "blob" {
			t.Errorf("ref[%d]: expected blob ref, got %q", i, blobRef.Kind)
		}

		// Verify blob content
		blobData, err := store.GetBlob(blobRef.ID)
		if err != nil {
			t.Errorf("ref[%d]: failed to get blob: %v", i, err)
			continue
		}

		// Verify content matches one of our test files
		found := false
		for _, expectedContent := range testFiles {
			if string(blobData) == expectedContent {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("ref[%d]: blob content doesn't match any test file", i)
		}
	}
}

func TestScanDirectory_WithSubdirectories(t *testing.T) {
	// Create temporary directory structure
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	// Create nested structure
	subDir := filepath.Join(tmpDir, "subdir")
	if err := os.MkdirAll(subDir, 0755); err != nil {
		t.Fatalf("failed to create subdir: %v", err)
	}

	// Create files
	if err := os.WriteFile(filepath.Join(tmpDir, "root.txt"), []byte("root file"), 0644); err != nil {
		t.Fatalf("failed to create root file: %v", err)
	}
	if err := os.WriteFile(filepath.Join(subDir, "nested.txt"), []byte("nested file"), 0644); err != nil {
		t.Fatalf("failed to create nested file: %v", err)
	}

	// Create scanner
	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan directory
	rootObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	// Verify root has 2 children (1 file + 1 folder)
	if len(rootObj.Refs) != 2 {
		t.Errorf("expected 2 refs in root, got %d", len(rootObj.Refs))
	}

	// Find the folder object
	var folderObj *objectstore.Store
	for _, ref := range rootObj.Refs {
		obj, err := store.GetObject(ref.ID)
		if err != nil {
			t.Errorf("failed to get object: %v", err)
			continue
		}

		if obj.Type == "folder" {
			// Verify folder has 1 file
			if len(obj.Refs) != 1 {
				t.Errorf("expected 1 ref in subfolder, got %d", len(obj.Refs))
			}
		}
	}

	if folderObj != nil {
		// Additional validation passed
	}
}

func TestScanDirectory_SkipsChemVCS(t *testing.T) {
	// Create temporary directory structure
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(repoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create test repo: %v", err)
	}

	// Create a file inside .chemvcs (should be ignored)
	if err := os.WriteFile(filepath.Join(repoPath, "test.txt"), []byte("should be ignored"), 0644); err != nil {
		t.Fatalf("failed to create test file: %v", err)
	}

	// Create a normal file
	if err := os.WriteFile(filepath.Join(tmpDir, "normal.txt"), []byte("normal file"), 0644); err != nil {
		t.Fatalf("failed to create normal file: %v", err)
	}

	// Create scanner
	store, err := objectstore.NewStore(tmpDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	// Scan directory
	rootObj, err := scanner.ScanDirectory(tmpDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	// Verify only 1 file is tracked (normal.txt, not .chemvcs)
	if len(rootObj.Refs) != 1 {
		t.Errorf("expected 1 ref (normal file only), got %d", len(rootObj.Refs))
	}
}

func TestRestoreDirectory_Simple(t *testing.T) {
	// Create source directory
	srcDir := t.TempDir()
	srcRepoPath := filepath.Join(srcDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(srcRepoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create source repo: %v", err)
	}

	// Create test files
	testContent := "test file content"
	if err := os.WriteFile(filepath.Join(srcDir, "test.txt"), []byte(testContent), 0644); err != nil {
		t.Fatalf("failed to create test file: %v", err)
	}

	// Scan directory
	store, err := objectstore.NewStore(srcDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	rootObj, err := scanner.ScanDirectory(srcDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	// Store root object
	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	// Create destination directory
	dstDir := t.TempDir()

	// Restore to destination
	if err := scanner.RestoreDirectory(rootHash, dstDir); err != nil {
		t.Fatalf("RestoreDirectory failed: %v", err)
	}

	// Verify file exists and has correct content
	restoredPath := filepath.Join(dstDir, "test.txt")
	data, err := os.ReadFile(restoredPath)
	if err != nil {
		t.Fatalf("failed to read restored file: %v", err)
	}

	if string(data) != testContent {
		t.Errorf("content mismatch: expected %q, got %q", testContent, string(data))
	}
}

func TestRestoreDirectory_WithSubdirectories(t *testing.T) {
	// Create source directory with nested structure
	srcDir := t.TempDir()
	srcRepoPath := filepath.Join(srcDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(srcRepoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create source repo: %v", err)
	}

	subDir := filepath.Join(srcDir, "subdir")
	if err := os.MkdirAll(subDir, 0755); err != nil {
		t.Fatalf("failed to create subdir: %v", err)
	}

	// Create files
	rootContent := "root content"
	nestedContent := "nested content"
	if err := os.WriteFile(filepath.Join(srcDir, "root.txt"), []byte(rootContent), 0644); err != nil {
		t.Fatalf("failed to create root file: %v", err)
	}
	if err := os.WriteFile(filepath.Join(subDir, "nested.txt"), []byte(nestedContent), 0644); err != nil {
		t.Fatalf("failed to create nested file: %v", err)
	}

	// Scan and store
	store, err := objectstore.NewStore(srcDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	rootObj, err := scanner.ScanDirectory(srcDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	// Restore to new location
	dstDir := t.TempDir()
	if err := scanner.RestoreDirectory(rootHash, dstDir); err != nil {
		t.Fatalf("RestoreDirectory failed: %v", err)
	}

	// Verify root file
	rootData, err := os.ReadFile(filepath.Join(dstDir, "root.txt"))
	if err != nil {
		t.Fatalf("failed to read restored root file: %v", err)
	}
	if string(rootData) != rootContent {
		t.Errorf("root file content mismatch")
	}

	// Verify nested file
	nestedData, err := os.ReadFile(filepath.Join(dstDir, "subdir", "nested.txt"))
	if err != nil {
		t.Fatalf("failed to read restored nested file: %v", err)
	}
	if string(nestedData) != nestedContent {
		t.Errorf("nested file content mismatch")
	}

	// Verify directory exists
	if info, err := os.Stat(filepath.Join(dstDir, "subdir")); err != nil {
		t.Errorf("subdir not restored: %v", err)
	} else if !info.IsDir() {
		t.Errorf("subdir is not a directory")
	}
}

func TestScanRestore_RoundTrip(t *testing.T) {
	// Create complex directory structure
	srcDir := t.TempDir()
	srcRepoPath := filepath.Join(srcDir, ".chemvcs")
	if err := os.MkdirAll(filepath.Join(srcRepoPath, "objects"), 0755); err != nil {
		t.Fatalf("failed to create source repo: %v", err)
	}

	// Create multiple levels
	dirs := []string{
		filepath.Join(srcDir, "dir1"),
		filepath.Join(srcDir, "dir2"),
		filepath.Join(srcDir, "dir1", "subdir"),
	}
	for _, dir := range dirs {
		if err := os.MkdirAll(dir, 0755); err != nil {
			t.Fatalf("failed to create dir: %v", err)
		}
	}

	// Create files
	files := map[string]string{
		"file1.txt":             "content 1",
		"dir1/file2.txt":        "content 2",
		"dir1/subdir/file3.txt": "content 3",
		"dir2/file4.txt":        "content 4",
	}
	for relPath, content := range files {
		path := filepath.Join(srcDir, relPath)
		if err := os.WriteFile(path, []byte(content), 0644); err != nil {
			t.Fatalf("failed to create file %s: %v", relPath, err)
		}
	}

	// Scan
	store, err := objectstore.NewStore(srcDir)
	if err != nil {
		t.Fatalf("failed to create store: %v", err)
	}
	scanner := NewScanner(store)

	rootObj, err := scanner.ScanDirectory(srcDir)
	if err != nil {
		t.Fatalf("ScanDirectory failed: %v", err)
	}

	rootHash, err := store.PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	// Restore
	dstDir := t.TempDir()
	if err := scanner.RestoreDirectory(rootHash, dstDir); err != nil {
		t.Fatalf("RestoreDirectory failed: %v", err)
	}

	// Verify all files
	for relPath, expectedContent := range files {
		restoredPath := filepath.Join(dstDir, relPath)
		data, err := os.ReadFile(restoredPath)
		if err != nil {
			t.Errorf("failed to read %s: %v", relPath, err)
			continue
		}
		if string(data) != expectedContent {
			t.Errorf("%s: content mismatch: expected %q, got %q",
				relPath, expectedContent, string(data))
		}
	}
}
