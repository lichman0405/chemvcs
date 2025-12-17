package objectstore

import (
	"crypto/sha256"
	"encoding/hex"
	"os"
	"path/filepath"
	"testing"
)

// Helper function to compute hash
func computeHash(data []byte) string {
	hash := sha256.Sum256(data)
	return hex.EncodeToString(hash[:])
}

// Helper function to put an object with a specific hash (for testing)
func putObjectWithHashFsck(store *Store, hash string, data []byte) error {
	shard := hash[:2]
	shardPath := filepath.Join(store.root, shard)
	if err := os.MkdirAll(shardPath, 0755); err != nil {
		return err
	}
	path := filepath.Join(shardPath, hash)
	return os.WriteFile(path, data, 0644)
}

func TestFsckBasic(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-fsck-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create valid objects
	hash1 := computeHash([]byte("test data 1"))
	hash2 := computeHash([]byte("test data 2"))

	putObjectWithHashFsck(store, hash1, []byte("test data 1"))
	putObjectWithHashFsck(store, hash2, []byte("test data 2"))

	refs := &MockRefsLister{refs: []string{hash1, hash2}}

	opts := FsckOptions{
		Full:    false,
		Verbose: false,
	}

	report, err := store.Fsck(refs, opts)
	if err != nil {
		t.Fatalf("Fsck failed: %v", err)
	}

	if report.ObjectsChecked != 2 {
		t.Errorf("Expected 2 objects checked, got %d", report.ObjectsChecked)
	}

	if len(report.Issues) != 0 {
		t.Errorf("Expected no issues, got %d", len(report.Issues))
	}

	if report.HasErrors {
		t.Errorf("Unexpected errors in clean repository")
	}
}

func TestFsckCorruptedObject(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-fsck-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create object with correct hash
	data := []byte("original data")
	hash := computeHash(data)
	putObjectWithHashFsck(store, hash, data)

	// Corrupt the object by writing wrong data
	path := store.objectPath(hash)
	os.WriteFile(path, []byte("corrupted data"), 0644)

	refs := &MockRefsLister{refs: []string{hash}}

	opts := FsckOptions{
		Full:    false,
		Verbose: false,
	}

	report, err := store.Fsck(refs, opts)
	if err != nil {
		t.Fatalf("Fsck failed: %v", err)
	}

	if !report.HasErrors {
		t.Errorf("Expected errors for corrupted object")
	}

	// Should have at least one error issue
	hasHashMismatch := false
	for _, issue := range report.Issues {
		if issue.Severity == "error" && issue.Object == hash {
			hasHashMismatch = true
			break
		}
	}

	if !hasHashMismatch {
		t.Errorf("Expected hash mismatch error for corrupted object")
	}
}

func TestFsckMissingRef(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-fsck-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create a ref pointing to non-existent object
	missingHash := "9999000000000000000000000000000000000000"
	refs := &MockRefsLister{refs: []string{missingHash}}

	opts := FsckOptions{
		Full:    false,
		Verbose: false,
	}

	report, err := store.Fsck(refs, opts)
	if err != nil {
		t.Fatalf("Fsck failed: %v", err)
	}

	if !report.HasErrors {
		t.Errorf("Expected error for missing ref")
	}

	// Should have error about missing object
	hasMissingRef := false
	for _, issue := range report.Issues {
		if issue.Severity == "error" && issue.Object == missingHash {
			hasMissingRef = true
			break
		}
	}

	if !hasMissingRef {
		t.Errorf("Expected missing ref error")
	}
}

func TestFsckWithPackfiles(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-fsck-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create a packfile
	pw, err := NewPackWriter(store.root, "test-pack")
	if err != nil {
		t.Fatalf("Failed to create pack writer: %v", err)
	}

	data1 := []byte("packed data 1")
	data2 := []byte("packed data 2")
	hash1 := computeHash(data1)
	hash2 := computeHash(data2)

	pw.AddObject(hash1, objTypeBlob, data1)
	pw.AddObject(hash2, objTypeBlob, data2)
	pw.Finalize()

	// Reload packs
	if err := store.ReloadPacks(); err != nil {
		t.Fatalf("Failed to reload packs: %v", err)
	}

	refs := &MockRefsLister{refs: []string{hash1, hash2}}

	opts := FsckOptions{
		Full:    true, // Enable pack checking
		Verbose: false,
	}

	report, err := store.Fsck(refs, opts)
	if err != nil {
		t.Fatalf("Fsck failed: %v", err)
	}

	if report.PacksChecked != 1 {
		t.Errorf("Expected 1 pack checked, got %d", report.PacksChecked)
	}

	if len(report.Issues) != 0 {
		t.Errorf("Expected no issues, got %d: %v", len(report.Issues), report.Issues)
	}
}

func TestFsckCorruptedPackfile(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-fsck-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create a packfile
	pw, err := NewPackWriter(store.root, "corrupt-pack")
	if err != nil {
		t.Fatalf("Failed to create pack writer: %v", err)
	}

	data := []byte("test data")
	hash := computeHash(data)
	pw.AddObject(hash, objTypeBlob, data)
	pw.Finalize()

	// Corrupt the packfile by appending garbage
	packPath := filepath.Join(store.root, "pack", "corrupt-pack.pack")
	f, err := os.OpenFile(packPath, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		t.Fatalf("Failed to open pack: %v", err)
	}
	f.WriteString("GARBAGE")
	f.Close()

	// Reload packs
	if err := store.ReloadPacks(); err != nil {
		t.Fatalf("Failed to reload packs: %v", err)
	}

	opts := FsckOptions{
		Full:    true,
		Verbose: false,
	}

	report, err := store.Fsck(nil, opts)
	if err != nil {
		t.Fatalf("Fsck failed: %v", err)
	}

	// Should detect checksum error
	if !report.HasErrors {
		t.Errorf("Expected error for corrupted packfile")
	}

	hasChecksumError := false
	for _, issue := range report.Issues {
		if issue.Severity == "error" && issue.Object == "corrupt-pack.pack" {
			hasChecksumError = true
			break
		}
	}

	if !hasChecksumError {
		t.Errorf("Expected checksum error for corrupted pack")
	}
}
