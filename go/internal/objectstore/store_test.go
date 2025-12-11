package objectstore

import (
	"io/ioutil"
	"os"
	"path/filepath"
	"sync"
	"testing"

	"github.com/lishi/chemvcs/internal/model"
)

func TestNewStore(t *testing.T) {
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, "test-repo")

	store, err := NewStore(repoPath)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	// Verify objects directory was created
	objectsPath := filepath.Join(repoPath, ".chemvcs", "objects")
	if _, err := os.Stat(objectsPath); os.IsNotExist(err) {
		t.Errorf("objects directory not created at %s", objectsPath)
	}

	if store.Root() != objectsPath {
		t.Errorf("expected root %s, got %s", objectsPath, store.Root())
	}
}

func TestPutGetBlob(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("hello world")

	// Put blob
	hash, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("PutBlob failed: %v", err)
	}

	if hash == "" {
		t.Error("hash should not be empty")
	}

	// Get blob
	retrieved, err := store.GetBlob(hash)
	if err != nil {
		t.Fatalf("GetBlob failed: %v", err)
	}

	if string(retrieved) != string(data) {
		t.Errorf("expected data %q, got %q", data, retrieved)
	}
}

func TestBlobIdempotency(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("test data")

	// Put blob twice
	hash1, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("first PutBlob failed: %v", err)
	}

	hash2, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("second PutBlob failed: %v", err)
	}

	if hash1 != hash2 {
		t.Errorf("hashes should be equal: %s vs %s", hash1, hash2)
	}
}

func TestBlobNotFound(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	// Try to get non-existent blob
	_, err = store.GetBlob("nonexistent1234567890")
	if err == nil {
		t.Error("expected error for non-existent blob")
	}
}

func TestBlobIntegrityCheck(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("original data")
	hash, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("PutBlob failed: %v", err)
	}

	// Corrupt the blob file
	path := store.objectPath(hash)
	if err := ioutil.WriteFile(path, []byte("corrupted"), 0644); err != nil {
		t.Fatalf("failed to corrupt blob: %v", err)
	}

	// GetBlob should detect corruption
	_, err = store.GetBlob(hash)
	if err == nil {
		t.Error("expected error for corrupted blob")
	}
}

func TestPutGetObject(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	obj := model.NewObject("test")
	obj.SetMeta("key", "value")
	obj.AddRef(model.NewBlobRef("abc123"))

	// Put object
	hash, err := store.PutObject(obj)
	if err != nil {
		t.Fatalf("PutObject failed: %v", err)
	}

	// Get object
	retrieved, err := store.GetObject(hash)
	if err != nil {
		t.Fatalf("GetObject failed: %v", err)
	}

	if retrieved.Type != "test" {
		t.Errorf("expected type 'test', got %q", retrieved.Type)
	}

	if val, ok := retrieved.GetMeta("key"); !ok || val != "value" {
		t.Errorf("expected meta key='value', got %v", val)
	}

	if len(retrieved.Refs) != 1 {
		t.Errorf("expected 1 ref, got %d", len(retrieved.Refs))
	}
}

func TestObjectIdempotency(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	obj := model.NewObject("test")
	obj.SetMeta("foo", "bar")

	// Put object twice
	hash1, err := store.PutObject(obj)
	if err != nil {
		t.Fatalf("first PutObject failed: %v", err)
	}

	hash2, err := store.PutObject(obj)
	if err != nil {
		t.Fatalf("second PutObject failed: %v", err)
	}

	if hash1 != hash2 {
		t.Errorf("hashes should be equal: %s vs %s", hash1, hash2)
	}
}

func TestPutGetSnapshot(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	snap := model.NewSnapshot("root123", "Alice <alice@example.com>", "Test commit", nil)
	snap.Timestamp = "2025-12-11T12:00:00Z" // Fixed for testing

	// Put snapshot
	hash, err := store.PutSnapshot(snap)
	if err != nil {
		t.Fatalf("PutSnapshot failed: %v", err)
	}

	// Get snapshot
	retrieved, err := store.GetSnapshot(hash)
	if err != nil {
		t.Fatalf("GetSnapshot failed: %v", err)
	}

	if retrieved.Root != "root123" {
		t.Errorf("expected root 'root123', got %q", retrieved.Root)
	}

	if retrieved.Message != "Test commit" {
		t.Errorf("expected message 'Test commit', got %q", retrieved.Message)
	}

	if retrieved.Author != "Alice <alice@example.com>" {
		t.Errorf("expected author 'Alice <alice@example.com>', got %q", retrieved.Author)
	}
}

func TestSnapshotWithParents(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	parents := []string{"parent1", "parent2"}
	snap := model.NewSnapshot("root456", "Bob <bob@example.com>", "Second commit", parents)

	hash, err := store.PutSnapshot(snap)
	if err != nil {
		t.Fatalf("PutSnapshot failed: %v", err)
	}

	retrieved, err := store.GetSnapshot(hash)
	if err != nil {
		t.Fatalf("GetSnapshot failed: %v", err)
	}

	if len(retrieved.Parents) != 2 {
		t.Errorf("expected 2 parents, got %d", len(retrieved.Parents))
	}

	if retrieved.Parents[0] != "parent1" || retrieved.Parents[1] != "parent2" {
		t.Errorf("parent hashes don't match")
	}
}

func TestShardedLayout(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("test sharding")
	hash, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("PutBlob failed: %v", err)
	}

	// Verify sharded path exists
	expectedPath := filepath.Join(store.Root(), hash[:2], hash)
	if _, err := os.Stat(expectedPath); os.IsNotExist(err) {
		t.Errorf("sharded path not created: %s", expectedPath)
	}
}

func TestConcurrentWrites(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("concurrent test")

	var wg sync.WaitGroup
	errors := make(chan error, 10)

	// Write same blob concurrently
	for i := 0; i < 10; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			_, err := store.PutBlob(data)
			if err != nil {
				errors <- err
			}
		}()
	}

	wg.Wait()
	close(errors)

	// Check for errors
	for err := range errors {
		t.Errorf("concurrent write error: %v", err)
	}

	// Verify blob exists and is correct
	hash := model.ComputeBlobHash(data)
	retrieved, err := store.GetBlob(hash)
	if err != nil {
		t.Fatalf("GetBlob after concurrent writes failed: %v", err)
	}

	if string(retrieved) != string(data) {
		t.Error("blob corrupted after concurrent writes")
	}
}

func TestHasBlob(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	data := []byte("test existence")
	hash, err := store.PutBlob(data)
	if err != nil {
		t.Fatalf("PutBlob failed: %v", err)
	}

	// Should exist
	exists, err := store.HasBlob(hash)
	if err != nil {
		t.Fatalf("HasBlob failed: %v", err)
	}
	if !exists {
		t.Error("blob should exist")
	}

	// Non-existent blob
	exists, err = store.HasBlob("nonexistent")
	if err != nil {
		t.Fatalf("HasBlob failed: %v", err)
	}
	if exists {
		t.Error("non-existent blob should not exist")
	}
}

func TestHasObject(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	obj := model.NewObject("test")
	hash, err := store.PutObject(obj)
	if err != nil {
		t.Fatalf("PutObject failed: %v", err)
	}

	exists, err := store.HasObject(hash)
	if err != nil {
		t.Fatalf("HasObject failed: %v", err)
	}
	if !exists {
		t.Error("object should exist")
	}
}

func TestHasSnapshot(t *testing.T) {
	tmpDir := t.TempDir()
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("NewStore failed: %v", err)
	}

	snap := model.NewSnapshot("root", "Author", "Message", nil)
	hash, err := store.PutSnapshot(snap)
	if err != nil {
		t.Fatalf("PutSnapshot failed: %v", err)
	}

	exists, err := store.HasSnapshot(hash)
	if err != nil {
		t.Fatalf("HasSnapshot failed: %v", err)
	}
	if !exists {
		t.Error("snapshot should exist")
	}
}
