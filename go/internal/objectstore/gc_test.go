package objectstore

import (
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/lishi/chemvcs/internal/model"
)

// MockRefsLister for testing
type MockRefsLister struct {
	refs []string
}

func (m *MockRefsLister) ListAllRefs() ([]string, error) {
	return m.refs, nil
}

// Helper function to put an object with a specific hash (for testing)
func putObjectWithHash(store *Store, hash string, data []byte) error {
	shard := hash[:2]
	shardPath := filepath.Join(store.root, shard)
	if err := os.MkdirAll(shardPath, 0755); err != nil {
		return err
	}
	path := filepath.Join(shardPath, hash)
	return os.WriteFile(path, data, 0644)
}

// Helper function to put a model.Object with a specific hash (for testing)
func putModelObjectWithHash(store *Store, hash string, obj *model.Object) error {
	data, err := obj.MarshalCanonical()
	if err != nil {
		return err
	}
	return putObjectWithHash(store, hash, data)
}

func TestGarbageCollectBasic(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-gc-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create some objects using real PutBlob (which computes correct hashes)
	reachableHash, _ := store.PutBlob([]byte("reachable blob"))
	unreachableHash, _ := store.PutBlob([]byte("unreachable blob"))

	// Only first object is reachable
	refs := &MockRefsLister{refs: []string{reachableHash}}

	// Run GC in dry-run mode
	opts := GCOptions{
		PruneAge: 0,
		DryRun:   true,
		Verbose:  false,
	}

	stats, err := store.GarbageCollect(refs, opts)
	if err != nil {
		t.Fatalf("GC failed: %v", err)
	}

	if stats.ReachableObjects != 1 {
		t.Errorf("Expected 1 reachable object, got %d", stats.ReachableObjects)
	}

	if stats.UnreachableObjects != 1 {
		t.Errorf("Expected 1 unreachable object, got %d", stats.UnreachableObjects)
	}

	// In dry-run mode, nothing should be deleted
	if _, err := store.GetBlob(unreachableHash); err != nil {
		t.Errorf("Unreachable object was deleted in dry-run mode")
	}

	// Run actual GC
	opts.DryRun = false
	stats, err = store.GarbageCollect(refs, opts)
	if err != nil {
		t.Fatalf("GC failed: %v", err)
	}

	if stats.RemovedObjects != 1 {
		t.Errorf("Expected 1 removed object, got %d", stats.RemovedObjects)
	}

	// Now unreachable object should be gone
	if _, err := store.GetBlob(unreachableHash); err == nil {
		t.Errorf("Unreachable object was not deleted")
	}

	// Reachable object should still exist
	if _, err := store.GetBlob(reachableHash); err != nil {
		t.Errorf("Reachable object was deleted: %v", err)
	}
}

func TestGarbageCollectGracePeriod(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-gc-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create unreachable object
	hash := "cccc000000000000000000000000000000000000"
	putObjectWithHash(store, hash, []byte("test blob"))

	refs := &MockRefsLister{refs: []string{}}

	// Run GC with 1 hour grace period - should not delete
	opts := GCOptions{
		PruneAge: 1 * time.Hour,
		DryRun:   false,
		Verbose:  false,
	}

	stats, err := store.GarbageCollect(refs, opts)
	if err != nil {
		t.Fatalf("GC failed: %v", err)
	}

	if stats.RemovedObjects != 0 {
		t.Errorf("Object was removed despite grace period")
	}

	// Object should still exist
	if _, err := store.GetBlob(hash); err != nil {
		t.Errorf("Object was deleted despite grace period")
	}

	// Run GC with immediate pruning
	opts.PruneAge = 0
	stats, err = store.GarbageCollect(refs, opts)
	if err != nil {
		t.Fatalf("GC failed: %v", err)
	}

	if stats.RemovedObjects != 1 {
		t.Errorf("Expected 1 removed object with prune=now")
	}
}

func TestGarbageCollectWithGraph(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-gc-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create object graph:
	// root -> child1 -> grandchild
	//      -> child2
	// orphan (unreachable)

	grandchildHash := "1111000000000000000000000000000000000000"
	child1Hash := "2222000000000000000000000000000000000000"
	child2Hash := "3333000000000000000000000000000000000000"
	rootHash := "4444000000000000000000000000000000000000"
	orphanHash := "5555000000000000000000000000000000000000"

	putObjectWithHash(store, grandchildHash, []byte("grandchild"))
	putModelObjectWithHash(store, child1Hash, &model.Object{
		Version: 1,
		Type:    "tree",
		Refs:    []model.Reference{{ID: grandchildHash, Kind: "blob"}},
	})
	putObjectWithHash(store, child2Hash, []byte("child2"))
	putModelObjectWithHash(store, rootHash, &model.Object{
		Version: 1,
		Type:    "tree",
		Refs: []model.Reference{
			{ID: child1Hash, Kind: "object"},
			{ID: child2Hash, Kind: "blob"},
		},
	})
	putObjectWithHash(store, orphanHash, []byte("orphan"))

	// Only root is directly referenced
	refs := &MockRefsLister{refs: []string{rootHash}}

	opts := GCOptions{
		PruneAge: 0,
		DryRun:   false,
		Verbose:  false,
	}

	stats, err := store.GarbageCollect(refs, opts)
	if err != nil {
		t.Fatalf("GC failed: %v", err)
	}

	// Should find 4 reachable (root + child1 + child2 + grandchild)
	if stats.ReachableObjects != 4 {
		t.Errorf("Expected 4 reachable objects, got %d", stats.ReachableObjects)
	}

	// Should find 1 unreachable (orphan)
	if stats.UnreachableObjects != 1 {
		t.Errorf("Expected 1 unreachable object, got %d", stats.UnreachableObjects)
	}

	// Should remove orphan
	if stats.RemovedObjects != 1 {
		t.Errorf("Expected 1 removed object, got %d", stats.RemovedObjects)
	}

	// Verify orphan is gone
	if _, err := store.GetBlob(orphanHash); err == nil {
		t.Errorf("Orphan was not deleted")
	}

	// Verify reachable objects still exist
	reachable := []string{rootHash, child1Hash, child2Hash, grandchildHash}
	for _, hash := range reachable {
		if _, err := store.GetBlob(hash); err != nil {
			if _, err := store.GetObject(hash); err != nil {
				t.Errorf("Reachable object %s was deleted", hash)
			}
		}
	}
}

func TestPackLooseObjects(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-pack-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create some loose objects
	hash1 := "dddd000000000000000000000000000000000000"
	hash2 := "eeee000000000000000000000000000000000000"
	hash3 := "ffff000000000000000000000000000000000000"

	putObjectWithHash(store, hash1, []byte("blob 1"))
	putObjectWithHash(store, hash2, []byte("blob 2"))
	putObjectWithHash(store, hash3, []byte("blob 3"))

	// All objects are reachable
	refs := &MockRefsLister{refs: []string{hash1, hash2, hash3}}

	opts := PackOptions{
		All:       false, // Only pack reachable
		KeepLoose: true,  // Keep loose for verification
		Verbose:   false,
	}

	packName, stats, err := store.PackLooseObjects(refs, opts)
	if err != nil {
		t.Fatalf("Pack failed: %v", err)
	}

	if stats.ObjectsPacked != 3 {
		t.Errorf("Expected 3 packed objects, got %d", stats.ObjectsPacked)
	}

	if packName == "" {
		t.Errorf("Pack name is empty")
	}

	// Verify pack file exists
	packPath := filepath.Join(tmpDir, "pack", packName+".pack")
	if _, err := os.Stat(packPath); os.IsNotExist(err) {
		t.Errorf("Pack file not created")
	}

	// Verify we can read from pack
	if !store.HasObjectInPack(hash1) {
		t.Errorf("Object not found in pack")
	}

	// If KeepLoose=false, loose objects should be removed
	opts.KeepLoose = false
	hash4 := "0000111111111111111111111111111111111111"
	putObjectWithHash(store, hash4, []byte("blob 4"))
	refs.refs = append(refs.refs, hash4)

	_, _, err = store.PackLooseObjects(refs, opts)
	if err != nil {
		t.Fatalf("Second pack failed: %v", err)
	}

	// Check that loose object was removed
	loosePath := store.objectPath(hash4)
	if _, err := os.Stat(loosePath); !os.IsNotExist(err) {
		t.Errorf("Loose object was not removed after packing")
	}

	// But we can still read from pack
	if _, err := store.GetBlob(hash4); err != nil {
		t.Errorf("Cannot read packed object: %v", err)
	}
}
