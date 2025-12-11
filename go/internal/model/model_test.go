package model

import (
	"testing"
)

func TestObjectHashDeterminism(t *testing.T) {
	// Create two objects with same content but different construction order
	obj1 := NewObject("test")
	obj1.SetMeta("foo", "bar")
	obj1.SetMeta("baz", 42)
	obj1.AddRef(NewBlobRef("abc123"))

	obj2 := NewObject("test")
	obj2.SetMeta("baz", 42)
	obj2.SetMeta("foo", "bar")
	obj2.AddRef(NewBlobRef("abc123"))

	hash1, err := obj1.Hash()
	if err != nil {
		t.Fatalf("failed to hash obj1: %v", err)
	}

	hash2, err := obj2.Hash()
	if err != nil {
		t.Fatalf("failed to hash obj2: %v", err)
	}

	if hash1 != hash2 {
		t.Errorf("hashes should be equal but got:\nobj1: %s\nobj2: %s", hash1, hash2)
	}
}

func TestObjectEmptyMeta(t *testing.T) {
	obj := NewObject("test")
	// No meta or refs added

	hash, err := obj.Hash()
	if err != nil {
		t.Fatalf("failed to hash object: %v", err)
	}

	if hash == "" {
		t.Error("hash should not be empty")
	}

	// Should be repeatable
	hash2, err := obj.Hash()
	if err != nil {
		t.Fatalf("failed to hash object second time: %v", err)
	}

	if hash != hash2 {
		t.Error("hash should be consistent")
	}
}

func TestObjectRefs(t *testing.T) {
	obj := NewObject("folder")

	obj.AddRef(NewObjectRef("hash1"))
	obj.AddRef(NewBlobRef("hash2"))

	if len(obj.Refs) != 2 {
		t.Errorf("expected 2 refs, got %d", len(obj.Refs))
	}

	if !obj.Refs[0].IsObject() {
		t.Error("first ref should be object type")
	}

	if !obj.Refs[1].IsBlob() {
		t.Error("second ref should be blob type")
	}
}

func TestSnapshotHash(t *testing.T) {
	snap1 := NewSnapshot("root123", "Alice <alice@example.com>", "Initial commit", nil)
	snap1.Timestamp = "2025-12-11T12:00:00Z" // Fixed timestamp for testing

	hash1, err := snap1.Hash()
	if err != nil {
		t.Fatalf("failed to hash snapshot: %v", err)
	}

	if hash1 == "" {
		t.Error("hash should not be empty")
	}

	// Create identical snapshot
	snap2 := NewSnapshot("root123", "Alice <alice@example.com>", "Initial commit", nil)
	snap2.Timestamp = "2025-12-11T12:00:00Z"

	hash2, err := snap2.Hash()
	if err != nil {
		t.Fatalf("failed to hash snapshot: %v", err)
	}

	if hash1 != hash2 {
		t.Error("identical snapshots should have same hash")
	}
}

func TestSnapshotWithParents(t *testing.T) {
	parents := []string{"parent1", "parent2"}
	snap := NewSnapshot("root456", "Bob <bob@example.com>", "Second commit", parents)

	if len(snap.Parents) != 2 {
		t.Errorf("expected 2 parents, got %d", len(snap.Parents))
	}

	if snap.IsInitial() {
		t.Error("snapshot with parents should not be initial")
	}
}

func TestSnapshotIsInitial(t *testing.T) {
	snap := NewSnapshot("root789", "Charlie <charlie@example.com>", "First", nil)

	if !snap.IsInitial() {
		t.Error("snapshot with no parents should be initial")
	}
}

func TestBlobHash(t *testing.T) {
	data := []byte("hello world")
	blob := NewBlob(data)

	hash1 := blob.Hash()
	if hash1 == "" {
		t.Error("hash should not be empty")
	}

	// Create another blob with same data
	blob2 := NewBlob(data)
	hash2 := blob2.Hash()

	if hash1 != hash2 {
		t.Error("blobs with same data should have same hash")
	}

	// Test utility function
	hash3 := ComputeBlobHash(data)
	if hash1 != hash3 {
		t.Error("ComputeBlobHash should match Blob.Hash()")
	}
}

func TestBlobSize(t *testing.T) {
	data := []byte("test data")
	blob := NewBlob(data)

	if blob.Size() != int64(len(data)) {
		t.Errorf("expected size %d, got %d", len(data), blob.Size())
	}
}

func TestReferenceTypes(t *testing.T) {
	objRef := NewObjectRef("abc123")
	if !objRef.IsObject() {
		t.Error("object ref should return true for IsObject()")
	}
	if objRef.IsBlob() {
		t.Error("object ref should return false for IsBlob()")
	}

	blobRef := NewBlobRef("def456")
	if blobRef.IsObject() {
		t.Error("blob ref should return false for IsObject()")
	}
	if !blobRef.IsBlob() {
		t.Error("blob ref should return true for IsBlob()")
	}
}
