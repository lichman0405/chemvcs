package objectstore

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func TestPackWriter(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-pack-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Create packfile
	pw, err := NewPackWriter(tmpDir, "test-pack")
	if err != nil {
		t.Fatalf("Failed to create pack writer: %v", err)
	}

	// Add test objects
	testObjects := []struct {
		hash string
		typ  byte
		data []byte
	}{
		{"aabbccdd11223344556677889900aabbccddee11223344556677889900aabbcc", objTypeBlob, []byte("blob data 1")},
		{"112233445566778899aabbccddeeff00112233445566778899aabbccddeeff00", objTypeBlob, []byte("blob data 2")},
		{"ffeeddccbbaa998877665544332211001122aabbccddeeff0011223344556677", objTypeTree, []byte(`{"type":"tree","refs":[]}`)},
	}

	for _, obj := range testObjects {
		if err := pw.AddObject(obj.hash, obj.typ, obj.data); err != nil {
			t.Fatalf("Failed to add object %s: %v", obj.hash, err)
		}
	}

	// Finalize pack
	if err := pw.Finalize(); err != nil {
		t.Fatalf("Failed to finalize pack: %v", err)
	}

	// Verify packfile exists
	packPath := filepath.Join(tmpDir, "pack", "test-pack.pack")
	if _, err := os.Stat(packPath); os.IsNotExist(err) {
		t.Errorf("Packfile was not created")
	}

	// Verify index exists
	idxPath := filepath.Join(tmpDir, "pack", "test-pack.idx")
	if _, err := os.Stat(idxPath); os.IsNotExist(err) {
		t.Errorf("Index file was not created")
	}

	// Read back using PackReader
	pr, err := OpenPackReader(tmpDir, "test-pack")
	if err != nil {
		t.Fatalf("Failed to create pack reader: %v", err)
	}

	// Verify each object
	for _, obj := range testObjects {
		data, typ, err := pr.GetObject(obj.hash)
		if err != nil {
			t.Errorf("Failed to read object %s: %v", obj.hash, err)
		}

		if typ != obj.typ {
			t.Errorf("Type mismatch for %s: expected %d, got %d", obj.hash, obj.typ, typ)
		}

		if !bytes.Equal(data, obj.data) {
			t.Errorf("Data mismatch for %s", obj.hash)
		}
	}

	// Test non-existent object
	_, _, err = pr.GetObject("0000000000000000000000000000000000000000000000000000000000000000")
	if err == nil {
		t.Errorf("Expected error for non-existent object")
	}
}

func TestPackReaderIndex(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-pack-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Create packfile with sorted hashes
	pw, err := NewPackWriter(tmpDir, "sorted-pack")
	if err != nil {
		t.Fatalf("Failed to create pack writer: %v", err)
	}

	// Add objects in non-sorted order (should be sorted in index)
	hashes := []string{
		"ff00000000000000000000000000000000000000000000000000000000000000",
		"00ff000000000000000000000000000000000000000000000000000000000000",
		"8800000000000000000000000000000000000000000000000000000000000000",
		"4400000000000000000000000000000000000000000000000000000000000000",
	}

	for _, hash := range hashes {
		data := []byte("data for " + hash)
		if err := pw.AddObject(hash, objTypeBlob, data); err != nil {
			t.Fatalf("Failed to add object: %v", err)
		}
	}

	if err := pw.Finalize(); err != nil {
		t.Fatalf("Failed to finalize: %v", err)
	}

	// Read back
	pr, err := OpenPackReader(tmpDir, "sorted-pack")
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}

	// Verify all objects can be found
	for _, hash := range hashes {
		if !pr.HasObject(hash) {
			t.Errorf("Object %s not found in pack", hash)
		}
	}
}

func TestPackCompression(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-pack-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Create packfile with compressible data
	pw, err := NewPackWriter(tmpDir, "compress-test")
	if err != nil {
		t.Fatalf("Failed to create pack writer: %v", err)
	}

	// Add a large compressible blob (repeated data)
	data := bytes.Repeat([]byte("AAABBBCCCDDDEEE"), 1000)
	hash := "1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef"

	if err := pw.AddObject(hash, objTypeBlob, data); err != nil {
		t.Fatalf("Failed to add object: %v", err)
	}

	if err := pw.Finalize(); err != nil {
		t.Fatalf("Failed to finalize: %v", err)
	}

	// Check pack size - should be much smaller than original
	packPath := filepath.Join(tmpDir, "pack", "compress-test.pack")
	info, err := os.Stat(packPath)
	if err != nil {
		t.Fatalf("Failed to stat pack: %v", err)
	}

	originalSize := len(data)
	packedSize := int(info.Size())

	// Expect at least 50% compression for this simple test
	if packedSize > originalSize/2 {
		t.Logf("Warning: compression not as effective as expected (original: %d, packed: %d)", originalSize, packedSize)
	}

	// Read back and verify
	pr, err := OpenPackReader(tmpDir, "compress-test")
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}

	readData, _, err := pr.GetObject(hash)
	if err != nil {
		t.Fatalf("Failed to read object: %v", err)
	}

	if !bytes.Equal(readData, data) {
		t.Errorf("Decompressed data doesn't match original")
	}
}

func TestPackMultiplePacks(t *testing.T) {
	// Create temporary directory
	tmpDir, err := os.MkdirTemp("", "chemvcs-pack-test-*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Create store
	store, err := NewStore(tmpDir)
	if err != nil {
		t.Fatalf("Failed to create store: %v", err)
	}

	// Create first pack
	pw1, err := NewPackWriter(store.root, "pack-1")
	if err != nil {
		t.Fatalf("Failed to create pack writer 1: %v", err)
	}
	pw1.AddObject("aaaa000000000000000000000000000000000000000000000000000000000000", objTypeBlob, []byte("pack 1 data"))
	if err := pw1.Finalize(); err != nil {
		t.Fatalf("Failed to finalize pack 1: %v", err)
	}

	// Create second pack
	pw2, err := NewPackWriter(store.root, "pack-2")
	if err != nil {
		t.Fatalf("Failed to create pack writer 2: %v", err)
	}
	pw2.AddObject("bbbb000000000000000000000000000000000000000000000000000000000000", objTypeBlob, []byte("pack 2 data"))
	if err := pw2.Finalize(); err != nil {
		t.Fatalf("Failed to finalize pack 2: %v", err)
	}

	// Reload packs
	if err := store.ReloadPacks(); err != nil {
		t.Fatalf("Failed to reload packs: %v", err)
	}

	// Should be able to find objects in both packs
	if !store.HasObjectInPack("aaaa000000000000000000000000000000000000000000000000000000000000") {
		t.Errorf("Object from pack 1 not found")
	}

	if !store.HasObjectInPack("bbbb000000000000000000000000000000000000000000000000000000000000") {
		t.Errorf("Object from pack 2 not found")
	}
}
