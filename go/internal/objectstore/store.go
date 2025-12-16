package objectstore

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/lishi/chemvcs/internal/model"
)

// Store provides content-addressable storage for blobs, objects, and snapshots.
// It implements a Git-style sharded directory layout under .chemvcs/objects/.
// It supports both loose objects and packfiles for efficient storage.
type Store struct {
	root  string        // Path to .chemvcs/objects directory
	mu    sync.RWMutex  // Protects concurrent writes
	packs []*PackReader // Loaded packfiles
}

// NewStore creates a new Store rooted at the given repository path.
func NewStore(repoPath string) (*Store, error) {
	objectsPath := filepath.Join(repoPath, ".chemvcs", "objects")

	// Ensure objects directory exists
	if err := os.MkdirAll(objectsPath, 0755); err != nil {
		return nil, fmt.Errorf("failed to create objects directory: %w", err)
	}

	s := &Store{
		root:  objectsPath,
		packs: []*PackReader{},
	}

	// Load existing packfiles
	if err := s.loadPacks(); err != nil {
		// Non-fatal: log but continue
		fmt.Fprintf(os.Stderr, "warning: failed to load some packfiles: %v\n", err)
	}

	return s, nil
}

// PutBlob stores a blob and returns its hash.
// If the blob already exists, it returns the hash without error (idempotent).
func (s *Store) PutBlob(data []byte) (string, error) {
	hash := model.ComputeBlobHash(data)
	path := s.objectPath(hash)

	// Check if already exists (idempotency)
	if exists, err := s.exists(path); err != nil {
		return "", err
	} else if exists {
		// Verify existing content matches
		existing, err := os.ReadFile(path)
		if err != nil {
			return "", fmt.Errorf("failed to read existing blob: %w", err)
		}
		if model.ComputeBlobHash(existing) != hash {
			return "", fmt.Errorf("hash collision detected for %s", hash)
		}
		return hash, nil
	}

	// Atomic write
	if err := s.atomicWrite(path, data); err != nil {
		return "", fmt.Errorf("failed to write blob: %w", err)
	}

	return hash, nil
}

// GetBlob retrieves a blob by its hash.
// It verifies the hash matches the content before returning.
// Tries loose objects first, then packfiles.
func (s *Store) GetBlob(hash string) ([]byte, error) {
	path := s.objectPath(hash)

	// Try loose object first
	data, err := os.ReadFile(path)
	if err == nil {
		// Verify hash integrity
		if model.ComputeBlobHash(data) != hash {
			return nil, fmt.Errorf("hash mismatch for blob %s: content corrupted", hash)
		}
		return data, nil
	}

	if !os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to read blob: %w", err)
	}

	// Try packfiles
	s.mu.RLock()
	defer s.mu.RUnlock()
	for _, pack := range s.packs {
		if packData, objType, err := pack.GetObject(hash); err == nil {
			if objType != objTypeBlob {
				return nil, fmt.Errorf("object %s is not a blob", hash)
			}
			// Verify hash
			if model.ComputeBlobHash(packData) != hash {
				return nil, fmt.Errorf("hash mismatch for packed blob %s", hash)
			}
			return packData, nil
		}
	}

	return nil, fmt.Errorf("blob not found: %s", hash)
}

// HasBlob checks if a blob exists in the store.
func (s *Store) HasBlob(hash string) (bool, error) {
	path := s.objectPath(hash)
	return s.exists(path)
}

// PutObject stores an object and returns its hash.
// If the object already exists, it returns the hash without error (idempotent).
func (s *Store) PutObject(obj *model.Object) (string, error) {
	hash, err := obj.Hash()
	if err != nil {
		return "", fmt.Errorf("failed to compute object hash: %w", err)
	}

	path := s.objectPath(hash)

	// Check if already exists (idempotency)
	if exists, err := s.exists(path); err != nil {
		return "", err
	} else if exists {
		// Verify existing content matches
		existing, err := s.GetObject(hash)
		if err != nil {
			return "", fmt.Errorf("failed to verify existing object: %w", err)
		}
		existingHash, err := existing.Hash()
		if err != nil || existingHash != hash {
			return "", fmt.Errorf("hash collision detected for %s", hash)
		}
		return hash, nil
	}

	// Serialise object to canonical JSON
	data, err := obj.MarshalCanonical()
	if err != nil {
		return "", fmt.Errorf("failed to marshal object: %w", err)
	}

	// Atomic write
	if err := s.atomicWrite(path, data); err != nil {
		return "", fmt.Errorf("failed to write object: %w", err)
	}

	return hash, nil
}

// GetObject retrieves an object by its hash.
// It verifies the hash matches the content before returning.
// Tries loose objects first, then packfiles.
func (s *Store) GetObject(hash string) (*model.Object, error) {
	path := s.objectPath(hash)

	// Try loose object first
	data, err := os.ReadFile(path)
	if err == nil {
		// Deserialise object
		var obj model.Object
		if err := json.Unmarshal(data, &obj); err != nil {
			return nil, fmt.Errorf("failed to unmarshal object: %w", err)
		}
		// Verify hash integrity
		objHash, err := obj.Hash()
		if err != nil {
			return nil, fmt.Errorf("failed to compute hash for verification: %w", err)
		}
		if objHash != hash {
			return nil, fmt.Errorf("hash mismatch for object %s: content corrupted", hash)
		}
		return &obj, nil
	}

	if !os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to read object: %w", err)
	}

	// Try packfiles
	s.mu.RLock()
	defer s.mu.RUnlock()
	for _, pack := range s.packs {
		if packData, objType, err := pack.GetObject(hash); err == nil {
			if objType != objTypeTree {
				continue // Not a tree object
			}
			var obj model.Object
			if err := json.Unmarshal(packData, &obj); err != nil {
				continue
			}
			// Verify hash
			objHash, err := obj.Hash()
			if err != nil || objHash != hash {
				return nil, fmt.Errorf("hash mismatch for packed object %s", hash)
			}
			return &obj, nil
		}
	}

	return nil, fmt.Errorf("object not found: %s", hash)
}

// getObjectRaw retrieves raw object data (for internal use).
func (s *Store) getObjectRaw(hash string) ([]byte, error) {
	path := s.objectPath(hash)

	// Try loose object first
	data, err := os.ReadFile(path)
	if err == nil {
		return data, nil
	}

	if !os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to read object: %w", err)
	}

	// Try packfiles
	s.mu.RLock()
	defer s.mu.RUnlock()
	for _, pack := range s.packs {
		if packData, _, err := pack.GetObject(hash); err == nil {
			return packData, nil
		}
	}

	return nil, fmt.Errorf("object not found: %s", hash)
}

// oldGetObject is the original GetObject code before pack support.
// Keep signature for testing but now delegates to new implementation.
func oldGetObject(s *Store, hash string) (*model.Object, error) {
	data, err := s.getObjectRaw(hash)
	if err != nil {
		return nil, err
	}

	// Deserialise object
	var obj model.Object
	if err := json.Unmarshal(data, &obj); err != nil {
		return nil, fmt.Errorf("failed to unmarshal object: %w", err)
	}

	// Verify hash integrity
	objHash, err := obj.Hash()
	if err != nil {
		return nil, fmt.Errorf("failed to compute hash for verification: %w", err)
	}
	if objHash != hash {
		return nil, fmt.Errorf("hash mismatch for object %s: content corrupted", hash)
	}

	return &obj, nil
}

// HasObject checks if an object exists in the store.
func (s *Store) HasObject(hash string) (bool, error) {
	path := s.objectPath(hash)
	return s.exists(path)
}

// PutSnapshot stores a snapshot and returns its hash.
// If the snapshot already exists, it returns the hash without error (idempotent).
func (s *Store) PutSnapshot(snap *model.Snapshot) (string, error) {
	hash, err := snap.Hash()
	if err != nil {
		return "", fmt.Errorf("failed to compute snapshot hash: %w", err)
	}

	path := s.objectPath(hash)

	// Check if already exists (idempotency)
	if exists, err := s.exists(path); err != nil {
		return "", err
	} else if exists {
		// Verify existing content matches
		existing, err := s.GetSnapshot(hash)
		if err != nil {
			return "", fmt.Errorf("failed to verify existing snapshot: %w", err)
		}
		existingHash, err := existing.Hash()
		if err != nil || existingHash != hash {
			return "", fmt.Errorf("hash collision detected for %s", hash)
		}
		return hash, nil
	}

	// Serialise snapshot to JSON
	data, err := snap.MarshalCanonical()
	if err != nil {
		return "", fmt.Errorf("failed to marshal snapshot: %w", err)
	}

	// Atomic write
	if err := s.atomicWrite(path, data); err != nil {
		return "", fmt.Errorf("failed to write snapshot: %w", err)
	}

	return hash, nil
}

// GetSnapshot retrieves a snapshot by its hash.
// It verifies the hash matches the content before returning.
// Tries loose objects first, then packfiles.
func (s *Store) GetSnapshot(hash string) (*model.Snapshot, error) {
	path := s.objectPath(hash)

	// Try loose object first
	data, err := os.ReadFile(path)
	if err == nil {
		// Deserialise snapshot
		snap, err := model.UnmarshalSnapshot(data)
		if err != nil {
			return nil, fmt.Errorf("failed to unmarshal snapshot: %w", err)
		}
		// Verify hash integrity
		snapHash, err := snap.Hash()
		if err != nil {
			return nil, fmt.Errorf("failed to compute hash for verification: %w", err)
		}
		if snapHash != hash {
			return nil, fmt.Errorf("hash mismatch for snapshot %s: content corrupted", hash)
		}
		return snap, nil
	}

	if !os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to read snapshot: %w", err)
	}

	// Try packfiles
	s.mu.RLock()
	defer s.mu.RUnlock()
	for _, pack := range s.packs {
		if packData, objType, err := pack.GetObject(hash); err == nil {
			if objType != objTypeSnapshot {
				continue // Not a snapshot
			}
			snap, err := model.UnmarshalSnapshot(packData)
			if err != nil {
				continue
			}
			// Verify hash
			snapHash, err := snap.Hash()
			if err != nil || snapHash != hash {
				return nil, fmt.Errorf("hash mismatch for packed snapshot %s", hash)
			}
			return snap, nil
		}
	}

	return nil, fmt.Errorf("snapshot not found: %s", hash)
}

// HasSnapshot checks if a snapshot exists in the store.
func (s *Store) HasSnapshot(hash string) (bool, error) {
	path := s.objectPath(hash)
	return s.exists(path)
}

// objectPath returns the filesystem path for a given hash.
// Uses Git-style sharding: objects/ab/cdef123456...
func (s *Store) objectPath(hash string) string {
	if len(hash) < 2 {
		// Should never happen with SHA-256, but be defensive
		return filepath.Join(s.root, hash)
	}
	// Shard by first two characters
	return filepath.Join(s.root, hash[:2], hash)
}

// exists checks if a file exists at the given path.
func (s *Store) exists(path string) (bool, error) {
	_, err := os.Stat(path)
	if err == nil {
		return true, nil
	}
	if os.IsNotExist(err) {
		return false, nil
	}
	return false, err
}

// atomicWrite writes data to a file atomically using a temporary file and rename.
// This ensures that the file either fully exists with correct content or doesn't exist at all.
func (s *Store) atomicWrite(path string, data []byte) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	// Ensure parent directory exists
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return fmt.Errorf("failed to create directory: %w", err)
	}

	// Create temporary file in the same directory
	tmpFile, err := os.CreateTemp(dir, ".tmp-*")
	if err != nil {
		return fmt.Errorf("failed to create temp file: %w", err)
	}
	tmpPath := tmpFile.Name()

	// Write data to temp file
	if _, err := tmpFile.Write(data); err != nil {
		tmpFile.Close()
		os.Remove(tmpPath)
		return fmt.Errorf("failed to write to temp file: %w", err)
	}

	// Sync to ensure data is on disk
	if err := tmpFile.Sync(); err != nil {
		tmpFile.Close()
		os.Remove(tmpPath)
		return fmt.Errorf("failed to sync temp file: %w", err)
	}

	tmpFile.Close()

	// Atomic rename
	if err := os.Rename(tmpPath, path); err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("failed to rename temp file: %w", err)
	}

	return nil
}

// Root returns the root path of the object store.
func (s *Store) Root() string {
	return s.root
}

// GetRaw retrieves raw bytes for a given object hash without validation.
// This is used for remote operations where validation happens elsewhere.
func (s *Store) GetRaw(hash string) ([]byte, error) {
	path := s.objectPath(hash)

	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("object not found: %s", hash)
		}
		return nil, fmt.Errorf("failed to read object: %w", err)
	}

	return data, nil
}

// PutRaw stores raw bytes with a given hash without validation.
// This is used for remote operations where validation happens elsewhere.
// The hash must match the content, but this method does not verify it.
func (s *Store) PutRaw(hash string, data []byte) error {
	path := s.objectPath(hash)

	// Check if already exists (idempotency)
	if exists, err := s.exists(path); err != nil {
		return err
	} else if exists {
		return nil
	}

	// Atomic write
	if err := s.atomicWrite(path, data); err != nil {
		return fmt.Errorf("failed to write object: %w", err)
	}

	return nil
}

// ReadFile reads a file relative to the repository root (parent of objects directory).
// This is a helper for reading config files, refs, etc.
func (s *Store) ReadFile(relativePath string) ([]byte, error) {
	// Go up one level from objects/ to .chemvcs/, then up again to repo root
	repoRoot := filepath.Dir(filepath.Dir(s.root))
	fullPath := filepath.Join(repoRoot, relativePath)

	data, err := os.ReadFile(fullPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("file not found: %s", relativePath)
		}
		return nil, fmt.Errorf("failed to read file: %w", err)
	}

	return data, nil
}

// WriteFile writes a file relative to the repository root.
// This is a helper for writing config files, refs, etc.
func (s *Store) WriteFile(relativePath string, data []byte) error {
	// Go up one level from objects/ to .chemvcs/, then up again to repo root
	repoRoot := filepath.Dir(filepath.Dir(s.root))
	fullPath := filepath.Join(repoRoot, relativePath)

	// Ensure parent directory exists
	dir := filepath.Dir(fullPath)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return fmt.Errorf("failed to create directory: %w", err)
	}

	// Write file atomically
	if err := s.atomicWrite(fullPath, data); err != nil {
		return fmt.Errorf("failed to write file: %w", err)
	}

	return nil
}

// ObjectInfo contains metadata about a stored object.
type ObjectInfo struct {
	Hash string
	Type string
}

// ListObjects returns information about all objects in the store.
// If typeFilter is non-empty, only objects of that type are returned.
func (s *Store) ListObjects(typeFilter string) ([]ObjectInfo, error) {
	s.mu.RLock()
	defer s.mu.RUnlock()

	var objects []ObjectInfo

	// Walk all shard directories (00-ff)
	for i := 0; i < 256; i++ {
		shard := fmt.Sprintf("%02x", i)
		shardPath := filepath.Join(s.root, shard)

		// Check if shard directory exists
		if _, err := os.Stat(shardPath); os.IsNotExist(err) {
			continue
		}

		// Read all files in this shard
		entries, err := os.ReadDir(shardPath)
		if err != nil {
			return nil, fmt.Errorf("failed to read shard %s: %w", shard, err)
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}

			// The filename is the full hash (not just the suffix)
			hash := entry.Name()

			// Read object to get its type
			obj, err := s.GetObject(hash)
			if err != nil {
				// Skip objects that can't be read (might be blobs or snapshots)
				continue
			}

			// Apply type filter if specified
			if typeFilter != "" && obj.Type != typeFilter {
				continue
			}

			objects = append(objects, ObjectInfo{
				Hash: hash,
				Type: obj.Type,
			})
		}
	}

	return objects, nil
}

// loadPacks scans the pack directory and loads all packfiles.
func (s *Store) loadPacks() error {
	packDir := filepath.Join(s.root, "pack")
	if _, err := os.Stat(packDir); os.IsNotExist(err) {
		return nil // No pack directory yet
	}

	entries, err := os.ReadDir(packDir)
	if err != nil {
		return fmt.Errorf("failed to read pack directory: %w", err)
	}

	for _, entry := range entries {
		if entry.IsDir() || !strings.HasSuffix(entry.Name(), ".idx") {
			continue
		}

		// Extract pack name (without .idx extension)
		packName := entry.Name()[:len(entry.Name())-4]

		pack, err := OpenPackReader(s.root, packName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "warning: failed to load pack %s: %v\n", packName, err)
			continue
		}

		s.packs = append(s.packs, pack)
	}

	return nil
}

// ReloadPacks reloads all packfiles from disk.
// Call this after creating new packfiles to make them available.
func (s *Store) ReloadPacks() error {
	s.mu.Lock()
	defer s.mu.Unlock()

	s.packs = []*PackReader{}
	return s.loadPacks()
}

// HasObjectInPack checks if any packfile contains the given hash.
func (s *Store) HasObjectInPack(hash string) bool {
	s.mu.RLock()
	defer s.mu.RUnlock()

	for _, pack := range s.packs {
		if pack.HasObject(hash) {
			return true
		}
	}
	return false
}
