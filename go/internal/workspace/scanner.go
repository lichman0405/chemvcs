package workspace

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/objectstore"
)

// Scanner handles working directory scanning and restoration.
type Scanner struct {
	store *objectstore.Store
}

// NewScanner creates a new Scanner instance.
func NewScanner(store *objectstore.Store) *Scanner {
	return &Scanner{store: store}
}

// ScanDirectory recursively scans a directory and creates object tree.
// Returns the root folder object (not yet stored).
func (s *Scanner) ScanDirectory(basePath string) (*model.Object, error) {
	// Convert to absolute path
	absPath, err := filepath.Abs(basePath)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve path: %w", err)
	}

	// Check path exists and is a directory
	info, err := os.Stat(absPath)
	if err != nil {
		return nil, fmt.Errorf("failed to stat path: %w", err)
	}
	if !info.IsDir() {
		return nil, fmt.Errorf("path is not a directory: %s", absPath)
	}

	// Scan from root
	return s.scanDir(absPath, "")
}

// scanDir recursively scans a directory and returns a folder object.
func (s *Scanner) scanDir(basePath, relPath string) (*model.Object, error) {
	fullPath := filepath.Join(basePath, relPath)

	// Read directory entries
	entries, err := os.ReadDir(fullPath)
	if err != nil {
		return nil, fmt.Errorf("failed to read directory %s: %w", fullPath, err)
	}

	// Create folder object
	folderName := filepath.Base(relPath)
	if relPath == "" {
		folderName = filepath.Base(basePath)
	}

	folderObj := &model.Object{
		Version: 1,
		Type:    "folder",
		Meta: map[string]interface{}{
			"name": folderName,
		},
		Refs: []model.Reference{},
	}

	// Process entries
	for _, entry := range entries {
		// Skip .chemvcs directory
		if entry.Name() == ".chemvcs" {
			continue
		}

		entryRelPath := filepath.Join(relPath, entry.Name())

		if entry.IsDir() {
			// Recurse into subdirectory
			childObj, err := s.scanDir(basePath, entryRelPath)
			if err != nil {
				return nil, err
			}

			// Store child object and add reference
			childHash, err := s.store.PutObject(childObj)
			if err != nil {
				return nil, fmt.Errorf("failed to store folder object: %w", err)
			}

			folderObj.Refs = append(folderObj.Refs, model.Reference{
				Kind: "object",
				ID:   childHash,
			})
		} else {
			// Scan file
			fileObj, err := s.scanFile(filepath.Join(basePath, entryRelPath), entry)
			if err != nil {
				return nil, err
			}

			// Store file object and add reference
			fileHash, err := s.store.PutObject(fileObj)
			if err != nil {
				return nil, fmt.Errorf("failed to store file object: %w", err)
			}

			folderObj.Refs = append(folderObj.Refs, model.Reference{
				Kind: "object",
				ID:   fileHash,
			})
		}
	}

	return folderObj, nil
}

// scanFile reads a file and creates a file object with blob reference.
func (s *Scanner) scanFile(path string, entry os.DirEntry) (*model.Object, error) {
	// Read file content
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("failed to read file %s: %w", path, err)
	}

	// Store blob
	blobHash, err := s.store.PutBlob(data)
	if err != nil {
		return nil, fmt.Errorf("failed to store blob: %w", err)
	}

	// Get file info
	info, err := entry.Info()
	if err != nil {
		return nil, fmt.Errorf("failed to get file info: %w", err)
	}

	// Create file object
	fileObj := &model.Object{
		Version: 1,
		Type:    "file",
		Meta: map[string]interface{}{
			"name": entry.Name(),
			"size": info.Size(),
			"mode": fmt.Sprintf("%o", info.Mode()),
		},
		Refs: []model.Reference{
			{Kind: "blob", ID: blobHash},
		},
	}

	return fileObj, nil
}

// RestoreDirectory restores a complete directory tree from an object hash.
func (s *Scanner) RestoreDirectory(rootObjHash, targetPath string) error {
	// Get root object
	rootObj, err := s.store.GetObject(rootObjHash)
	if err != nil {
		return fmt.Errorf("failed to get root object: %w", err)
	}

	// Ensure target path is absolute
	absTarget, err := filepath.Abs(targetPath)
	if err != nil {
		return fmt.Errorf("failed to resolve target path: %w", err)
	}

	// Restore recursively
	return s.restoreObject(rootObj, absTarget)
}

// restoreObject restores a single object (folder or file).
func (s *Scanner) restoreObject(obj *model.Object, path string) error {
	switch obj.Type {
	case "folder":
		return s.restoreFolder(obj, path)
	case "file":
		return s.restoreFile(obj, path)
	default:
		return fmt.Errorf("unknown object type: %s", obj.Type)
	}
}

// restoreFolder restores a folder and its contents.
func (s *Scanner) restoreFolder(obj *model.Object, path string) error {
	// Create directory
	if err := os.MkdirAll(path, 0755); err != nil {
		return fmt.Errorf("failed to create directory %s: %w", path, err)
	}

	// Restore children
	for _, ref := range obj.Refs {
		if ref.Kind != "object" {
			continue
		}

		// Get child object
		childObj, err := s.store.GetObject(ref.ID)
		if err != nil {
			return fmt.Errorf("failed to get child object %s: %w", ref.ID, err)
		}

		// Get child name from metadata
		childName, ok := childObj.Meta["name"].(string)
		if !ok {
			return fmt.Errorf("child object missing 'name' in metadata")
		}

		childPath := filepath.Join(path, childName)

		// Restore child
		if err := s.restoreObject(childObj, childPath); err != nil {
			return err
		}
	}

	return nil
}

// restoreFile restores a file from its blob reference.
func (s *Scanner) restoreFile(obj *model.Object, path string) error {
	// Get blob reference
	if len(obj.Refs) == 0 {
		return fmt.Errorf("file object has no blob reference")
	}

	blobRef := obj.Refs[0]
	if blobRef.Kind != "blob" {
		return fmt.Errorf("first reference is not a blob: %s", blobRef.Kind)
	}

	// Get blob data
	data, err := s.store.GetBlob(blobRef.ID)
	if err != nil {
		return fmt.Errorf("failed to get blob %s: %w", blobRef.ID, err)
	}

	// Parse file mode
	mode := os.FileMode(0644)
	if modeStr, ok := obj.Meta["mode"].(string); ok {
		if modeInt, err := strconv.ParseUint(modeStr, 8, 32); err == nil {
			mode = os.FileMode(modeInt)
		}
	}

	// Write file
	if err := os.WriteFile(path, data, mode); err != nil {
		return fmt.Errorf("failed to write file %s: %w", path, err)
	}

	return nil
}
