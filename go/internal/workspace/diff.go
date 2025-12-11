package workspace

import (
	"fmt"
	"sort"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/objectstore"
)

// Change represents a single file change.
type Change struct {
	Type    string // "added", "modified", "deleted"
	Path    string
	OldHash string
	NewHash string
}

// Changes represents all changes in a directory tree.
type Changes struct {
	Added    []string
	Modified []string
	Deleted  []string
}

// ComputeChanges compares two object trees and returns the differences.
// Both oldRootHash and newRootHash should be folder object hashes.
func ComputeChanges(oldRootHash, newRootHash string, store *objectstore.Store) (*Changes, error) {
	changes := &Changes{
		Added:    []string{},
		Modified: []string{},
		Deleted:  []string{},
	}

	// Get both root objects
	oldObj, err := store.GetObject(oldRootHash)
	if err != nil {
		return nil, fmt.Errorf("failed to get old root object: %w", err)
	}

	newObj, err := store.GetObject(newRootHash)
	if err != nil {
		return nil, fmt.Errorf("failed to get new root object: %w", err)
	}

	// Compare recursively
	if err := compareObjects(oldObj, newObj, "", changes, store); err != nil {
		return nil, err
	}

	// Sort results for consistent output
	sort.Strings(changes.Added)
	sort.Strings(changes.Modified)
	sort.Strings(changes.Deleted)

	return changes, nil
}

// compareObjects recursively compares two objects.
func compareObjects(oldObj, newObj *model.Object, currentPath string, changes *Changes, store *objectstore.Store) error {
	// Both must be folders for recursive comparison
	if oldObj.Type != "folder" || newObj.Type != "folder" {
		return fmt.Errorf("can only compare folder objects")
	}

	// Build maps of children by name
	oldChildren := make(map[string]string) // name -> object hash
	newChildren := make(map[string]string)

	// Map old children
	for _, ref := range oldObj.Refs {
		if ref.Kind != "object" {
			continue
		}

		childObj, err := store.GetObject(ref.ID)
		if err != nil {
			return fmt.Errorf("failed to get old child object: %w", err)
		}

		name, ok := childObj.Meta["name"].(string)
		if !ok {
			return fmt.Errorf("child object missing name metadata")
		}

		oldChildren[name] = ref.ID
	}

	// Map new children
	for _, ref := range newObj.Refs {
		if ref.Kind != "object" {
			continue
		}

		childObj, err := store.GetObject(ref.ID)
		if err != nil {
			return fmt.Errorf("failed to get new child object: %w", err)
		}

		name, ok := childObj.Meta["name"].(string)
		if !ok {
			return fmt.Errorf("child object missing name metadata")
		}

		newChildren[name] = ref.ID
	}

	// Find added, deleted, and potentially modified
	allNames := make(map[string]bool)
	for name := range oldChildren {
		allNames[name] = true
	}
	for name := range newChildren {
		allNames[name] = true
	}

	for name := range allNames {
		oldHash, inOld := oldChildren[name]
		newHash, inNew := newChildren[name]

		// Compute full path
		var fullPath string
		if currentPath == "" {
			fullPath = name
		} else {
			fullPath = currentPath + "/" + name
		}

		if !inOld && inNew {
			// Added
			if err := markAsAdded(newHash, fullPath, changes, store); err != nil {
				return err
			}
		} else if inOld && !inNew {
			// Deleted
			if err := markAsDeleted(oldHash, fullPath, changes, store); err != nil {
				return err
			}
		} else {
			// Both exist - check if modified
			if oldHash != newHash {
				// Hashes differ - need to check type
				oldChild, err := store.GetObject(oldHash)
				if err != nil {
					return fmt.Errorf("failed to get old child: %w", err)
				}

				newChild, err := store.GetObject(newHash)
				if err != nil {
					return fmt.Errorf("failed to get new child: %w", err)
				}

				if oldChild.Type != newChild.Type {
					// Type changed - treat as delete + add
					if err := markAsDeleted(oldHash, fullPath, changes, store); err != nil {
						return err
					}
					if err := markAsAdded(newHash, fullPath, changes, store); err != nil {
						return err
					}
				} else if oldChild.Type == "file" {
					// File modified
					changes.Modified = append(changes.Modified, fullPath)
				} else if oldChild.Type == "folder" {
					// Recurse into folder
					if err := compareObjects(oldChild, newChild, fullPath, changes, store); err != nil {
						return err
					}
				}
			}
			// If hashes are same, no change
		}
	}

	return nil
}

// markAsAdded recursively marks all files under an object as added.
func markAsAdded(objHash, basePath string, changes *Changes, store *objectstore.Store) error {
	obj, err := store.GetObject(objHash)
	if err != nil {
		return fmt.Errorf("failed to get object: %w", err)
	}

	if obj.Type == "file" {
		changes.Added = append(changes.Added, basePath)
	} else if obj.Type == "folder" {
		// Recurse into children
		for _, ref := range obj.Refs {
			if ref.Kind != "object" {
				continue
			}

			childObj, err := store.GetObject(ref.ID)
			if err != nil {
				return fmt.Errorf("failed to get child object: %w", err)
			}

			name, ok := childObj.Meta["name"].(string)
			if !ok {
				return fmt.Errorf("child object missing name metadata")
			}

			childPath := basePath + "/" + name
			if err := markAsAdded(ref.ID, childPath, changes, store); err != nil {
				return err
			}
		}
	}

	return nil
}

// markAsDeleted recursively marks all files under an object as deleted.
func markAsDeleted(objHash, basePath string, changes *Changes, store *objectstore.Store) error {
	obj, err := store.GetObject(objHash)
	if err != nil {
		return fmt.Errorf("failed to get object: %w", err)
	}

	if obj.Type == "file" {
		changes.Deleted = append(changes.Deleted, basePath)
	} else if obj.Type == "folder" {
		// Recurse into children
		for _, ref := range obj.Refs {
			if ref.Kind != "object" {
				continue
			}

			childObj, err := store.GetObject(ref.ID)
			if err != nil {
				return fmt.Errorf("failed to get child object: %w", err)
			}

			name, ok := childObj.Meta["name"].(string)
			if !ok {
				return fmt.Errorf("child object missing name metadata")
			}

			childPath := basePath + "/" + name
			if err := markAsDeleted(ref.ID, childPath, changes, store); err != nil {
				return err
			}
		}
	}

	return nil
}
