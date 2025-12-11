package repo

import (
	"fmt"

	"github.com/lishi/chemvcs/internal/model"
)

// MergeResult represents the outcome of a merge operation.
type MergeResult struct {
	FastForward bool     // True if merge was a fast-forward
	Conflicts   []string // List of conflicting file paths
}

// FindCommonAncestor finds the common ancestor of two snapshots.
// Returns empty string if no common ancestor exists.
func (r *Repository) FindCommonAncestor(hash1, hash2 string) (string, error) {
	// Get all ancestors of hash1
	ancestors1, err := r.collectAncestors(hash1)
	if err != nil {
		return "", fmt.Errorf("failed to collect ancestors of %s: %w", hash1, err)
	}

	// Walk ancestors of hash2 until we find one in ancestors1
	return r.findFirstCommonAncestor(hash2, ancestors1)
}

// collectAncestors collects all ancestors of a snapshot in a map.
func (r *Repository) collectAncestors(snapshotHash string) (map[string]bool, error) {
	ancestors := make(map[string]bool)
	queue := []string{snapshotHash}

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]

		if ancestors[current] {
			continue // Already visited
		}
		ancestors[current] = true

		// Get snapshot
		snapshot, err := r.store.GetSnapshot(current)
		if err != nil {
			return nil, fmt.Errorf("failed to get snapshot %s: %w", current, err)
		}

		// Add parents to queue
		for _, parentHash := range snapshot.Parents {
			if !ancestors[parentHash] {
				queue = append(queue, parentHash)
			}
		}
	}

	return ancestors, nil
}

// findFirstCommonAncestor walks ancestors of startHash until finding one in ancestorSet.
func (r *Repository) findFirstCommonAncestor(startHash string, ancestorSet map[string]bool) (string, error) {
	queue := []string{startHash}
	visited := make(map[string]bool)

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]

		if visited[current] {
			continue
		}
		visited[current] = true

		// Check if current is in ancestor set
		if ancestorSet[current] {
			return current, nil
		}

		// Get snapshot and add parents to queue
		snapshot, err := r.store.GetSnapshot(current)
		if err != nil {
			return "", fmt.Errorf("failed to get snapshot %s: %w", current, err)
		}

		for _, parentHash := range snapshot.Parents {
			if !visited[parentHash] {
				queue = append(queue, parentHash)
			}
		}
	}

	return "", nil // No common ancestor
}

// ThreeWayMerge performs a three-way merge between current and target using base as common ancestor.
// Returns the merged root object hash and any conflicts.
func (r *Repository) ThreeWayMerge(baseHash, currentHash, targetHash string) (*model.Object, []string, error) {
	// Get root objects
	baseObj, err := r.store.GetObject(baseHash)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to get base object: %w", err)
	}

	currentObj, err := r.store.GetObject(currentHash)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to get current object: %w", err)
	}

	targetObj, err := r.store.GetObject(targetHash)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to get target object: %w", err)
	}

	// Perform recursive merge
	conflicts := []string{}
	mergedObj, mergeConflicts, err := r.mergeObjects(baseObj, currentObj, targetObj, "")
	if err != nil {
		return nil, nil, err
	}

	conflicts = append(conflicts, mergeConflicts...)

	return mergedObj, conflicts, nil
}

// mergeObjects recursively merges three objects (base, current, target).
// relPath is the relative path for conflict reporting.
func (r *Repository) mergeObjects(base, current, target *model.Object, relPath string) (*model.Object, []string, error) {
	conflicts := []string{}

	// Handle different object types
	if base.Type != current.Type || base.Type != target.Type {
		// Type changed - this is a conflict
		conflicts = append(conflicts, relPath)
		// For now, prefer current
		return current, conflicts, nil
	}

	switch base.Type {
	case "folder":
		return r.mergeFolders(base, current, target, relPath)
	case "file":
		return r.mergeFiles(base, current, target, relPath)
	default:
		// Unknown type, treat as conflict
		conflicts = append(conflicts, relPath)
		return current, conflicts, nil
	}
}

// mergeFolders merges three folder objects.
func (r *Repository) mergeFolders(base, current, target *model.Object, relPath string) (*model.Object, []string, error) {
	conflicts := []string{}

	// Build maps of entries by name
	baseEntries := make(map[string]model.Reference)
	currentEntries := make(map[string]model.Reference)
	targetEntries := make(map[string]model.Reference)

	for _, ref := range base.Refs {
		name := getRefName(r, ref)
		baseEntries[name] = ref
	}

	for _, ref := range current.Refs {
		name := getRefName(r, ref)
		currentEntries[name] = ref
	}

	for _, ref := range target.Refs {
		name := getRefName(r, ref)
		targetEntries[name] = ref
	}

	// Merge entries
	mergedRefs := []model.Reference{}
	allNames := make(map[string]bool)

	// Collect all names
	for name := range baseEntries {
		allNames[name] = true
	}
	for name := range currentEntries {
		allNames[name] = true
	}
	for name := range targetEntries {
		allNames[name] = true
	}

	// Process each name
	for name := range allNames {
		baseRef, inBase := baseEntries[name]
		currentRef, inCurrent := currentEntries[name]
		targetRef, inTarget := targetEntries[name]

		var entryPath string
		if relPath == "" {
			entryPath = name
		} else {
			entryPath = relPath + "/" + name
		}

		if !inBase {
			// Added in current or target or both
			if inCurrent && inTarget {
				if currentRef.ID == targetRef.ID {
					// Both added the same thing
					mergedRefs = append(mergedRefs, currentRef)
				} else {
					// Both added different things - conflict
					conflicts = append(conflicts, entryPath)
					mergedRefs = append(mergedRefs, currentRef) // Prefer current
				}
			} else if inCurrent {
				// Added only in current
				mergedRefs = append(mergedRefs, currentRef)
			} else {
				// Added only in target
				mergedRefs = append(mergedRefs, targetRef)
			}
		} else {
			// Existed in base
			if !inCurrent && !inTarget {
				// Deleted in both - OK, don't include
			} else if !inCurrent {
				// Deleted in current
				if targetRef.ID == baseRef.ID {
					// Target unchanged, respect deletion
				} else {
					// Modified in target, deleted in current - conflict
					conflicts = append(conflicts, entryPath)
					// Deletion wins for now
				}
			} else if !inTarget {
				// Deleted in target
				if currentRef.ID == baseRef.ID {
					// Current unchanged, respect deletion
				} else {
					// Modified in current, deleted in target - conflict
					conflicts = append(conflicts, entryPath)
					mergedRefs = append(mergedRefs, currentRef) // Keep current
				}
			} else {
				// Modified in current or target or both
				if currentRef.ID == targetRef.ID {
					// Same modification or both unchanged
					mergedRefs = append(mergedRefs, currentRef)
				} else if currentRef.ID == baseRef.ID {
					// Only target modified
					mergedRefs = append(mergedRefs, targetRef)
				} else if targetRef.ID == baseRef.ID {
					// Only current modified
					mergedRefs = append(mergedRefs, currentRef)
				} else {
					// Both modified differently - recurse if same type
					baseObj, err := r.store.GetObject(baseRef.ID)
					if err != nil {
						return nil, nil, fmt.Errorf("failed to get base object: %w", err)
					}

					currentObj, err := r.store.GetObject(currentRef.ID)
					if err != nil {
						return nil, nil, fmt.Errorf("failed to get current object: %w", err)
					}

					targetObj, err := r.store.GetObject(targetRef.ID)
					if err != nil {
						return nil, nil, fmt.Errorf("failed to get target object: %w", err)
					}

					mergedObj, subConflicts, err := r.mergeObjects(baseObj, currentObj, targetObj, entryPath)
					if err != nil {
						return nil, nil, err
					}

					conflicts = append(conflicts, subConflicts...)

					// Store merged object
					mergedHash, err := r.store.PutObject(mergedObj)
					if err != nil {
						return nil, nil, fmt.Errorf("failed to store merged object: %w", err)
					}

					mergedRefs = append(mergedRefs, model.Reference{
						Kind: "object",
						ID:   mergedHash,
					})
				}
			}
		}
	}

	// Create merged folder object
	mergedFolder := &model.Object{
		Version: 1,
		Type:    "folder",
		Meta:    current.Meta, // Use current metadata
		Refs:    mergedRefs,
	}

	return mergedFolder, conflicts, nil
}

// mergeFiles merges three file objects.
// For now, files are treated atomically - if both changed, it's a conflict.
func (r *Repository) mergeFiles(base, current, target *model.Object, relPath string) (*model.Object, []string, error) {
	conflicts := []string{}

	// Get blob references
	if len(base.Refs) == 0 || len(current.Refs) == 0 || len(target.Refs) == 0 {
		return nil, nil, fmt.Errorf("file object missing blob reference")
	}

	baseBlob := base.Refs[0].ID
	currentBlob := current.Refs[0].ID
	targetBlob := target.Refs[0].ID

	// Check if changed
	if currentBlob == targetBlob {
		// Same content or both unchanged
		return current, conflicts, nil
	} else if currentBlob == baseBlob {
		// Only target changed
		return target, conflicts, nil
	} else if targetBlob == baseBlob {
		// Only current changed
		return current, conflicts, nil
	} else {
		// Both changed - conflict
		conflicts = append(conflicts, relPath)
		return current, conflicts, nil // Prefer current
	}
}

// getRefName gets the name of a reference by looking up its object.
func getRefName(r *Repository, ref model.Reference) string {
	if ref.Kind != "object" {
		return ""
	}

	obj, err := r.store.GetObject(ref.ID)
	if err != nil {
		return ""
	}

	if name, ok := obj.Meta["name"].(string); ok {
		return name
	}

	return ""
}
