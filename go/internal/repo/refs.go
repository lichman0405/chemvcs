package repo

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// RefManager handles reading and writing references (branches, HEAD).
type RefManager struct {
	repoPath string
}

// NewRefManager creates a new RefManager for the given repository path.
func NewRefManager(repoPath string) *RefManager {
	return &RefManager{
		repoPath: repoPath,
	}
}

// ResolveHEAD resolves HEAD to a snapshot hash.
// HEAD can be either:
// - A symbolic reference: "ref: refs/heads/main"
// - A direct hash: "abc123..."
// Returns empty string if HEAD points to a non-existent ref (e.g., initial repo).
func (rm *RefManager) ResolveHEAD() (string, error) {
	headPath := filepath.Join(rm.repoPath, ".chemvcs", "HEAD")

	data, err := os.ReadFile(headPath)
	if err != nil {
		if os.IsNotExist(err) {
			return "", nil // HEAD doesn't exist yet
		}
		return "", fmt.Errorf("failed to read HEAD: %w", err)
	}

	content := strings.TrimSpace(string(data))

	// Check if symbolic reference
	if strings.HasPrefix(content, "ref:") {
		refName := strings.TrimSpace(content[4:])
		return rm.ResolveRef(refName)
	}

	// Direct hash
	return content, nil
}

// ResolveRef resolves a reference name to a snapshot hash.
// Returns empty string if the ref doesn't exist yet.
func (rm *RefManager) ResolveRef(refName string) (string, error) {
	refPath := filepath.Join(rm.repoPath, ".chemvcs", refName)

	data, err := os.ReadFile(refPath)
	if err != nil {
		if os.IsNotExist(err) {
			return "", nil // Ref doesn't exist yet
		}
		return "", fmt.Errorf("failed to read ref %s: %w", refName, err)
	}

	return strings.TrimSpace(string(data)), nil
}

// CurrentBranch returns the current branch name (e.g., "refs/heads/main").
// Returns empty string if HEAD is detached.
func (rm *RefManager) CurrentBranch() (string, error) {
	headPath := filepath.Join(rm.repoPath, ".chemvcs", "HEAD")

	data, err := os.ReadFile(headPath)
	if err != nil {
		if os.IsNotExist(err) {
			return "", nil
		}
		return "", fmt.Errorf("failed to read HEAD: %w", err)
	}

	content := strings.TrimSpace(string(data))

	// Check if symbolic reference
	if strings.HasPrefix(content, "ref:") {
		return strings.TrimSpace(content[4:]), nil
	}

	// Detached HEAD
	return "", nil
}

// UpdateRef atomically updates a reference to point to a new snapshot hash.
func (rm *RefManager) UpdateRef(refName, snapshotHash string) error {
	refPath := filepath.Join(rm.repoPath, ".chemvcs", refName)

	// Ensure parent directory exists
	if err := os.MkdirAll(filepath.Dir(refPath), 0755); err != nil {
		return fmt.Errorf("failed to create ref directory: %w", err)
	}

	// Atomic write using temp file + rename
	return atomicWriteFile(refPath, []byte(snapshotHash+"\n"))
}

// SetHEAD sets HEAD to point to a reference or direct hash.
// For symbolic refs, use "ref: refs/heads/main".
// For direct hash, just pass the hash.
func (rm *RefManager) SetHEAD(target string) error {
	headPath := filepath.Join(rm.repoPath, ".chemvcs", "HEAD")

	content := target
	if !strings.HasPrefix(target, "ref:") && len(target) > 0 {
		// If not already a ref: prefix and not empty, ensure it ends with newline
		content = target + "\n"
	} else if strings.HasPrefix(target, "ref:") {
		// Ensure ref: format has newline
		if !strings.HasSuffix(target, "\n") {
			content = target + "\n"
		}
	}

	return atomicWriteFile(headPath, []byte(content))
}

// ListBranches returns a list of all branch names (without refs/heads/ prefix).
func (rm *RefManager) ListBranches() ([]string, error) {
	headsPath := filepath.Join(rm.repoPath, ".chemvcs", "refs", "heads")

	entries, err := os.ReadDir(headsPath)
	if err != nil {
		if os.IsNotExist(err) {
			return []string{}, nil
		}
		return nil, fmt.Errorf("failed to read branches: %w", err)
	}

	branches := []string{}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		branches = append(branches, entry.Name())
	}

	return branches, nil
}

// CreateBranch creates a new branch pointing to the given snapshot hash.
func (rm *RefManager) CreateBranch(branchName, snapshotHash string) error {
	refName := "refs/heads/" + branchName
	return rm.UpdateRef(refName, snapshotHash)
}

// atomicWriteFile writes data to a file atomically using temp file + rename.
func atomicWriteFile(path string, data []byte) error {
	dir := filepath.Dir(path)

	// Create temp file in same directory
	tmpFile, err := os.CreateTemp(dir, ".tmp-*")
	if err != nil {
		return fmt.Errorf("failed to create temp file: %w", err)
	}
	tmpPath := tmpFile.Name()

	// Write data
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
