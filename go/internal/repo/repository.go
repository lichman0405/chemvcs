package repo

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/objectstore"
)

// Repository represents a ChemVCS repository.
type Repository struct {
	path  string             // Root path of the repository (contains .chemvcs)
	store *objectstore.Store // Object storage
	refs  *RefManager        // Reference manager
}

// Init initialises a new ChemVCS repository in the specified directory.
// Creates the .chemvcs directory structure and initial HEAD pointing to main branch.
func Init(path string) (*Repository, error) {
	absPath, err := filepath.Abs(path)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve absolute path: %w", err)
	}

	chemvcsPath := filepath.Join(absPath, ".chemvcs")

	// Check if repository already exists
	if _, err := os.Stat(chemvcsPath); err == nil {
		return nil, fmt.Errorf("repository already exists at %s", absPath)
	}

	// Create directory structure
	dirs := []string{
		chemvcsPath,
		filepath.Join(chemvcsPath, "objects"),
		filepath.Join(chemvcsPath, "refs", "heads"),
	}

	for _, dir := range dirs {
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("failed to create directory %s: %w", dir, err)
		}
	}

	// Create initial HEAD pointing to main branch
	headPath := filepath.Join(chemvcsPath, "HEAD")
	initialHEAD := "ref: refs/heads/main"
	if err := os.WriteFile(headPath, []byte(initialHEAD), 0644); err != nil {
		return nil, fmt.Errorf("failed to create HEAD: %w", err)
	}

	// Create optional config file placeholder
	configPath := filepath.Join(chemvcsPath, "config")
	configContent := "# ChemVCS repository configuration\n"
	if err := os.WriteFile(configPath, []byte(configContent), 0644); err != nil {
		return nil, fmt.Errorf("failed to create config: %w", err)
	}

	return Open(absPath)
}

// Open opens an existing ChemVCS repository.
// It searches upwards from the given path to find a .chemvcs directory.
func Open(path string) (*Repository, error) {
	absPath, err := filepath.Abs(path)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve absolute path: %w", err)
	}

	repoPath := findRepository(absPath)
	if repoPath == "" {
		return nil, fmt.Errorf("not a chemvcs repository (or any of the parent directories)")
	}

	store, err := objectstore.NewStore(repoPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open object store: %w", err)
	}

	return &Repository{
		path:  repoPath,
		store: store,
		refs:  NewRefManager(repoPath),
	}, nil
}

// findRepository searches upwards from startPath for a .chemvcs directory.
// Returns the repository root path, or empty string if not found.
func findRepository(startPath string) string {
	current := startPath

	for {
		chemvcsPath := filepath.Join(current, ".chemvcs")
		if info, err := os.Stat(chemvcsPath); err == nil && info.IsDir() {
			return current
		}

		parent := filepath.Dir(current)
		if parent == current {
			// Reached filesystem root
			return ""
		}
		current = parent
	}
}

// Path returns the root path of the repository.
func (r *Repository) Path() string {
	return r.path
}

// Store returns the object store.
func (r *Repository) Store() *objectstore.Store {
	return r.store
}

// Refs returns the reference manager.
func (r *Repository) Refs() *RefManager {
	return r.refs
}

// CreateSnapshot creates a new snapshot with the given root object hash, message, and author.
// It automatically determines the parent snapshot from the current HEAD.
// Returns the hash of the created snapshot.
func (r *Repository) CreateSnapshot(rootObjectHash, message, author string) (string, error) {
	// Get parent snapshot
	parentHash, err := r.refs.ResolveHEAD()
	if err != nil {
		return "", fmt.Errorf("failed to resolve HEAD: %w", err)
	}

	var parents []string
	if parentHash != "" {
		parents = []string{parentHash}
	}

	// Create snapshot
	snapshot := model.NewSnapshot(rootObjectHash, author, message, parents)

	// Store snapshot in object store
	// CRITICAL: Objects are written first, refs updated last
	hash, err := r.store.PutSnapshot(snapshot)
	if err != nil {
		return "", fmt.Errorf("failed to store snapshot: %w", err)
	}

	// Update current branch ref
	currentBranch, err := r.refs.CurrentBranch()
	if err != nil {
		return "", fmt.Errorf("failed to get current branch: %w", err)
	}

	if currentBranch != "" {
		if err := r.refs.UpdateRef(currentBranch, hash); err != nil {
			return "", fmt.Errorf("failed to update branch ref: %w", err)
		}
	} else {
		// Detached HEAD - update HEAD directly
		if err := r.refs.SetHEAD(hash); err != nil {
			return "", fmt.Errorf("failed to update HEAD: %w", err)
		}
	}

	return hash, nil
}

// GetSnapshot retrieves a snapshot by its hash.
func (r *Repository) GetSnapshot(hash string) (*model.Snapshot, error) {
	return r.store.GetSnapshot(hash)
}

// Log returns the history of snapshots starting from HEAD, up to the specified limit.
// If limit is 0 or negative, returns all snapshots in the history.
func (r *Repository) Log(limit int) ([]*model.Snapshot, error) {
	currentHash, err := r.refs.ResolveHEAD()
	if err != nil {
		return nil, fmt.Errorf("failed to resolve HEAD: %w", err)
	}

	if currentHash == "" {
		// No commits yet
		return []*model.Snapshot{}, nil
	}

	snapshots := []*model.Snapshot{}

	for currentHash != "" {
		// Check limit
		if limit > 0 && len(snapshots) >= limit {
			break
		}

		snap, err := r.store.GetSnapshot(currentHash)
		if err != nil {
			return nil, fmt.Errorf("failed to get snapshot %s: %w", currentHash, err)
		}

		snapshots = append(snapshots, snap)

		// Follow parent chain (MVP: linear history only)
		if len(snap.Parents) > 0 {
			currentHash = snap.Parents[0]
		} else {
			currentHash = ""
		}
	}

	return snapshots, nil
}

// CreateBranch creates a new branch at the current HEAD.
func (r *Repository) CreateBranch(branchName string) error {
	currentHash, err := r.refs.ResolveHEAD()
	if err != nil {
		return fmt.Errorf("failed to resolve HEAD: %w", err)
	}

	if currentHash == "" {
		return fmt.Errorf("cannot create branch: no commits yet")
	}

	return r.refs.CreateBranch(branchName, currentHash)
}

// ListBranches returns a list of all branches.
func (r *Repository) ListBranches() ([]string, error) {
	return r.refs.ListBranches()
}

// CheckoutBranch switches to the specified branch.
// This updates HEAD to point to the branch.
func (r *Repository) CheckoutBranch(branchName string) error {
	// Verify branch exists
	refName := "refs/heads/" + branchName
	hash, err := r.refs.ResolveRef(refName)
	if err != nil {
		return err
	}

	if hash == "" {
		return fmt.Errorf("branch '%s' not found", branchName)
	}

	// Update HEAD to point to the branch
	return r.refs.SetHEAD("ref: " + refName)
}

// CheckoutSnapshot switches to a specific snapshot (detached HEAD).
func (r *Repository) CheckoutSnapshot(snapshotHash string) error {
	// Verify snapshot exists
	_, err := r.store.GetSnapshot(snapshotHash)
	if err != nil {
		return fmt.Errorf("snapshot not found: %s", snapshotHash)
	}

	// Update HEAD to point directly to the snapshot
	return r.refs.SetHEAD(snapshotHash)
}

// CurrentBranch returns the current branch name (without refs/heads/ prefix).
// Returns empty string if HEAD is detached.
func (r *Repository) CurrentBranch() (string, error) {
	branch, err := r.refs.CurrentBranch()
	if err != nil {
		return "", err
	}

	// Remove refs/heads/ prefix if present
	if branch != "" && filepath.Base(branch) != branch {
		return filepath.Base(branch), nil
	}

	return branch, nil
}

// IsDetached returns true if HEAD is detached (not on a branch).
func (r *Repository) IsDetached() (bool, error) {
	branch, err := r.refs.CurrentBranch()
	if err != nil {
		return false, err
	}
	return branch == "", nil
}

// IsAncestor checks if ancestorHash is an ancestor of descendantHash.
// This walks the parent chain from descendant to see if ancestor is found.
func (r *Repository) IsAncestor(ancestorHash, descendantHash string) (bool, error) {
	if ancestorHash == descendantHash {
		return true, nil
	}

	current := descendantHash
	visited := make(map[string]bool) // Prevent infinite loops

	for current != "" {
		if visited[current] {
			// Cycle detected (shouldn't happen in normal operation)
			return false, nil
		}
		visited[current] = true

		if current == ancestorHash {
			return true, nil
		}

		snap, err := r.store.GetSnapshot(current)
		if err != nil {
			return false, fmt.Errorf("failed to get snapshot %s: %w", current, err)
		}

		// Follow first parent (linear history for now)
		if len(snap.Parents) > 0 {
			current = snap.Parents[0]
		} else {
			// Reached root
			return false, nil
		}
	}

	return false, nil
}

// Merge performs a merge of the specified branch into the current branch.
// Currently only supports fast-forward merges.
func (r *Repository) Merge(branchName string) (bool, error) {
	// Check if on a branch (can't merge into detached HEAD)
	currentBranch, err := r.refs.CurrentBranch()
	if err != nil {
		return false, err
	}
	if currentBranch == "" {
		return false, fmt.Errorf("cannot merge: HEAD is detached")
	}

	// Resolve current HEAD
	currentHash, err := r.refs.ResolveHEAD()
	if err != nil {
		return false, fmt.Errorf("failed to resolve HEAD: %w", err)
	}
	if currentHash == "" {
		return false, fmt.Errorf("cannot merge: no commits yet")
	}

	// Resolve target branch
	targetRef := "refs/heads/" + branchName
	targetHash, err := r.refs.ResolveRef(targetRef)
	if err != nil {
		return false, fmt.Errorf("failed to resolve branch: %w", err)
	}
	if targetHash == "" {
		return false, fmt.Errorf("branch '%s' not found", branchName)
	}

	// Check if already up to date
	if currentHash == targetHash {
		return false, nil // Already up to date
	}

	// Check if fast-forward is possible
	// Target must be an ancestor of current, OR current must be an ancestor of target
	targetIsAncestor, err := r.IsAncestor(targetHash, currentHash)
	if err != nil {
		return false, fmt.Errorf("failed to check ancestry: %w", err)
	}

	if targetIsAncestor {
		// Target is behind current - already up to date
		return false, nil
	}

	// Check if current is ancestor of target (fast-forward possible)
	currentIsAncestor, err := r.IsAncestor(currentHash, targetHash)
	if err != nil {
		return false, fmt.Errorf("failed to check ancestry: %w", err)
	}

	if !currentIsAncestor {
		// Diverged branches - need three-way merge (not supported yet)
		return false, fmt.Errorf("cannot fast-forward: branches have diverged (three-way merge not yet implemented)")
	}

	// Perform fast-forward merge: update current branch to point to target
	if err := r.refs.UpdateRef(currentBranch, targetHash); err != nil {
		return false, fmt.Errorf("failed to update branch: %w", err)
	}

	return true, nil // Fast-forward merge successful
}

// ResolveRef resolves a ref name to a snapshot hash.
// Returns empty string if ref doesn't exist (not an error).
func (r *Repository) ResolveRef(refName string) (string, error) {
	return r.refs.ResolveRef(refName)
}

// UpdateRef updates a ref to point to a new snapshot hash.
func (r *Repository) UpdateRef(refName, snapshotHash string) error {
	return r.refs.UpdateRef(refName, snapshotHash)
}

// GitDir returns the path to the .chemvcs directory.
func (r *Repository) GitDir() string {
	return filepath.Join(r.path, ".chemvcs")
}
