package remote

import (
	"fmt"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/repo"
)

// PushOptions contains options for a push operation.
type PushOptions struct {
	RepoID    string // Remote repository ID
	LocalRef  string // Local ref to push (e.g., "refs/heads/main")
	RemoteRef string // Remote ref to update (e.g., "refs/heads/main")
	Force     bool   // Force non-fast-forward updates (not implemented in MVP)
}

// PullOptions contains options for a pull operation.
type PullOptions struct {
	RepoID    string // Remote repository ID
	RemoteRef string // Remote ref to pull (e.g., "refs/heads/main")
	LocalRef  string // Local ref to update (e.g., "refs/heads/main")
}

// Push pushes a local branch to a remote repository.
func (c *Client) Push(r *repo.Repository, opts PushOptions) error {
	// 1. Get local snapshot
	localSnap, err := r.ResolveRef(opts.LocalRef)
	if err != nil {
		return fmt.Errorf("resolving local ref: %w", err)
	}
	if localSnap == "" {
		return fmt.Errorf("local ref %s does not exist", opts.LocalRef)
	}

	// 2. Collect all objects reachable from local snapshot
	objectIDs, err := collectReachableObjects(r, localSnap)
	if err != nil {
		return fmt.Errorf("collecting reachable objects: %w", err)
	}

	// 3. Check which objects exist on remote
	existsResp, err := c.CheckObjectsExist(opts.RepoID, objectIDs)
	if err != nil {
		return fmt.Errorf("checking remote objects: %w", err)
	}

	// 4. Upload missing objects
	for _, objID := range existsResp.Missing {
		data, err := r.Store().GetRaw(objID)
		if err != nil {
			return fmt.Errorf("reading object %s: %w", objID, err)
		}

		if err := c.UploadObject(opts.RepoID, objID, data); err != nil {
			return fmt.Errorf("uploading object %s: %w", objID, err)
		}
	}

	// 5. Get current remote ref (for concurrency control)
	var oldTarget string
	remoteRef, err := c.GetRef(opts.RepoID, opts.RemoteRef)
	if err != nil {
		// Ref might not exist yet (first push), which is OK
		oldTarget = ""
	} else {
		oldTarget = remoteRef.Target

		// Check if this is a fast-forward (MVP: only allow fast-forward)
		if oldTarget != "" && !opts.Force {
			isAncestor, err := r.IsAncestor(oldTarget, localSnap)
			if err != nil {
				return fmt.Errorf("checking ancestry: %w", err)
			}
			if !isAncestor {
				return fmt.Errorf("push rejected: not a fast-forward (remote has diverged)")
			}
		}
	}

	// 6. Update remote ref
	if err := c.UpdateRef(opts.RepoID, opts.RemoteRef, oldTarget, localSnap); err != nil {
		return fmt.Errorf("updating remote ref: %w", err)
	}

	return nil
}

// Pull fetches objects from a remote repository and updates the local ref.
func (c *Client) Pull(r *repo.Repository, opts PullOptions) error {
	// 1. Get remote ref
	remoteRef, err := c.GetRef(opts.RepoID, opts.RemoteRef)
	if err != nil {
		return fmt.Errorf("fetching remote ref: %w", err)
	}

	remoteSnap := remoteRef.Target

	// 2. Check if we already have this snapshot
	if _, err := r.GetSnapshot(remoteSnap); err == nil {
		// We already have it, just update local ref
		return r.UpdateRef(opts.LocalRef, remoteSnap)
	}

	// 3. Collect all objects we need to download
	objectIDs, err := c.collectRemoteObjects(opts.RepoID, remoteSnap)
	if err != nil {
		return fmt.Errorf("collecting remote objects: %w", err)
	}

	// 4. Check which objects we already have locally
	var missing []string
	for _, objID := range objectIDs {
		if _, err := r.Store().GetRaw(objID); err != nil {
			missing = append(missing, objID)
		}
	}

	// 5. Download missing objects
	for _, objID := range missing {
		data, err := c.DownloadObject(opts.RepoID, objID)
		if err != nil {
			return fmt.Errorf("downloading object %s: %w", objID, err)
		}

		// Store directly using PutRaw (assumes data is already properly formatted)
		if err := r.Store().PutRaw(objID, data); err != nil {
			return fmt.Errorf("storing object %s: %w", objID, err)
		}
	}

	// 6. Check if pull is fast-forward
	localSnap, err := r.ResolveRef(opts.LocalRef)
	if err == nil && localSnap != "" {
		// Local ref exists, check if this is fast-forward
		isAncestor, err := r.IsAncestor(localSnap, remoteSnap)
		if err != nil {
			return fmt.Errorf("checking ancestry: %w", err)
		}
		if !isAncestor {
			return fmt.Errorf("pull rejected: not a fast-forward (local branch has diverged)")
		}
	}

	// 7. Update local ref
	if err := r.UpdateRef(opts.LocalRef, remoteSnap); err != nil {
		return fmt.Errorf("updating local ref: %w", err)
	}

	return nil
}

// collectReachableObjects collects all object IDs reachable from a snapshot.
func collectReachableObjects(r *repo.Repository, snapID string) ([]string, error) {
	visited := make(map[string]bool)
	var result []string

	var visit func(id string) error
	visit = func(id string) error {
		if visited[id] {
			return nil
		}
		visited[id] = true
		result = append(result, id)

		// Try to load as snapshot
		snap, err := r.GetSnapshot(id)
		if err == nil {
			// Visit root object
			if err := visit(snap.Root); err != nil {
				return err
			}
			// Visit parents
			for _, parent := range snap.Parents {
				if err := visit(parent); err != nil {
					return err
				}
			}
			return nil
		}

		// Try to load as object (folder)
		obj, err := r.Store().GetObject(id)
		if err == nil {
			for _, ref := range obj.Refs {
				if err := visit(ref.ID); err != nil {
					return err
				}
			}
			return nil
		}

		// Otherwise, it's a blob (no children to visit)
		return nil
	}

	if err := visit(snapID); err != nil {
		return nil, err
	}

	return result, nil
}

// collectRemoteObjects collects object IDs from a remote snapshot by recursively downloading.
func (c *Client) collectRemoteObjects(repoID, snapID string) ([]string, error) {
	visited := make(map[string]bool)
	var result []string

	var visit func(id string) error
	visit = func(id string) error {
		if visited[id] {
			return nil
		}
		visited[id] = true
		result = append(result, id)

		// Download the object
		data, err := c.DownloadObject(repoID, id)
		if err != nil {
			return fmt.Errorf("downloading %s: %w", id, err)
		}

		// Try to parse as snapshot
		var snap model.Snapshot
		if err := model.DecodeJSON(data, &snap); err == nil {
			// Visit root and parents
			if err := visit(snap.Root); err != nil {
				return err
			}
			for _, parent := range snap.Parents {
				if err := visit(parent); err != nil {
					return err
				}
			}
			return nil
		}

		// Try to parse as object (folder)
		var obj model.Object
		if err := model.DecodeJSON(data, &obj); err == nil {
			for _, ref := range obj.Refs {
				if err := visit(ref.ID); err != nil {
					return err
				}
			}
			return nil
		}

		// Otherwise, it's a blob (no children)
		return nil
	}

	if err := visit(snapID); err != nil {
		return nil, err
	}

	return result, nil
}

// FetchOptions contains options for a fetch operation.
type FetchOptions struct {
	RepoID    string // Remote repository ID
	RemoteRef string // Remote ref to fetch (e.g., "refs/heads/main")
}

// Fetch downloads objects from a remote repository without updating any local refs.
// Returns the snapshot ID of the fetched ref.
func (c *Client) Fetch(r *repo.Repository, opts FetchOptions) (string, error) {
	// Get remote ref
	remoteRef, err := c.GetRef(opts.RepoID, opts.RemoteRef)
	if err != nil {
		return "", fmt.Errorf("fetching remote ref: %w", err)
	}

	remoteSnap := remoteRef.Target

	// Check if we already have this snapshot
	if _, err := r.GetSnapshot(remoteSnap); err == nil {
		return remoteSnap, nil
	}

	// Collect and download all objects
	objectIDs, err := c.collectRemoteObjects(opts.RepoID, remoteSnap)
	if err != nil {
		return "", fmt.Errorf("collecting remote objects: %w", err)
	}

	// Check which objects we need
	var missing []string
	for _, objID := range objectIDs {
		if _, err := r.Store().GetRaw(objID); err != nil {
			missing = append(missing, objID)
		}
	}

	// Download missing objects
	for _, objID := range missing {
		data, err := c.DownloadObject(opts.RepoID, objID)
		if err != nil {
			return "", fmt.Errorf("downloading object %s: %w", objID, err)
		}

		if err := r.Store().PutRaw(objID, data); err != nil {
			return "", fmt.Errorf("storing object %s: %w", objID, err)
		}
	}

	return remoteSnap, nil
}

// AddRemoteConfig adds a remote configuration to the repository.
func AddRemoteConfig(r *repo.Repository, name, url string) error {
	// For MVP, we'll store remotes in a simple config file format
	// This is a simplified implementation; a full version would use repo.Config
	remotesPath := r.GitDir() + "/remotes/" + name

	config := fmt.Sprintf("url=%s\n", url)

	return r.Store().WriteFile(remotesPath, []byte(config))
}

// GetRemoteURL retrieves the URL for a named remote.
func GetRemoteURL(r *repo.Repository, name string) (string, error) {
	remotesPath := r.GitDir() + "/remotes/" + name

	data, err := r.Store().ReadFile(remotesPath)
	if err != nil {
		return "", fmt.Errorf("remote %q not found", name)
	}

	// Parse simple config format
	lines := string(data)
	if len(lines) > 4 && lines[:4] == "url=" {
		url := lines[4:]
		// Remove trailing newline if present
		if len(url) > 0 && url[len(url)-1] == '\n' {
			url = url[:len(url)-1]
		}
		return url, nil
	}

	return "", fmt.Errorf("invalid remote config for %q", name)
}
