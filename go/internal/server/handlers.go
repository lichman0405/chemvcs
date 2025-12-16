package server

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"strconv"
	"strings"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/repo"
)

// handleRepoInfo handles GET /repos/{repoId} - get repository information.
func (s *Server) handleRepoInfo(w http.ResponseWriter, r *http.Request, repoID string) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// Get default branch
	currentBranch, _ := repo.CurrentBranch()
	if currentBranch == "" {
		currentBranch = "main"
	}

	s.sendJSON(w, http.StatusOK, map[string]string{
		"id":             repoID,
		"description":    "ChemVCS repository",
		"default_branch": currentBranch,
	})
}

// handleObjectsExist handles POST /repos/{repoId}/objects/exists.
func (s *Server) handleObjectsExist(w http.ResponseWriter, r *http.Request, repoID string) {
	if r.Method != http.MethodPost {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// Parse request
	var req struct {
		IDs []string `json:"ids"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid JSON")
		return
	}

	// Check which objects exist
	var present, missing []string
	for _, id := range req.IDs {
		if _, err := repo.Store().GetRaw(id); err == nil {
			present = append(present, id)
		} else {
			missing = append(missing, id)
		}
	}

	s.sendJSON(w, http.StatusOK, map[string][]string{
		"present": present,
		"missing": missing,
	})
}

// handleObjectsUpload handles POST /repos/{repoId}/objects/upload.
func (s *Server) handleObjectsUpload(w http.ResponseWriter, r *http.Request, repoID string) {
	if r.Method != http.MethodPost {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// Read body with simple framing: hash\nlength\nbytes
	reader := bufio.NewReader(r.Body)
	var stored []string

	for {
		// Read hash
		hashLine, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			s.sendError(w, http.StatusBadRequest, "bad_request", "Failed to read hash")
			return
		}
		hash := strings.TrimSpace(hashLine)

		// Read length
		lengthLine, err := reader.ReadString('\n')
		if err != nil {
			s.sendError(w, http.StatusBadRequest, "bad_request", "Failed to read length")
			return
		}
		length, err := strconv.Atoi(strings.TrimSpace(lengthLine))
		if err != nil {
			s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid length")
			return
		}

		// Read data
		data := make([]byte, length)
		if _, err := io.ReadFull(reader, data); err != nil {
			s.sendError(w, http.StatusBadRequest, "bad_request", "Failed to read data")
			return
		}

		// Verify hash
		computedHash := model.ComputeBlobHash(data)
		if computedHash != hash {
			s.sendError(w, http.StatusBadRequest, "bad_request",
				fmt.Sprintf("Hash mismatch: expected %s, got %s", hash, computedHash))
			return
		}

		// Store object
		if err := repo.Store().PutRaw(hash, data); err != nil {
			s.sendError(w, http.StatusInternalServerError, "internal_error",
				fmt.Sprintf("Failed to store object: %v", err))
			return
		}

		stored = append(stored, hash)
	}

	s.sendJSON(w, http.StatusOK, map[string][]string{
		"stored": stored,
	})
}

// handleObjectDownload handles GET /repos/{repoId}/objects/{hash}.
func (s *Server) handleObjectDownload(w http.ResponseWriter, r *http.Request, repoID, hash string) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// Get object data
	data, err := repo.Store().GetRaw(hash)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", fmt.Sprintf("Object not found: %s", hash))
		return
	}

	// Send raw bytes
	w.Header().Set("Content-Type", "application/octet-stream")
	w.Header().Set("Content-Length", strconv.Itoa(len(data)))
	w.WriteHeader(http.StatusOK)
	w.Write(data)
}

// handleRefs handles GET /repos/{repoId}/refs - list all refs.
func (s *Server) handleRefs(w http.ResponseWriter, r *http.Request, repoID string) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// List all branches
	branches, err := repo.ListBranches()
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error",
			fmt.Sprintf("Failed to list branches: %v", err))
		return
	}

	// Build ref list
	var refs []map[string]string
	for _, branch := range branches {
		refName := "refs/heads/" + branch
		target, err := repo.ResolveRef(refName)
		if err != nil || target == "" {
			continue
		}

		refs = append(refs, map[string]string{
			"name":   refName,
			"target": target,
		})
	}

	s.sendJSON(w, http.StatusOK, map[string]interface{}{
		"refs": refs,
	})
}

// handleRef handles GET/POST /repos/{repoId}/refs/{refName}.
func (s *Server) handleRef(w http.ResponseWriter, r *http.Request, repoID, refName string) {
	repo, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	// Prepend "refs/" if not present
	fullRefName := refName
	if !strings.HasPrefix(refName, "refs/") {
		fullRefName = "refs/" + refName
	}

	switch r.Method {
	case http.MethodGet:
		s.handleRefGet(w, r, repo, fullRefName)
	case http.MethodPost:
		s.handleRefUpdate(w, r, repo, fullRefName)
	default:
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
	}
}

// handleRefGet handles GET /repos/{repoId}/refs/{refName}.
func (s *Server) handleRefGet(w http.ResponseWriter, r *http.Request, repo *repo.Repository, refName string) {
	target, err := repo.ResolveRef(refName)
	if err != nil || target == "" {
		s.sendError(w, http.StatusNotFound, "not_found", fmt.Sprintf("Ref not found: %s", refName))
		return
	}

	s.sendJSON(w, http.StatusOK, map[string]string{
		"name":   refName,
		"target": target,
	})
}

// handleRefUpdate handles POST /repos/{repoId}/refs/{refName}.
func (s *Server) handleRefUpdate(w http.ResponseWriter, r *http.Request, repo *repo.Repository, refName string) {
	// Parse request
	var req struct {
		OldTarget string `json:"old_target"`
		NewTarget string `json:"new_target"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid JSON")
		return
	}

	// Check if new target exists
	if _, err := repo.GetSnapshot(req.NewTarget); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request",
			fmt.Sprintf("Target snapshot not found: %s", req.NewTarget))
		return
	}

	// Optimistic concurrency control (always enforced).
	// If OldTarget is omitted/empty, the update is only allowed when the ref does not exist.
	currentTarget, _ := repo.ResolveRef(refName)
	if currentTarget != req.OldTarget {
		s.sendError(w, http.StatusConflict, "conflict",
			fmt.Sprintf("Ref has changed: expected %s, got %s", req.OldTarget, currentTarget))
		return
	}

	// Update ref
	if err := repo.UpdateRef(refName, req.NewTarget); err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error",
			fmt.Sprintf("Failed to update ref: %v", err))
		return
	}

	s.sendJSON(w, http.StatusOK, map[string]string{
		"name":   refName,
		"target": req.NewTarget,
	})
}
