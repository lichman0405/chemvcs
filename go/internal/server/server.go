// Package server implements the ChemVCS HTTP server.
package server

import (
	"encoding/json"
	"fmt"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"strings"

	"github.com/lishi/chemvcs/internal/repo"
)

// Server represents a ChemVCS HTTP server.
type Server struct {
	RepoRoot string                      // Root directory containing repositories
	Repos    map[string]*repo.Repository // Cached repository instances
	mux      *http.ServeMux              // HTTP request multiplexer
}

// Config holds server configuration.
type Config struct {
	RepoRoot string // Root directory for repositories
	Port     int    // Port to listen on
}

// NewServer creates a new ChemVCS server.
func NewServer(config Config) (*Server, error) {
	// Ensure repository root exists
	if err := os.MkdirAll(config.RepoRoot, 0755); err != nil {
		return nil, fmt.Errorf("failed to create repo root: %w", err)
	}

	s := &Server{
		RepoRoot: config.RepoRoot,
		Repos:    make(map[string]*repo.Repository),
		mux:      http.NewServeMux(),
	}

	s.setupRoutes()
	return s, nil
}

// setupRoutes configures HTTP routes.
func (s *Server) setupRoutes() {
	// Repository endpoints
	s.mux.HandleFunc("/chemvcs/v1/repos", s.handleRepos)
	s.mux.HandleFunc("/chemvcs/v1/repos/", s.handleRepoRoutes)
}

// ServeHTTP implements http.Handler.
func (s *Server) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	log.Printf("%s %s", r.Method, r.URL.Path)
	s.mux.ServeHTTP(w, r)
}

// Start starts the HTTP server on the configured port.
func (s *Server) Start(port int) error {
	addr := fmt.Sprintf(":%d", port)
	log.Printf("ChemVCS server starting on %s", addr)
	log.Printf("Repository root: %s", s.RepoRoot)
	return http.ListenAndServe(addr, s)
}

// handleRepos handles GET /repos (list all repositories).
func (s *Server) handleRepos(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	repos, err := s.listRepositories()
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", err.Error())
		return
	}

	s.sendJSON(w, http.StatusOK, map[string]interface{}{
		"repositories": repos,
	})
}

// handleRepoRoutes routes requests under /repos/{repoId}/...
func (s *Server) handleRepoRoutes(w http.ResponseWriter, r *http.Request) {
	// Parse path: /chemvcs/v1/repos/{owner}/{repo}/{resource}...
	path := strings.TrimPrefix(r.URL.Path, "/chemvcs/v1/repos/")
	parts := strings.Split(path, "/")

	if len(parts) < 2 {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid repository path")
		return
	}

	repoID := parts[0] + "/" + parts[1]
	resource := ""
	if len(parts) > 2 {
		resource = strings.Join(parts[2:], "/")
	}

	// Route based on resource
	switch {
	case resource == "":
		// GET /repos/{repoId} - repository info
		s.handleRepoInfo(w, r, repoID)
	case strings.HasPrefix(resource, "objects/exists"):
		s.handleObjectsExist(w, r, repoID)
	case strings.HasPrefix(resource, "objects/upload"):
		s.handleObjectsUpload(w, r, repoID)
	case strings.HasPrefix(resource, "objects/"):
		// GET /repos/{repoId}/objects/{hash}
		hash := strings.TrimPrefix(resource, "objects/")
		s.handleObjectDownload(w, r, repoID, hash)
	case resource == "refs":
		s.handleRefs(w, r, repoID)
	case strings.HasPrefix(resource, "refs/"):
		refName := strings.TrimPrefix(resource, "refs/")
		s.handleRef(w, r, repoID, refName)
	default:
		s.sendError(w, http.StatusNotFound, "not_found", "Resource not found")
	}
}

// getRepository retrieves a repository instance, opening it if necessary.
func (s *Server) getRepository(repoID string) (*repo.Repository, error) {
	// Check cache
	if r, ok := s.Repos[repoID]; ok {
		return r, nil
	}

	// Open repository
	repoPath := filepath.Join(s.RepoRoot, repoID)
	if _, err := os.Stat(filepath.Join(repoPath, ".chemvcs")); os.IsNotExist(err) {
		return nil, fmt.Errorf("repository not found: %s", repoID)
	}

	r, err := repo.Open(repoPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open repository: %w", err)
	}

	// Cache
	s.Repos[repoID] = r
	return r, nil
}

// listRepositories returns a list of all repositories in the root directory.
func (s *Server) listRepositories() ([]map[string]string, error) {
	var repos []map[string]string

	// Walk the repository root directory
	err := filepath.Walk(s.RepoRoot, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}

		// Check if this is a .chemvcs directory
		if info.IsDir() && info.Name() == ".chemvcs" {
			// Get repository ID (relative path from root)
			repoPath := filepath.Dir(path)
			relPath, err := filepath.Rel(s.RepoRoot, repoPath)
			if err != nil {
				return err
			}

			// Normalize to forward slashes for repository ID
			repoID := strings.ReplaceAll(relPath, string(filepath.Separator), "/")

			repos = append(repos, map[string]string{
				"id":          repoID,
				"description": "ChemVCS repository",
			})

			// Don't descend into .chemvcs
			return filepath.SkipDir
		}

		return nil
	})

	if err != nil {
		return nil, err
	}

	return repos, nil
}

// sendJSON sends a JSON response.
func (s *Server) sendJSON(w http.ResponseWriter, status int, data interface{}) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(status)
	json.NewEncoder(w).Encode(data)
}

// sendError sends an error response.
func (s *Server) sendError(w http.ResponseWriter, status int, errorCode, message string) {
	s.sendJSON(w, status, map[string]string{
		"error":   errorCode,
		"message": message,
	})
}
