// Package server implements the ChemVCS HTTP server.
package server

import (
	"crypto/rand"
	"crypto/subtle"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/lishi/chemvcs/internal/repo"
)

// Server represents a ChemVCS HTTP server.
type Server struct {
	RepoRoot string                      // Root directory containing repositories
	Repos    map[string]*repo.Repository // Cached repository instances
	mux      *http.ServeMux              // HTTP request multiplexer

	AuthToken    string
	AuthRepos    map[string]struct{}
	AuthAllowAll bool
	AdminToken   string
	lockByRepoID sync.Map // map[string]*sync.Mutex
}

// Config holds server configuration.
type Config struct {
	RepoRoot string // Root directory for repositories
	Port     int    // Port to listen on

	// Authentication and authorization (optional).
	// If AuthToken is empty, the server runs without authentication.
	AuthToken  string
	AuthRepos  []string // Allowed repos for AuthToken (e.g. []{"owner/repo"} or []{"*"})
	AdminToken string   // Optional token that bypasses repo scoping and allows listing repos
}

// NewServer creates a new ChemVCS server.
func NewServer(config Config) (*Server, error) {
	// Ensure repository root exists
	if err := os.MkdirAll(config.RepoRoot, 0755); err != nil {
		return nil, fmt.Errorf("failed to create repo root: %w", err)
	}

	s := &Server{
		RepoRoot:   config.RepoRoot,
		Repos:      make(map[string]*repo.Repository),
		mux:        http.NewServeMux(),
		AuthToken:  strings.TrimSpace(config.AuthToken),
		AdminToken: strings.TrimSpace(config.AdminToken),
	}

	// Parse allowed repos.
	if s.AuthToken != "" {
		s.AuthRepos = make(map[string]struct{})
		for _, r := range config.AuthRepos {
			r = strings.TrimSpace(r)
			if r == "" {
				continue
			}
			if r == "*" {
				s.AuthAllowAll = true
				continue
			}
			s.AuthRepos[r] = struct{}{}
		}
		if !s.AuthAllowAll && len(s.AuthRepos) == 0 {
			// Keep behavior predictable: if auth is enabled but no repos are configured,
			// deny all repo-scoped access (admin token can still be used).
			log.Printf("Auth enabled but no allowed repos configured; denying repo access unless admin token is used")
		}
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
	start := time.Now()

	requestID := strings.TrimSpace(r.Header.Get("X-Request-Id"))
	if requestID == "" {
		buf := make([]byte, 16)
		if _, err := rand.Read(buf); err == nil {
			requestID = hex.EncodeToString(buf)
		} else {
			requestID = fmt.Sprintf("req-%d", time.Now().UnixNano())
		}
	}
	w.Header().Set("X-Request-Id", requestID)

	sw := &statusWriter{ResponseWriter: w, status: http.StatusOK}

	repoID, repoScoped, isRepoList := parseRepoFromPath(r.URL.Path)
	authPrincipal := "anonymous"

	if s.AuthToken != "" || s.AdminToken != "" {
		tok, ok := parseBearerToken(r.Header.Get("Authorization"))
		if !ok {
			s.sendError(sw, http.StatusUnauthorized, "unauthorized", "missing or invalid Authorization header")
			s.auditLog(r, repoID, requestID, authPrincipal, sw.status, time.Since(start))
			return
		}

		isAdmin := s.AdminToken != "" && subtle.ConstantTimeCompare([]byte(tok), []byte(s.AdminToken)) == 1
		isUser := s.AuthToken != "" && subtle.ConstantTimeCompare([]byte(tok), []byte(s.AuthToken)) == 1
		if !isAdmin && !isUser {
			s.sendError(sw, http.StatusUnauthorized, "unauthorized", "invalid token")
			s.auditLog(r, repoID, requestID, authPrincipal, sw.status, time.Since(start))
			return
		}

		if isAdmin {
			authPrincipal = "admin"
		} else {
			authPrincipal = "token"
		}

		// Authorization: repo scoping and repo list.
		if isRepoList && !isAdmin {
			s.sendError(sw, http.StatusForbidden, "forbidden", "repo listing requires admin token")
			s.auditLog(r, repoID, requestID, authPrincipal, sw.status, time.Since(start))
			return
		}
		if repoScoped && !isAdmin {
			if !s.isRepoAllowed(repoID) {
				s.sendError(sw, http.StatusForbidden, "forbidden", "token not authorized for this repo")
				s.auditLog(r, repoID, requestID, authPrincipal, sw.status, time.Since(start))
				return
			}
		}
	}

	// Concurrency control: serialize mutations per repo.
	if repoScoped && shouldLockRequest(r) {
		mu := s.getRepoLock(repoID)
		mu.Lock()
		defer mu.Unlock()
	}

	log.Printf("request_id=%s method=%s path=%s repo=%s", requestID, r.Method, r.URL.Path, repoID)
	s.mux.ServeHTTP(sw, r)
	s.auditLog(r, repoID, requestID, authPrincipal, sw.status, time.Since(start))
}

type statusWriter struct {
	http.ResponseWriter
	status int
}

func (w *statusWriter) WriteHeader(statusCode int) {
	w.status = statusCode
	w.ResponseWriter.WriteHeader(statusCode)
}

func (s *Server) auditLog(r *http.Request, repoID, requestID, principal string, status int, dur time.Duration) {
	log.Printf(
		"request_id=%s principal=%s method=%s path=%s repo=%s status=%d duration_ms=%d",
		requestID,
		principal,
		r.Method,
		r.URL.Path,
		repoID,
		status,
		dur.Milliseconds(),
	)
}

func parseRepoFromPath(path string) (repoID string, repoScoped bool, isRepoList bool) {
	if path == "/chemvcs/v1/repos" || path == "/chemvcs/v1/repos/" {
		return "", false, true
	}
	prefix := "/chemvcs/v1/repos/"
	if !strings.HasPrefix(path, prefix) {
		return "", false, false
	}
	rest := strings.TrimPrefix(path, prefix)
	rest = strings.Trim(rest, "/")
	parts := strings.Split(rest, "/")
	if len(parts) < 2 {
		return "", false, false
	}
	return parts[0] + "/" + parts[1], true, false
}

func parseBearerToken(authHeader string) (string, bool) {
	authHeader = strings.TrimSpace(authHeader)
	if authHeader == "" {
		return "", false
	}
	const prefix = "Bearer "
	if !strings.HasPrefix(authHeader, prefix) {
		return "", false
	}
	tok := strings.TrimSpace(strings.TrimPrefix(authHeader, prefix))
	if tok == "" {
		return "", false
	}
	return tok, true
}

func (s *Server) isRepoAllowed(repoID string) bool {
	if s.AuthAllowAll {
		return true
	}
	_, ok := s.AuthRepos[repoID]
	return ok
}

func (s *Server) getRepoLock(repoID string) *sync.Mutex {
	if repoID == "" {
		// Should not happen for scoped requests, but keep safe.
		return &sync.Mutex{}
	}
	if v, ok := s.lockByRepoID.Load(repoID); ok {
		return v.(*sync.Mutex)
	}
	mu := &sync.Mutex{}
	actual, _ := s.lockByRepoID.LoadOrStore(repoID, mu)
	return actual.(*sync.Mutex)
}

func shouldLockRequest(r *http.Request) bool {
	if r.Method != http.MethodGet && r.Method != http.MethodHead {
		return true
	}
	// Some GETs have side effects (e.g., refresh=1 updates HPC job records).
	if r.URL.Query().Get("refresh") == "1" || strings.EqualFold(r.URL.Query().Get("refresh"), "true") {
		return true
	}
	return false
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
	case strings.HasPrefix(resource, "hpc/") || resource == "hpc":
		s.handleHPC(w, r, repoID, resource)
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
