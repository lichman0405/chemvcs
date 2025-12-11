package server

import (
	"bytes"
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/repo"
)

func TestNewServer(t *testing.T) {
	tmpDir := t.TempDir()

	config := Config{
		RepoRoot: tmpDir,
		Port:     8080,
	}

	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	if srv.RepoRoot != tmpDir {
		t.Errorf("expected RepoRoot %s, got %s", tmpDir, srv.RepoRoot)
	}
}

func TestListRepositories(t *testing.T) {
	tmpDir := t.TempDir()

	// Create test repositories
	repo1Path := filepath.Join(tmpDir, "owner1", "repo1")
	repo2Path := filepath.Join(tmpDir, "owner2", "repo2")

	for _, path := range []string{repo1Path, repo2Path} {
		os.MkdirAll(path, 0755)
		if _, err := repo.Init(path); err != nil {
			t.Fatalf("Failed to init repo: %v", err)
		}
	}

	config := Config{RepoRoot: tmpDir, Port: 8080}
	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	// Test list repositories endpoint
	req := httptest.NewRequest("GET", "/chemvcs/v1/repos", nil)
	w := httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	var resp struct {
		Repositories []map[string]string `json:"repositories"`
	}
	if err := json.NewDecoder(w.Body).Decode(&resp); err != nil {
		t.Fatalf("Failed to decode response: %v", err)
	}

	if len(resp.Repositories) != 2 {
		t.Errorf("expected 2 repositories, got %d", len(resp.Repositories))
	}
}

func TestRepoInfo(t *testing.T) {
	tmpDir := t.TempDir()

	// Create test repository
	repoPath := filepath.Join(tmpDir, "owner", "repo")
	os.MkdirAll(repoPath, 0755)
	if _, err := repo.Init(repoPath); err != nil {
		t.Fatalf("Failed to init repo: %v", err)
	}

	config := Config{RepoRoot: tmpDir, Port: 8080}
	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	// Test repository info endpoint
	req := httptest.NewRequest("GET", "/chemvcs/v1/repos/owner/repo", nil)
	w := httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	var resp map[string]string
	if err := json.NewDecoder(w.Body).Decode(&resp); err != nil {
		t.Fatalf("Failed to decode response: %v", err)
	}

	if resp["id"] != "owner/repo" {
		t.Errorf("expected id 'owner/repo', got %s", resp["id"])
	}
}

func TestObjectsExist(t *testing.T) {
	tmpDir := t.TempDir()

	// Create test repository
	repoPath := filepath.Join(tmpDir, "test", "repo")
	os.MkdirAll(repoPath, 0755)
	r, err := repo.Init(repoPath)
	if err != nil {
		t.Fatalf("Failed to init repo: %v", err)
	}

	// Store a test blob
	hash1, err := r.Store().PutBlob([]byte("test data 1"))
	if err != nil {
		t.Fatalf("Failed to store blob: %v", err)
	}

	config := Config{RepoRoot: tmpDir, Port: 8080}
	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	// Test objects/exists endpoint
	reqBody, _ := json.Marshal(map[string][]string{
		"ids": {hash1, "nonexistent"},
	})

	req := httptest.NewRequest("POST", "/chemvcs/v1/repos/test/repo/objects/exists", bytes.NewReader(reqBody))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	var resp map[string][]string
	if err := json.NewDecoder(w.Body).Decode(&resp); err != nil {
		t.Fatalf("Failed to decode response: %v", err)
	}

	if len(resp["present"]) != 1 || resp["present"][0] != hash1 {
		t.Errorf("expected present=[%s], got %v", hash1, resp["present"])
	}
	if len(resp["missing"]) != 1 || resp["missing"][0] != "nonexistent" {
		t.Errorf("expected missing=['nonexistent'], got %v", resp["missing"])
	}
}

func TestObjectUploadDownload(t *testing.T) {
	tmpDir := t.TempDir()

	// Create test repository
	repoPath := filepath.Join(tmpDir, "test", "repo")
	os.MkdirAll(repoPath, 0755)
	if _, err := repo.Init(repoPath); err != nil {
		t.Fatalf("Failed to init repo: %v", err)
	}

	config := Config{RepoRoot: tmpDir, Port: 8080}
	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	// Upload an object
	testData := []byte("test object data")
	testHash := "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" // SHA-256 of empty string - will be computed

	// Actually compute the correct hash
	r, _ := repo.Open(repoPath)
	testHash, _ = r.Store().PutBlob(testData)

	var uploadBody bytes.Buffer
	uploadBody.WriteString(testHash + "\n")
	uploadBody.WriteString("16\n") // len("test object data")
	uploadBody.Write(testData)

	req := httptest.NewRequest("POST", "/chemvcs/v1/repos/test/repo/objects/upload", &uploadBody)
	req.Header.Set("Content-Type", "application/octet-stream")
	w := httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d: %s", w.Code, w.Body.String())
	}

	// Download the object
	req = httptest.NewRequest("GET", "/chemvcs/v1/repos/test/repo/objects/"+testHash, nil)
	w = httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	downloadedData, _ := io.ReadAll(w.Body)
	if !bytes.Equal(downloadedData, testData) {
		t.Errorf("downloaded data mismatch")
	}
}

func TestRefOperations(t *testing.T) {
	tmpDir := t.TempDir()

	// Create test repository with a commit
	repoPath := filepath.Join(tmpDir, "test", "repo")
	os.MkdirAll(repoPath, 0755)
	r, err := repo.Init(repoPath)
	if err != nil {
		t.Fatalf("Failed to init repo: %v", err)
	}

	// Create a commit
	blobHash, _ := r.Store().PutBlob([]byte("test"))
	obj := model.NewObject("folder")
	obj.AddRef(model.NewBlobRef(blobHash))
	obj.SetMeta("name", "root")
	objHash, _ := r.Store().PutObject(obj)
	snapHash, _ := r.CreateSnapshot(objHash, "Test commit", "Test <test@example.com>")

	config := Config{RepoRoot: tmpDir, Port: 8080}
	srv, err := NewServer(config)
	if err != nil {
		t.Fatalf("NewServer failed: %v", err)
	}

	// Test list refs
	req := httptest.NewRequest("GET", "/chemvcs/v1/repos/test/repo/refs", nil)
	w := httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	var listResp struct {
		Refs []map[string]string `json:"refs"`
	}
	if err := json.NewDecoder(w.Body).Decode(&listResp); err != nil {
		t.Fatalf("Failed to decode response: %v", err)
	}

	if len(listResp.Refs) != 1 {
		t.Errorf("expected 1 ref, got %d", len(listResp.Refs))
	}

	// Test get ref
	req = httptest.NewRequest("GET", "/chemvcs/v1/repos/test/repo/refs/heads/main", nil)
	w = httptest.NewRecorder()

	srv.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("expected status 200, got %d", w.Code)
	}

	var getResp map[string]string
	if err := json.NewDecoder(w.Body).Decode(&getResp); err != nil {
		t.Fatalf("Failed to decode response: %v", err)
	}

	if getResp["target"] != snapHash {
		t.Errorf("expected target %s, got %s", snapHash, getResp["target"])
	}
}
