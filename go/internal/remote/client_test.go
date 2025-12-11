package remote

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"testing"
)

func TestNewClient(t *testing.T) {
	client := NewClient("https://example.com/api", "test-token")

	if client.BaseURL != "https://example.com/api/" {
		t.Errorf("expected BaseURL to have trailing slash, got %q", client.BaseURL)
	}
	if client.Token != "test-token" {
		t.Errorf("expected token %q, got %q", "test-token", client.Token)
	}
}

func TestGetRepositoryInfo(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path != "/repos/test/repo" {
			t.Errorf("unexpected path: %s", r.URL.Path)
		}
		if r.Method != "GET" {
			t.Errorf("unexpected method: %s", r.Method)
		}

		info := RepositoryInfo{
			ID:            "test/repo",
			Description:   "Test repository",
			DefaultBranch: "main",
		}
		json.NewEncoder(w).Encode(info)
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	info, err := client.GetRepositoryInfo("test/repo")
	if err != nil {
		t.Fatalf("GetRepositoryInfo failed: %v", err)
	}

	if info.ID != "test/repo" {
		t.Errorf("expected ID %q, got %q", "test/repo", info.ID)
	}
	if info.Description != "Test repository" {
		t.Errorf("expected description %q, got %q", "Test repository", info.Description)
	}
	if info.DefaultBranch != "main" {
		t.Errorf("expected default branch %q, got %q", "main", info.DefaultBranch)
	}
}

func TestListRefs(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path != "/repos/test/repo/refs" {
			t.Errorf("unexpected path: %s", r.URL.Path)
		}

		result := struct {
			Refs []RefInfo `json:"refs"`
		}{
			Refs: []RefInfo{
				{Name: "refs/heads/main", Target: "abc123"},
				{Name: "refs/heads/feature", Target: "def456"},
			},
		}
		json.NewEncoder(w).Encode(result)
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	refs, err := client.ListRefs("test/repo")
	if err != nil {
		t.Fatalf("ListRefs failed: %v", err)
	}

	if len(refs) != 2 {
		t.Fatalf("expected 2 refs, got %d", len(refs))
	}
	if refs[0].Name != "refs/heads/main" || refs[0].Target != "abc123" {
		t.Errorf("unexpected first ref: %+v", refs[0])
	}
}

func TestGetRef(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path != "/repos/test/repo/refs/refs/heads/main" {
			t.Errorf("unexpected path: %s", r.URL.Path)
		}

		ref := RefInfo{
			Name:   "refs/heads/main",
			Target: "abc123",
		}
		json.NewEncoder(w).Encode(ref)
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	ref, err := client.GetRef("test/repo", "refs/heads/main")
	if err != nil {
		t.Fatalf("GetRef failed: %v", err)
	}

	if ref.Name != "refs/heads/main" {
		t.Errorf("expected name %q, got %q", "refs/heads/main", ref.Name)
	}
	if ref.Target != "abc123" {
		t.Errorf("expected target %q, got %q", "abc123", ref.Target)
	}
}

func TestUpdateRef(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method != "POST" {
			t.Errorf("expected POST, got %s", r.Method)
		}
		if r.URL.Path != "/repos/test/repo/refs/refs/heads/main" {
			t.Errorf("unexpected path: %s", r.URL.Path)
		}

		var req RefUpdateRequest
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			t.Errorf("decoding request: %v", err)
		}

		if req.OldTarget != "old123" {
			t.Errorf("expected old_target %q, got %q", "old123", req.OldTarget)
		}
		if req.NewTarget != "new456" {
			t.Errorf("expected new_target %q, got %q", "new456", req.NewTarget)
		}

		w.WriteHeader(http.StatusOK)
		json.NewEncoder(w).Encode(RefInfo{Name: "refs/heads/main", Target: "new456"})
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	err := client.UpdateRef("test/repo", "refs/heads/main", "old123", "new456")
	if err != nil {
		t.Fatalf("UpdateRef failed: %v", err)
	}
}

func TestCheckObjectsExist(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method != "POST" {
			t.Errorf("expected POST, got %s", r.Method)
		}

		var req ObjectExistsRequest
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			t.Errorf("decoding request: %v", err)
		}

		if len(req.IDs) != 3 {
			t.Errorf("expected 3 IDs, got %d", len(req.IDs))
		}

		resp := ObjectExistsResponse{
			Present: []string{"abc123"},
			Missing: []string{"def456", "ghi789"},
		}
		json.NewEncoder(w).Encode(resp)
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	result, err := client.CheckObjectsExist("test/repo", []string{"abc123", "def456", "ghi789"})
	if err != nil {
		t.Fatalf("CheckObjectsExist failed: %v", err)
	}

	if len(result.Present) != 1 || result.Present[0] != "abc123" {
		t.Errorf("unexpected present list: %v", result.Present)
	}
	if len(result.Missing) != 2 {
		t.Errorf("expected 2 missing objects, got %d", len(result.Missing))
	}
}

func TestUploadObject(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method != "POST" {
			t.Errorf("expected POST, got %s", r.Method)
		}
		if r.Header.Get("Content-Type") != "application/octet-stream" {
			t.Errorf("expected Content-Type application/octet-stream")
		}

		// Read and validate framing
		// Simple validation: should have hash\nlength\ndata format
		w.WriteHeader(http.StatusOK)
		json.NewEncoder(w).Encode(map[string]interface{}{"stored": []string{"abc123"}})
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	err := client.UploadObject("test/repo", "abc123", []byte("test data"))
	if err != nil {
		t.Fatalf("UploadObject failed: %v", err)
	}
}

func TestDownloadObject(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method != "GET" {
			t.Errorf("expected GET, got %s", r.Method)
		}
		if r.URL.Path != "/repos/test/repo/objects/abc123" {
			t.Errorf("unexpected path: %s", r.URL.Path)
		}

		w.Header().Set("Content-Type", "application/octet-stream")
		w.Write([]byte("test data"))
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	data, err := client.DownloadObject("test/repo", "abc123")
	if err != nil {
		t.Fatalf("DownloadObject failed: %v", err)
	}

	if string(data) != "test data" {
		t.Errorf("expected %q, got %q", "test data", string(data))
	}
}

func TestErrorHandling(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusNotFound)
		json.NewEncoder(w).Encode(ErrorResponse{
			Error:   "not_found",
			Message: "Repository not found",
		})
	}))
	defer server.Close()

	client := NewClient(server.URL, "")
	_, err := client.GetRepositoryInfo("nonexistent/repo")
	if err == nil {
		t.Fatal("expected error, got nil")
	}

	if !contains(err.Error(), "not_found") || !contains(err.Error(), "Repository not found") {
		t.Errorf("unexpected error message: %v", err)
	}
}

func contains(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr || len(s) > len(substr) && containsAt(s, substr))
}

func containsAt(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
