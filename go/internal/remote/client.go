// Package remote provides client and server functionality for remote ChemVCS repositories.
package remote

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"strings"
)

// Client represents a client for interacting with a remote ChemVCS server.
type Client struct {
	BaseURL string       // Base URL of the remote server (e.g., "https://example.com/chemvcs/v1")
	Token   string       // Authentication token (optional)
	HTTP    *http.Client // HTTP client for requests
}

// NewClient creates a new remote client with the given base URL and token.
func NewClient(baseURL, token string) *Client {
	if !strings.HasSuffix(baseURL, "/") {
		baseURL += "/"
	}
	return &Client{
		BaseURL: baseURL,
		Token:   token,
		HTTP:    &http.Client{},
	}
}

// RepositoryInfo represents basic information about a remote repository.
type RepositoryInfo struct {
	ID            string `json:"id"`
	Description   string `json:"description"`
	DefaultBranch string `json:"default_branch"`
}

// RefInfo represents a remote ref (branch).
type RefInfo struct {
	Name   string `json:"name"`
	Target string `json:"target"` // Snapshot hash
}

// ObjectExistsRequest represents a request to check object existence.
type ObjectExistsRequest struct {
	IDs []string `json:"ids"`
}

// ObjectExistsResponse represents the response from object existence check.
type ObjectExistsResponse struct {
	Present []string `json:"present"`
	Missing []string `json:"missing"`
}

// RefUpdateRequest represents a request to update a ref.
type RefUpdateRequest struct {
	OldTarget string `json:"old_target,omitempty"` // Optional, for concurrency control
	NewTarget string `json:"new_target"`
}

// ErrorResponse represents an error response from the server.
type ErrorResponse struct {
	Error   string `json:"error"`
	Message string `json:"message"`
}

// doRequest performs an HTTP request and handles common error cases.
func (c *Client) doRequest(method, path string, body io.Reader) (*http.Response, error) {
	url := c.BaseURL + strings.TrimPrefix(path, "/")

	req, err := http.NewRequest(method, url, body)
	if err != nil {
		return nil, fmt.Errorf("creating request: %w", err)
	}

	if c.Token != "" {
		req.Header.Set("Authorization", "Bearer "+c.Token)
	}

	if body != nil {
		req.Header.Set("Content-Type", "application/json")
	}

	resp, err := c.HTTP.Do(req)
	if err != nil {
		return nil, fmt.Errorf("executing request: %w", err)
	}

	// Handle HTTP error status codes
	if resp.StatusCode >= 400 {
		defer resp.Body.Close()

		var errResp ErrorResponse
		if err := json.NewDecoder(resp.Body).Decode(&errResp); err != nil {
			return nil, fmt.Errorf("HTTP %d: %s", resp.StatusCode, resp.Status)
		}
		return nil, fmt.Errorf("server error (%s): %s", errResp.Error, errResp.Message)
	}

	return resp, nil
}

// GetRepositoryInfo retrieves information about a specific repository.
func (c *Client) GetRepositoryInfo(repoID string) (*RepositoryInfo, error) {
	resp, err := c.doRequest("GET", "repos/"+repoID, nil)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var info RepositoryInfo
	if err := json.NewDecoder(resp.Body).Decode(&info); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}

	return &info, nil
}

// ListRefs retrieves all refs in a repository.
func (c *Client) ListRefs(repoID string) ([]RefInfo, error) {
	resp, err := c.doRequest("GET", "repos/"+repoID+"/refs", nil)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var result struct {
		Refs []RefInfo `json:"refs"`
	}
	if err := json.NewDecoder(resp.Body).Decode(&result); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}

	return result.Refs, nil
}

// GetRef retrieves a specific ref from a repository.
func (c *Client) GetRef(repoID, refName string) (*RefInfo, error) {
	resp, err := c.doRequest("GET", "repos/"+repoID+"/refs/"+refName, nil)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var ref RefInfo
	if err := json.NewDecoder(resp.Body).Decode(&ref); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}

	return &ref, nil
}

// UpdateRef updates a ref on the remote server.
func (c *Client) UpdateRef(repoID, refName string, oldTarget, newTarget string) error {
	reqBody := RefUpdateRequest{
		OldTarget: oldTarget,
		NewTarget: newTarget,
	}

	data, err := json.Marshal(reqBody)
	if err != nil {
		return fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/refs/"+refName, bytes.NewReader(data))
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	return nil
}

// CheckObjectsExist checks which objects exist on the remote server.
func (c *Client) CheckObjectsExist(repoID string, objectIDs []string) (*ObjectExistsResponse, error) {
	reqBody := ObjectExistsRequest{IDs: objectIDs}

	data, err := json.Marshal(reqBody)
	if err != nil {
		return nil, fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/objects/exists", bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var result ObjectExistsResponse
	if err := json.NewDecoder(resp.Body).Decode(&result); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}

	return &result, nil
}

// UploadObject uploads a single object to the remote server.
func (c *Client) UploadObject(repoID, objectID string, data []byte) error {
	// Simple framing: hash\nlength\nbytes
	var buf bytes.Buffer
	buf.WriteString(objectID + "\n")
	buf.WriteString(fmt.Sprintf("%d\n", len(data)))
	buf.Write(data)

	url := c.BaseURL + "repos/" + repoID + "/objects/upload"

	req, err := http.NewRequest("POST", url, &buf)
	if err != nil {
		return fmt.Errorf("creating request: %w", err)
	}

	if c.Token != "" {
		req.Header.Set("Authorization", "Bearer "+c.Token)
	}
	req.Header.Set("Content-Type", "application/octet-stream")

	resp, err := c.HTTP.Do(req)
	if err != nil {
		return fmt.Errorf("executing request: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode >= 400 {
		var errResp ErrorResponse
		if err := json.NewDecoder(resp.Body).Decode(&errResp); err != nil {
			return fmt.Errorf("HTTP %d: %s", resp.StatusCode, resp.Status)
		}
		return fmt.Errorf("server error (%s): %s", errResp.Error, errResp.Message)
	}

	return nil
}

// DownloadObject downloads a single object from the remote server.
func (c *Client) DownloadObject(repoID, objectID string) ([]byte, error) {
	url := c.BaseURL + "repos/" + repoID + "/objects/" + objectID

	req, err := http.NewRequest("GET", url, nil)
	if err != nil {
		return nil, fmt.Errorf("creating request: %w", err)
	}

	if c.Token != "" {
		req.Header.Set("Authorization", "Bearer "+c.Token)
	}

	resp, err := c.HTTP.Do(req)
	if err != nil {
		return nil, fmt.Errorf("executing request: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode >= 400 {
		var errResp ErrorResponse
		if err := json.NewDecoder(resp.Body).Decode(&errResp); err != nil {
			return nil, fmt.Errorf("HTTP %d: %s", resp.StatusCode, resp.Status)
		}
		return nil, fmt.Errorf("server error (%s): %s", errResp.Error, errResp.Message)
	}

	data, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("reading response: %w", err)
	}

	return data, nil
}
