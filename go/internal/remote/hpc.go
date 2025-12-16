package remote

import (
	"archive/zip"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"net/url"
	"os"
	"path/filepath"
	"strings"
)

// HPCJobInfo mirrors the server-side job representation.
// Keep this decoupled from internal/hpc to avoid package dependency tangles.
type HPCJobInfo struct {
	JobID      string `json:"job_id"`
	RunHash    string `json:"run_hash"`
	Status     string `json:"status"`
	JobSystem  string `json:"job_system"`
	Queue      string `json:"queue_name"`
	Submitted  string `json:"submitted_at"`
	Updated    string `json:"updated_at"`
	WorkingDir string `json:"working_dir"`
}

type SubmitHPCRequest struct {
	RunHash    string `json:"run_hash"`
	Script     string `json:"script"`
	CaptureEnv bool   `json:"capture_env"`
}

type SubmitHPCResponse struct {
	JobID     string `json:"job_id"`
	RunHash   string `json:"run_hash"`
	JobSystem string `json:"job_system"`
}

type CancelHPCRequest struct {
	Identifier string `json:"identifier"`
}

type CancelHPCResponse struct {
	JobID  string `json:"job_id"`
	Status string `json:"status"`
}

type RetrieveHPCRequest struct {
	Patterns []string `json:"patterns"`
}

// SubmitHPC submits a job on the remote server.
func (c *Client) SubmitHPC(repoID string, req SubmitHPCRequest) (*SubmitHPCResponse, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return nil, err
	}

	data, err := json.Marshal(req)
	if err != nil {
		return nil, fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/hpc/submit", bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var out SubmitHPCResponse
	if err := json.NewDecoder(resp.Body).Decode(&out); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}
	return &out, nil
}

// ListHPCJobs lists tracked jobs for the repository.
func (c *Client) ListHPCJobs(repoID, status string, refresh bool) ([]HPCJobInfo, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return nil, err
	}

	q := ""
	if status != "" {
		q += "status=" + url.QueryEscape(status)
	}
	if refresh {
		if q != "" {
			q += "&"
		}
		q += "refresh=1"
	}

	path := "repos/" + repoID + "/hpc/jobs"
	if q != "" {
		path += "?" + q
	}

	resp, err := c.doRequest("GET", path, nil)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var result struct {
		Jobs []HPCJobInfo `json:"jobs"`
	}
	if err := json.NewDecoder(resp.Body).Decode(&result); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}
	return result.Jobs, nil
}

// GetHPCJob retrieves a single job record by run hash or job id.
func (c *Client) GetHPCJob(repoID, identifier string, refresh bool) (*HPCJobInfo, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return nil, err
	}

	path := "repos/" + repoID + "/hpc/jobs/" + url.PathEscape(identifier)
	if refresh {
		path += "?refresh=1"
	}

	resp, err := c.doRequest("GET", path, nil)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var out HPCJobInfo
	if err := json.NewDecoder(resp.Body).Decode(&out); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}
	return &out, nil
}

// CancelHPC cancels a job by run hash or job id.
func (c *Client) CancelHPC(repoID, identifier string) (*CancelHPCResponse, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return nil, err
	}

	data, err := json.Marshal(CancelHPCRequest{Identifier: identifier})
	if err != nil {
		return nil, fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/hpc/cancel", bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var out CancelHPCResponse
	if err := json.NewDecoder(resp.Body).Decode(&out); err != nil {
		return nil, fmt.Errorf("decoding response: %w", err)
	}
	return &out, nil
}

// RetrieveHPC downloads a zip of files matching patterns.
func (c *Client) RetrieveHPC(repoID string, patterns []string) ([]byte, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return nil, err
	}

	data, err := json.Marshal(RetrieveHPCRequest{Patterns: patterns})
	if err != nil {
		return nil, fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/hpc/retrieve", bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	return io.ReadAll(resp.Body)
}

// RetrieveHPCTo streams a zip of files matching patterns into w.
// This avoids holding the entire zip in memory.
func (c *Client) RetrieveHPCTo(repoID string, patterns []string, w io.Writer) (int64, error) {
	repoID, err := c.normalizeRepoID(repoID)
	if err != nil {
		return 0, err
	}

	data, err := json.Marshal(RetrieveHPCRequest{Patterns: patterns})
	if err != nil {
		return 0, fmt.Errorf("encoding request: %w", err)
	}

	resp, err := c.doRequest("POST", "repos/"+repoID+"/hpc/retrieve", bytes.NewReader(data))
	if err != nil {
		return 0, err
	}
	defer resp.Body.Close()

	return io.Copy(w, resp.Body)
}

// ExtractZipTo extracts zip data into destination directory.
func ExtractZipTo(zipData []byte, destDir string) ([]string, error) {
	if err := os.MkdirAll(destDir, 0755); err != nil {
		return nil, fmt.Errorf("creating destination: %w", err)
	}

	zr, err := zip.NewReader(bytes.NewReader(zipData), int64(len(zipData)))
	if err != nil {
		return nil, fmt.Errorf("reading zip: %w", err)
	}

	var written []string
	for _, f := range zr.File {
		name := filepath.Clean(f.Name)
		name = strings.TrimPrefix(name, string(filepath.Separator))
		if name == "." || name == "" {
			continue
		}
		// Prevent zip slip.
		outPath := filepath.Join(destDir, name)
		if !strings.HasPrefix(filepath.Clean(outPath)+string(filepath.Separator), filepath.Clean(destDir)+string(filepath.Separator)) {
			continue
		}

		if f.FileInfo().IsDir() {
			_ = os.MkdirAll(outPath, 0755)
			continue
		}

		if err := os.MkdirAll(filepath.Dir(outPath), 0755); err != nil {
			return written, fmt.Errorf("creating parent dir: %w", err)
		}

		rc, err := f.Open()
		if err != nil {
			return written, fmt.Errorf("opening zip entry: %w", err)
		}

		outFile, err := os.Create(outPath)
		if err != nil {
			_ = rc.Close()
			return written, fmt.Errorf("creating file: %w", err)
		}

		_, copyErr := io.Copy(outFile, rc)
		_ = outFile.Close()
		_ = rc.Close()
		if copyErr != nil {
			return written, fmt.Errorf("writing file: %w", copyErr)
		}

		written = append(written, outPath)
	}

	return written, nil
}

// ExtractZipFileTo extracts a zip file on disk into destination directory.
func ExtractZipFileTo(zipPath string, destDir string) ([]string, error) {
	if err := os.MkdirAll(destDir, 0755); err != nil {
		return nil, fmt.Errorf("creating destination: %w", err)
	}

	zr, err := zip.OpenReader(zipPath)
	if err != nil {
		return nil, fmt.Errorf("reading zip: %w", err)
	}
	defer zr.Close()

	var written []string
	for _, f := range zr.File {
		name := filepath.Clean(f.Name)
		name = strings.TrimPrefix(name, string(filepath.Separator))
		if name == "." || name == "" {
			continue
		}
		// Prevent zip slip.
		outPath := filepath.Join(destDir, name)
		if !strings.HasPrefix(filepath.Clean(outPath)+string(filepath.Separator), filepath.Clean(destDir)+string(filepath.Separator)) {
			continue
		}

		if f.FileInfo().IsDir() {
			_ = os.MkdirAll(outPath, 0755)
			continue
		}

		if err := os.MkdirAll(filepath.Dir(outPath), 0755); err != nil {
			return written, fmt.Errorf("creating parent dir: %w", err)
		}

		rc, err := f.Open()
		if err != nil {
			return written, fmt.Errorf("opening zip entry: %w", err)
		}

		outFile, err := os.Create(outPath)
		if err != nil {
			_ = rc.Close()
			return written, fmt.Errorf("creating file: %w", err)
		}

		_, copyErr := io.Copy(outFile, rc)
		_ = outFile.Close()
		_ = rc.Close()
		if copyErr != nil {
			return written, fmt.Errorf("writing file: %w", copyErr)
		}

		written = append(written, outPath)
	}

	return written, nil
}
