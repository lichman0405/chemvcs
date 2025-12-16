package server

import (
	"archive/zip"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"time"

	"github.com/lishi/chemvcs/internal/hpc"
	"github.com/lishi/chemvcs/internal/repo"
)

type hpcJobRecord struct {
	RunHash    string        `json:"run_hash"`
	JobID      string        `json:"job_id"`
	JobSystem  string        `json:"job_system"`
	Queue      string        `json:"queue_name"`
	Status     hpc.JobStatus `json:"status"`
	Submitted  string        `json:"submitted_at"`
	Updated    string        `json:"updated_at"`
	WorkingDir string        `json:"working_dir"`
}

var (
	sbatchJobIDRe = regexp.MustCompile(`^([0-9]+)(?:;.*)?$`)
)

const (
	maxHPCRequestBodyBytes = 1 << 20 // 1 MiB
	maxHPCScriptBytes      = 512 << 10
	maxHPCPatterns         = 64
	maxHPCPatternLen       = 256
	maxRetrieveFiles       = 2000
	maxRetrieveBytes       = int64(512 << 20) // 512 MiB

	sbatchTimeout  = 30 * time.Second
	squeueTimeout  = 10 * time.Second
	sacctTimeout   = 10 * time.Second
	scancelTimeout = 10 * time.Second
)

func runCommand(timeout time.Duration, name string, args ...string) ([]byte, error) {
	ctx, cancel := context.WithTimeout(context.Background(), timeout)
	defer cancel()
	cmd := exec.CommandContext(ctx, name, args...)
	out, err := cmd.CombinedOutput()
	if ctx.Err() == context.DeadlineExceeded {
		return out, fmt.Errorf("command timed out: %s", name)
	}
	return out, err
}

func (s *Server) handleHPC(w http.ResponseWriter, req *http.Request, repoID, resource string) {
	repository, err := s.getRepository(repoID)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	sub := strings.TrimPrefix(resource, "hpc")
	sub = strings.TrimPrefix(sub, "/")

	switch {
	case sub == "submit":
		s.handleHPCSubmit(w, req, repository)
	case sub == "jobs":
		s.handleHPCJobs(w, req, repository)
	case strings.HasPrefix(sub, "jobs/"):
		identifier := strings.TrimPrefix(sub, "jobs/")
		s.handleHPCJobGet(w, req, repository, identifier)
	case sub == "cancel":
		s.handleHPCCancel(w, req, repository)
	case sub == "retrieve":
		s.handleHPCRetrieve(w, req, repository)
	default:
		s.sendError(w, http.StatusNotFound, "not_found", "Resource not found")
	}
}

func (s *Server) handleHPCSubmit(w http.ResponseWriter, r *http.Request, repo *repo.Repository) {
	if r.Method != http.MethodPost {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	r.Body = http.MaxBytesReader(w, r.Body, maxHPCRequestBodyBytes)

	var req struct {
		RunHash    string `json:"run_hash"`
		Script     string `json:"script"`
		CaptureEnv bool   `json:"capture_env"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid JSON")
		return
	}
	if strings.TrimSpace(req.RunHash) == "" {
		s.sendError(w, http.StatusBadRequest, "bad_request", "run_hash is required")
		return
	}
	if req.Script == "" {
		s.sendError(w, http.StatusBadRequest, "bad_request", "script is required")
		return
	}
	if int64(len(req.Script)) > maxHPCScriptBytes {
		s.sendError(w, http.StatusBadRequest, "bad_request", fmt.Sprintf("script is too large (max %d bytes)", maxHPCScriptBytes))
		return
	}

	repoPath := repo.Path()
	tmp, err := os.CreateTemp(repoPath, "chemvcs-sbatch-*.sh")
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to create temp script: %v", err))
		return
	}
	defer os.Remove(tmp.Name())

	if err := tmp.Chmod(0755); err != nil {
		_ = tmp.Close()
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to chmod temp script: %v", err))
		return
	}
	if _, err := tmp.WriteString(req.Script); err != nil {
		_ = tmp.Close()
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to write temp script: %v", err))
		return
	}
	_ = tmp.Close()

	// Prefer --parsable for stable job id parsing.
	args := []string{"--parsable", "--chdir", repoPath}
	if req.CaptureEnv {
		args = append(args, "--export=ALL")
	}
	args = append(args, tmp.Name())

	out, err := runCommand(sbatchTimeout, "sbatch", args...)
	if err != nil {
		s.sendError(w, http.StatusBadRequest, "submit_failed", strings.TrimSpace(string(out)))
		return
	}

	jobID := strings.TrimSpace(string(out))
	m := sbatchJobIDRe.FindStringSubmatch(jobID)
	if m == nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("unexpected sbatch output: %q", jobID))
		return
	}
	jobID = m[1]

	now := time.Now().UTC().Format(time.RFC3339)
	rec := hpcJobRecord{
		RunHash:    req.RunHash,
		JobID:      jobID,
		JobSystem:  "slurm",
		Status:     hpc.JobStatusPending,
		Submitted:  now,
		Updated:    now,
		WorkingDir: repoPath,
	}

	data, err := json.Marshal(rec)
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to encode job record: %v", err))
		return
	}
	blobHash, err := repo.Store().PutBlob(data)
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to store job record: %v", err))
		return
	}

	if err := repo.Refs().UpdateRef(hpcRunRef(req.RunHash), blobHash); err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to update run ref: %v", err))
		return
	}
	if err := repo.Refs().UpdateRef(hpcJobRef(jobID), blobHash); err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", fmt.Sprintf("failed to update job ref: %v", err))
		return
	}

	s.sendJSON(w, http.StatusOK, map[string]string{
		"job_id":     jobID,
		"run_hash":   req.RunHash,
		"job_system": "slurm",
	})
}

func (s *Server) handleHPCJobs(w http.ResponseWriter, r *http.Request, repo *repo.Repository) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	statusFilter := strings.TrimSpace(r.URL.Query().Get("status"))
	refresh := r.URL.Query().Get("refresh") == "1" || strings.EqualFold(r.URL.Query().Get("refresh"), "true")

	recs, err := s.loadAllHPCJobRecords(repo)
	if err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", err.Error())
		return
	}

	if refresh {
		for i := range recs {
			updated, _ := s.refreshHPCJobStatus(repo, &recs[i])
			recs[i] = updated
		}
	}

	var jobs []hpc.JobInfo
	for _, rec := range recs {
		if statusFilter != "" && !strings.EqualFold(string(rec.Status), statusFilter) {
			continue
		}
		jobs = append(jobs, hpc.JobInfo{
			JobID:      rec.JobID,
			RunHash:    rec.RunHash,
			Status:     rec.Status,
			JobSystem:  rec.JobSystem,
			Queue:      rec.Queue,
			Submitted:  rec.Submitted,
			Updated:    rec.Updated,
			WorkingDir: rec.WorkingDir,
		})
	}

	// Deterministic order for clients.
	sort.Slice(jobs, func(i, j int) bool {
		if jobs[i].Submitted == jobs[j].Submitted {
			return jobs[i].JobID < jobs[j].JobID
		}
		return jobs[i].Submitted < jobs[j].Submitted
	})

	s.sendJSON(w, http.StatusOK, map[string]interface{}{
		"jobs": jobs,
	})
}

func (s *Server) handleHPCJobGet(w http.ResponseWriter, r *http.Request, repo *repo.Repository, identifier string) {
	if r.Method != http.MethodGet {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	refresh := r.URL.Query().Get("refresh") == "1" || strings.EqualFold(r.URL.Query().Get("refresh"), "true")

	rec, err := s.loadHPCJobRecordByIdentifier(repo, identifier)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}
	if refresh {
		updated, _ := s.refreshHPCJobStatus(repo, &rec)
		rec = updated
	}

	s.sendJSON(w, http.StatusOK, rec)
}

func (s *Server) handleHPCCancel(w http.ResponseWriter, r *http.Request, repo *repo.Repository) {
	if r.Method != http.MethodPost {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	r.Body = http.MaxBytesReader(w, r.Body, maxHPCRequestBodyBytes)

	var req struct {
		Identifier string `json:"identifier"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid JSON")
		return
	}
	if strings.TrimSpace(req.Identifier) == "" {
		s.sendError(w, http.StatusBadRequest, "bad_request", "identifier is required")
		return
	}

	jobID, rec, err := s.resolveHPCJobID(repo, req.Identifier)
	if err != nil {
		s.sendError(w, http.StatusNotFound, "not_found", err.Error())
		return
	}

	out, err := runCommand(scancelTimeout, "scancel", jobID)
	if err != nil {
		s.sendError(w, http.StatusBadRequest, "cancel_failed", strings.TrimSpace(string(out)))
		return
	}

	now := time.Now().UTC().Format(time.RFC3339)
	rec.Status = hpc.JobStatusCancelled
	rec.Updated = now

	if err := s.persistHPCJobRecord(repo, rec); err != nil {
		s.sendError(w, http.StatusInternalServerError, "internal_error", err.Error())
		return
	}

	s.sendJSON(w, http.StatusOK, map[string]string{
		"job_id": jobID,
		"status": string(hpc.JobStatusCancelled),
	})
}

func (s *Server) handleHPCRetrieve(w http.ResponseWriter, r *http.Request, repo *repo.Repository) {
	if r.Method != http.MethodPost {
		s.sendError(w, http.StatusMethodNotAllowed, "method_not_allowed", "Method not allowed")
		return
	}

	r.Body = http.MaxBytesReader(w, r.Body, maxHPCRequestBodyBytes)

	var req struct {
		Patterns []string `json:"patterns"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", "Invalid JSON")
		return
	}

	patterns := req.Patterns
	if len(patterns) == 0 {
		patterns = []string{"*"}
	}
	if len(patterns) > maxHPCPatterns {
		s.sendError(w, http.StatusBadRequest, "bad_request", fmt.Sprintf("too many patterns (max %d)", maxHPCPatterns))
		return
	}
	for i := range patterns {
		patterns[i] = strings.TrimSpace(patterns[i])
		if patterns[i] == "" {
			continue
		}
		if len(patterns[i]) > maxHPCPatternLen {
			s.sendError(w, http.StatusBadRequest, "bad_request", fmt.Sprintf("pattern too long (max %d)", maxHPCPatternLen))
			return
		}
	}

	repoRoot := repo.Path()

	files, err := collectFiles(repoRoot, patterns)
	if err != nil {
		s.sendError(w, http.StatusBadRequest, "bad_request", err.Error())
		return
	}
	if len(files) > maxRetrieveFiles {
		s.sendError(w, http.StatusBadRequest, "bad_request", fmt.Sprintf("too many files matched (max %d)", maxRetrieveFiles))
		return
	}

	// Preflight size limits before streaming a zip.
	var total int64
	for _, absPath := range files {
		info, err := os.Stat(absPath)
		if err != nil || info.IsDir() {
			continue
		}
		total += info.Size()
		if total > maxRetrieveBytes {
			s.sendError(w, http.StatusBadRequest, "bad_request", fmt.Sprintf("matched files too large (max %d bytes)", maxRetrieveBytes))
			return
		}
	}

	w.Header().Set("Content-Type", "application/zip")
	w.Header().Set("Content-Disposition", "attachment; filename=chemvcs-results.zip")
	w.WriteHeader(http.StatusOK)

	zw := zip.NewWriter(w)
	for _, absPath := range files {
		rel, err := filepath.Rel(repoRoot, absPath)
		if err != nil {
			continue
		}
		// Safety: never include paths outside repoRoot.
		if rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
			continue
		}
		rel = filepath.ToSlash(rel)
		if rel == "." || rel == "" {
			continue
		}
		if rel == ".." || strings.HasPrefix(rel, "../") {
			continue
		}
		if strings.HasPrefix(rel, ".chemvcs/") || rel == ".chemvcs" {
			continue
		}

		info, err := os.Stat(absPath)
		if err != nil || info.IsDir() {
			continue
		}

		h, err := zip.FileInfoHeader(info)
		if err != nil {
			continue
		}
		h.Name = rel
		h.Method = zip.Deflate

		fw, err := zw.CreateHeader(h)
		if err != nil {
			continue
		}

		f, err := os.Open(absPath)
		if err != nil {
			continue
		}
		_, _ = io.Copy(fw, f)
		_ = f.Close()
	}
	_ = zw.Close()
}

func collectFiles(repoRoot string, patterns []string) ([]string, error) {
	seen := map[string]struct{}{}
	var all []string

	for _, pat := range patterns {
		pat = strings.TrimSpace(pat)
		if pat == "" {
			continue
		}
		if filepath.IsAbs(pat) {
			return nil, fmt.Errorf("invalid pattern %q: absolute paths are not allowed", pat)
		}
		// Prevent escaping the repository via ../ or Windows drive-letter tricks.
		if pat == ".." || strings.HasPrefix(pat, "../") || strings.HasPrefix(pat, "..\\") {
			return nil, fmt.Errorf("invalid pattern %q: path traversal is not allowed", pat)
		}
		if strings.Contains(pat, ":") {
			return nil, fmt.Errorf("invalid pattern %q: drive letters are not allowed", pat)
		}
		// Restrict to repository root.
		glob := filepath.Join(repoRoot, pat)
		matches, err := filepath.Glob(glob)
		if err != nil {
			return nil, fmt.Errorf("invalid pattern %q", pat)
		}
		for _, m := range matches {
			abs := m
			rel, err := filepath.Rel(repoRoot, abs)
			if err != nil {
				continue
			}
			if rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
				continue
			}
			if _, ok := seen[abs]; ok {
				continue
			}
			seen[abs] = struct{}{}
			all = append(all, abs)
		}
	}

	sort.Strings(all)
	return all, nil
}

func hpcRunRef(runHash string) string {
	return "refs/hpc/runs/" + runHash
}

func hpcJobRef(jobID string) string {
	return "refs/hpc/jobs/" + jobID
}

func (s *Server) loadHPCJobRecordByIdentifier(repo *repo.Repository, identifier string) (hpcJobRecord, error) {
	ref := hpcRunRef(identifier)
	if isDigits(identifier) {
		ref = hpcJobRef(identifier)
	}

	hash, err := repo.Refs().ResolveRef(ref)
	if err != nil || hash == "" {
		return hpcJobRecord{}, fmt.Errorf("job not found: %s", identifier)
	}
	data, err := repo.Store().GetBlob(hash)
	if err != nil {
		return hpcJobRecord{}, fmt.Errorf("job record not found: %s", identifier)
	}

	var rec hpcJobRecord
	if err := json.Unmarshal(data, &rec); err != nil {
		return hpcJobRecord{}, fmt.Errorf("invalid job record")
	}
	return rec, nil
}

func (s *Server) loadAllHPCJobRecords(repo *repo.Repository) ([]hpcJobRecord, error) {
	runsDir := filepath.Join(repo.Path(), ".chemvcs", "refs", "hpc", "runs")
	entries, err := os.ReadDir(runsDir)
	if err != nil {
		if os.IsNotExist(err) {
			return []hpcJobRecord{}, nil
		}
		return nil, fmt.Errorf("failed to read hpc refs: %v", err)
	}

	var recs []hpcJobRecord
	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		runHash := e.Name()
		rec, err := s.loadHPCJobRecordByIdentifier(repo, runHash)
		if err != nil {
			continue
		}
		recs = append(recs, rec)
	}
	return recs, nil
}

func (s *Server) resolveHPCJobID(repo *repo.Repository, identifier string) (string, hpcJobRecord, error) {
	rec, err := s.loadHPCJobRecordByIdentifier(repo, identifier)
	if err == nil {
		return rec.JobID, rec, nil
	}

	// If identifier is job id but no record exists, still allow cancellation.
	if isDigits(identifier) {
		return identifier, hpcJobRecord{JobID: identifier, JobSystem: "slurm"}, nil
	}

	return "", hpcJobRecord{}, err
}

func (s *Server) refreshHPCJobStatus(repo *repo.Repository, rec *hpcJobRecord) (hpcJobRecord, error) {
	status, queue := slurmStatus(rec.JobID)
	if status == "" {
		return *rec, nil
	}

	now := time.Now().UTC().Format(time.RFC3339)
	updated := *rec
	updated.Status = status
	if queue != "" {
		updated.Queue = queue
	}
	updated.Updated = now

	if err := s.persistHPCJobRecord(repo, updated); err != nil {
		return *rec, err
	}
	return updated, nil
}

func (s *Server) persistHPCJobRecord(repo *repo.Repository, rec hpcJobRecord) error {
	data, err := json.Marshal(rec)
	if err != nil {
		return fmt.Errorf("failed to encode job record: %v", err)
	}
	blobHash, err := repo.Store().PutBlob(data)
	if err != nil {
		return fmt.Errorf("failed to store job record: %v", err)
	}
	if rec.RunHash != "" {
		if err := repo.Refs().UpdateRef(hpcRunRef(rec.RunHash), blobHash); err != nil {
			return fmt.Errorf("failed to update run ref: %v", err)
		}
	}
	if rec.JobID != "" {
		if err := repo.Refs().UpdateRef(hpcJobRef(rec.JobID), blobHash); err != nil {
			return fmt.Errorf("failed to update job ref: %v", err)
		}
	}
	return nil
}

func isDigits(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return false
		}
	}
	return true
}

func slurmStatus(jobID string) (hpc.JobStatus, string) {
	// Try squeue first.
	{
		out, err := runCommand(squeueTimeout, "squeue", "-j", jobID, "-h", "-o", "%T|%P")
		if err == nil {
			line := strings.TrimSpace(string(out))
			if line != "" {
				parts := strings.SplitN(line, "|", 2)
				state := strings.TrimSpace(parts[0])
				queue := ""
				if len(parts) == 2 {
					queue = strings.TrimSpace(parts[1])
				}
				return mapSlurmState(state), queue
			}
		}
	}

	// Fallback to sacct.
	{
		out, err := runCommand(sacctTimeout, "sacct", "-j", jobID, "-n", "-o", "State", "-P")
		if err == nil {
			line := strings.TrimSpace(string(out))
			if line != "" {
				// sacct may output multiple lines; take first non-empty.
				lines := strings.Split(line, "\n")
				for _, l := range lines {
					l = strings.TrimSpace(l)
					if l == "" {
						continue
					}
					// -P uses '|' separator; state is first field.
					state := strings.SplitN(l, "|", 2)[0]
					return mapSlurmState(state), ""
				}
			}
		}
	}

	return "", ""
}

func mapSlurmState(state string) hpc.JobStatus {
	s := strings.ToUpper(strings.TrimSpace(state))
	s = strings.SplitN(s, "+", 2)[0]
	s = strings.SplitN(s, ":", 2)[0]

	switch {
	case s == "PENDING" || s == "CONFIGURING" || s == "SUSPENDED":
		return hpc.JobStatusPending
	case s == "RUNNING" || s == "COMPLETING":
		return hpc.JobStatusRunning
	case s == "COMPLETED":
		return hpc.JobStatusCompleted
	case strings.Contains(s, "CANCEL"):
		return hpc.JobStatusCancelled
	case s == "FAILED" || s == "TIMEOUT" || s == "NODE_FAIL" || s == "OUT_OF_MEMORY":
		return hpc.JobStatusFailed
	default:
		return hpc.JobStatusUnknown
	}
}
