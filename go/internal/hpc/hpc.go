package hpc

import (
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
)

// JobStatus represents the status of an HPC job
type JobStatus string

const (
	JobStatusPending   JobStatus = "PENDING"
	JobStatusRunning   JobStatus = "RUNNING"
	JobStatusCompleted JobStatus = "COMPLETED"
	JobStatusFailed    JobStatus = "FAILED"
	JobStatusCancelled JobStatus = "CANCELLED"
	JobStatusUnknown   JobStatus = "UNKNOWN"
)

// JobInfo contains information about an HPC job
type JobInfo struct {
	JobID     string    `json:"job_id"`
	RunHash   string    `json:"run_hash"`
	Status    JobStatus `json:"status"`
	JobSystem string    `json:"job_system"`
	Queue     string    `json:"queue_name"`
	Submitted string    `json:"submitted_at"`
}

// SubmitOptions contains options for job submission
type SubmitOptions struct {
	RunHash    string
	ScriptPath string
	CaptureEnv bool
	WorkingDir string
}

// SubmitResult contains the result of job submission
type SubmitResult struct {
	JobID     string `json:"job_id"`
	RunHash   string `json:"run_hash"`
	JobSystem string `json:"job_system"`
}

// RetrieveOptions contains options for result retrieval
type RetrieveOptions struct {
	RunHash     string
	Patterns    []string
	Destination string
}

// findPythonExecutable finds the Python interpreter
func findPythonExecutable() (string, error) {
	// Try python3 first, then python
	for _, name := range []string{"python3", "python"} {
		path, err := exec.LookPath(name)
		if err == nil {
			return path, nil
		}
	}
	return "", fmt.Errorf("python executable not found")
}

// findChemVCSPyModule finds the chemvcs_py module path
func findChemVCSPyModule() (string, error) {
	// Try to find the module in the repository
	// Assume this Go code is in go/internal/hpc, so chemvcs_py is in ../../python
	execPath, err := os.Executable()
	if err != nil {
		return "", fmt.Errorf("failed to get executable path: %w", err)
	}

	// Navigate from executable to repository root
	repoRoot := filepath.Dir(filepath.Dir(filepath.Dir(execPath)))
	pythonPath := filepath.Join(repoRoot, "python")

	// Check if the path exists
	if _, err := os.Stat(pythonPath); err == nil {
		return pythonPath, nil
	}

	// If not found relative to executable, try relative to current directory
	cwd, err := os.Getwd()
	if err != nil {
		return "", fmt.Errorf("failed to get current directory: %w", err)
	}

	// Try to find python directory
	for dir := cwd; dir != filepath.Dir(dir); dir = filepath.Dir(dir) {
		pythonPath := filepath.Join(dir, "python")
		if _, err := os.Stat(pythonPath); err == nil {
			return pythonPath, nil
		}
	}

	return "", fmt.Errorf("chemvcs_py module not found")
}

// runPythonScript executes a Python script with chemvcs_py in PYTHONPATH
func runPythonScript(script string, args ...string) ([]byte, error) {
	pythonExec, err := findPythonExecutable()
	if err != nil {
		return nil, err
	}

	pythonPath, err := findChemVCSPyModule()
	if err != nil {
		return nil, err
	}

	// Build command
	cmdArgs := append([]string{"-c", script}, args...)
	cmd := exec.Command(pythonExec, cmdArgs...)

	// Set PYTHONPATH to include chemvcs_py
	env := os.Environ()
	pythonPathEnv := fmt.Sprintf("PYTHONPATH=%s", pythonPath)
	cmd.Env = append(env, pythonPathEnv)

	// Run command and capture output
	output, err := cmd.CombinedOutput()
	if err != nil {
		return output, fmt.Errorf("python script failed: %w\nOutput: %s", err, string(output))
	}

	return output, nil
}

// SubmitJob submits an HPC job for a run
func SubmitJob(repoPath string, opts SubmitOptions) (*SubmitResult, error) {
	script := `
import sys
import json
from chemvcs_py.api import Repo
from chemvcs_py.hpc import SlurmAdapter, JobSubmitter

repo_path = sys.argv[1]
run_hash = sys.argv[2]
script_path = sys.argv[3]
capture_env = sys.argv[4] == "true"

try:
    repo = Repo(repo_path)
    adapter = SlurmAdapter()
    submitter = JobSubmitter(repo, adapter)
    
    job_id = submitter.submit_run(
        run_hash=run_hash,
        script_path=script_path,
        capture_env=capture_env
    )
    
    # Get updated run info
    run = repo.get_run(run_hash)
    result = {
        "job_id": job_id,
        "run_hash": run_hash,
        "job_system": run.job_system
    }
    
    print(json.dumps(result))
except Exception as e:
    print(json.dumps({"error": str(e)}), file=sys.stderr)
    sys.exit(1)
`

	captureEnvStr := "false"
	if opts.CaptureEnv {
		captureEnvStr = "true"
	}

	output, err := runPythonScript(script, repoPath, opts.RunHash, opts.ScriptPath, captureEnvStr)
	if err != nil {
		return nil, err
	}

	var result SubmitResult
	if err := json.Unmarshal(output, &result); err != nil {
		return nil, fmt.Errorf("failed to parse submission result: %w", err)
	}

	return &result, nil
}

// ListJobs lists all tracked HPC jobs
func ListJobs(repoPath string, statusFilter string) ([]JobInfo, error) {
	script := `
import sys
import json
from chemvcs_py.api import Repo

repo_path = sys.argv[1]
status_filter = sys.argv[2] if len(sys.argv) > 2 else ""

try:
    repo = Repo(repo_path)
    
    # Get all runs with job_id
    runs = []
    for obj_hash in repo.list_objects("run"):
        run = repo.get_run(obj_hash)
        if run.job_id:
            # Filter by status if specified
            if not status_filter or run.status.value.upper() == status_filter.upper():
                runs.append({
                    "job_id": run.job_id,
                    "run_hash": obj_hash,
                    "status": run.status.value.upper(),
                    "job_system": run.job_system or "",
                    "queue_name": run.queue_name or "",
                    "submitted_at": run.submitted_at.isoformat() if run.submitted_at else ""
                })
    
    print(json.dumps(runs))
except Exception as e:
    print(json.dumps({"error": str(e)}), file=sys.stderr)
    sys.exit(1)
`

	args := []string{repoPath}
	if statusFilter != "" {
		args = append(args, statusFilter)
	}

	output, err := runPythonScript(script, args...)
	if err != nil {
		return nil, err
	}

	var jobs []JobInfo
	if err := json.Unmarshal(output, &jobs); err != nil {
		return nil, fmt.Errorf("failed to parse jobs list: %w", err)
	}

	return jobs, nil
}

// CheckJobStatus checks the status of a specific job
func CheckJobStatus(repoPath string, runHash string) (*JobInfo, error) {
	script := `
import sys
import json
from chemvcs_py.api import Repo
from chemvcs_py.hpc import SlurmAdapter, JobTracker

repo_path = sys.argv[1]
run_hash = sys.argv[2]

try:
    repo = Repo(repo_path)
    adapter = SlurmAdapter()
    tracker = JobTracker(repo, adapter)
    
    # Update status
    status = tracker.check_status(run_hash)
    
    # Get updated run info
    run = repo.get_run(run_hash)
    result = {
        "job_id": run.job_id or "",
        "run_hash": run_hash,
        "status": run.status.value.upper(),
        "job_system": run.job_system or "",
        "queue_name": run.queue_name or "",
        "submitted_at": run.submitted_at.isoformat() if run.submitted_at else ""
    }
    
    print(json.dumps(result))
except Exception as e:
    print(json.dumps({"error": str(e)}), file=sys.stderr)
    sys.exit(1)
`

	output, err := runPythonScript(script, repoPath, runHash)
	if err != nil {
		return nil, err
	}

	var info JobInfo
	if err := json.Unmarshal(output, &info); err != nil {
		return nil, fmt.Errorf("failed to parse job info: %w", err)
	}

	return &info, nil
}

// RetrieveResults retrieves results for a completed job
func RetrieveResults(repoPath string, opts RetrieveOptions) ([]string, error) {
	script := `
import sys
import json
from pathlib import Path
from chemvcs_py.api import Repo
from chemvcs_py.hpc import SlurmAdapter, JobRetriever

repo_path = sys.argv[1]
run_hash = sys.argv[2]
patterns = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] else []
destination = sys.argv[4] if len(sys.argv) > 4 else "."

try:
    repo = Repo(repo_path)
    adapter = SlurmAdapter()
    retriever = JobRetriever(repo, adapter)
    
    # Retrieve results
    files = retriever.retrieve_results(
        run_hash=run_hash,
        output_patterns=patterns if patterns else None,
        destination=Path(destination)
    )
    
    result = {
        "files": [str(f) for f in files],
        "count": len(files)
    }
    
    print(json.dumps(result))
except Exception as e:
    print(json.dumps({"error": str(e)}), file=sys.stderr)
    sys.exit(1)
`

	patternsStr := ""
	if len(opts.Patterns) > 0 {
		patternsStr = filepath.Join(opts.Patterns...)
	}

	destination := opts.Destination
	if destination == "" {
		destination = "."
	}

	output, err := runPythonScript(script, repoPath, opts.RunHash, patternsStr, destination)
	if err != nil {
		return nil, err
	}

	var result struct {
		Files []string `json:"files"`
		Count int      `json:"count"`
	}
	if err := json.Unmarshal(output, &result); err != nil {
		return nil, fmt.Errorf("failed to parse retrieval result: %w", err)
	}

	return result.Files, nil
}
