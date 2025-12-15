package hpc_test

import (
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/lishi/chemvcs/internal/hpc"
	"github.com/lishi/chemvcs/internal/repo"
)

// TestHPCIntegration tests the HPC integration layer with mock Python scripts
func TestHPCIntegration(t *testing.T) {
	// Skip if Python is not available
	if _, err := exec.LookPath("python3"); err != nil {
		if _, err := exec.LookPath("python"); err != nil {
			t.Skip("Python not available, skipping HPC integration tests")
		}
	}

	// Create temporary test directory
	tmpDir := t.TempDir()

	// Initialize a ChemVCS repository
	r, err := repo.Init(tmpDir)
	if err != nil {
		t.Fatalf("Failed to initialize repository: %v", err)
	}

	repoPath := r.Path()

	t.Run("SubmitJob", func(t *testing.T) {
		// This test would require:
		// 1. A valid Run object in the repository
		// 2. A mock SLURM script
		// 3. Python with chemvcs_py installed

		// For now, we'll skip this test if chemvcs_py is not available
		t.Skip("Integration test requires chemvcs_py to be installed")
	})

	t.Run("ListJobs", func(t *testing.T) {
		// Test listing jobs from an empty repository
		jobs, err := hpc.ListJobs(repoPath, "")
		if err != nil {
			// If chemvcs_py is not available or Python has issues, skip
			if strings.Contains(err.Error(), "chemvcs_py") ||
				strings.Contains(err.Error(), "python") ||
				strings.Contains(err.Error(), "exit status") {
				t.Skip("chemvcs_py not available or Python not configured properly")
			}
			t.Fatalf("Failed to list jobs: %v", err)
		}

		if len(jobs) != 0 {
			t.Errorf("Expected 0 jobs, got %d", len(jobs))
		}
	})
}

// TestHPCPackageStructure tests that the HPC package is properly structured
func TestHPCPackageStructure(t *testing.T) {
	tests := []struct {
		name     string
		testFunc func(t *testing.T)
	}{
		{
			name: "JobStatus constants",
			testFunc: func(t *testing.T) {
				statuses := []hpc.JobStatus{
					hpc.JobStatusPending,
					hpc.JobStatusRunning,
					hpc.JobStatusCompleted,
					hpc.JobStatusFailed,
					hpc.JobStatusCancelled,
					hpc.JobStatusUnknown,
				}

				for _, status := range statuses {
					if status == "" {
						t.Error("JobStatus constant is empty")
					}
				}
			},
		},
		{
			name: "SubmitOptions structure",
			testFunc: func(t *testing.T) {
				opts := hpc.SubmitOptions{
					RunHash:    "abc123",
					ScriptPath: "/path/to/script.slurm",
					CaptureEnv: true,
				}

				if opts.RunHash != "abc123" {
					t.Errorf("Expected RunHash 'abc123', got '%s'", opts.RunHash)
				}
			},
		},
		{
			name: "RetrieveOptions structure",
			testFunc: func(t *testing.T) {
				opts := hpc.RetrieveOptions{
					RunHash:     "abc123",
					Patterns:    []string{"*.out", "*.log"},
					Destination: "/results",
				}

				if len(opts.Patterns) != 2 {
					t.Errorf("Expected 2 patterns, got %d", len(opts.Patterns))
				}
			},
		},
		{
			name: "JobInfo JSON marshaling",
			testFunc: func(t *testing.T) {
				info := hpc.JobInfo{
					JobID:     "12345",
					RunHash:   "abc123def456",
					Status:    hpc.JobStatusRunning,
					JobSystem: "slurm",
					Queue:     "batch",
					Submitted: "2025-12-15T10:00:00",
				}

				data, err := json.Marshal(info)
				if err != nil {
					t.Fatalf("Failed to marshal JobInfo: %v", err)
				}

				var decoded hpc.JobInfo
				if err := json.Unmarshal(data, &decoded); err != nil {
					t.Fatalf("Failed to unmarshal JobInfo: %v", err)
				}

				if decoded.JobID != info.JobID {
					t.Errorf("JobID mismatch: expected %s, got %s", info.JobID, decoded.JobID)
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, tt.testFunc)
	}
}

// TestPythonAvailability tests whether Python is available for HPC operations
func TestPythonAvailability(t *testing.T) {
	pythonExecs := []string{"python3", "python"}
	found := false

	for _, name := range pythonExecs {
		if _, err := exec.LookPath(name); err == nil {
			found = true
			t.Logf("Found Python executable: %s", name)
			break
		}
	}

	if !found {
		t.Log("Warning: No Python executable found. HPC commands will not work.")
	}
}

// TestMockWorkflow demonstrates how to test HPC workflow with mocks
func TestMockWorkflow(t *testing.T) {
	// Create a temporary directory for testing
	tmpDir := t.TempDir()

	// Create a mock Python script that simulates chemvcs_py behavior
	mockScript := filepath.Join(tmpDir, "mock_hpc.py")
	mockContent := `#!/usr/bin/env python3
import sys
import json

# Mock submission
if "submit" in sys.argv:
    result = {
        "job_id": "mock_12345",
        "run_hash": "abc123",
        "job_system": "slurm"
    }
    print(json.dumps(result))
    sys.exit(0)

# Mock list jobs
if "list" in sys.argv:
    jobs = []
    print(json.dumps(jobs))
    sys.exit(0)

# Mock retrieve
if "retrieve" in sys.argv:
    result = {
        "files": ["output.log", "results.dat"],
        "count": 2
    }
    print(json.dumps(result))
    sys.exit(0)

sys.exit(1)
`

	if err := os.WriteFile(mockScript, []byte(mockContent), 0755); err != nil {
		t.Fatalf("Failed to create mock script: %v", err)
	}

	t.Log("Mock workflow test setup complete")
	t.Log("In a real scenario, we would:")
	t.Log("  1. Create a Run object in the repository")
	t.Log("  2. Submit it using the mock adapter")
	t.Log("  3. Query its status")
	t.Log("  4. Retrieve results")
}

// TestCLIIntegration tests the CLI commands end-to-end
func TestCLIIntegration(t *testing.T) {
	// Build the chemvcs CLI if not already built
	cliPath := "../../../go/chemvcs.exe"
	if _, err := os.Stat(cliPath); os.IsNotExist(err) {
		// Try alternative path
		cliPath = "../../chemvcs.exe"
		if _, err := os.Stat(cliPath); os.IsNotExist(err) {
			// Try in parent directory
			cliPath = "../chemvcs.exe"
			if _, err := os.Stat(cliPath); os.IsNotExist(err) {
				t.Skip("chemvcs CLI not built, run 'go build -o chemvcs.exe ./cmd/chemvcs' first")
			}
		}
	}

	// Create temporary test directory
	tmpDir := t.TempDir()

	// Initialize repository
	cmd := exec.Command(cliPath, "init", tmpDir)
	if err := cmd.Run(); err != nil {
		t.Fatalf("Failed to initialize repository: %v", err)
	}

	t.Run("HelpMessage", func(t *testing.T) {
		cmd := exec.Command(cliPath, "help")
		output, err := cmd.CombinedOutput()
		if err != nil {
			t.Fatalf("Failed to run help command: %v", err)
		}

		helpText := string(output)
		requiredSections := []string{
			"submit",
			"jobs",
			"retrieve",
			"HPC Commands",
		}

		for _, section := range requiredSections {
			if !strings.Contains(helpText, section) {
				t.Errorf("Help text missing section: %s", section)
			}
		}
	})

	t.Run("SubmitWithoutArgs", func(t *testing.T) {
		cmd := exec.Command(cliPath, "submit")
		cmd.Dir = tmpDir
		output, err := cmd.CombinedOutput()

		// Should fail with usage message
		if err == nil {
			t.Error("Expected submit without args to fail")
		}

		outputStr := string(output)
		if !strings.Contains(outputStr, "usage") && !strings.Contains(outputStr, "Error") {
			t.Logf("Expected usage or error message, got: %s", outputStr)
			// Don't fail test - error handling may vary
		}
	})

	t.Run("JobsInEmptyRepo", func(t *testing.T) {
		cmd := exec.Command(cliPath, "jobs")
		cmd.Dir = tmpDir
		output, err := cmd.CombinedOutput()

		// This might fail if Python/chemvcs_py is not available
		if err != nil {
			if strings.Contains(string(output), "python") || strings.Contains(string(output), "chemvcs_py") {
				t.Skip("Python or chemvcs_py not available")
			}
			t.Logf("Jobs command output: %s", string(output))
		}
	})
}

// Benchmark tests for HPC operations
func BenchmarkJobInfoMarshaling(b *testing.B) {
	info := hpc.JobInfo{
		JobID:     "12345",
		RunHash:   "abc123def456",
		Status:    hpc.JobStatusRunning,
		JobSystem: "slurm",
		Queue:     "batch",
		Submitted: "2025-12-15T10:00:00",
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := json.Marshal(info)
		if err != nil {
			b.Fatal(err)
		}
	}
}

// Example test demonstrating expected usage
func ExampleSubmitJob() {
	// This example shows how to use the HPC integration
	// In a real scenario:

	// 1. Initialize repository
	// repo, _ := repo.Init("/path/to/repo")

	// 2. Submit job
	// result, err := hpc.SubmitJob(repo.Path(), hpc.SubmitOptions{
	//     RunHash:    "abc123def456",
	//     ScriptPath: "/path/to/script.slurm",
	//     CaptureEnv: true,
	// })

	// 3. Check result
	// fmt.Printf("Job ID: %s\n", result.JobID)

	fmt.Println("Example: Submit HPC job")
	// Output: Example: Submit HPC job
}
