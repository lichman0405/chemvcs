package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/lishi/chemvcs/internal/hpc"
	"github.com/lishi/chemvcs/internal/objectstore"
	"github.com/lishi/chemvcs/internal/remote"
	"github.com/lishi/chemvcs/internal/repo"
	"github.com/lishi/chemvcs/internal/workspace"
)

const version = "0.1.0-dev"

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(2)
	}

	cmd := os.Args[1]
	args := os.Args[2:]

	var err error
	switch cmd {
	case "init":
		err = handleInit(args)
	case "commit":
		err = handleCommit(args)
	case "log":
		err = handleLog(args)
	case "branch":
		err = handleBranch(args)
	case "checkout":
		err = handleCheckout(args)
	case "status":
		err = handleStatus(args)
	case "merge":
		err = handleMerge(args)
	case "remote":
		err = handleRemote(args)
	case "push":
		err = handlePush(args)
	case "pull":
		err = handlePull(args)
	case "fetch":
		err = handleFetch(args)
	case "inspect-object":
		err = handleInspectObject(args)
	case "list-objects":
		err = handleListObjects(args)
	case "submit":
		err = handleSubmit(args)
	case "jobs":
		err = handleJobs(args)
	case "retrieve":
		err = handleRetrieve(args)
	case "cancel":
		err = handleCancel(args)
	case "watch":
		err = handleWatch(args)
	case "pack":
		err = handlePack(args)
	case "gc":
		err = handleGC(args)
	case "fsck":
		err = handleFsck(args)
	case "version":
		handleVersion()
	case "help", "--help", "-h":
		printUsage()
	default:
		fmt.Fprintf(os.Stderr, "Unknown command: %s\n\n", cmd)
		printUsage()
		os.Exit(2)
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func printUsage() {
	fmt.Println("ChemVCS - A domain-aware version control system for computational chemistry")
	fmt.Println()
	fmt.Println("Usage: chemvcs <command> [options]")
	fmt.Println()
	fmt.Println("Commands:")
	fmt.Println("  init [path]           Initialise a new repository")
	fmt.Println("  commit -m <message>   Create a new snapshot")
	fmt.Println("  log [-n <count>]      Show snapshot history")
	fmt.Println("  branch [name]         List or create branches")
	fmt.Println("  checkout <target>     Switch branches or snapshots")
	fmt.Println("  status                Show working directory changes")
	fmt.Println("  merge <branch>        Merge a branch into current branch")
	fmt.Println("  remote add <name> <url>  Add a remote repository")
	fmt.Println("  push <remote> <branch>   Push a branch to remote")
	fmt.Println("  pull <remote> <branch>   Pull a branch from remote")
	fmt.Println("  fetch <remote> <branch>  Fetch objects from remote")
	fmt.Println("  inspect-object <hash> [--format=json]  Inspect an object")
	fmt.Println("  list-objects [--type=<type>]           List all objects")
	fmt.Println()
	fmt.Println("HPC Commands:")
	fmt.Println("  submit <run-hash> <script>        Submit HPC job for a run")
	fmt.Println("  jobs [--status=<status>] [<run-hash|job-id>]  List tracked HPC jobs")
	fmt.Println("  retrieve <run-hash> [--patterns=<patterns>] [--dest=<path>] [--commit] [--commit-message=<msg>]")
	fmt.Println("                                    Retrieve job results")
	fmt.Println("  cancel <run-hash|job-id>          Cancel a running job")
	fmt.Println("  watch <run-hash|job-id> [--interval=<sec>] [--timeout=<sec>]")
	fmt.Println("                                    Monitor job until completion")
	fmt.Println()
	fmt.Println("Maintenance Commands:")
	fmt.Println("  pack [--all] [--keep-loose]       Pack loose objects into packfile")
	fmt.Println("  gc [--prune=<age>] [--dry-run]    Garbage collect unreachable objects")
	fmt.Println("  fsck [--full]                     Verify repository integrity")
	fmt.Println()
	fmt.Println("Other Commands:")
	fmt.Println("  version               Show version information")
	fmt.Println("  help                  Show this help message")
	fmt.Println()
	fmt.Println("Environment variables:")
	fmt.Println("  CHEMVCS_AUTHOR_NAME   Default author name")
	fmt.Println("  CHEMVCS_AUTHOR_EMAIL  Default author email")
}

func handleVersion() {
	fmt.Printf("ChemVCS version %s\n", version)
}

func handleInit(args []string) error {
	var path string
	if len(args) > 0 {
		path = args[0]
	} else {
		path = "."
	}

	absPath, err := filepath.Abs(path)
	if err != nil {
		return fmt.Errorf("failed to resolve path: %w", err)
	}

	// Check if directory exists
	if _, err := os.Stat(absPath); os.IsNotExist(err) {
		return fmt.Errorf("directory does not exist: %s", absPath)
	}

	_, err = repo.Init(absPath)
	if err != nil {
		return err
	}

	fmt.Printf("Initialised empty ChemVCS repository in %s\n",
		filepath.Join(absPath, ".chemvcs"))
	return nil
}

func handleCommit(args []string) error {
	fs := flag.NewFlagSet("commit", flag.ExitOnError)
	message := fs.String("m", "", "Commit message (required)")
	author := fs.String("author", "", "Author (overrides environment)")
	fs.Parse(args)

	if *message == "" {
		return fmt.Errorf("commit message required (use -m)")
	}

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Get author
	authorStr := getAuthor(*author)

	snapHash, err := commitWorkingDirectory(r, *message, authorStr)
	if err != nil {
		return err
	}

	printCommitResult(r, snapHash, *message)
	return nil
}

func commitWorkingDirectory(r *repo.Repository, message string, authorStr string) (string, error) {
	// Scan working directory to create root object
	scanner := workspace.NewScanner(r.Store())
	rootObj, err := scanner.ScanDirectory(".")
	if err != nil {
		return "", fmt.Errorf("failed to scan working directory: %w", err)
	}

	// Store root object
	rootHash, err := r.Store().PutObject(rootObj)
	if err != nil {
		return "", fmt.Errorf("failed to store root object: %w", err)
	}

	// Create snapshot
	snapHash, err := r.CreateSnapshot(rootHash, message, authorStr)
	if err != nil {
		return "", err
	}

	return snapHash, nil
}

func printCommitResult(r *repo.Repository, snapHash string, message string) {
	shortHash := snapHash
	if len(snapHash) > 8 {
		shortHash = snapHash[:8]
	}

	branch, _ := r.CurrentBranch()
	if branch != "" {
		fmt.Printf("[%s %s] %s\n", branch, shortHash, message)
		return
	}
	fmt.Printf("[detached %s] %s\n", shortHash, message)
}

func getAuthor(flagValue string) string {
	if flagValue != "" {
		return flagValue
	}

	// Try environment variables
	name := os.Getenv("CHEMVCS_AUTHOR_NAME")
	email := os.Getenv("CHEMVCS_AUTHOR_EMAIL")

	if name != "" && email != "" {
		return fmt.Sprintf("%s <%s>", name, email)
	}

	if name != "" {
		return name
	}

	// Default
	return "Unknown <unknown@localhost>"
}

func handleLog(args []string) error {
	fs := flag.NewFlagSet("log", flag.ExitOnError)
	limit := fs.Int("n", 20, "Number of snapshots to show")
	oneline := fs.Bool("oneline", false, "Show condensed one-line format")
	fs.Parse(args)

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	snapshots, err := r.Log(*limit)
	if err != nil {
		return err
	}

	if len(snapshots) == 0 {
		fmt.Println("No commits yet")
		return nil
	}

	// Display snapshots
	for _, snap := range snapshots {
		shortHash := snap.ShortHash()

		if *oneline {
			// Condensed format
			fmt.Printf("%s %s\n", shortHash, snap.Message)
		} else {
			// Full format
			fmt.Printf("snapshot %s\n", shortHash)
			fmt.Printf("Author: %s\n", snap.Author)
			fmt.Printf("Date:   %s\n", snap.Timestamp)
			fmt.Printf("\n    %s\n\n", snap.Message)
		}
	}

	return nil
}

func handleBranch(args []string) error {
	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// No arguments: list branches
	if len(args) == 0 {
		return listBranches(r)
	}

	// Create new branch
	branchName := args[0]

	// Validate branch name
	if strings.Contains(branchName, " ") || strings.Contains(branchName, "/") {
		return fmt.Errorf("invalid branch name: %s", branchName)
	}

	err = r.CreateBranch(branchName)
	if err != nil {
		return err
	}

	// Get snapshot hash for display
	hash, _ := r.Refs().ResolveHEAD()
	shortHash := "unknown"
	if hash != "" && len(hash) > 8 {
		shortHash = hash[:8]
	}

	fmt.Printf("Created branch '%s' at %s\n", branchName, shortHash)
	return nil
}

func listBranches(r *repo.Repository) error {
	branches, err := r.ListBranches()
	if err != nil {
		return err
	}

	if len(branches) == 0 {
		fmt.Println("No branches yet (create a commit first)")
		return nil
	}

	currentBranch, _ := r.CurrentBranch()

	for _, branch := range branches {
		marker := " "
		if branch == currentBranch {
			marker = "*"
		}
		fmt.Printf("%s %s\n", marker, branch)
	}

	return nil
}

func handleCheckout(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("usage: chemvcs checkout <branch|snapshot>")
	}

	target := args[0]

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Try as branch name first
	branches, _ := r.ListBranches()
	isBranch := false
	for _, b := range branches {
		if b == target {
			isBranch = true
			break
		}
	}

	if isBranch {
		err = r.CheckoutBranch(target)
		if err != nil {
			return err
		}
		fmt.Printf("Switched to branch '%s'\n", target)

		// Restore working directory
		err = restoreWorkingDirectory(r)
		if err != nil {
			return fmt.Errorf("failed to restore working directory: %w", err)
		}
	} else {
		// Try as snapshot hash
		err = r.CheckoutSnapshot(target)
		if err != nil {
			return fmt.Errorf("reference not found: %s", target)
		}
		shortHash := target
		if len(target) > 8 {
			shortHash = target[:8]
		}
		fmt.Printf("Switched to snapshot %s (detached HEAD)\n", shortHash)

		// Restore working directory
		err = restoreWorkingDirectory(r)
		if err != nil {
			return fmt.Errorf("failed to restore working directory: %w", err)
		}
	}

	return nil
}

func handleStatus(args []string) error {
	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Get current HEAD snapshot
	currentHash, err := r.Refs().ResolveHEAD()
	if err != nil {
		return err
	}

	if currentHash == "" {
		fmt.Println("No commits yet")
		return nil
	}

	// Get current snapshot
	currentSnap, err := r.GetSnapshot(currentHash)
	if err != nil {
		return fmt.Errorf("failed to get current snapshot: %w", err)
	}

	// Scan working directory
	scanner := workspace.NewScanner(r.Store())
	workingRoot, err := scanner.ScanDirectory(".")
	if err != nil {
		return fmt.Errorf("failed to scan working directory: %w", err)
	}

	// Store working root temporarily (for comparison)
	workingRootHash, err := r.Store().PutObject(workingRoot)
	if err != nil {
		return fmt.Errorf("failed to store working root: %w", err)
	}

	// Compute differences
	changes, err := workspace.ComputeChanges(currentSnap.Root, workingRootHash, r.Store())
	if err != nil {
		return fmt.Errorf("failed to compute changes: %w", err)
	}

	// Display current branch
	currentBranch, _ := r.CurrentBranch()
	if currentBranch != "" {
		fmt.Printf("On branch %s\n", currentBranch)
	} else {
		shortHash := currentHash
		if len(currentHash) > 8 {
			shortHash = currentHash[:8]
		}
		fmt.Printf("HEAD detached at %s\n", shortHash)
	}

	// Check if there are any changes
	if len(changes.Added) == 0 && len(changes.Modified) == 0 && len(changes.Deleted) == 0 {
		fmt.Println("\nNo changes (working directory clean)")
		return nil
	}

	// Print changes
	fmt.Println("\nChanges in working directory:")
	fmt.Println()

	if len(changes.Added) > 0 {
		for _, path := range changes.Added {
			fmt.Printf("  \x1b[32mA\x1b[0m  %s\n", path)
		}
	}

	if len(changes.Modified) > 0 {
		for _, path := range changes.Modified {
			fmt.Printf("  \x1b[33mM\x1b[0m  %s\n", path)
		}
	}

	if len(changes.Deleted) > 0 {
		for _, path := range changes.Deleted {
			fmt.Printf("  \x1b[31mD\x1b[0m  %s\n", path)
		}
	}

	return nil
}

func restoreWorkingDirectory(r *repo.Repository) error {
	// Get current HEAD snapshot
	currentHash, err := r.Refs().ResolveHEAD()
	if err != nil {
		return err
	}

	if currentHash == "" {
		// No snapshot to restore
		return nil
	}

	// Get snapshot
	snap, err := r.GetSnapshot(currentHash)
	if err != nil {
		return fmt.Errorf("failed to get snapshot: %w", err)
	}

	// Restore working directory
	scanner := workspace.NewScanner(r.Store())
	if err := scanner.RestoreDirectory(snap.Root, "."); err != nil {
		return err
	}

	return nil
}

func handleMerge(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("usage: chemvcs merge <branch>")
	}

	branchName := args[0]

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Check if HEAD is detached
	detached, err := r.IsDetached()
	if err != nil {
		return err
	}
	if detached {
		return fmt.Errorf("cannot merge: HEAD is detached (checkout a branch first)")
	}

	// Get current branch for display
	currentBranch, err := r.CurrentBranch()
	if err != nil {
		return err
	}

	// Perform merge
	result, err := r.Merge(branchName)
	if err != nil {
		if result != nil && len(result.Conflicts) > 0 {
			fmt.Fprintf(os.Stderr, "Merge conflicts detected in the following files:\n")
			for _, conflict := range result.Conflicts {
				fmt.Fprintf(os.Stderr, "  %s\n", conflict)
			}
			return fmt.Errorf("merge aborted due to conflicts")
		}
		return err
	}

	// Get new HEAD hash for display
	newHash, err := r.Refs().ResolveHEAD()
	if err != nil {
		return err
	}

	shortHash := newHash
	if len(newHash) > 8 {
		shortHash = newHash[:8]
	}

	if result.FastForward {
		fmt.Printf("Fast-forward merge '%s' into '%s'\n", branchName, currentBranch)
	} else if len(result.Conflicts) == 0 {
		// Three-way merge succeeded
		fmt.Printf("Merge branch '%s' into '%s'\n", branchName, currentBranch)
	} else {
		// Already up to date
		fmt.Printf("Already up-to-date with '%s'\n", branchName)
		return nil
	}

	fmt.Printf("Updated to %s\n", shortHash)

	// Restore working directory to match new HEAD
	err = restoreWorkingDirectory(r)
	if err != nil {
		return fmt.Errorf("merge succeeded but failed to update working directory: %w", err)
	}

	return nil
}

func handleRemote(args []string) error {
	if len(args) == 0 {
		return fmt.Errorf("usage: chemvcs remote add <name> <url>")
	}

	subCmd := args[0]
	subArgs := args[1:]

	switch subCmd {
	case "add":
		return handleRemoteAdd(subArgs)
	default:
		return fmt.Errorf("unknown remote subcommand: %s", subCmd)
	}
}

func handleRemoteAdd(args []string) error {
	if len(args) < 2 {
		return fmt.Errorf("usage: chemvcs remote add <name> <url>")
	}

	name := args[0]
	url := args[1]

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	if err := remote.AddRemoteConfig(r, name, url); err != nil {
		return err
	}

	fmt.Printf("Added remote '%s' at %s\n", name, url)
	return nil
}

func handlePush(args []string) error {
	if len(args) < 2 {
		return fmt.Errorf("usage: chemvcs push <remote> <branch>")
	}

	remoteName := args[0]
	branchName := args[1]

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Get remote URL
	remoteURL, err := remote.GetRemoteURL(r, remoteName)
	if err != nil {
		return err
	}

	// Create client
	client := remote.NewClient(remoteURL, remoteTokenFor(remoteName))

	// Parse repository ID from URL (simplified: assume URL ends with repo ID)
	// In a full implementation, this would be more sophisticated
	repoID := extractRepoID(remoteURL)

	// Set up push options
	localRef := "refs/heads/" + branchName
	opts := remote.PushOptions{
		RepoID:    repoID,
		LocalRef:  localRef,
		RemoteRef: localRef,
		Force:     false,
	}

	// Perform push
	fmt.Printf("Pushing %s to %s (%s)...\n", branchName, remoteName, remoteURL)
	if err := client.Push(r, opts); err != nil {
		return err
	}

	fmt.Printf("Successfully pushed %s to %s\n", branchName, remoteName)
	return nil
}

func handlePull(args []string) error {
	if len(args) < 2 {
		return fmt.Errorf("usage: chemvcs pull <remote> <branch>")
	}

	remoteName := args[0]
	branchName := args[1]

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Get remote URL
	remoteURL, err := remote.GetRemoteURL(r, remoteName)
	if err != nil {
		return err
	}

	// Create client
	client := remote.NewClient(remoteURL, remoteTokenFor(remoteName))

	repoID := extractRepoID(remoteURL)

	// Set up pull options
	remoteRef := "refs/heads/" + branchName
	localRef := "refs/heads/" + branchName
	opts := remote.PullOptions{
		RepoID:    repoID,
		RemoteRef: remoteRef,
		LocalRef:  localRef,
	}

	// Perform pull
	fmt.Printf("Pulling %s from %s (%s)...\n", branchName, remoteName, remoteURL)
	if err := client.Pull(r, opts); err != nil {
		return err
	}

	// Restore working directory
	newHash, err := r.Refs().ResolveHEAD()
	if err != nil {
		return err
	}

	scanner := workspace.NewScanner(r.Store())
	if err := scanner.RestoreDirectory(newHash, r.Path()); err != nil {
		return fmt.Errorf("failed to restore working directory: %w", err)
	}

	fmt.Printf("Successfully pulled %s from %s\n", branchName, remoteName)
	return nil
}

func handleFetch(args []string) error {
	if len(args) < 2 {
		return fmt.Errorf("usage: chemvcs fetch <remote> <branch>")
	}

	remoteName := args[0]
	branchName := args[1]

	r, err := repo.Open(".")
	if err != nil {
		return err
	}

	// Get remote URL
	remoteURL, err := remote.GetRemoteURL(r, remoteName)
	if err != nil {
		return err
	}

	// Create client
	client := remote.NewClient(remoteURL, remoteTokenFor(remoteName))

	repoID := extractRepoID(remoteURL)

	// Set up fetch options
	remoteRef := "refs/heads/" + branchName
	opts := remote.FetchOptions{
		RepoID:    repoID,
		RemoteRef: remoteRef,
	}

	// Perform fetch
	fmt.Printf("Fetching %s from %s (%s)...\n", branchName, remoteName, remoteURL)
	snapID, err := client.Fetch(r, opts)
	if err != nil {
		return err
	}

	shortHash := snapID
	if len(snapID) > 8 {
		shortHash = snapID[:8]
	}

	fmt.Printf("Successfully fetched %s from %s (snapshot %s)\n", branchName, remoteName, shortHash)
	return nil
}

// extractRepoID extracts repository ID from a remote URL.
// This is a simplified implementation; a full version would handle various URL formats.
func extractRepoID(url string) string {
	// Remove trailing slash
	url = strings.TrimSuffix(url, "/")

	// For URLs like "http://server/chemvcs/v1/repos/owner/repo", extract "owner/repo"
	parts := strings.Split(url, "/repos/")
	if len(parts) > 1 {
		return parts[1]
	}

	// Fallback: use last two path components
	parts = strings.Split(url, "/")
	if len(parts) >= 2 {
		return parts[len(parts)-2] + "/" + parts[len(parts)-1]
	}

	return "default/repo"
}

// remoteTokenFor resolves the auth token to use for a given remote.
//
// Priority:
// 1) CHEMVCS_REMOTE_TOKEN_<REMOTE_NAME>
// 2) CHEMVCS_REMOTE_TOKEN
func remoteTokenFor(remoteName string) string {
	if remoteName != "" {
		key := "CHEMVCS_REMOTE_TOKEN_" + envKeySuffix(remoteName)
		if v := os.Getenv(key); v != "" {
			return v
		}
	}
	return os.Getenv("CHEMVCS_REMOTE_TOKEN")
}

func envKeySuffix(s string) string {
	s = strings.TrimSpace(s)
	if s == "" {
		return ""
	}
	var b strings.Builder
	b.Grow(len(s))
	for _, r := range s {
		if r >= 'a' && r <= 'z' {
			b.WriteRune(r - ('a' - 'A'))
			continue
		}
		if (r >= 'A' && r <= 'Z') || (r >= '0' && r <= '9') {
			b.WriteRune(r)
			continue
		}
		b.WriteByte('_')
	}
	return b.String()
}

func handleInspectObject(args []string) error {
	// Parse flags
	fs := flag.NewFlagSet("inspect-object", flag.ExitOnError)
	format := fs.String("format", "text", "Output format: text or json")
	if err := fs.Parse(args); err != nil {
		return err
	}

	if fs.NArg() < 1 {
		return fmt.Errorf("usage: chemvcs inspect-object <hash> [--format=json]")
	}

	hash := fs.Arg(0)

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	// Get object
	obj, err := r.Store().GetObject(hash)
	if err != nil {
		return fmt.Errorf("failed to get object: %w", err)
	}

	// Output based on format
	if *format == "json" {
		// Pretty-print JSON
		data, err := obj.MarshalCanonical()
		if err != nil {
			return fmt.Errorf("failed to marshal object: %w", err)
		}
		fmt.Println(string(data))
	} else {
		// Human-readable text format
		fmt.Printf("Object: %s\n", hash)
		fmt.Printf("Type: %s\n", obj.Type)
		fmt.Printf("Version: %d\n", obj.Version)

		if len(obj.Meta) > 0 {
			fmt.Println("\nMetadata:")
			for k, v := range obj.Meta {
				fmt.Printf("  %s: %v\n", k, v)
			}
		}

		if len(obj.Refs) > 0 {
			fmt.Println("\nReferences:")
			for i, ref := range obj.Refs {
				fmt.Printf("  [%d] kind=%s id=%s\n", i, ref.Kind, ref.ID)
			}
		}
	}

	return nil
}

func handleListObjects(args []string) error {
	// Parse flags
	fs := flag.NewFlagSet("list-objects", flag.ExitOnError)
	typeFilter := fs.String("type", "", "Filter by object type")
	format := fs.String("format", "text", "Output format: text or json")
	if err := fs.Parse(args); err != nil {
		return err
	}

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	// List objects
	objects, err := r.Store().ListObjects(*typeFilter)
	if err != nil {
		return fmt.Errorf("failed to list objects: %w", err)
	}

	// Output based on format
	if *format == "json" {
		// Ensure we output [] instead of null for empty results
		if objects == nil {
			objects = []objectstore.ObjectInfo{}
		}
		// Output as JSON array
		data, err := json.Marshal(objects)
		if err != nil {
			return fmt.Errorf("failed to marshal objects: %w", err)
		}
		fmt.Println(string(data))
	} else {
		// Human-readable text format
		if len(objects) == 0 {
			fmt.Println("No objects found")
			return nil
		}

		fmt.Printf("Found %d object(s):\n", len(objects))
		for _, obj := range objects {
			shortHash := obj.Hash
			if len(obj.Hash) > 12 {
				shortHash = obj.Hash[:12]
			}
			fmt.Printf("  %s  %s\n", shortHash, obj.Type)
		}
	}

	return nil
}

func handleSubmit(args []string) error {
	// Parse flags
	fs := flag.NewFlagSet("submit", flag.ExitOnError)
	captureEnv := fs.Bool("capture-env", true, "Capture environment (modules, env vars)")
	remoteName := fs.String("remote", "", "Remote name (submit via chemvcs-server HPC gateway)")
	if err := fs.Parse(args); err != nil {
		return err
	}

	// Check arguments
	remainingArgs := fs.Args()
	if len(remainingArgs) < 2 {
		return fmt.Errorf("usage: chemvcs submit <run-hash> <script> [--capture-env]")
	}

	runHash := remainingArgs[0]
	scriptPath := remainingArgs[1]

	// Verify script exists
	if _, err := os.Stat(scriptPath); os.IsNotExist(err) {
		return fmt.Errorf("script file not found: %s", scriptPath)
	}

	// Get absolute script path
	absScriptPath, err := filepath.Abs(scriptPath)
	if err != nil {
		return fmt.Errorf("failed to resolve script path: %w", err)
	}

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	repoPath := r.Path()

	// Submit job
	fmt.Printf("Submitting HPC job for run %s...\n", runHash[:8])

	if *remoteName != "" {
		remoteURL, err := remote.GetRemoteURL(r, *remoteName)
		if err != nil {
			return err
		}
		client := remote.NewClient(remoteURL, remoteTokenFor(*remoteName))
		repoID := extractRepoID(remoteURL)

		scriptBytes, err := os.ReadFile(absScriptPath)
		if err != nil {
			return fmt.Errorf("failed to read script: %w", err)
		}

		resp, err := client.SubmitHPC(repoID, remote.SubmitHPCRequest{
			RunHash:    runHash,
			Script:     string(scriptBytes),
			CaptureEnv: *captureEnv,
		})
		if err != nil {
			return fmt.Errorf("failed to submit job: %w", err)
		}

		fmt.Printf("Job submitted successfully!\n")
		fmt.Printf("  Job ID: %s\n", resp.JobID)
		fmt.Printf("  Job System: %s\n", resp.JobSystem)
		fmt.Printf("  Run: %s\n", runHash[:8])
		return nil
	}

	result, err := hpc.SubmitJob(repoPath, hpc.SubmitOptions{
		RunHash:    runHash,
		ScriptPath: absScriptPath,
		CaptureEnv: *captureEnv,
	})
	if err != nil {
		return fmt.Errorf("failed to submit job: %w", err)
	}

	fmt.Printf("Job submitted successfully!\n")
	fmt.Printf("  Job ID: %s\n", result.JobID)
	fmt.Printf("  Job System: %s\n", result.JobSystem)
	fmt.Printf("  Run: %s\n", runHash[:8])

	return nil
}

func handleJobs(args []string) error {
	// Parse flags
	fs := flag.NewFlagSet("jobs", flag.ExitOnError)
	statusFilter := fs.String("status", "", "Filter by status (e.g., RUNNING, COMPLETED)")
	verbose := fs.Bool("v", false, "Verbose output")
	remoteName := fs.String("remote", "", "Remote name (list via chemvcs-server HPC gateway)")
	if err := fs.Parse(args); err != nil {
		return err
	}

	identifier := ""
	if fs.NArg() >= 1 {
		identifier = fs.Arg(0)
	}

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	repoPath := r.Path()

	if *remoteName != "" {
		remoteURL, err := remote.GetRemoteURL(r, *remoteName)
		if err != nil {
			return err
		}
		client := remote.NewClient(remoteURL, remoteTokenFor(*remoteName))
		repoID := extractRepoID(remoteURL)

		jobs, err := client.ListHPCJobs(repoID, *statusFilter, true)
		if err != nil {
			return fmt.Errorf("failed to list jobs: %w", err)
		}
		if identifier != "" {
			jobs, err = filterRemoteJobs(jobs, identifier)
			if err != nil {
				return err
			}
		}

		if len(jobs) == 0 {
			if *statusFilter != "" {
				fmt.Printf("No jobs found with status: %s\n", *statusFilter)
			} else {
				fmt.Println("No tracked jobs found")
			}
			return nil
		}

		fmt.Printf("Found %d job(s):\n\n", len(jobs))
		for _, job := range jobs {
			shortHash := job.RunHash
			if len(job.RunHash) > 8 {
				shortHash = job.RunHash[:8]
			}

			fmt.Printf("Job ID: %s\n", job.JobID)
			fmt.Printf("  Run:        %s\n", shortHash)
			fmt.Printf("  Status:     %s\n", job.Status)
			fmt.Printf("  System:     %s\n", job.JobSystem)
			if job.Queue != "" {
				fmt.Printf("  Queue:      %s\n", job.Queue)
			}
			if *verbose {
				if job.Submitted != "" {
					fmt.Printf("  Submitted:  %s\n", job.Submitted)
				}
				if job.Updated != "" {
					fmt.Printf("  Updated:    %s\n", job.Updated)
				}
				if job.WorkingDir != "" {
					fmt.Printf("  Workdir:    %s\n", job.WorkingDir)
				}
				if job.Started != "" {
					fmt.Printf("  Started:    %s\n", job.Started)
				}
				if job.Ended != "" {
					fmt.Printf("  Ended:      %s\n", job.Ended)
				}
				if job.ElapsedSec != "" {
					fmt.Printf("  Elapsed(s): %s\n", job.ElapsedSec)
				}
				if job.ExitCode != "" {
					fmt.Printf("  ExitCode:   %s\n", job.ExitCode)
				}
				if job.Reason != "" {
					fmt.Printf("  Reason:     %s\n", job.Reason)
				}
			}
			fmt.Println()
		}

		return nil
	}

	// List jobs
	jobs, err := hpc.ListJobs(repoPath, *statusFilter)
	if err != nil {
		return fmt.Errorf("failed to list jobs: %w", err)
	}
	if identifier != "" {
		jobs, err = filterLocalJobs(jobs, identifier)
		if err != nil {
			return err
		}
	}

	if len(jobs) == 0 {
		if *statusFilter != "" {
			fmt.Printf("No jobs found with status: %s\n", *statusFilter)
		} else {
			fmt.Println("No tracked jobs found")
		}
		return nil
	}

	// Display jobs
	fmt.Printf("Found %d job(s):\n\n", len(jobs))
	for _, job := range jobs {
		shortHash := job.RunHash
		if len(job.RunHash) > 8 {
			shortHash = job.RunHash[:8]
		}

		fmt.Printf("Job ID: %s\n", job.JobID)
		fmt.Printf("  Run:        %s\n", shortHash)
		fmt.Printf("  Status:     %s\n", job.Status)
		fmt.Printf("  System:     %s\n", job.JobSystem)
		if job.Queue != "" {
			fmt.Printf("  Queue:      %s\n", job.Queue)
		}
		if *verbose {
			if job.Submitted != "" {
				fmt.Printf("  Submitted:  %s\n", job.Submitted)
			}
			if job.Updated != "" {
				fmt.Printf("  Updated:    %s\n", job.Updated)
			}
			if job.WorkingDir != "" {
				fmt.Printf("  Workdir:    %s\n", job.WorkingDir)
			}
			if job.Started != "" {
				fmt.Printf("  Started:    %s\n", job.Started)
			}
			if job.Ended != "" {
				fmt.Printf("  Ended:      %s\n", job.Ended)
			}
			if job.ElapsedSec != "" {
				fmt.Printf("  Elapsed(s): %s\n", job.ElapsedSec)
			}
			if job.ExitCode != "" {
				fmt.Printf("  ExitCode:   %s\n", job.ExitCode)
			}
			if job.Reason != "" {
				fmt.Printf("  Reason:     %s\n", job.Reason)
			}
		}
		fmt.Println()
	}

	return nil
}

func handleRetrieve(args []string) error {
	// Parse flags
	fs := flag.NewFlagSet("retrieve", flag.ExitOnError)
	patterns := fs.String("patterns", "", "Comma-separated file patterns (e.g., *.out,*.log)")
	destination := fs.String("dest", ".", "Destination directory for retrieved files")
	commitAfter := fs.Bool("commit", false, "Commit retrieved files as a new snapshot")
	commitMessage := fs.String("commit-message", "", "Commit message (optional; defaults to a generated message)")
	remoteName := fs.String("remote", "", "Remote name (retrieve via chemvcs-server HPC gateway)")
	if err := fs.Parse(args); err != nil {
		return err
	}

	// Check arguments
	remainingArgs := fs.Args()
	if len(remainingArgs) < 1 {
		return fmt.Errorf("usage: chemvcs retrieve <run-hash> [--patterns=<patterns>] [--dest=<path>] [--commit] [--commit-message=<msg>]")
	}

	runHash := remainingArgs[0]

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	repoPath := r.Path()

	// Parse patterns
	var patternList []string
	if *patterns != "" {
		patternList = strings.Split(*patterns, ",")
	}

	// Retrieve results
	fmt.Printf("Retrieving results for run %s...\n", runHash[:8])
	if *remoteName != "" {
		remoteURL, err := remote.GetRemoteURL(r, *remoteName)
		if err != nil {
			return err
		}
		client := remote.NewClient(remoteURL, remoteTokenFor(*remoteName))
		repoID := extractRepoID(remoteURL)

		tmp, err := os.CreateTemp("", "chemvcs-results-*.zip")
		if err != nil {
			return fmt.Errorf("failed to create temp file: %w", err)
		}
		tmpPath := tmp.Name()
		defer func() { _ = os.Remove(tmpPath) }()

		if _, err := client.RetrieveHPCTo(repoID, patternList, tmp); err != nil {
			_ = tmp.Close()
			return fmt.Errorf("failed to retrieve results: %w", err)
		}
		if err := tmp.Close(); err != nil {
			return fmt.Errorf("failed to finalize download: %w", err)
		}

		files, err := remote.ExtractZipFileTo(tmpPath, *destination)
		if err != nil {
			return fmt.Errorf("failed to extract results: %w", err)
		}

		fmt.Printf("Retrieved %d file(s):\n", len(files))
		for _, file := range files {
			fmt.Printf("  %s\n", file)
		}

		if *commitAfter {
			if err := ensureDestinationInRepo(*destination); err != nil {
				return err
			}
			msg := *commitMessage
			if msg == "" {
				msg = fmt.Sprintf("Retrieve results for run %s", runHash)
			}
			authorStr := getAuthor("")
			snapHash, err := commitWorkingDirectory(r, msg, authorStr)
			if err != nil {
				return err
			}
			printCommitResult(r, snapHash, msg)
		}
		return nil
	}

	files, err := hpc.RetrieveResults(repoPath, hpc.RetrieveOptions{
		RunHash:     runHash,
		Patterns:    patternList,
		Destination: *destination,
	})
	if err != nil {
		return fmt.Errorf("failed to retrieve results: %w", err)
	}

	fmt.Printf("Retrieved %d file(s):\n", len(files))
	for _, file := range files {
		fmt.Printf("  %s\n", file)
	}

	if *commitAfter {
		if err := ensureDestinationInRepo(*destination); err != nil {
			return err
		}
		msg := *commitMessage
		if msg == "" {
			msg = fmt.Sprintf("Retrieve results for run %s", runHash)
		}
		authorStr := getAuthor("")
		snapHash, err := commitWorkingDirectory(r, msg, authorStr)
		if err != nil {
			return err
		}
		printCommitResult(r, snapHash, msg)
	}

	return nil
}

func ensureDestinationInRepo(destination string) error {
	absDest, err := filepath.Abs(destination)
	if err != nil {
		return fmt.Errorf("failed to resolve destination path: %w", err)
	}
	absRepo, err := filepath.Abs(".")
	if err != nil {
		return fmt.Errorf("failed to resolve repo path: %w", err)
	}

	rel, err := filepath.Rel(absRepo, absDest)
	if err != nil {
		return fmt.Errorf("failed to resolve destination relative path: %w", err)
	}
	if rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
		return fmt.Errorf("cannot use --commit when --dest is outside the repository (dest=%s)", destination)
	}
	return nil
}

func filterRemoteJobs(jobs []remote.HPCJobInfo, identifier string) ([]remote.HPCJobInfo, error) {
	identifier = strings.TrimSpace(identifier)
	if identifier == "" {
		return jobs, nil
	}

	var out []remote.HPCJobInfo
	if isDigits(identifier) {
		for _, j := range jobs {
			if j.JobID == identifier {
				out = append(out, j)
			}
		}
		if len(out) == 0 {
			return nil, fmt.Errorf("job not found: %s", identifier)
		}
		return out, nil
	}

	for _, j := range jobs {
		if strings.HasPrefix(j.RunHash, identifier) {
			out = append(out, j)
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("job not found: %s", identifier)
	}
	if len(out) > 1 {
		return nil, fmt.Errorf("ambiguous run hash %q (matched %d jobs); use a longer prefix", identifier, len(out))
	}
	return out, nil
}

func filterLocalJobs(jobs []hpc.JobInfo, identifier string) ([]hpc.JobInfo, error) {
	identifier = strings.TrimSpace(identifier)
	if identifier == "" {
		return jobs, nil
	}

	var out []hpc.JobInfo
	if isDigits(identifier) {
		for _, j := range jobs {
			if j.JobID == identifier {
				out = append(out, j)
			}
		}
		if len(out) == 0 {
			return nil, fmt.Errorf("job not found: %s", identifier)
		}
		return out, nil
	}

	for _, j := range jobs {
		if strings.HasPrefix(j.RunHash, identifier) {
			out = append(out, j)
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("job not found: %s", identifier)
	}
	if len(out) > 1 {
		return nil, fmt.Errorf("ambiguous run hash %q (matched %d jobs); use a longer prefix", identifier, len(out))
	}
	return out, nil
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

// handleCancel handles the cancel command
func handleCancel(args []string) error {
	fs := flag.NewFlagSet("cancel", flag.ExitOnError)
	remoteName := fs.String("remote", "", "Remote name (cancel via chemvcs-server HPC gateway)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() < 1 {
		return fmt.Errorf("usage: chemvcs cancel <run-hash|job-id> [--remote=<name>]")
	}

	identifier := fs.Arg(0)

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	repoPath := r.Path()

	if *remoteName != "" {
		remoteURL, err := remote.GetRemoteURL(r, *remoteName)
		if err != nil {
			return err
		}
		client := remote.NewClient(remoteURL, remoteTokenFor(*remoteName))
		repoID := extractRepoID(remoteURL)

		fmt.Printf("Cancelling job/run %s...\n", identifier)
		_, err = client.CancelHPC(repoID, identifier)
		if err != nil {
			return fmt.Errorf("failed to cancel job: %w", err)
		}
		fmt.Println("Job cancelled successfully!")
		return nil
	}

	// Cancel job
	fmt.Printf("Cancelling job/run %s...\n", identifier)
	err = hpc.CancelJob(repoPath, identifier)
	if err != nil {
		return fmt.Errorf("failed to cancel job: %w", err)
	}

	fmt.Println("Job cancelled successfully!")

	return nil
}

// handleWatch handles the watch command
func handleWatch(args []string) error {
	fs := flag.NewFlagSet("watch", flag.ExitOnError)
	interval := fs.Int("interval", 30, "Polling interval in seconds")
	timeout := fs.Int("timeout", 0, "Maximum watch time in seconds (0 for indefinite)")
	remoteName := fs.String("remote", "", "Remote name (watch via chemvcs-server HPC gateway)")

	if err := fs.Parse(args); err != nil {
		return err
	}

	// Check arguments
	remainingArgs := fs.Args()
	if len(remainingArgs) < 1 {
		return fmt.Errorf("usage: chemvcs watch <run-hash|job-id> [--interval=<sec>] [--timeout=<sec>]")
	}

	identifier := remainingArgs[0]

	// Open repository
	r, err := repo.Open(".")
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	repoPath := r.Path()

	if *remoteName != "" {
		remoteURL, err := remote.GetRemoteURL(r, *remoteName)
		if err != nil {
			return err
		}
		client := remote.NewClient(remoteURL, remoteTokenFor(*remoteName))
		repoID := extractRepoID(remoteURL)

		start := time.Now()
		lastStatus := ""
		for {
			job, err := client.GetHPCJob(repoID, identifier, true)
			if err != nil {
				return fmt.Errorf("failed to watch job: %w", err)
			}
			if job.Status != lastStatus {
				fmt.Printf("Status: %s\n", job.Status)
				lastStatus = job.Status
			}

			upper := strings.ToUpper(job.Status)
			if upper == "COMPLETED" || upper == "FAILED" || upper == "CANCELLED" {
				return nil
			}

			if *timeout > 0 && time.Since(start) > time.Duration(*timeout)*time.Second {
				return fmt.Errorf("watch timeout reached")
			}

			time.Sleep(time.Duration(*interval) * time.Second)
		}
	}

	// Watch job
	err = hpc.WatchJob(repoPath, identifier, *interval, *timeout)
	if err != nil {
		return fmt.Errorf("failed to watch job: %w", err)
	}

	return nil
}

// handlePack packs loose objects into a packfile.
func handlePack(args []string) error {
	fs := flag.NewFlagSet("pack", flag.ExitOnError)
	all := fs.Bool("all", false, "Pack all loose objects (not just reachable)")
	keepLoose := fs.Bool("keep-loose", false, "Keep loose objects after packing")
	verbose := fs.Bool("v", false, "Verbose output")
	fs.Parse(args)

	repoPath, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get working directory: %w", err)
	}

	repoInstance, err := repo.Open(repoPath)
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	opts := objectstore.PackOptions{
		All:       *all,
		KeepLoose: *keepLoose,
		Verbose:   *verbose,
	}

	packName, stats, err := repoInstance.Store().PackLooseObjects(repoInstance, opts)
	if err != nil {
		return fmt.Errorf("pack failed: %w", err)
	}

	if stats.ObjectsPacked == 0 {
		fmt.Println("Nothing to pack")
		return nil
	}

	if !*verbose {
		ratio := 100.0 * (1.0 - float64(stats.BytesAfterPack)/float64(stats.BytesBeforePack))
		fmt.Printf("Created pack %s with %d objects (%.1f%% size reduction)\n",
			packName, stats.ObjectsPacked, ratio)
	}

	return nil
}

// handleGC performs garbage collection.
func handleGC(args []string) error {
	fs := flag.NewFlagSet("gc", flag.ExitOnError)
	prune := fs.String("prune", "2w", "Prune objects older than <time> (e.g., 1h, 2w, now)")
	dryRun := fs.Bool("dry-run", false, "Show what would be deleted without actually deleting")
	verbose := fs.Bool("v", false, "Verbose output")
	fs.Parse(args)

	repoPath, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get working directory: %w", err)
	}

	repoInstance, err := repo.Open(repoPath)
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	// Parse prune age
	pruneAge, err := parseDuration(*prune)
	if err != nil {
		return fmt.Errorf("invalid prune duration: %w", err)
	}

	opts := objectstore.GCOptions{
		PruneAge: pruneAge,
		DryRun:   *dryRun,
		Verbose:  *verbose,
	}

	stats, err := repoInstance.Store().GarbageCollect(repoInstance, opts)
	if err != nil {
		return fmt.Errorf("garbage collection failed: %w", err)
	}

	if !*verbose {
		if *dryRun {
			fmt.Printf("Would remove %d unreachable objects\n", stats.UnreachableObjects)
		} else {
			fmt.Printf("Removed %d objects, reclaimed %.2f MB\n",
				stats.RemovedObjects, float64(stats.ReclaimedBytes)/(1024*1024))
		}
	}

	return nil
}

// handleFsck performs filesystem check.
func handleFsck(args []string) error {
	fs := flag.NewFlagSet("fsck", flag.ExitOnError)
	full := fs.Bool("full", false, "Full check including packfile verification")
	verbose := fs.Bool("v", false, "Verbose output")
	fs.Parse(args)

	repoPath, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get working directory: %w", err)
	}

	repoInstance, err := repo.Open(repoPath)
	if err != nil {
		return fmt.Errorf("failed to open repository: %w", err)
	}

	opts := objectstore.FsckOptions{
		Full:    *full,
		Verbose: *verbose,
	}

	report, err := repoInstance.Store().Fsck(repoInstance, opts)
	if err != nil {
		return fmt.Errorf("fsck failed: %w", err)
	}

	// Display issues
	if len(report.Issues) > 0 {
		fmt.Printf("\nIssues found:\n")
		for _, issue := range report.Issues {
			prefix := ""
			switch issue.Severity {
			case "error":
				prefix = "ERROR"
			case "warning":
				prefix = "WARN"
			case "info":
				prefix = "INFO"
			}

			if issue.Object != "" {
				fmt.Printf("[%s] %s: %s\n", prefix, issue.Object[:8], issue.Message)
			} else {
				fmt.Printf("[%s] %s\n", prefix, issue.Message)
			}
		}
	}

	if !*verbose {
		if report.HasErrors {
			fmt.Printf("\nFilesystem check FAILED with %d errors\n", len(report.Issues))
			os.Exit(1)
		} else if len(report.Issues) > 0 {
			fmt.Printf("\nFilesystem check completed with %d warnings\n", len(report.Issues))
		} else {
			fmt.Println("\nFilesystem check OK")
		}
	}

	return nil
}

// parseDuration parses a duration string like "1h", "2w", "now".
func parseDuration(s string) (time.Duration, error) {
	if s == "now" {
		return 0, nil
	}

	// Parse format like "1h", "2w", "30m"
	if len(s) < 2 {
		return 0, fmt.Errorf("invalid duration format")
	}

	valueStr := s[:len(s)-1]
	unit := s[len(s)-1:]

	var value int
	_, err := fmt.Sscanf(valueStr, "%d", &value)
	if err != nil {
		return 0, fmt.Errorf("invalid duration value: %w", err)
	}

	switch unit {
	case "s":
		return time.Duration(value) * time.Second, nil
	case "m":
		return time.Duration(value) * time.Minute, nil
	case "h":
		return time.Duration(value) * time.Hour, nil
	case "d":
		return time.Duration(value) * 24 * time.Hour, nil
	case "w":
		return time.Duration(value) * 7 * 24 * time.Hour, nil
	default:
		return 0, fmt.Errorf("unknown duration unit: %s", unit)
	}
}
