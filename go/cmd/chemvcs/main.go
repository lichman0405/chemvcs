package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

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

	// Scan working directory to create root object
	scanner := workspace.NewScanner(r.Store())
	rootObj, err := scanner.ScanDirectory(".")
	if err != nil {
		return fmt.Errorf("failed to scan working directory: %w", err)
	}

	// Store root object
	rootHash, err := r.Store().PutObject(rootObj)
	if err != nil {
		return fmt.Errorf("failed to store root object: %w", err)
	}

	// Create snapshot
	snapHash, err := r.CreateSnapshot(rootHash, *message, authorStr)
	if err != nil {
		return err
	}

	// Get short hash for display
	shortHash := snapHash
	if len(snapHash) > 8 {
		shortHash = snapHash[:8]
	}

	// Get current branch for display
	branch, _ := r.CurrentBranch()
	if branch != "" {
		fmt.Printf("[%s %s] %s\n", branch, shortHash, *message)
	} else {
		fmt.Printf("[detached %s] %s\n", shortHash, *message)
	}

	return nil
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
	client := remote.NewClient(remoteURL, "")

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
	client := remote.NewClient(remoteURL, "")

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
	client := remote.NewClient(remoteURL, "")

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

