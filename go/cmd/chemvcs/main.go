package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

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
	merged, err := r.Merge(branchName)
	if err != nil {
		return err
	}

	if !merged {
		// Already up to date
		fmt.Printf("Already up-to-date with '%s'\n", branchName)
		return nil
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

	fmt.Printf("Fast-forward merge '%s' into '%s'\n", branchName, currentBranch)
	fmt.Printf("Updated to %s\n", shortHash)

	// Restore working directory to match new HEAD
	err = restoreWorkingDirectory(r)
	if err != nil {
		return fmt.Errorf("merge succeeded but failed to update working directory: %w", err)
	}

	return nil
}
