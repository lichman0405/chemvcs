package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/lishi/chemvcs/internal/model"
	"github.com/lishi/chemvcs/internal/repo"
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

	// MVP: Create a simple root object
	// Full workspace scanning will be implemented in Milestone 2
	rootHash, err := createSimpleRoot(r)
	if err != nil {
		return fmt.Errorf("failed to create root object: %w", err)
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

// createSimpleRoot creates a minimal root object for MVP.
// This will be replaced with full workspace scanning in Milestone 2.
func createSimpleRoot(r *repo.Repository) (string, error) {
	obj := model.NewObject("generic")
	obj.SetMeta("mvp", true)
	obj.SetMeta("type", "simple_root")

	return r.Store().PutObject(obj)
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
	}

	// Note: Actual working directory restoration will be in Milestone 2
	fmt.Println("Note: Working directory update not yet implemented (Milestone 2)")

	return nil
}
