package repo

import (
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"

	"github.com/lishi/chemvcs/internal/model"
)

func TestInitRepository(t *testing.T) {
	tmpDir := t.TempDir()
	repoPath := filepath.Join(tmpDir, "test-repo")

	if err := os.MkdirAll(repoPath, 0755); err != nil {
		t.Fatalf("failed to create test directory: %v", err)
	}

	repo, err := Init(repoPath)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Verify directory structure
	chemvcsPath := filepath.Join(repoPath, ".chemvcs")
	if _, err := os.Stat(chemvcsPath); os.IsNotExist(err) {
		t.Error(".chemvcs directory not created")
	}

	objectsPath := filepath.Join(chemvcsPath, "objects")
	if _, err := os.Stat(objectsPath); os.IsNotExist(err) {
		t.Error("objects directory not created")
	}

	headsPath := filepath.Join(chemvcsPath, "refs", "heads")
	if _, err := os.Stat(headsPath); os.IsNotExist(err) {
		t.Error("refs/heads directory not created")
	}

	// Verify HEAD
	headPath := filepath.Join(chemvcsPath, "HEAD")
	data, err := ioutil.ReadFile(headPath)
	if err != nil {
		t.Fatalf("failed to read HEAD: %v", err)
	}

	if string(data) != "ref: refs/heads/main" {
		t.Errorf("expected HEAD to be 'ref: refs/heads/main', got %q", string(data))
	}

	// Verify repo path
	if repo.Path() != repoPath {
		t.Errorf("expected path %s, got %s", repoPath, repo.Path())
	}
}

func TestInitAlreadyExists(t *testing.T) {
	tmpDir := t.TempDir()

	// Init once
	_, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("first Init failed: %v", err)
	}

	// Init again should fail
	_, err = Init(tmpDir)
	if err == nil {
		t.Error("expected error when initialising existing repository")
	}
}

func TestOpenRepository(t *testing.T) {
	tmpDir := t.TempDir()

	// Init first
	_, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Open
	repo, err := Open(tmpDir)
	if err != nil {
		t.Fatalf("Open failed: %v", err)
	}

	if repo.Path() != tmpDir {
		t.Errorf("expected path %s, got %s", tmpDir, repo.Path())
	}
}

func TestOpenFromSubdirectory(t *testing.T) {
	tmpDir := t.TempDir()

	// Init in root
	_, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create subdirectory
	subDir := filepath.Join(tmpDir, "subdir", "nested")
	if err := os.MkdirAll(subDir, 0755); err != nil {
		t.Fatalf("failed to create subdirectory: %v", err)
	}

	// Open from subdirectory should find parent repo
	repo, err := Open(subDir)
	if err != nil {
		t.Fatalf("Open from subdirectory failed: %v", err)
	}

	if repo.Path() != tmpDir {
		t.Errorf("expected path %s, got %s", tmpDir, repo.Path())
	}
}

func TestOpenNonExistent(t *testing.T) {
	tmpDir := t.TempDir()

	_, err := Open(tmpDir)
	if err == nil {
		t.Error("expected error when opening non-existent repository")
	}
}

func TestCreateSnapshot(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create a dummy root object
	rootObj := model.NewObject("generic")
	rootObj.SetMeta("test", true)
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	// Create snapshot
	snapHash, err := repo.CreateSnapshot(rootHash, "Initial commit", "Alice <alice@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	if snapHash == "" {
		t.Error("snapshot hash should not be empty")
	}

	// Verify snapshot was stored
	snap, err := repo.GetSnapshot(snapHash)
	if err != nil {
		t.Fatalf("failed to retrieve snapshot: %v", err)
	}

	if snap.Root != rootHash {
		t.Errorf("expected root %s, got %s", rootHash, snap.Root)
	}

	if snap.Message != "Initial commit" {
		t.Errorf("expected message 'Initial commit', got %q", snap.Message)
	}

	if len(snap.Parents) != 0 {
		t.Errorf("initial commit should have no parents, got %d", len(snap.Parents))
	}

	// Verify main branch was updated
	mainHash, err := repo.Refs().ResolveRef("refs/heads/main")
	if err != nil {
		t.Fatalf("failed to resolve main: %v", err)
	}

	if mainHash != snapHash {
		t.Errorf("main branch should point to %s, got %s", snapHash, mainHash)
	}
}

func TestCreateSnapshotWithParent(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create first snapshot
	rootObj1 := model.NewObject("generic")
	rootHash1, err := repo.Store().PutObject(rootObj1)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	snap1Hash, err := repo.CreateSnapshot(rootHash1, "First", "Alice <alice@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Create second snapshot
	rootObj2 := model.NewObject("generic")
	rootObj2.SetMeta("version", 2)
	rootHash2, err := repo.Store().PutObject(rootObj2)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	snap2Hash, err := repo.CreateSnapshot(rootHash2, "Second", "Bob <bob@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Verify second snapshot has first as parent
	snap2, err := repo.GetSnapshot(snap2Hash)
	if err != nil {
		t.Fatalf("failed to retrieve snapshot: %v", err)
	}

	if len(snap2.Parents) != 1 {
		t.Fatalf("expected 1 parent, got %d", len(snap2.Parents))
	}

	if snap2.Parents[0] != snap1Hash {
		t.Errorf("expected parent %s, got %s", snap1Hash, snap2.Parents[0])
	}
}

func TestLog(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create chain of snapshots
	messages := []string{"First", "Second", "Third"}
	for _, msg := range messages {
		rootObj := model.NewObject("generic")
		rootObj.SetMeta("message", msg)
		rootHash, err := repo.Store().PutObject(rootObj)
		if err != nil {
			t.Fatalf("failed to store root object: %v", err)
		}

		_, err = repo.CreateSnapshot(rootHash, msg, "Author <author@example.com>")
		if err != nil {
			t.Fatalf("CreateSnapshot failed: %v", err)
		}
	}

	// Get log
	snapshots, err := repo.Log(10)
	if err != nil {
		t.Fatalf("Log failed: %v", err)
	}

	if len(snapshots) != 3 {
		t.Fatalf("expected 3 snapshots, got %d", len(snapshots))
	}

	// Verify order (most recent first)
	for i, snap := range snapshots {
		expectedMsg := messages[len(messages)-1-i]
		if snap.Message != expectedMsg {
			t.Errorf("snapshot %d: expected message %q, got %q", i, expectedMsg, snap.Message)
		}
	}
}

func TestLogLimit(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create 5 snapshots
	for i := 1; i <= 5; i++ {
		rootObj := model.NewObject("generic")
		rootHash, err := repo.Store().PutObject(rootObj)
		if err != nil {
			t.Fatalf("failed to store root object: %v", err)
		}

		_, err = repo.CreateSnapshot(rootHash, "Commit", "Author <author@example.com>")
		if err != nil {
			t.Fatalf("CreateSnapshot failed: %v", err)
		}
	}

	// Get log with limit
	snapshots, err := repo.Log(2)
	if err != nil {
		t.Fatalf("Log failed: %v", err)
	}

	if len(snapshots) != 2 {
		t.Errorf("expected 2 snapshots with limit, got %d", len(snapshots))
	}
}

func TestLogEmpty(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// No commits yet
	snapshots, err := repo.Log(10)
	if err != nil {
		t.Fatalf("Log failed: %v", err)
	}

	if len(snapshots) != 0 {
		t.Errorf("expected 0 snapshots, got %d", len(snapshots))
	}
}

func TestCreateBranch(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create initial snapshot
	rootObj := model.NewObject("generic")
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	snapHash, err := repo.CreateSnapshot(rootHash, "Initial", "Author <author@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Create branch
	err = repo.CreateBranch("develop")
	if err != nil {
		t.Fatalf("CreateBranch failed: %v", err)
	}

	// Verify branch points to same snapshot
	developHash, err := repo.Refs().ResolveRef("refs/heads/develop")
	if err != nil {
		t.Fatalf("failed to resolve develop: %v", err)
	}

	if developHash != snapHash {
		t.Errorf("develop should point to %s, got %s", snapHash, developHash)
	}
}

func TestListBranches(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Initially no branches (main doesn't exist until first commit)
	branches, err := repo.ListBranches()
	if err != nil {
		t.Fatalf("ListBranches failed: %v", err)
	}

	if len(branches) != 0 {
		t.Errorf("expected 0 branches initially, got %d", len(branches))
	}

	// Create initial commit
	rootObj := model.NewObject("generic")
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	_, err = repo.CreateSnapshot(rootHash, "Initial", "Author <author@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Now main should exist
	branches, err = repo.ListBranches()
	if err != nil {
		t.Fatalf("ListBranches failed: %v", err)
	}

	if len(branches) != 1 || branches[0] != "main" {
		t.Errorf("expected [main], got %v", branches)
	}

	// Create more branches
	repo.CreateBranch("develop")
	repo.CreateBranch("feature")

	branches, err = repo.ListBranches()
	if err != nil {
		t.Fatalf("ListBranches failed: %v", err)
	}

	if len(branches) != 3 {
		t.Errorf("expected 3 branches, got %d", len(branches))
	}
}

func TestCheckoutBranch(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create initial commit
	rootObj := model.NewObject("generic")
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	_, err = repo.CreateSnapshot(rootHash, "Initial", "Author <author@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Create and checkout branch
	repo.CreateBranch("develop")

	err = repo.CheckoutBranch("develop")
	if err != nil {
		t.Fatalf("CheckoutBranch failed: %v", err)
	}

	// Verify current branch
	current, err := repo.CurrentBranch()
	if err != nil {
		t.Fatalf("CurrentBranch failed: %v", err)
	}

	if current != "develop" {
		t.Errorf("expected current branch 'develop', got %q", current)
	}
}

func TestCheckoutSnapshot(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create snapshot
	rootObj := model.NewObject("generic")
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	snapHash, err := repo.CreateSnapshot(rootHash, "Initial", "Author <author@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Checkout snapshot directly (detached HEAD)
	err = repo.CheckoutSnapshot(snapHash)
	if err != nil {
		t.Fatalf("CheckoutSnapshot failed: %v", err)
	}

	// Verify detached HEAD
	isDetached, err := repo.IsDetached()
	if err != nil {
		t.Fatalf("IsDetached failed: %v", err)
	}

	if !isDetached {
		t.Error("HEAD should be detached")
	}

	current, err := repo.CurrentBranch()
	if err != nil {
		t.Fatalf("CurrentBranch failed: %v", err)
	}

	if current != "" {
		t.Errorf("detached HEAD should have no current branch, got %q", current)
	}
}

func TestRefsResolveHEAD(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	refs := repo.Refs()

	// Initially HEAD points to non-existent main branch
	hash, err := refs.ResolveHEAD()
	if err != nil {
		t.Fatalf("ResolveHEAD failed: %v", err)
	}

	if hash != "" {
		t.Errorf("expected empty hash for non-existent branch, got %q", hash)
	}

	// Create commit
	rootObj := model.NewObject("generic")
	rootHash, err := repo.Store().PutObject(rootObj)
	if err != nil {
		t.Fatalf("failed to store root object: %v", err)
	}

	snapHash, err := repo.CreateSnapshot(rootHash, "Initial", "Author <author@example.com>")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Now HEAD should resolve
	hash, err = refs.ResolveHEAD()
	if err != nil {
		t.Fatalf("ResolveHEAD failed: %v", err)
	}

	if hash != snapHash {
		t.Errorf("expected HEAD to resolve to %s, got %s", snapHash, hash)
	}
}
