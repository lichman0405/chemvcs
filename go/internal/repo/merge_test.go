package repo

import (
	"testing"

	"github.com/lishi/chemvcs/internal/model"
)

func TestFindCommonAncestor_Linear(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create linear history: commit1 -> commit2 -> commit3
	obj1 := model.NewObject("generic")
	root1, _ := repo.Store().PutObject(obj1)
	hash1, _ := repo.CreateSnapshot(root1, "Commit 1", "Author")

	obj2 := model.NewObject("generic")
	root2, _ := repo.Store().PutObject(obj2)
	hash2, _ := repo.CreateSnapshot(root2, "Commit 2", "Author")

	obj3 := model.NewObject("generic")
	root3, _ := repo.Store().PutObject(obj3)
	hash3, _ := repo.CreateSnapshot(root3, "Commit 3", "Author")

	// Common ancestor of hash2 and hash3 should be hash2
	ancestor, err := repo.FindCommonAncestor(hash2, hash3)
	if err != nil {
		t.Fatalf("FindCommonAncestor failed: %v", err)
	}

	if ancestor != hash2 {
		t.Errorf("expected ancestor %s, got %s", hash2, ancestor)
	}

	// Common ancestor of hash1 and hash3 should be hash1
	ancestor, err = repo.FindCommonAncestor(hash1, hash3)
	if err != nil {
		t.Fatalf("FindCommonAncestor failed: %v", err)
	}

	if ancestor != hash1 {
		t.Errorf("expected ancestor %s, got %s", hash1, ancestor)
	}
}

func TestFindCommonAncestor_Diverged(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create base commit
	obj1 := model.NewObject("generic")
	root1, _ := repo.Store().PutObject(obj1)
	baseHash, _ := repo.CreateSnapshot(root1, "Base commit", "Author")

	// Create branch and make commit on it
	repo.CreateBranch("feature")
	repo.CheckoutBranch("feature")

	obj2 := model.NewObject("generic")
	root2, _ := repo.Store().PutObject(obj2)
	featureHash, _ := repo.CreateSnapshot(root2, "Feature commit", "Author")

	// Switch back to main and make commit
	repo.CheckoutBranch("main")

	obj3 := model.NewObject("generic")
	root3, _ := repo.Store().PutObject(obj3)
	mainHash, _ := repo.CreateSnapshot(root3, "Main commit", "Author")

	// Common ancestor should be base commit
	ancestor, err := repo.FindCommonAncestor(mainHash, featureHash)
	if err != nil {
		t.Fatalf("FindCommonAncestor failed: %v", err)
	}

	if ancestor != baseHash {
		t.Errorf("expected ancestor %s, got %s", baseHash, ancestor)
	}
}

func TestThreeWayMerge_NoConflicts(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create base: folder with file1.txt
	baseFolder := model.NewObject("folder")
	baseFolder.SetMeta("name", "root")

	file1 := model.NewObject("file")
	file1.SetMeta("name", "file1.txt")
	blob1Hash, _ := repo.Store().PutBlob([]byte("original content"))
	file1.Refs = []model.Reference{{Kind: "blob", ID: blob1Hash}}
	file1Hash, _ := repo.Store().PutObject(file1)

	baseFolder.Refs = []model.Reference{{Kind: "object", ID: file1Hash}}
	baseFolderHash, _ := repo.Store().PutObject(baseFolder)

	// Current: added file2.txt
	currentFolder := model.NewObject("folder")
	currentFolder.SetMeta("name", "root")

	file2 := model.NewObject("file")
	file2.SetMeta("name", "file2.txt")
	blob2Hash, _ := repo.Store().PutBlob([]byte("current content"))
	file2.Refs = []model.Reference{{Kind: "blob", ID: blob2Hash}}
	file2Hash, _ := repo.Store().PutObject(file2)

	currentFolder.Refs = []model.Reference{
		{Kind: "object", ID: file1Hash},
		{Kind: "object", ID: file2Hash},
	}
	currentFolderHash, _ := repo.Store().PutObject(currentFolder)

	// Target: added file3.txt
	targetFolder := model.NewObject("folder")
	targetFolder.SetMeta("name", "root")

	file3 := model.NewObject("file")
	file3.SetMeta("name", "file3.txt")
	blob3Hash, _ := repo.Store().PutBlob([]byte("target content"))
	file3.Refs = []model.Reference{{Kind: "blob", ID: blob3Hash}}
	file3Hash, _ := repo.Store().PutObject(file3)

	targetFolder.Refs = []model.Reference{
		{Kind: "object", ID: file1Hash},
		{Kind: "object", ID: file3Hash},
	}
	targetFolderHash, _ := repo.Store().PutObject(targetFolder)

	// Merge
	mergedObj, conflicts, err := repo.ThreeWayMerge(baseFolderHash, currentFolderHash, targetFolderHash)
	if err != nil {
		t.Fatalf("ThreeWayMerge failed: %v", err)
	}

	if len(conflicts) > 0 {
		t.Errorf("unexpected conflicts: %v", conflicts)
	}

	// Verify merged folder has all three files
	if len(mergedObj.Refs) != 3 {
		t.Errorf("expected 3 files in merged folder, got %d", len(mergedObj.Refs))
	}
}

func TestThreeWayMerge_WithConflicts(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create base: file1.txt with original content
	baseFolder := model.NewObject("folder")
	baseFolder.SetMeta("name", "root")

	file1 := model.NewObject("file")
	file1.SetMeta("name", "file1.txt")
	blob1Hash, _ := repo.Store().PutBlob([]byte("original content"))
	file1.Refs = []model.Reference{{Kind: "blob", ID: blob1Hash}}
	file1Hash, _ := repo.Store().PutObject(file1)

	baseFolder.Refs = []model.Reference{{Kind: "object", ID: file1Hash}}
	baseFolderHash, _ := repo.Store().PutObject(baseFolder)

	// Current: file1.txt modified to "current content"
	currentFolder := model.NewObject("folder")
	currentFolder.SetMeta("name", "root")

	currentFile1 := model.NewObject("file")
	currentFile1.SetMeta("name", "file1.txt")
	currentBlob1Hash, _ := repo.Store().PutBlob([]byte("current content"))
	currentFile1.Refs = []model.Reference{{Kind: "blob", ID: currentBlob1Hash}}
	currentFile1Hash, _ := repo.Store().PutObject(currentFile1)

	currentFolder.Refs = []model.Reference{{Kind: "object", ID: currentFile1Hash}}
	currentFolderHash, _ := repo.Store().PutObject(currentFolder)

	// Target: file1.txt modified to "target content"
	targetFolder := model.NewObject("folder")
	targetFolder.SetMeta("name", "root")

	targetFile1 := model.NewObject("file")
	targetFile1.SetMeta("name", "file1.txt")
	targetBlob1Hash, _ := repo.Store().PutBlob([]byte("target content"))
	targetFile1.Refs = []model.Reference{{Kind: "blob", ID: targetBlob1Hash}}
	targetFile1Hash, _ := repo.Store().PutObject(targetFile1)

	targetFolder.Refs = []model.Reference{{Kind: "object", ID: targetFile1Hash}}
	targetFolderHash, _ := repo.Store().PutObject(targetFolder)

	// Merge - should detect conflict
	_, conflicts, err := repo.ThreeWayMerge(baseFolderHash, currentFolderHash, targetFolderHash)
	if err != nil {
		t.Fatalf("ThreeWayMerge failed: %v", err)
	}

	if len(conflicts) == 0 {
		t.Error("expected conflicts, got none")
	}

	// Verify conflict is for file1.txt
	found := false
	for _, conflict := range conflicts {
		if conflict == "file1.txt" {
			found = true
			break
		}
	}

	if !found {
		t.Errorf("expected conflict for file1.txt, got %v", conflicts)
	}
}

func TestMerge_ThreeWayIntegration(t *testing.T) {
	tmpDir := t.TempDir()
	repo, err := Init(tmpDir)
	if err != nil {
		t.Fatalf("Init failed: %v", err)
	}

	// Create base commit with file1.txt
	baseFolder := model.NewObject("folder")
	baseFolder.SetMeta("name", "root")

	file1 := model.NewObject("file")
	file1.SetMeta("name", "file1.txt")
	blob1Hash, _ := repo.Store().PutBlob([]byte("base content"))
	file1.Refs = []model.Reference{{Kind: "blob", ID: blob1Hash}}
	file1Hash, _ := repo.Store().PutObject(file1)

	baseFolder.Refs = []model.Reference{{Kind: "object", ID: file1Hash}}
	baseFolderHash, _ := repo.Store().PutObject(baseFolder)

	_, err = repo.CreateSnapshot(baseFolderHash, "Base commit", "Author")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Create feature branch and add file2.txt
	if err := repo.CreateBranch("feature"); err != nil {
		t.Fatalf("CreateBranch failed: %v", err)
	}
	if err := repo.CheckoutBranch("feature"); err != nil {
		t.Fatalf("CheckoutBranch failed: %v", err)
	}

	featureFolder := model.NewObject("folder")
	featureFolder.SetMeta("name", "root")

	file2 := model.NewObject("file")
	file2.SetMeta("name", "file2.txt")
	blob2Hash, _ := repo.Store().PutBlob([]byte("feature content"))
	file2.Refs = []model.Reference{{Kind: "blob", ID: blob2Hash}}
	file2Hash, _ := repo.Store().PutObject(file2)

	featureFolder.Refs = []model.Reference{
		{Kind: "object", ID: file1Hash},
		{Kind: "object", ID: file2Hash},
	}
	featureFolderHash, _ := repo.Store().PutObject(featureFolder)

	_, err = repo.CreateSnapshot(featureFolderHash, "Feature commit", "Author")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Switch to main and add file3.txt
	if err := repo.CheckoutBranch("main"); err != nil {
		t.Fatalf("CheckoutBranch failed: %v", err)
	}

	mainFolder := model.NewObject("folder")
	mainFolder.SetMeta("name", "root")

	file3 := model.NewObject("file")
	file3.SetMeta("name", "file3.txt")
	blob3Hash, _ := repo.Store().PutBlob([]byte("main content"))
	file3.Refs = []model.Reference{{Kind: "blob", ID: blob3Hash}}
	file3Hash, _ := repo.Store().PutObject(file3)

	mainFolder.Refs = []model.Reference{
		{Kind: "object", ID: file1Hash},
		{Kind: "object", ID: file3Hash},
	}
	mainFolderHash, _ := repo.Store().PutObject(mainFolder)

	_, err = repo.CreateSnapshot(mainFolderHash, "Main commit", "Author")
	if err != nil {
		t.Fatalf("CreateSnapshot failed: %v", err)
	}

	// Merge feature into main
	result, err := repo.Merge("feature")
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if result.FastForward {
		t.Error("expected three-way merge, not fast-forward")
	}

	if len(result.Conflicts) > 0 {
		t.Errorf("unexpected conflicts: %v", result.Conflicts)
	}

	// Verify merge commit has both parents
	mainHash, _ := repo.Refs().ResolveRef("refs/heads/main")
	mergeSnapshot, _ := repo.Store().GetSnapshot(mainHash)

	if len(mergeSnapshot.Parents) != 2 {
		t.Errorf("expected 2 parents, got %d", len(mergeSnapshot.Parents))
	}

	// Verify merged tree has all three files
	mergedRoot, _ := repo.Store().GetObject(mergeSnapshot.Root)
	if len(mergedRoot.Refs) != 3 {
		t.Errorf("expected 3 files in merged root, got %d", len(mergedRoot.Refs))
	}
}
