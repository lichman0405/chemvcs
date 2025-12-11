# ChemVCS Execution Plan (v0.1)

## Overview

This document provides a concrete, actionable execution plan for implementing ChemVCS Milestones 1-3. It translates the high-level development plan into specific tasks with implementation details, time estimates, and testing strategies.

---

## Phase 0: Project Foundation (1 day)

### Tasks

1. **Initialise Go module**
   ```bash
   cd go
   go mod init github.com/lishi/chemvcs
   ```

2. **Create directory structure**
   ```
   go/
   ├── cmd/
   │   └── chemvcs/
   │       └── main.go
   ├── internal/
   │   ├── model/
   │   │   ├── blob.go
   │   │   ├── object.go
   │   │   ├── snapshot.go
   │   │   └── ref.go
   │   ├── objectstore/
   │   │   ├── store.go
   │   │   └── store_test.go
   │   ├── repo/
   │   │   ├── repository.go
   │   │   ├── refs.go
   │   │   └── repository_test.go
   │   └── workspace/
   │       ├── scanner.go
   │       ├── diff.go
   │       └── scanner_test.go
   └── go.mod
   ```

3. **Set up CI/CD**
   - Create `.github/workflows/go.yml`
   - Configure: `go test ./...`, `go vet`, `gofmt -l`

---

## Milestone 1: Local VCS Core MVP (2-3 weeks)

### Step 1: Model Layer (3-4 days)

#### Objective
Implement core data structures with deterministic hashing.

#### Implementation Details

**1.1 Blob Model (`internal/model/blob.go`)**
```go
type Blob struct {
    data []byte
}

func (b *Blob) Hash() string {
    return sha256Hex(b.data)
}

func (b *Blob) Size() int64 {
    return int64(len(b.data))
}
```

**1.2 Object Model (`internal/model/object.go`)**
```go
type Reference struct {
    Kind string `json:"kind"` // "object" or "blob"
    ID   string `json:"id"`   // SHA-256 hash
}

type Object struct {
    Version int                    `json:"version"`
    Type    string                 `json:"type"`
    Meta    map[string]interface{} `json:"meta,omitempty"`
    Refs    []Reference            `json:"refs,omitempty"`
}

// Critical: canonical JSON serialisation for deterministic hashing
func (o *Object) Hash() string {
    // 1. Sort map keys lexicographically
    // 2. Use compact JSON (no whitespace)
    // 3. Return SHA-256 hex
}
```

**Key Requirements:**
- Map key ordering must be stable
- JSON output must be deterministic
- Handle edge cases: empty meta, nil refs, Unicode

**1.3 Snapshot Model (`internal/model/snapshot.go`)**
```go
type Snapshot struct {
    Version   int      `json:"version"`
    Root      string   `json:"root"`      // Hash of root object
    Parents   []string `json:"parents"`   // Parent snapshot hashes
    Author    string   `json:"author"`
    Timestamp string   `json:"timestamp"` // RFC3339 format
    Message   string   `json:"message"`
}

func (s *Snapshot) Hash() string {
    // Canonical JSON → SHA-256
}
```

**1.4 Ref Model (`internal/model/ref.go`)**
```go
type Ref struct {
    Name string // e.g., "refs/heads/main"
    Target string // Snapshot hash or symbolic ref
}

func (r *Ref) IsSymbolic() bool {
    return strings.HasPrefix(r.Target, "ref:")
}
```

#### Testing Strategy
```go
// model_test.go
func TestObjectHashDeterminism(t *testing.T) {
    // Same content, different construction order → same hash
}

func TestObjectRoundTrip(t *testing.T) {
    // Serialise → deserialise → should be equal
}

func TestCanonicalJSONOrdering(t *testing.T) {
    // Verify map keys are sorted
}
```

#### Success Criteria
- All model types serialise/deserialise correctly
- Hashes are deterministic across runs
- Unit test coverage > 90%

---

### Step 2: ObjectStore Layer (4-5 days)

#### Objective
Implement content-addressable storage with atomic writes and integrity verification.

#### Implementation Details

**2.1 Store Structure (`internal/objectstore/store.go`)**
```go
type Store struct {
    root string       // Path to .chemvcs/objects
    mu   sync.RWMutex // Protect concurrent writes
}

func NewStore(repoPath string) (*Store, error) {
    objectsPath := filepath.Join(repoPath, ".chemvcs", "objects")
    return &Store{root: objectsPath}, nil
}
```

**2.2 Blob Storage**
```go
func (s *Store) PutBlob(data []byte) (string, error) {
    hash := sha256Hex(data)
    path := s.blobPath(hash)
    
    // Check if already exists (idempotency)
    if exists, _ := s.hasBlob(hash); exists {
        return hash, nil
    }
    
    // Atomic write: temp file + rename
    return hash, s.atomicWrite(path, data)
}

func (s *Store) GetBlob(hash string) ([]byte, error) {
    data, err := ioutil.ReadFile(s.blobPath(hash))
    if err != nil {
        return nil, err
    }
    
    // Verify hash integrity
    if sha256Hex(data) != hash {
        return nil, ErrHashMismatch
    }
    
    return data, nil
}

// Sharded layout: objects/ab/cdef123456...
func (s *Store) blobPath(hash string) string {
    return filepath.Join(s.root, hash[:2], hash)
}
```

**2.3 Object Storage**
```go
func (s *Store) PutObject(obj *model.Object) (string, error) {
    jsonData, err := obj.MarshalCanonical()
    if err != nil {
        return "", err
    }
    
    hash := obj.Hash()
    path := s.objectPath(hash)
    
    // Idempotent write
    if exists, _ := s.hasObject(hash); exists {
        return hash, nil
    }
    
    return hash, s.atomicWrite(path, jsonData)
}

func (s *Store) GetObject(hash string) (*model.Object, error) {
    data, err := ioutil.ReadFile(s.objectPath(hash))
    if err != nil {
        return nil, err
    }
    
    obj := &model.Object{}
    if err := json.Unmarshal(data, obj); err != nil {
        return nil, err
    }
    
    // Verify hash
    if obj.Hash() != hash {
        return nil, ErrHashMismatch
    }
    
    return obj, nil
}
```

**2.4 Atomic Write Helper**
```go
func (s *Store) atomicWrite(path string, data []byte) error {
    dir := filepath.Dir(path)
    if err := os.MkdirAll(dir, 0755); err != nil {
        return err
    }
    
    // Create temp file in same directory
    tmp, err := ioutil.TempFile(dir, ".tmp-*")
    if err != nil {
        return err
    }
    tmpPath := tmp.Name()
    
    // Write data
    if _, err := tmp.Write(data); err != nil {
        tmp.Close()
        os.Remove(tmpPath)
        return err
    }
    tmp.Close()
    
    // Atomic rename
    return os.Rename(tmpPath, path)
}
```

#### Testing Strategy
```go
func TestPutGetBlob(t *testing.T) {
    // Write → read → verify content match
}

func TestBlobIntegrityCheck(t *testing.T) {
    // Corrupt file → GetBlob should fail
}

func TestConcurrentWrites(t *testing.T) {
    // Multiple goroutines writing same blob → should be safe
}

func TestIdempotency(t *testing.T) {
    // Writing same blob twice → no error, same hash
}
```

#### Success Criteria
- Atomic writes prevent corruption
- Concurrent writes are safe
- Hash verification catches corruption
- Sharded directory layout works correctly

---

### Step 3: Repository Layer (5-6 days)

#### Objective
Implement repository initialisation, snapshot management, and ref handling.

#### Implementation Details

**3.1 Repository Structure (`internal/repo/repository.go`)**
```go
type Repository struct {
    path  string
    store *objectstore.Store
    refs  *RefManager
}

func Init(path string) (*Repository, error) {
    repoPath := filepath.Join(path, ".chemvcs")
    
    // Create directory structure
    dirs := []string{
        repoPath,
        filepath.Join(repoPath, "objects"),
        filepath.Join(repoPath, "refs", "heads"),
    }
    
    for _, dir := range dirs {
        if err := os.MkdirAll(dir, 0755); err != nil {
            return nil, err
        }
    }
    
    // Create initial HEAD
    headPath := filepath.Join(repoPath, "HEAD")
    initialHEAD := "ref: refs/heads/main"
    if err := ioutil.WriteFile(headPath, []byte(initialHEAD), 0644); err != nil {
        return nil, err
    }
    
    return Open(path)
}

func Open(path string) (*Repository, error) {
    repoPath := findRepository(path)
    if repoPath == "" {
        return nil, ErrNotARepository
    }
    
    store, err := objectstore.NewStore(repoPath)
    if err != nil {
        return nil, err
    }
    
    return &Repository{
        path:  repoPath,
        store: store,
        refs:  NewRefManager(repoPath),
    }, nil
}

// Search upwards for .chemvcs directory
func findRepository(startPath string) string {
    current := startPath
    for {
        chemvcsPath := filepath.Join(current, ".chemvcs")
        if info, err := os.Stat(chemvcsPath); err == nil && info.IsDir() {
            return current
        }
        
        parent := filepath.Dir(current)
        if parent == current {
            return "" // Reached filesystem root
        }
        current = parent
    }
}
```

**3.2 Ref Manager (`internal/repo/refs.go`)**
```go
type RefManager struct {
    repoPath string
}

func (rm *RefManager) ResolveHEAD() (string, error) {
    headPath := filepath.Join(rm.repoPath, ".chemvcs", "HEAD")
    data, err := ioutil.ReadFile(headPath)
    if err != nil {
        return "", err
    }
    
    content := strings.TrimSpace(string(data))
    
    // Check if symbolic ref
    if strings.HasPrefix(content, "ref:") {
        refName := strings.TrimSpace(content[4:])
        return rm.ResolveRef(refName)
    }
    
    // Direct hash
    return content, nil
}

func (rm *RefManager) ResolveRef(refName string) (string, error) {
    refPath := filepath.Join(rm.repoPath, ".chemvcs", refName)
    data, err := ioutil.ReadFile(refPath)
    if err != nil {
        if os.IsNotExist(err) {
            return "", nil // Ref doesn't exist yet (e.g., initial commit)
        }
        return "", err
    }
    
    return strings.TrimSpace(string(data)), nil
}

func (rm *RefManager) UpdateRef(refName, snapshotHash string) error {
    refPath := filepath.Join(rm.repoPath, ".chemvcs", refName)
    
    // Ensure parent directory exists
    if err := os.MkdirAll(filepath.Dir(refPath), 0755); err != nil {
        return err
    }
    
    // Atomic write
    return atomicWrite(refPath, []byte(snapshotHash+"\n"))
}

func (rm *RefManager) CurrentBranch() (string, error) {
    headPath := filepath.Join(rm.repoPath, ".chemvcs", "HEAD")
    data, err := ioutil.ReadFile(headPath)
    if err != nil {
        return "", err
    }
    
    content := strings.TrimSpace(string(data))
    if strings.HasPrefix(content, "ref:") {
        return strings.TrimSpace(content[4:]), nil
    }
    
    return "", nil // Detached HEAD
}
```

**3.3 Snapshot Operations**
```go
func (r *Repository) CreateSnapshot(rootObjectHash, message, author string) (string, error) {
    // Get parent snapshot
    parentHash, err := r.refs.ResolveHEAD()
    if err != nil {
        return "", err
    }
    
    var parents []string
    if parentHash != "" {
        parents = []string{parentHash}
    }
    
    // Create snapshot
    snapshot := &model.Snapshot{
        Version:   1,
        Root:      rootObjectHash,
        Parents:   parents,
        Author:    author,
        Timestamp: time.Now().UTC().Format(time.RFC3339),
        Message:   message,
    }
    
    // Store snapshot (writes to objectstore)
    hash, err := r.store.PutSnapshot(snapshot)
    if err != nil {
        return "", err
    }
    
    // Update current branch ref
    // CRITICAL: Objects written first, refs updated last
    currentBranch, err := r.refs.CurrentBranch()
    if err != nil {
        return "", err
    }
    
    if currentBranch != "" {
        if err := r.refs.UpdateRef(currentBranch, hash); err != nil {
            return "", err
        }
    }
    
    return hash, nil
}

func (r *Repository) GetSnapshot(hash string) (*model.Snapshot, error) {
    return r.store.GetSnapshot(hash)
}

func (r *Repository) Log(limit int) ([]*model.Snapshot, error) {
    currentHash, err := r.refs.ResolveHEAD()
    if err != nil || currentHash == "" {
        return nil, err
    }
    
    snapshots := []*model.Snapshot{}
    
    for currentHash != "" && len(snapshots) < limit {
        snap, err := r.GetSnapshot(currentHash)
        if err != nil {
            return nil, err
        }
        
        snapshots = append(snapshots, snap)
        
        // Follow parent chain (MVP: linear history only)
        if len(snap.Parents) > 0 {
            currentHash = snap.Parents[0]
        } else {
            currentHash = ""
        }
    }
    
    return snapshots, nil
}
```

#### Testing Strategy
```go
func TestInitRepository(t *testing.T) {
    tmpDir := t.TempDir()
    repo, err := Init(tmpDir)
    assert.NoError(t, err)
    assert.DirExists(t, filepath.Join(tmpDir, ".chemvcs"))
}

func TestCreateSnapshot(t *testing.T) {
    // Create dummy root object
    // Call CreateSnapshot
    // Verify snapshot stored and ref updated
}

func TestLogTraversal(t *testing.T) {
    // Create chain of snapshots
    // Verify Log returns them in correct order
}
```

#### Success Criteria
- Repository can be initialised and opened
- Snapshots can be created and retrieved
- Refs are updated atomically
- Log traverses history correctly

---

### Step 4: CLI Implementation (3-4 days)

#### Objective
Provide user-facing commands with clear error messages.

#### Implementation Details

**4.1 Main Entry Point (`cmd/chemvcs/main.go`)**
```go
package main

import (
    "fmt"
    "os"
    
    "github.com/lishi/chemvcs/internal/repo"
)

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
    case "version":
        handleVersion()
    default:
        fmt.Fprintf(os.Stderr, "Unknown command: %s\n", cmd)
        printUsage()
        os.Exit(2)
    }
    
    if err != nil {
        fmt.Fprintf(os.Stderr, "Error: %v\n", err)
        os.Exit(1)
    }
}

func printUsage() {
    fmt.Println("Usage: chemvcs <command> [options]")
    fmt.Println()
    fmt.Println("Commands:")
    fmt.Println("  init              Initialise a new repository")
    fmt.Println("  commit -m <msg>   Create a new snapshot")
    fmt.Println("  log               Show snapshot history")
    fmt.Println("  version           Show version information")
}
```

**4.2 Init Command**
```go
func handleInit(args []string) error {
    var path string
    if len(args) > 0 {
        path = args[0]
    } else {
        path = "."
    }
    
    absPath, err := filepath.Abs(path)
    if err != nil {
        return err
    }
    
    _, err = repo.Init(absPath)
    if err != nil {
        return err
    }
    
    fmt.Printf("Initialised empty ChemVCS repository in %s\n", 
               filepath.Join(absPath, ".chemvcs"))
    return nil
}
```

**4.3 Commit Command (MVP: simplified)**
```go
func handleCommit(args []string) error {
    // Parse flags
    fs := flag.NewFlagSet("commit", flag.ExitOnError)
    message := fs.String("m", "", "Commit message")
    author := fs.String("author", "", "Author (defaults to config)")
    fs.Parse(args)
    
    if *message == "" {
        return fmt.Errorf("commit message required (-m)")
    }
    
    // Open repository
    r, err := repo.Open(".")
    if err != nil {
        return err
    }
    
    // Get author (from flag or config or environment)
    authorStr := getAuthor(*author)
    
    // MVP: Create a simple root object
    // (Full workspace scanning comes in Milestone 2)
    rootHash, err := createSimpleRoot(r)
    if err != nil {
        return err
    }
    
    // Create snapshot
    snapHash, err := r.CreateSnapshot(rootHash, *message, authorStr)
    if err != nil {
        return err
    }
    
    fmt.Printf("[%s] %s\n", snapHash[:8], *message)
    return nil
}

// MVP: Creates a minimal root object for testing
func createSimpleRoot(r *repo.Repository) (string, error) {
    obj := &model.Object{
        Version: 1,
        Type:    "generic",
        Meta:    map[string]interface{}{"mvp": true},
    }
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
    
    // Default
    return "Unknown <unknown@localhost>"
}
```

**4.4 Log Command**
```go
func handleLog(args []string) error {
    fs := flag.NewFlagSet("log", flag.ExitOnError)
    limit := fs.Int("n", 20, "Number of snapshots to show")
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
        fmt.Println("No snapshots yet")
        return nil
    }
    
    for _, snap := range snapshots {
        fmt.Printf("snapshot %s\n", snap.Hash()[:8])
        fmt.Printf("Author: %s\n", snap.Author)
        fmt.Printf("Date:   %s\n", snap.Timestamp)
        fmt.Printf("\n    %s\n\n", snap.Message)
    }
    
    return nil
}
```

#### Testing Strategy
```go
func TestCLI_InitCommitLog(t *testing.T) {
    tmpDir := t.TempDir()
    
    // Run init
    cmd := exec.Command("chemvcs", "init", tmpDir)
    output, err := cmd.CombinedOutput()
    assert.NoError(t, err)
    assert.Contains(t, string(output), "Initialised")
    
    // Run commit
    cmd = exec.Command("chemvcs", "commit", "-m", "test commit")
    cmd.Dir = tmpDir
    output, err = cmd.CombinedOutput()
    assert.NoError(t, err)
    
    // Run log
    cmd = exec.Command("chemvcs", "log")
    cmd.Dir = tmpDir
    output, err = cmd.CombinedOutput()
    assert.NoError(t, err)
    assert.Contains(t, string(output), "test commit")
}
```

#### Success Criteria
- Commands execute without panics
- Error messages are clear and actionable
- Exit codes follow conventions
- Help text is informative

---

## Milestone 2: Working Directory and Status (2-3 weeks)

### Step 5: Workspace Scanner (5-6 days)

#### Objective
Implement bidirectional mapping between filesystem and object graph.

#### Implementation Details

**5.1 Scanner (`internal/workspace/scanner.go`)**
```go
type Scanner struct {
    store *objectstore.Store
}

func (s *Scanner) ScanDirectory(root string) (*model.Object, error) {
    return s.scanDir(root, "")
}

func (s *Scanner) scanDir(basePath, relPath string) (*model.Object, error) {
    fullPath := filepath.Join(basePath, relPath)
    entries, err := ioutil.ReadDir(fullPath)
    if err != nil {
        return nil, err
    }
    
    folderObj := &model.Object{
        Version: 1,
        Type:    "folder",
        Meta:    map[string]interface{}{
            "name": filepath.Base(relPath),
        },
        Refs: []model.Reference{},
    }
    
    for _, entry := range entries {
        // Skip .chemvcs directory
        if entry.Name() == ".chemvcs" {
            continue
        }
        
        entryPath := filepath.Join(relPath, entry.Name())
        
        if entry.IsDir() {
            // Recurse into subdirectory
            childObj, err := s.scanDir(basePath, entryPath)
            if err != nil {
                return nil, err
            }
            childHash, err := s.store.PutObject(childObj)
            if err != nil {
                return nil, err
            }
            folderObj.Refs = append(folderObj.Refs, model.Reference{
                Kind: "object",
                ID:   childHash,
            })
        } else {
            // Create file object
            fileObj, err := s.scanFile(filepath.Join(basePath, entryPath), entry)
            if err != nil {
                return nil, err
            }
            fileHash, err := s.store.PutObject(fileObj)
            if err != nil {
                return nil, err
            }
            folderObj.Refs = append(folderObj.Refs, model.Reference{
                Kind: "object",
                ID:   fileHash,
            })
        }
    }
    
    return folderObj, nil
}

func (s *Scanner) scanFile(path string, info os.FileInfo) (*model.Object, error) {
    // Read file content
    data, err := ioutil.ReadFile(path)
    if err != nil {
        return nil, err
    }
    
    // Store blob
    blobHash, err := s.store.PutBlob(data)
    if err != nil {
        return nil, err
    }
    
    // Create file object
    fileObj := &model.Object{
        Version: 1,
        Type:    "file",
        Meta: map[string]interface{}{
            "name": info.Name(),
            "size": info.Size(),
            "mode": fmt.Sprintf("%o", info.Mode()),
        },
        Refs: []model.Reference{
            {Kind: "blob", ID: blobHash},
        },
    }
    
    return fileObj, nil
}
```

**5.2 Restorer**
```go
func (s *Scanner) RestoreDirectory(rootObjHash, targetPath string) error {
    rootObj, err := s.store.GetObject(rootObjHash)
    if err != nil {
        return err
    }
    
    return s.restoreObject(rootObj, targetPath)
}

func (s *Scanner) restoreObject(obj *model.Object, path string) error {
    switch obj.Type {
    case "folder":
        return s.restoreFolder(obj, path)
    case "file":
        return s.restoreFile(obj, path)
    default:
        return fmt.Errorf("unknown object type: %s", obj.Type)
    }
}

func (s *Scanner) restoreFolder(obj *model.Object, path string) error {
    if err := os.MkdirAll(path, 0755); err != nil {
        return err
    }
    
    for _, ref := range obj.Refs {
        if ref.Kind != "object" {
            continue
        }
        
        childObj, err := s.store.GetObject(ref.ID)
        if err != nil {
            return err
        }
        
        childName := childObj.Meta["name"].(string)
        childPath := filepath.Join(path, childName)
        
        if err := s.restoreObject(childObj, childPath); err != nil {
            return err
        }
    }
    
    return nil
}

func (s *Scanner) restoreFile(obj *model.Object, path string) error {
    if len(obj.Refs) == 0 {
        return fmt.Errorf("file object has no blob reference")
    }
    
    blobHash := obj.Refs[0].ID
    data, err := s.store.GetBlob(blobHash)
    if err != nil {
        return err
    }
    
    // Parse mode if present
    mode := os.FileMode(0644)
    if modeStr, ok := obj.Meta["mode"].(string); ok {
        if modeInt, err := strconv.ParseUint(modeStr, 8, 32); err == nil {
            mode = os.FileMode(modeInt)
        }
    }
    
    return ioutil.WriteFile(path, data, mode)
}
```

### Step 6: Status Command (3-4 days)

**6.1 Diff Engine (`internal/workspace/diff.go`)**
```go
type Change struct {
    Type string // "added", "modified", "deleted"
    Path string
    OldHash string
    NewHash string
}

type Changes struct {
    Added    []string
    Modified []string
    Deleted  []string
}

func ComputeChanges(oldRootHash, newRootHash string, store *objectstore.Store) (*Changes, error) {
    changes := &Changes{}
    
    oldObj, err := store.GetObject(oldRootHash)
    if err != nil {
        return nil, err
    }
    
    newObj, err := store.GetObject(newRootHash)
    if err != nil {
        return nil, err
    }
    
    // Recursive comparison
    compareObjects(oldObj, newObj, "", changes, store)
    
    return changes, nil
}
```

**6.2 CLI Integration**
```go
func handleStatus(args []string) error {
    r, err := repo.Open(".")
    if err != nil {
        return err
    }
    
    // Get current snapshot
    currentHash, err := r.Refs().ResolveHEAD()
    if err != nil || currentHash == "" {
        fmt.Println("No commits yet")
        return nil
    }
    
    currentSnap, err := r.GetSnapshot(currentHash)
    if err != nil {
        return err
    }
    
    // Scan working directory
    scanner := workspace.NewScanner(r.Store())
    workingRoot, err := scanner.ScanDirectory(".")
    if err != nil {
        return err
    }
    
    // Compute differences
    changes, err := workspace.ComputeChanges(
        currentSnap.Root,
        workingRoot.Hash(),
        r.Store(),
    )
    if err != nil {
        return err
    }
    
    // Print changes
    if len(changes.Added) == 0 && len(changes.Modified) == 0 && len(changes.Deleted) == 0 {
        fmt.Println("No changes")
        return nil
    }
    
    for _, path := range changes.Added {
        fmt.Printf("  A  %s\n", path)
    }
    for _, path := range changes.Modified {
        fmt.Printf("  M  %s\n", path)
    }
    for _, path := range changes.Deleted {
        fmt.Printf("  D  %s\n", path)
    }
    
    return nil
}
```

### Step 7: Checkout Command (2-3 days)

**7.1 CLI Integration**
```go
func handleCheckout(args []string) error {
    if len(args) < 1 {
        return fmt.Errorf("usage: chemvcs checkout <branch|snapshot>")
    }
    
    target := args[0]
    
    r, err := repo.Open(".")
    if err != nil {
        return err
    }
    
    // Check for uncommitted changes
    // ... (safety check)
    
    // Resolve target
    var targetHash string
    if strings.HasPrefix(target, "refs/") {
        targetHash, err = r.Refs().ResolveRef(target)
    } else {
        // Try as branch name
        targetHash, err = r.Refs().ResolveRef("refs/heads/" + target)
        if err != nil || targetHash == "" {
            // Try as direct snapshot hash
            targetHash = target
        }
    }
    
    // Get snapshot
    snap, err := r.GetSnapshot(targetHash)
    if err != nil {
        return err
    }
    
    // Restore working directory
    scanner := workspace.NewScanner(r.Store())
    if err := scanner.RestoreDirectory(snap.Root, "."); err != nil {
        return err
    }
    
    fmt.Printf("Checked out snapshot %s\n", targetHash[:8])
    return nil
}
```

---

## Milestone 3: Branches and Merge (1-2 weeks)

### Step 8: Branch Management (4-5 days)

**8.1 Branch Command**
```go
func handleBranch(args []string) error {
    r, err := repo.Open(".")
    if err != nil {
        return err
    }
    
    if len(args) == 0 {
        // List branches
        return listBranches(r)
    }
    
    // Create new branch
    branchName := args[0]
    currentHash, err := r.Refs().ResolveHEAD()
    if err != nil || currentHash == "" {
        return fmt.Errorf("cannot create branch: no commits yet")
    }
    
    refName := "refs/heads/" + branchName
    if err := r.Refs().UpdateRef(refName, currentHash); err != nil {
        return err
    }
    
    fmt.Printf("Created branch '%s' at %s\n", branchName, currentHash[:8])
    return nil
}

func listBranches(r *repo.Repository) error {
    headsPath := filepath.Join(r.Path(), ".chemvcs", "refs", "heads")
    entries, err := ioutil.ReadDir(headsPath)
    if err != nil {
        return err
    }
    
    currentBranch, _ := r.Refs().CurrentBranch()
    
    for _, entry := range entries {
        if entry.IsDir() {
            continue
        }
        
        marker := " "
        branchName := "refs/heads/" + entry.Name()
        if branchName == currentBranch {
            marker = "*"
        }
        
        fmt.Printf("%s %s\n", marker, entry.Name())
    }
    
    return nil
}
```

### Step 9: Fast-Forward Merge (3-4 days)

**9.1 Merge Logic**
```go
func (r *Repository) Merge(branchName string) error {
    // Resolve current HEAD
    currentHash, err := r.refs.ResolveHEAD()
    if err != nil || currentHash == "" {
        return fmt.Errorf("cannot merge: no current snapshot")
    }
    
    // Resolve target branch
    targetRef := "refs/heads/" + branchName
    targetHash, err := r.refs.ResolveRef(targetRef)
    if err != nil || targetHash == "" {
        return fmt.Errorf("branch '%s' not found", branchName)
    }
    
    // Check if already up to date
    if currentHash == targetHash {
        return fmt.Errorf("already up to date")
    }
    
    // Check if fast-forward is possible
    if !r.isAncestor(targetHash, currentHash) {
        return fmt.Errorf("cannot fast-forward merge")
    }
    
    // Update current branch
    currentBranch, err := r.refs.CurrentBranch()
    if err != nil {
        return err
    }
    
    if err := r.refs.UpdateRef(currentBranch, targetHash); err != nil {
        return err
    }
    
    fmt.Printf("Fast-forwarded to %s\n", targetHash[:8])
    return nil
}

func (r *Repository) isAncestor(ancestor, descendant string) bool {
    current := descendant
    for current != "" {
        if current == ancestor {
            return true
        }
        
        snap, err := r.GetSnapshot(current)
        if err != nil || len(snap.Parents) == 0 {
            return false
        }
        
        current = snap.Parents[0]
    }
    return false
}
```

---

## Testing and Quality Assurance

### Unit Test Coverage Goals
- `model`: > 95%
- `objectstore`: > 90%
- `repo`: > 85%
- `workspace`: > 80%

### Integration Test Scenarios
1. Full workflow: init → create files → commit → modify → status → commit → log
2. Branch workflow: create branch → switch → commit → merge
3. Corruption detection: modify .chemvcs files → verify errors
4. Concurrent access: multiple processes accessing same repo

### Performance Benchmarks
- Commit 1000 small files: < 5 seconds
- Log 1000 snapshots: < 1 second
- Status on 1000 files: < 3 seconds

---

## Timeline Summary

| Phase | Duration | Key Deliverables |
|-------|----------|------------------|
| Phase 0 | 1 day | Project structure, CI setup |
| M1 Step 1-2 | 1 week | Model + ObjectStore with tests |
| M1 Step 3-4 | 1.5 weeks | Repository + CLI MVP |
| M2 Step 5-7 | 2 weeks | Workspace scanning, status, checkout |
| M3 Step 8-9 | 1 week | Branch management, merge |
| **Total** | **6-7 weeks** | **Fully functional local VCS** |

---

## Risk Mitigation

### Risk 1: Hash Non-Determinism
**Mitigation:** Extensive testing of JSON canonicalisation. Create test suite with various data types, edge cases, and Unicode.

### Risk 2: Windows Path Handling
**Mitigation:** Use `filepath` package exclusively. Internally normalise to forward slashes. Test on Windows CI.

### Risk 3: Large File Performance
**Mitigation:** Add file size warnings in commit. Document limitations. Plan streaming for future.

### Risk 4: Concurrent Write Conflicts
**Mitigation:** Use atomic file operations. Consider adding lockfiles if needed.

---

## Success Criteria for Each Milestone

### Milestone 1
- [ ] Can initialise a repository
- [ ] Can create snapshots with message and author
- [ ] Can view history with `log`
- [ ] All objects are content-addressable and verifiable
- [ ] Unit tests pass with >85% coverage

### Milestone 2
- [ ] Can detect changes in working directory
- [ ] Can restore working directory from any snapshot
- [ ] Status accurately reports added/modified/deleted files
- [ ] Integration tests pass for full workflows

### Milestone 3
- [ ] Can create and list branches
- [ ] Can switch between branches
- [ ] Can fast-forward merge branches
- [ ] Non-fast-forward merges are rejected gracefully

---

## Next Steps After Milestone 3

1. **Milestone 4:** Remote protocol and server
2. **Milestone 5:** Python domain layer
3. **Milestone 6:** HPC scheduler integration

Each subsequent milestone builds on the solid foundation of Milestones 1-3.
