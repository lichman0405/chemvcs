# ChemVCS vs Git: Detailed Comparison (v0.1)

## 1. Introduction

This document provides a comprehensive comparison between ChemVCS and Git. While ChemVCS draws inspiration from Git's elegant design, it makes deliberate departures to better serve computational chemistry and scientific computing workflows.

Understanding these differences helps clarify:
- Why ChemVCS exists as a separate system
- When to use ChemVCS vs Git
- What trade-offs were made and why

---

## 2. Shared Foundations

Both systems share these core concepts:

### 2.1 Content-Addressable Storage
- Objects identified by SHA-256 hash of their content
- Immutable objects (write once, reference forever)
- Deduplication through hash-based storage

### 2.2 Merkle DAG Structure
- Snapshots/commits form a directed acyclic graph
- Each snapshot references parent(s)
- History is cryptographically verifiable

### 2.3 Distributed Architecture
- Clone entire repository with full history
- Push, pull, fetch operations for synchronisation
- Work offline, sync later

### 2.4 Branch-Based Workflow
- Create lightweight branches
- Switch between branches
- Merge branches (fast-forward or three-way)

---

## 3. Critical Architectural Differences

### 3.1 Object Model (⭐ Most Important Difference)

**Git: Fixed Four-Type System**
```
blob    → Raw file content (opaque bytes)
tree    → Directory listing (name → hash mapping)
commit  → Snapshot with parents, author, message
tag     → Named reference to a commit
```

**ChemVCS: Extensible Object System**
```go
type Object struct {
    Version int                    // Schema version
    Type    string                 // Extensible: "file", "folder", or "structure", "calculation"
    Meta    map[string]interface{} // Arbitrary domain metadata
    Refs    []Reference            // Typed references (object/blob)
}
```

**Implications:**
1. **Git**: To add chemistry semantics, you must layer external metadata systems on top
2. **ChemVCS**: Chemistry semantics can live **inside** the VCS as first-class objects

**Example Future Usage (M5+):**
```go
// A molecular structure as a VCS object
{
    Type: "structure",
    Meta: {
        "formula": "C6H6",
        "atoms": 12,
        "charge": 0,
        "multiplicity": 1
    },
    Refs: [
        {Kind: "blob", ID: "abc123..."}  // XYZ coordinates
    ]
}

// A calculation as a VCS object
{
    Type: "calculation",
    Meta: {
        "method": "B3LYP/6-31G*",
        "energy": -232.123456,
        "converged": true,
        "software": "Gaussian16"
    },
    Refs: [
        {Kind: "object", ID: "def456..."},  // Input structure
        {Kind: "blob", ID: "ghi789..."}     // Output file
    ]
}
```

This is **architecturally impossible** in Git without external databases.

---

### 3.2 Staging Area / Index

**Git: Three-Stage Workflow**
```
Working Directory → (git add) → Staging Area → (git commit) → Repository
```
- Fine-grained control over what goes into each commit
- Can stage parts of files (hunks)
- Complex but powerful

**ChemVCS: Two-Stage Workflow**
```
Working Directory → (chemvcs commit) → Repository
```
- No staging area
- Commit captures entire working directory state
- Simpler mental model

**Rationale for ChemVCS Design:**
1. **Scientific calculations are atomic**: A calculation either completed or it didn't. Partial staging doesn't map to scientific practice.
2. **Snapshot semantics**: Each commit should represent a complete, reproducible computational state.
3. **Reduced complexity**: One fewer concept to learn and maintain.

**Trade-off:**
- ❌ Less flexibility for crafting "clean" commits
- ✅ Simpler workflow aligned with scientific thinking
- ✅ Impossible to commit incomplete states

---

### 3.3 History Integrity Constraints

**Git Philosophy: Flexibility**
- `git rebase` - Rewrite commit history for cleaner lineage
- `git cherry-pick` - Copy commits across branches
- `git commit --amend` - Modify the last commit
- `git reset --hard` - Discard commits and changes
- **Use case**: Clean up messy development history before sharing

**ChemVCS Philosophy: Authenticity**
- **Rebase**: ❌ Not implemented, not planned
- **Cherry-pick**: ❌ Not implemented, not planned
- **Amend**: ❌ Not supported
- **Reset**: ❌ No staging area, so most modes don't apply
- **Use case**: Preserve authentic computational history

**Rationale:**
1. **Scientific reproducibility**: The order of calculations matters. Failures teach us as much as successes.
2. **Provenance tracking**: "I tried method A, it failed, then tried method B" is valuable context.
3. **Regulatory compliance**: Some scientific domains require unmodified audit trails.
4. **Simplicity**: Fewer dangerous commands to learn and avoid misusing.

**Trade-off:**
- ❌ Can't beautify history for presentations
- ✅ History is always truthful and complete
- ✅ No risk of accidentally destroying work

---

### 3.4 Terminology and Conceptual Framing

| Concept | Git | ChemVCS | Reason for Difference |
|---------|-----|---------|----------------------|
| Version point | Commit | Snapshot | Emphasizes capturing complete state |
| Directory tree | Tree object | Folder object | Unified object model |
| File content | Blob | Blob | Same (no difference needed) |
| Working changes | Staged/unstaged | Modified | No staging area |

**Philosophy:**
- Git terminology reflects source code development
- ChemVCS terminology reflects scientific data capture

---

## 4. Feature Comparison Table

| Feature | Git | ChemVCS | Notes |
|---------|-----|---------|-------|
| Content-addressable storage | ✅ | ✅ | Core shared concept |
| Merkle DAG | ✅ | ✅ | Core shared concept |
| Branching | ✅ | ✅ | Both support |
| Fast-forward merge | ✅ | ✅ | Both support |
| Three-way merge | ✅ | ✅ | Both support |
| Conflict detection | ✅ | ✅ | Both support |
| Remote sync (push/pull/fetch) | ✅ | ✅ | Both support |
| Staging area | ✅ | ❌ | Deliberately removed |
| Rebase | ✅ | ❌ | Violates provenance |
| Cherry-pick | ✅ | ❌ | Violates provenance |
| Commit --amend | ✅ | ❌ | Violates provenance |
| Interactive rebase | ✅ | ❌ | Violates provenance |
| Reset --hard | ✅ | ⚠️ | Limited (no staging) |
| Stash | ✅ | 📋 | Planned (low priority) |
| Submodules | ✅ | 📋 | Planned |
| Hooks | ✅ | 📋 | Planned (high priority) |
| Reflog | ✅ | 📋 | Planned (low priority) |
| Bisect | ✅ | 📋 | Planned (low priority) |
| **Domain object types** | ❌ | ✅ | Python layer: Structure, Run, Workflow (stored as typed objects) |
| **Extensible metadata** | ⚠️ | ✅ | Notes vs native Meta |
| **Chemistry file diff** | ❌ | 📋 | M5: XYZ, CIF awareness |
| **HPC job tracking** | ❌ | ✅ | SLURM/PBS/LSF adapters; local + remote gateway workflows |
| **Provenance capture** | ❌ | ✅ | Scheduler directives + environment/modules captured on submit |
| **Computational queries** | ❌ | 📋 | M5: search by properties |

Legend:
- ✅ Fully supported
- ⚠️ Partial support or different approach
- ❌ Not supported / deliberately excluded
- 📋 Planned in future milestones

---

## 5. Use Case Comparison

### 5.1 When to Use Git

**Ideal for:**
- Source code development
- Collaborative software projects
- Fine-grained commit control needed
- History cleanup/beautification desired
- Mature ecosystem of tools required

**Example:**
```bash
# Typical software development
git checkout -b feature
# ... make changes ...
git add src/module.py          # Stage specific file
git commit -m "Add feature X"
git rebase main                # Clean up history
git push origin feature
```

### 5.2 When to Use ChemVCS (Current State)

**Ideal for:**
- Tracking computational chemistry projects
- Capturing complete calculation states
- Simpler workflow (no staging)
- Authentic history required

**Example:**
```bash
# Typical computational workflow
chemvcs init
# ... run calculation, generate outputs ...
chemvcs commit -m "B3LYP/6-31G* optimization of benzene"
chemvcs push origin main
```

### 5.3 When to Use ChemVCS (Future M5+)

**Will be ideal for:**
- Domain-aware provenance tracking
- Querying calculations by properties
- HPC job correlation
- Chemistry file format intelligence

**Example (Future):**
```python
# M5+ Python API
import chemvcs

repo = chemvcs.Repository.open('.')

# Commit with chemistry semantics
repo.commit_calculation(
    structure=mol,
    method='B3LYP/6-31G*',
    energy=-232.123,
    output_file='benzene.out'
)

# Query by properties
results = repo.find_calculations(
    method='B3LYP',
    energy_range=(-240, -230),
    converged=True
)

# Chemistry-aware diff
diff = repo.diff_structure('HEAD', 'HEAD~1')
print(f"RMSD: {diff.rmsd} Å")
```

---

## 6. Migration and Interoperability

### 6.1 Git → ChemVCS

**Not a simple import** because:
- Git commits don't map cleanly to ChemVCS snapshots (no staging semantics to preserve)
- Git metadata (committer vs author, GPG signatures) not fully modeled
- History rewriting in Git contradicts ChemVCS philosophy

**Possible approach:**
- Export Git commits as ChemVCS snapshots (one-way conversion)
- Preserve timestamps and messages
- Lose staging information and merge details

### 6.2 ChemVCS → Git

**Even harder** because:
- ChemVCS extensible objects don't map to Git's fixed types
- Domain metadata would be lost or externalized

### 6.3 Hybrid Workflows

**Recommended:**
- Use Git for source code (scripts, analysis tools)
- Use ChemVCS for computational data and results
- Keep them as separate repositories

**Example Project Structure:**
```
research-project/
├── code/               # Git repository (Python scripts, analysis)
│   └── .git/
└── data/               # ChemVCS repository (calculations, structures)
    └── .chemvcs/
```

---

## 7. Performance Characteristics

### 7.1 Storage

**Git:**
- Packfiles for efficient storage
- Delta compression for similar objects
- Garbage collection to remove unreachable objects

**ChemVCS (Current):**
- Individual objects in sharded directories
- Optional packfiles for efficient storage (`chemvcs pack`)
- Garbage collection (`chemvcs gc`) and integrity verification (`chemvcs fsck`)

**Implication:** ChemVCS now has the core maintenance primitives, but still does not (yet) do Git-style delta compression.

### 7.2 Network Transfer

**Git:**
- Smart protocol with server-side delta computation
- Thin packs for efficient transfer

**ChemVCS:**
- Simple HTTP protocol
- Full object transfer (no deltas yet)

**Implication:** ChemVCS remote operations less efficient for large histories, but adequate for scientific data scales.

---

## 8. Ecosystem and Tooling

**Git Advantages:**
- Decades of tools (GitHub, GitLab, Gitea, etc.)
- IDE integrations (VS Code, IntelliJ, etc.)
- CI/CD pipelines designed around Git
- Massive community and documentation

**ChemVCS Current State:**
- Basic CLI and HTTP server
- No GUI yet
- No hosting platforms
- Small/experimental community

**ChemVCS Future Vision:**
- Expand the Python API for workflow integration (M5)
- Extend HPC scheduler integration (M6)
- Chemistry-aware visualisation tools
- Domain-specific hosting platforms

---

## 9. Design Philosophy Summary

### Git Philosophy
> A distributed version control system for coordinating work among programmers. Emphasizes flexibility, powerful tools, and a rich command set for managing source code history.

**Values:**
1. Flexibility (rebase, cherry-pick, etc.)
2. Efficiency (packfiles, smart protocol)
3. Universality (works for any project)

### ChemVCS Philosophy
> A domain-native version control system for computational chemistry. Emphasizes scientific integrity, provenance, and extensibility for domain-specific semantics.

**Values:**
1. **Scientific integrity** (authentic history)
2. **Domain semantics** (chemistry as first-class)
3. **Simplicity** (fewer concepts, clearer model)
4. **Provenance** (reproducibility first)

---

## 10. Conclusion

### Current Reality (M1-M4 Complete)

ChemVCS today is best described as:
> "A simplified Git-like core with real HPC tracking and a domain layer, built to make scientific provenance first-class"

The differentiation is both **architectural** and **partly functional**:
- ✅ Extensible object model (core)
- ✅ Simplified workflow (no staging)
- ✅ History integrity (no rewrites)
- ✅ Repository maintenance (pack/gc/fsck)
- ✅ HPC job tracking and provenance capture (Python layer + CLI integration)
- ⚠️ Chemistry-aware diff/merge and rich domain queries are still early / incomplete

### Future Vision (M5-M6)

ChemVCS will become:
> "A chemistry-native version control system with computational provenance built in"

True differentiation emerges when:
- 🔬 Chemistry objects are first-class VCS entities
- 🖥️ HPC workflows are natively integrated
- 📊 Queries can search by scientific properties
- 🔗 Provenance is automatic, not bolted on

### Bottom Line

**Use Git if:** You're managing source code and want mature tooling

**Use ChemVCS if:** You're managing computational chemistry data and want:
- Authentic computational history
- Simpler workflow aligned with scientific practice
- Architecture ready for domain-specific extensions (once M5/M6 land)

**The real test:** ChemVCS must deliver on M5/M6 to justify its existence as a separate system. Until then, it's an interesting experiment with a promising architecture.

---

## References

- Git Official Documentation: https://git-scm.com/doc
- ChemVCS Vision and Scope: [01-vision-and-scope.md](01-vision-and-scope.md)
- ChemVCS Architecture: [02-architecture-overview.md](02-architecture-overview.md)
- ChemVCS Feature Status / Roadmap Notes: [FEATURE_STATUS.md](../FEATURE_STATUS.md)
