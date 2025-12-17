# ChemVCS – Object Model and Storage Specification (v0.1)

## 1. Introduction

This document defines the core object model and on-disk storage format for ChemVCS.
It specifies:

- The types of core entities (Blob, Object, Snapshot, Ref);
- How these entities are serialised and addressed;
- The layout of the `.chemvcs/` directory in a repository;
- Integrity constraints and invariants that must always hold;
- Examples of JSON encodings for objects and snapshots.

The goal is to provide a stable foundation for the Go core, the Python domain layer,
and any external tools that wish to read or write ChemVCS repositories.

This specification focuses on the **MVP format** for Milestones 1–3. Future versions
may introduce additional fields or optimisations. Packfiles were later introduced
as an optional storage optimisation (P3); loose objects remain supported, and the
core on-disk invariants described here still apply.

---

## 2. Object Model Overview

ChemVCS uses a small set of core entity types:

1. **Blob** – raw binary data (e.g. a file’s contents).  
2. **Object** – typed node with metadata and references to other blobs/objects.  
3. **Snapshot** – a versioned “commit” representing a project state at a point in time.  
4. **Ref** – a named pointer to a snapshot.

All of these are addressed by a **content hash** (SHA-256, hex-encoded), making the
storage **content-addressable**. The logical object model forms a **Merkle DAG**;
snapshots reference a root object, which in turn references other objects and blobs.

The core model is intentionally domain-agnostic. Computational chemistry semantics
(e.g. “structure”, “run”, “workflow”) are expressed as specific `Object.type` values
and `Object.meta` contents defined in the domain layer, not in the core.

---

## 3. Identifiers and Hashing

### 3.1 Hash Function

- All core entities are addressed by a **SHA-256 hash** of their serialised content.
- Hashes are represented as lower-case hexadecimal strings, length 64 characters.

Example:

```text
a3f4b1c9d7e8... (64 hex characters total)
```

### 3.2 Hash Input Rules

The exact bytes being hashed depend on the entity type:

1. **Blob**  
   - Hash input: the raw binary data of the blob.  
   - No prefix or header is added in the MVP (future versions may define an envelope).

2. **Object**  
   - Hash input: the UTF-8 bytes of the object’s JSON representation, in a canonical
     form (see §4.2).

3. **Snapshot**  
   - Hash input: the UTF-8 bytes of the snapshot’s JSON representation, in a canonical
     form (see §5.2).

Canonical JSON encoding is used to ensure deterministic hashes for logically equal
objects and snapshots.

### 3.3 Hash Collisions

ChemVCS assumes SHA-256 collisions are practically infeasible. If a collision is
detected (two different payloads produce the same hash), the implementation must
refuse to overwrite existing content and report an error.

Collision handling beyond error reporting is out of scope for the MVP.

---

## 4. Object: Structure and Encoding

### 4.1 Object Semantics

An **Object** is a typed, immutable node with:

- A **type** string (e.g. `generic`, `folder`, `structure`, `run`, `workflow`);
- A **metadata** map (`meta`) containing key-value pairs;
- A list of **references** (`refs`) to other core entities (typically other Objects
  or Blobs).

This design allows domain layers to express arbitrary graphs of entities while
keeping the core type system simple.

### 4.2 JSON Schema (MVP)

In the MVP, Objects are stored on disk as JSON documents with the following schema
(before hashing/compression):

```json
{
  "version": 1,
  "type": "generic",
  "meta": { /* arbitrary JSON object, optional but recommended */ },
  "refs": [
    {
      "kind": "object",
      "id": "..." 
    },
    {
      "kind": "blob",
      "id": "..."
    }
  ]
}
```

Fields:

- `version` (integer, required):
  - Version of the Object encoding schema. For the MVP, this must be `1`.
  - Reserved for future evolution of the on-disk format.

- `type` (string, required):
  - Logical type identifier, e.g.:
    - `generic` – untyped or generic node;
    - `folder` – used by the generic workspace mapping for directory trees;
    - `file` – used by the generic workspace mapping for files;
    - `structure`, `run`, `workflow` – domain-specific types defined by the Python layer.

- `meta` (object, optional, default `{}`):
  - Arbitrary JSON object that stores attributes for this object.
  - Domain layers are responsible for defining and documenting their own metadata keys.

- `refs` (array of reference descriptors, optional, default `[]`):
  - Each element has the shape:
    - `kind` (string): `"object"` or `"blob"`;
    - `id` (string): SHA-256 hash (hex) of the referenced entity.

#### Canonicalisation Rules

To ensure deterministic hashes:

- Keys in `meta` and in each `refs` object should be serialised in a stable order
  (e.g. lexicographic) by the implementation;
- The top-level keys in the Object JSON (`version`, `type`, `meta`, `refs`) are
  serialised in a fixed order;
- No insignificant whitespace (pretty-printing) is included when computing the hash.

Exact canonicalisation details are an implementation concern, but the effect must
be: equal logical Objects → equal hashes.

### 4.3 Example: Folder Object (Generic Workspace Mapping)

```json
{
  "version": 1,
  "type": "folder",
  "meta": {
    "name": "src"
  },
  "refs": [
    {
      "kind": "object",
      "id": "8f17a1c3..."  // child file object
    },
    {
      "kind": "object",
      "id": "c2b4d9e0..."  // child folder object
    }
  ]
}
```

### 4.4 Example: File Object (Generic Workspace Mapping)

```json
{
  "version": 1,
  "type": "file",
  "meta": {
    "name": "main.py",
    "mode": "0644"
  },
  "refs": [
    {
      "kind": "blob",
      "id": "a9d4f8..."  // blob containing file contents
    }
  ]
}
```

### 4.5 Example: Domain Object (Run)

```json
{
  "version": 1,
  "type": "run",
  "meta": {
    "code": "vasp",
    "code_version": "6.4.1",
    "parameters": {
      "encut": 520,
      "kpoints": "4x4x4",
      "functional": "PBE"
    },
    "resources": {
      "nnodes": 2,
      "ntasks_per_node": 16,
      "walltime": "02:00:00"
    },
    "status": "finished"
  },
  "refs": [
    {
      "kind": "object",
      "id": "b1c2d3..."    // structure object
    },
    {
      "kind": "blob",
      "id": "c4d5e6..."    // INCAR blob
    },
    {
      "kind": "blob",
      "id": "d7e8f9..."    // POSCAR blob
    },
    {
      "kind": "blob",
      "id": "abcd12..."    // OUTCAR blob
    }
  ]
}
```

---

## 5. Snapshot: Structure and Encoding

### 5.1 Snapshot Semantics

A **Snapshot** is an immutable record of the repository state at a point in time,
similar to a git commit. It references:

- A root Object representing the top-level project state (e.g. a folder tree);
- Zero or more parent snapshots (for history and branching);
- Metadata about the commit (author, timestamp, message).

### 5.2 JSON Schema (MVP)

Snapshots are stored as JSON documents with the following schema:

```json
{
  "version": 1,
  "root": "...",
  "parents": ["...", "..."],
  "author": "Name <email@example.com>",
  "timestamp": "2025-01-23T10:15:30Z",
  "message": "Commit message"
}
```

Fields:

- `version` (integer, required):
  - Schema version for snapshot encoding. MVP uses `1`.

- `root` (string, required):
  - SHA-256 hash (hex) of the root Object for this snapshot.

- `parents` (array of strings, required, possibly empty):
  - List of parent snapshot hashes (hex). For:
    - Initial commits: `[]`;
    - Normal commits: a single parent;
    - Merge commits (later phases): multiple parents.

- `author` (string, required):
  - Free-form author identifier, typically `"Name <email>"`.

- `timestamp` (string, required):
  - UTC timestamp in ISO 8601 format (`YYYY-MM-DDTHH:MM:SSZ`).

- `message` (string, required):
  - Commit message.

Canonicalisation rules analogous to Objects apply: keys are serialised in a stable
order, minimal whitespace is used for hashing.

### 5.3 Example Snapshot

```json
{
  "version": 1,
  "root": "f3a1b2c3d4e5f6...",
  "parents": [
    "ab12cd34ef56..."
  ],
  "author": "Alice Researcher <alice@example.edu>",
  "timestamp": "2025-01-23T10:15:30Z",
  "message": "Add initial VASP calculation for system A"
}
```

---

## 6. Blob: Structure and Encoding

### 6.1 Blob Semantics

A **Blob** represents raw binary content, typically corresponding to a file’s
contents (input files, output logs, binary checkpoints, etc.). Blobs are opaque
to the core VCS; interpretation is handled by higher layers.

### 6.2 Blob Storage Format

In the MVP:

- Blob content is stored as raw bytes in a file under `.chemvcs/objects/`;
- The blob’s hash is computed directly from these bytes (no JSON wrapper).

Future enhancements may introduce chunking or further pack/index improvements, but these must preserve
the logical semantics that a blob is addressed by its content hash.

---

## 7. Refs and HEAD

### 7.1 Ref Semantics

A **Ref** is a named pointer to a Snapshot. The core supports:

- Branch refs under `refs/heads/` (e.g. `refs/heads/main`);
- HEAD, which typically points to a branch ref.

Refs are updated as part of commit and branch operations.

### 7.2 On-Disk Representation

Ref files are plain text files containing either:

- A snapshot hash, e.g.:

  ```text
  a1b2c3d4e5f6...
  ```

- Or (for `HEAD`) a symbolic reference:

  ```text
  ref: refs/heads/main
  ```

There is no JSON wrapper for ref files in the MVP.

### 7.3 HEAD Resolution

- `HEAD` is a file in `.chemvcs/` that usually contains `ref: refs/heads/<branch>`.
- To resolve `HEAD` to a snapshot:
  1. Read `HEAD`;
  2. If it starts with `ref: `, interpret the remainder as a ref path and read that file;
  3. The resulting contents are the current snapshot hash.

Detached HEAD states (HEAD pointing directly to a snapshot hash) may be supported
in later milestones; the on-disk format is the same as for branch refs.

---

## 8. Repository Layout

### 8.1 Top-Level Layout

Within a project root directory, ChemVCS stores its data in a `.chemvcs/` directory.
MVP layout:

```text
.chemvcs/
├─ objects/
│  ├─ ab/
│  │  └─ cdef1234...        # object or blob file, named by full hash
│  └─ ...
├─ refs/
│  └─ heads/
│     └─ main               # branch ref
├─ HEAD                     # current HEAD
└─ config                   # repository configuration (optional)
```

### 8.2 Objects Directory

The `objects/` directory contains both blobs and objects, named by their hash.

- The directory is sharded to avoid too many files in a single directory:
  - The first two hex characters of the hash define the shard directory, e.g. `ab/`;
  - The remaining 62 characters form the filename.

Example:

- Hash: `ab1234...xyz`
- Path: `.chemvcs/objects/ab/1234...xyz`

The file content is either:

- Raw blob bytes; or
- JSON representation of an Object or Snapshot.

The implementation MUST be able to distinguish between blob and non-blob files,
for example by:

- Keeping a small header or prefix (implementation detail), or
- Using a simple magic marker for JSON-encoded entities.

This detail is left to the implementation, as long as the logical model remains
consistent.

### 8.3 Refs Directory

The `refs/` directory holds named references to snapshots. For the MVP:

- `refs/heads/` contains branch refs;
- Other namespaces (e.g. `refs/tags/`) may be introduced later.

Ref paths are always relative to `.chemvcs/` and are plain text files containing
a snapshot hash.

### 8.4 Config File

The `config` file in `.chemvcs/` stores repository-level configuration. Its format
is intentionally simple in the MVP (e.g. key=value lines or a small JSON object).
The specific schema is defined in the CLI and repo implementation and may include:

- User name and email defaults;
- Preferred compression options;
- Future flags or feature toggles.

---

## 9. Integrity Constraints and Invariants

ChemVCS relies on several invariants to guarantee repository correctness.

### 9.1 Object and Snapshot Validity

For every stored Object and Snapshot:

- Its hash must be the SHA-256 of its canonical serialisation;
- `version` MUST be a supported version (currently `1`);
- `type` MUST be a non-empty string for Objects;
- `root` MUST be a valid hash string for Snapshots;
- `parents` MUST be an array of valid hash strings (can be empty);
- `refs[*].id` MUST be valid hash strings when present.

### 9.2 Ref Validity

For every ref under `refs/`:

- The file contents must either be:
  - A valid snapshot hash; or
  - For HEAD only, a `ref: <path>` string.
- If a snapshot hash is present, the referenced Snapshot SHOULD exist on disk.
  Implementations may accept temporarily dangling refs (e.g. during fetch), but
  must repair or report such states.

### 9.3 Hash and Content Consistency

When writing objects or snapshots:

- The implementation must compute the hash from the exact bytes that are written.
- It must refuse to overwrite existing files if the hash already exists with
  different content (potential collision or corruption).

When reading:

- The implementation should verify that the hash of the read content matches the
  expected hash (optional but recommended for robustness).

### 9.4 Order of Operations

To avoid inconsistent states:

- Objects and snapshots must be written **before** any refs that point to them.
- Ref updates should be atomic at the filesystem level where possible (e.g. write
  to a temporary file and rename).

---

## 10. Examples

### 10.1 Minimal Repository After Initial Commit

Assume a new repository with a single snapshot and branch `main`.

- `.chemvcs/HEAD`:

  ```text
  ref: refs/heads/main
  ```

- `.chemvcs/refs/heads/main`:

  ```text
  7f3c2b1a...
  ```

  Where `7f3c2b1a...` is the hash of the initial Snapshot.

- `.chemvcs/objects/7f/3c2b1a...`:

  ```json
  {
    "version": 1,
    "root": "a1b2c3d4...",
    "parents": [],
    "author": "Alice <alice@example.edu>",
    "timestamp": "2025-01-23T10:15:30Z",
    "message": "Initial commit"
  }
  ```

- `.chemvcs/objects/a1/b2c3d4...` (root Object, e.g. a folder):

  ```json
  {
    "version": 1,
    "type": "folder",
    "meta": {
      "name": "root"
    },
    "refs": [
      {
        "kind": "object",
        "id": "bcdef123..."
      }
    ]
  }
  ```

- `.chemvcs/objects/bc/def123...` (child Object, e.g. a file):

  ```json
  {
    "version": 1,
    "type": "file",
    "meta": {
      "name": "README.md",
      "mode": "0644"
    },
    "refs": [
      {
        "kind": "blob",
        "id": "e9f8a7b6..."
      }
    ]
  }
  ```

- `.chemvcs/objects/e9/f8a7b6...` (Blob, raw bytes of `README.md`):

  ```text
  # My Project
  ```

  (stored as raw bytes, not JSON).

### 10.2 Evolving the Model

Future versions of ChemVCS may:

- Add optional fields to `Object.meta` or snapshot JSON;
- Introduce additional `kind` types in `refs` (e.g. for remote references);
- Introduce and evolve packfile formats that group multiple objects in a single file.

Such changes must preserve the basic invariants and seek to remain compatible
with repositories created according to this specification.

---

## 11. Versioning and Compatibility

### 11.1 Schema Version Fields

Both Object and Snapshot encodings include a `version` field. This field is used to:

- Distinguish between different on-disk schemas;
- Allow the implementation to evolve while maintaining backward compatibility.

For the MVP:

- `version = 1` is the only allowed value for Objects and Snapshots.
- Future versions may introduce `version = 2`, etc., but must preserve the ability
  to read `version = 1` entities.

### 11.2 Repository Format Version

The repository-level format version (e.g. written to `config`) may also track
breaking changes such as packfile introduction or layout modifications. The exact
mechanism will be defined in later documents.

---

This specification defines the core object model and storage format for ChemVCS.
It is the primary reference for implementing the `model`, `objectstore`, and
`repo` packages in Go, as well as for any external tools that need to inspect or
manipulate `.chemvcs` contents directly.
