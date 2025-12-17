# ChemVCS – CLI Specification (MVP, v0.1)

## 1. Introduction

This document specifies the behaviour of the ChemVCS command-line interface (CLI)
for the initial milestones (1–3). It covers:

- Global CLI behaviour and conventions;
- Core commands:
  - `chemvcs init`
  - `chemvcs commit`
  - `chemvcs log`
  - `chemvcs status`
  - `chemvcs branch`
  - `chemvcs checkout`
- Command arguments and options;
- Exit codes and error handling expectations;
- Examples for typical usage patterns.

The CLI is the primary user-facing interface to the Go-based VCS core. It is also
intended to be a stable integration point for scripts and tools that wish to drive
ChemVCS programmatically (e.g. the Python domain layer).

This specification focuses on behaviour for local repositories. Remote-related
commands (`remote`, `push`, `pull`, `fetch`) will be defined in the remote
protocol and server specification.

Note: ChemVCS has grown beyond this MVP scope (remote sync, HPC commands, and
repository maintenance). This document remains the source of truth for the
Milestone 1–3 command semantics; newer commands are documented in:
- Remote: `docs/05-remote-protocol-and-server.md`
- HPC: `docs/09-hpc-integration-design.md` and `docs/10-hpc-user-guide.md`
- Current feature surface: `FEATURE_STATUS.md`

---

## 2. General CLI Conventions

### 2.1 Command Name

The CLI executable is named:

```bash
chemvcs
```

Subcommands are invoked as:

```bash
chemvcs <command> [options] [arguments]
```

### 2.2 Working Directory and Repository Discovery

Unless otherwise specified by an option, commands operate on the repository
containing the current working directory. Repository discovery follows these rules:

1. Starting from the current working directory, search upwards for a `.chemvcs`
   directory.
2. If found, that directory is considered the repository root.
3. If no `.chemvcs` is found up to the filesystem root, the command behaves as if
   no repository exists and may fail (except for `init`, which can create one).

Some commands may accept an optional path argument to override the default working
directory (see command-specific sections).

### 2.3 Output Modes

By default, commands produce **human-oriented text output**. Where appropriate,
a `--format` or `--json` option may be added in later iterations to support
machine-readable output. For the MVP:

- Primary output is simple text;
- Exit codes are used to indicate success or failure;
- Error messages are printed to `stderr`.

### 2.4 Exit Codes

Unless otherwise specified, commands follow these exit code conventions:

- `0` – successful execution;
- `1` – general error (invalid arguments, repository not found, etc.);
- `2` – usage error (e.g. missing required options, invalid flags);
- Other codes may be introduced for specific error classes if needed.

### 2.5 Environment Variables

The following environment variables may affect CLI behaviour (MVP suggestions):

- `CHEMVCS_AUTHOR_NAME` – default author name if not specified in config;
- `CHEMVCS_AUTHOR_EMAIL` – default author email;
- `CHEMVCS_REPO` – optional explicit path to a repository root (overrides discovery).

Exact precedence between CLI flags, config, and environment variables is an
implementation detail, but the typical priority is:

1. CLI flags;
2. Environment variables;
3. Repository config defaults.

---

## 3. `chemvcs init`

### 3.1 Synopsis

```bash
chemvcs init [<directory>]
```

### 3.2 Description

Initialises a new ChemVCS repository in the specified directory or in the current
working directory if none is given.

- Creates a `.chemvcs/` directory with the standard layout:
  - `objects/`
  - `refs/heads/`
  - `HEAD`
  - `config` (optional)
- Sets `HEAD` to point to `refs/heads/main` by default.
- Does **not** create any initial snapshot or commit.

### 3.3 Arguments

- `<directory>` (optional): path to initialise the repository in. If omitted,
  uses the current working directory.

### 3.4 Options

None in the MVP. Future options might include:

- `--bare` – create a repository without a working directory (not in MVP).

### 3.5 Behaviour

- If `<directory>` does not exist, the command MUST fail with an error.
- If `.chemvcs/` already exists in the target directory, the command SHOULD
  fail with a clear error indicating that a repository already exists.
- No changes are made to other files in the directory.

### 3.6 Output

On success, prints a short confirmation message, for example:

```text
Initialised empty ChemVCS repository in /path/to/project/.chemvcs
```

On failure, prints a descriptive error message to `stderr`.

### 3.7 Exit Codes

- `0` – repository initialised successfully;
- `1` – failure (e.g. directory not found, `.chemvcs` already exists).

---

## 4. `chemvcs commit`

### 4.1 Synopsis

```bash
chemvcs commit -m <message> [--author <author>] [--allow-empty]
```

### 4.2 Description

Creates a new snapshot (commit) that represents the current state of the working
directory, according to the generic file-tree mapping for the MVP.

The command:

1. Determines the current repository via discovery or `CHEMVCS_REPO`.
2. Resolves `HEAD` to find the current branch and parent snapshot (if any).
3. Uses the `workspace` logic to:
   - Scan the working directory;
   - Construct an object graph representing directories and files;
   - Compute hashes and store any new blobs/objects.
4. Builds a new Snapshot object with:
   - `root` = root object hash;
   - `parents` = [current snapshot hash] (or empty if initial commit);
   - `author`, `timestamp`, `message` provided or derived from defaults.
5. Stores the Snapshot and updates the current branch ref to point to it.

### 4.3 Arguments

None (for the MVP). Future versions might support a pathspec to limit the set of
tracked files, but MVP treats the working directory as a whole.

### 4.4 Options

- `-m, --message <message>` (required):
  - Commit message to use for the snapshot.

- `--author <author>` (optional):
  - Overrides the author string, which otherwise defaults to some combination of
    config and environment variables (e.g. `Name <email>`).

- `--allow-empty` (optional):
  - If present, allows creating a commit even if the resulting snapshot is
    identical to the parent snapshot (no changes).
  - Without this flag, attempting to commit with no changes SHOULD result in an error.

Future options (not in MVP) might include `--amend`, `--no-verify`, etc.

### 4.5 Behaviour

- If no repository is found, the command MUST fail with an error.
- If `-m/--message` is missing, the command MUST fail with a usage error.
- If there are no differences compared to the parent snapshot and
  `--allow-empty` is not set, the command SHOULD:
  - Print a message indicating there is nothing to commit;
  - Exit with a non-zero status (usage or general error).

### 4.6 Output

On success, prints at least the new snapshot hash, for example:

```text
[main 7f3c2b1a...] Commit message here
```

Exact formatting is implementation-specific but should include:

- The branch name (if known);
- The snapshot hash (abbreviated or full);
- The commit message.

On failure, prints a descriptive error message to `stderr`.

### 4.7 Exit Codes

- `0` – snapshot created successfully;
- `1` – general failure (e.g. repository not found, internal error);
- `2` – usage error (e.g. missing `-m`, invalid option).

---

## 5. `chemvcs log`

### 5.1 Synopsis

```bash
chemvcs log [<ref>] [--max-count <n>] [--oneline]
```

### 5.2 Description

Displays the snapshot history starting from a given ref or from `HEAD` if no ref
is specified. Traverses parent links backwards (a simple linear walk for MVP).

### 5.3 Arguments

- `<ref>` (optional):
  - A ref name (`main`, `feature-x`, `refs/heads/main`) or a snapshot hash.
  - If omitted, `HEAD` is used.

### 5.4 Options

- `--max-count <n>` (optional):
  - Limit the number of snapshots printed.
  - Default: no limit (or a sensible implementation default, e.g. 100).

- `--oneline` (optional):
  - Print each commit on a single line (hash and message), suitable for quick scanning.

Future options (not in MVP) might include formatting templates or filtering by author.

### 5.5 Behaviour

- Resolves `<ref>` or `HEAD` to a snapshot hash.
- Follows the `parents` chain from that snapshot backwards until:
  - No more parents exist; or
  - `--max-count` is reached.
- Prints snapshots in reverse chronological order (most recent first).

### 5.6 Output

Default (multi-line) format per snapshot, e.g.:

```text
commit 7f3c2b1a...
Author: Alice Researcher <alice@example.edu>
Date:   2025-01-23T10:15:30Z

    Add initial VASP calculation for system A
```

With `--oneline`, e.g.:

```text
7f3c2b1a... Add initial VASP calculation for system A
ab12cd34... Initialise repository
```

Exact formatting is flexible as long as it is stable and documented.

### 5.7 Exit Codes

- `0` – log printed successfully;
- `1` – error (e.g. repository not found, ref not found).

---

## 6. `chemvcs status`

### 6.1 Synopsis

```bash
chemvcs status
```

### 6.2 Description

Shows the difference between the current working directory and the snapshot
pointed to by `HEAD` (or indicates that no commits exist yet). This is based on
the generic file-tree mapping in the `workspace` component.

### 6.3 Arguments and Options

No arguments or options in the MVP. Future extensions may include:

- `--short` / `--porcelain` for machine-readable output;
- Pathspecs to restrict status to certain subdirectories.

### 6.4 Behaviour

- If no repository is found, fails with an error.
- If no snapshot is pointed to by the current branch (no commits yet):
  - Treats all tracked files as “added” (or simply reports that the branch has no commits).
- Otherwise:
  1. Loads the snapshot referenced by `HEAD`;
  2. Constructs an object graph representing the snapshot’s file tree;
  3. Constructs an object graph representing the current working directory;
  4. Compares them to detect:
     - Added files;
     - Modified files;
     - Deleted files.

The exact diffing algorithm is an implementation detail; at minimum it should
compare file paths and content hashes.

### 6.5 Output

A simple, human-friendly format, for example:

```text
On branch main
Changes to be committed:
  (use "chemvcs commit" to record these changes)

    modified:   src/main.py
    added:      data/input_02.in
    deleted:    old/unused.txt
```

For the MVP, there is no staging area; status compares the working directory
directly against the latest snapshot.

### 6.6 Exit Codes

- `0` – status printed successfully;
- `1` – error (e.g. no repository);
- Future: non-zero codes may be used to indicate “has changes” vs “clean”, if useful
  for scripting.

---

## 7. `chemvcs branch`

### 7.1 Synopsis

```bash
chemvcs branch            # list branches
chemvcs branch <name>     # create branch at current HEAD
```

(Deletion or renaming of branches is not part of the MVP.)

### 7.2 Description

Manages branches (refs under `refs/heads/`). In the MVP:

- Without arguments, lists all branches, indicating the current one.
- With a `<name>` argument, creates a new branch that points to the current
  snapshot (as determined by `HEAD`).

### 7.3 Arguments

- `<name>` (optional):
  - Name of the branch to create (e.g. `feature-x`).
  - If omitted, the command lists existing branches.

### 7.4 Options

None in the MVP. Future options might include:

- `-d` / `--delete` for branch deletion;
- `-m` / `--move` for renaming.

### 7.5 Behaviour

- Listing:
  - Reads all ref files under `refs/heads/`;
  - Resolves `HEAD` to identify the current branch;
  - Prints one line per branch, marking the current one.

- Creation:
  - Ensures `HEAD` is currently pointing to a valid snapshot;
  - Creates `refs/heads/<name>` and writes the current snapshot hash into it;
  - Does not switch the working directory or `HEAD` to the new branch (MVP may
    choose either behaviour; it must be documented. A common choice is not to
    switch automatically).

### 7.6 Output

- Listing example:

  ```text
* main
  feature-x
  test-branch
  ```

  Here `*` marks the current branch.

- Creation example:

  ```text
Created branch feature-x at 7f3c2b1a...
  ```

### 7.7 Exit Codes

- `0` – success (listing or creation);
- `1` – failure (e.g. invalid repository, branch already exists, HEAD not pointing
  to a valid snapshot).

---

## 8. `chemvcs checkout`

### 8.1 Synopsis

```bash
chemvcs checkout <branch-or-hash>
```

### 8.2 Description

Switches the current branch or moves `HEAD` to a given snapshot, and updates the
working directory to reflect the selected snapshot (in later milestones).

For the MVP, checkout focuses on:

- Changing the `HEAD` symbolic ref to point to another branch; and
- Optionally materialising the snapshot into the working directory using the
  `workspace` component.

Complex scenarios (uncommitted changes, merge conflicts) will be handled in a
conservative way for the MVP (e.g. refusing to overwrite changes without explicit
confirmation or requiring a clean working tree).

### 8.3 Arguments

- `<branch-or-hash>` (required):
  - A branch name (e.g. `feature-x`) or a snapshot hash.

### 8.4 Options

No options in the MVP. Future options may include:

- `--detach` – explicitly enter detached HEAD state;
- `--force` – force checkout discarding local changes.

### 8.5 Behaviour

- If `<branch-or-hash>` resolves to a branch:
  - `HEAD` is updated to `ref: refs/heads/<branch>`;
  - The current branch is now `<branch>`.
- If `<branch-or-hash>` resolves to a snapshot hash:
  - MVP may either:
    - Enter detached HEAD state (HEAD points directly to the hash), or
    - Refuse and require a branch name.
  - The chosen semantics must be documented and consistent.

- Working directory update:
  - If the implementation has `workspace` materialisation ready:
    - It replaces the working directory contents to match the snapshot;
    - It must protect or refuse checkout if there are uncommitted modifications,
      unless a force flag is introduced.
  - In early MVP implementations, checkout may update only HEAD and print a warning
    that working directory synchronisation is not yet supported.

### 8.6 Output

Example:

```text
Switched to branch feature-x
```

or, in case of detached HEAD:

```text
Note: checking out '7f3c2b1a...'. You are in a detached HEAD state.
```

On failure (e.g. branch not found, uncommitted changes), prints a descriptive
error message.

### 8.7 Exit Codes

- `0` – checkout successful;
- `1` – failure (invalid ref, conflicts with uncommitted changes, repository issues).

---

## 9. Error Handling and Messages

ChemVCS CLI should use clear, actionable error messages. Typical error classes:

- **Repository errors:**
  - No `.chemvcs` directory found;
  - `.chemvcs` is corrupt or partially initialised;
  - Refs or snapshots pointing to missing objects.

- **Usage errors:**
  - Missing required options (`-m` for `commit`);
  - Unknown commands or options;
  - Invalid arguments (e.g. illegal branch names).

- **Operational errors:**
  - Filesystem errors (permissions, out-of-space);
  - Hash mismatches or integrity failures (rare, serious).

In all cases:

- Error messages MUST be written to `stderr`;
- Exit codes MUST be non-zero;
- Where possible, messages SHOULD suggest a remediation step or at least indicate
  the failing component (e.g. “failed to open repository”, “failed to write object”).

---

## 10. Future Extensions (Beyond MVP)

The CLI is expected to grow with additional commands and options as ChemVCS evolves:

- Remote operations: `remote`, `push`, `pull`, `fetch`;
- Tagging: `tag` for lightweight tags;
- More advanced branching features: `branch -d`, `merge`, etc.;
- Domain-aware commands surfaced from the Python layer (`chemvcs run ...`,
  `chemvcs query ...`).

These extensions will follow the same conventions as specified in this document
and will be captured in updated CLI specifications as the project matures.
