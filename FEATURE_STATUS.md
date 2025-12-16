# ChemVCS Feature Status

**Last updated**: 2025-12-16

This document is a short, source-of-truth snapshot of what ChemVCS can do *today* and what is still missing. It intentionally avoids volatile line counts, test counts, and commit hashes.

---

## ✅ Implemented

### Core VCS (Go)
- Content-addressable object storage (SHA-256) with snapshots and refs
- Local repo workflow: `init`, `status`, `commit`, `log`, `branch`, `checkout`, `merge`

### Remote repository (Go)
- HTTP server (`chemvcs-server`) and remote client
- Repo-scoped endpoints for object transfer and refs
- CLI remotes: `remote add`, `remote list` and sync commands (push/pull/fetch/clone)

### Python domain layer
- Domain objects: Structure, Run, Workflow
- IO: XYZ and POSCAR
- JSON interop with Go core

### HPC integration (local execution)
- Python HPC adapters (SLURM + PBS + LSF) with provenance capture
- CLI commands that can submit/track/retrieve when scheduler tools are available locally

### Remote HPC gateway (no SSH)
- Server-side HPC endpoints under the repo-scoped API (submit/jobs/get/cancel/retrieve)
- CLI `--remote=<name>` support for HPC commands (laptop talks HTTP only)
- Job record persistence in the repo (blobs + `refs/hpc/...`)
- Hardened retrieve: repo-boundary enforcement and path traversal protections
- Workflow closure: `chemvcs retrieve --commit` to snapshot retrieved outputs
- Remote retrieve streams zip to disk (avoids large in-memory buffers)
- `chemvcs jobs [<run-hash|job-id>]` filtering

---

## ⚠️ Implemented but still limited (remaining gaps)

### Multi-user / identity
- Server auth is token-based (MVP) and repo-scoped; no per-user OS identity mapping
- All HPC jobs run under the server service account (by design for MVP)

### Advanced repo operations
- No pack/gc/fsck maintenance tooling yet
- No rich query/indexing for domain objects

---

## ⏭️ Next priorities (recommended)

1. Enrich HPC job status details (exit code/reason/timestamps) and expose them via API/CLI
2. Make retrieve/commit workflows more explicit (e.g., commit only retrieved dest, not whole tree)
3. Add repo maintenance commands (pack/gc/fsck) for large workspaces

