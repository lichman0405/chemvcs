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

---

## ⚠️ Implemented but not production-hardened (gaps)

### Security
- No authentication/authorization on `chemvcs-server` (repo operations and HPC gateway)
- No audit log / request tracing suitable for multi-user environments

### Concurrency and correctness
- Limited conflict handling for concurrent ref updates (multiple clients)
- HPC job record/ref updates need stronger locking or conflict strategy

### Operations
- No deployment guide for systemd/reverse proxy/TLS
- No quotas / rate limiting / timeouts policy for HPC endpoints

---

## ⏭️ Next priorities (recommended)

1. Add server authentication (token-based for MVP)
2. Add server-side authorization rules (repo scope + HPC scope)
3. Add structured audit logging for HPC actions (submit/cancel/retrieve)
4. Add concurrency control for ref writes (optimistic concurrency or server-side locking)
5. Add “retrieve → commit” loop closure (optional auto-commit of retrieved outputs)

