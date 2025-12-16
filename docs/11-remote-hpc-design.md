# ChemVCS – Remote HPC Gateway Design (Scheme C)

**Status**: Draft (implemented MVP)

## 1. Goal
Enable ChemVCS users to run the CLI on a laptop and interact with a remote SLURM environment **without SSH**.

The laptop talks HTTP to `chemvcs-server`. The server runs `sbatch/squeue/sacct/scancel` locally on the SLURM host (single-user service account for MVP).

## 2. Key Constraints (from current ChemVCS architecture)
- Current Go CLI HPC flows call Python (`chemvcs_py`) which assumes SLURM commands exist on the same machine.
- Running Python domain/HPC logic inside the server would significantly increase deployment complexity.
- ChemVCS already has an HTTP server (`chemvcs-server`) and a Go HTTP client (`internal/remote`).

## 3. Architecture
### 3.1 High-level
- Laptop: `chemvcs` CLI
- Server (SLURM host): `chemvcs-server` + SLURM CLI tools

The server is the **HPC gateway**.

### 3.2 Single-user MVP
- `chemvcs-server` runs under a single Unix account (service account).
- All jobs are submitted as that account.
- Future multi-user support is preserved by keeping the API and persistence format identity-agnostic.

## 4. API (implemented)
All routes are scoped to a repository:

Base prefix: `/chemvcs/v1/repos/{owner}/{repo}`

### 4.1 Submit
- `POST /hpc/submit`
- Request:
  - `run_hash` (string)
  - `script` (string, full job script content)
  - `capture_env` (bool)
- Response:
  - `job_id`, `run_hash`, `job_system`

Server behavior:
- Writes `script` to a temp file.
- Calls `sbatch --parsable --chdir <repoRoot> [--export=ALL] <tempScript>`.
- Persists a job record (see Section 5).

### 4.2 List jobs
- `GET /hpc/jobs?status=<STATUS>&refresh=1`
- Response: `{ "jobs": [ ... ] }`

`refresh=1` triggers server-side status refresh using `squeue` and fallback to `sacct`.

### 4.3 Get a single job
- `GET /hpc/jobs/{identifier}?refresh=1`

`identifier` can be either a `run_hash` or a numeric `job_id`.

### 4.4 Cancel
- `POST /hpc/cancel`
- Request: `{ "identifier": "<run-hash|job-id>" }`

Server calls `scancel <job_id>`.

### 4.5 Retrieve
- `POST /hpc/retrieve`
- Request: `{ "patterns": ["*.out", "*.log"] }`
- Response: `application/zip`

The server zips matching files from the repository working directory (excluding `.chemvcs`).

## 5. Persistence Model (implemented)
ChemVCS HPC state is persisted inside the repository as:
- Blob object: JSON `hpcJobRecord`
- Refs:
  - `refs/hpc/runs/<run_hash>` → blob hash
  - `refs/hpc/jobs/<job_id>` → blob hash

This keeps HPC state separate from snapshot history and avoids rewriting existing domain objects.

## 6. CLI UX (implemented)
HPC commands accept `--remote=<name>` to use the HTTP gateway:
- `chemvcs submit --remote=<name> <run-hash> <script> [--capture-env]`
- `chemvcs jobs --remote=<name> [--status=...] [<run-hash|job-id>]`
- `chemvcs watch --remote=<name> <run-hash|job-id> [--interval=...] [--timeout=...]`
- `chemvcs cancel --remote=<name> <run-hash|job-id>`
- `chemvcs retrieve --remote=<name> <run-hash> [--patterns=...] [--dest=...] [--commit] [--commit-message=...]`

Remote configuration uses existing:
- `chemvcs remote add <name> <url>`

**Recommended URL format** (repo-scoped):
- `http://<host>:<port>/chemvcs/v1/repos/<owner>/<repo>`

## 7. Deployment Notes (SLURM host)
- `chemvcs-server` must be installed on the SLURM host (or a host that can run `sbatch/squeue/sacct/scancel`).
- Ensure PATH/environment for the service account can invoke SLURM commands.

### 7.1 Authentication (recommended)
`chemvcs-server` supports a simple bearer token.

- Configure a token:
  - Flag: `--auth-token <token>`
  - Or env: `CHEMVCS_SERVER_AUTH_TOKEN=<token>`
- Configure repo scope for that token:
  - Flag: `--auth-repos owner/repo,owner2/repo2` (or `*`)
  - Or env: `CHEMVCS_SERVER_AUTH_REPOS=owner/repo,owner2/repo2`
- Optional admin token (allows listing repos and bypasses repo scoping):
  - Flag: `--admin-token <token>`
  - Or env: `CHEMVCS_SERVER_ADMIN_TOKEN=<token>`

On the client side, set one of:
- `CHEMVCS_REMOTE_TOKEN=<token>` (applies to all remotes)
- `CHEMVCS_REMOTE_TOKEN_<REMOTE_NAME>=<token>` (per remote, name is uppercased and non-alnum becomes `_`)

### 7.2 Network and TLS
For non-toy deployments:
- Bind `chemvcs-server` to localhost and put it behind a reverse proxy (nginx/caddy) with TLS.
- Restrict ingress to trusted networks (VPN, security group, firewall rules).
- Consider rate limits and request size limits at the proxy layer.

### 7.3 Example: systemd + reverse proxy

#### systemd unit

Create `/etc/systemd/system/chemvcs-server.service`:

```ini
[Unit]
Description=ChemVCS Server
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
User=chemvcs
Group=chemvcs
WorkingDirectory=/srv/chemvcs
EnvironmentFile=-/etc/chemvcs-server.env
ExecStart=/usr/local/bin/chemvcs-server \
  --listen 127.0.0.1:8080 \
  --root /srv/chemvcs
Restart=on-failure
RestartSec=2

# Optional hardening
NoNewPrivileges=true
PrivateTmp=true
ProtectSystem=strict
ProtectHome=true
ReadWritePaths=/srv/chemvcs

[Install]
WantedBy=multi-user.target
```

Then create `/etc/chemvcs-server.env` (permissions `0600` recommended):

```bash
CHEMVCS_SERVER_AUTH_TOKEN=change-me
CHEMVCS_SERVER_AUTH_REPOS=owner/repo
# CHEMVCS_SERVER_ADMIN_TOKEN=change-me-too
```

Reload and start:

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now chemvcs-server
sudo systemctl status chemvcs-server
```

#### Reverse proxy (nginx)

Example `server` block (TLS config omitted):

```nginx
server {
  listen 443 ssl;
  server_name chemvcs.example.com;

  # Add your ssl_certificate / ssl_certificate_key here

  client_max_body_size 10m;

  location / {
    proxy_pass http://127.0.0.1:8080;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
  }
}
```

## 8. Future Work (multi-user)
- Add authentication and request identity propagation.
- Map identity → Unix user, then submit jobs as that user (e.g., via `sudo -u`, setuid helper, or a dedicated job-submission service).
- Audit logging and per-user job visibility.

