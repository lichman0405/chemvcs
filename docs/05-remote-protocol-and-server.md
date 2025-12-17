# ChemVCS – Remote Protocol and Server Design (v0.1)

## 1. Introduction

This document specifies the design of the ChemVCS remote protocol and the
corresponding server component. It defines:

- The responsibilities and architecture of the remote server (`chemvcs-server`);
- How local repositories interact with remote repositories;
- The HTTP-based protocol for listing repositories, pushing, pulling, and
  updating refs;
- Authentication, access control, and basic security considerations;
- Expected error handling and deployment guidelines.

The goal is a protocol that is:

- Simple to implement and debug;
- Sufficient for ChemVCS push / pull workflows;
- Robust against partial failures and corrupted data;
- Flexible enough to allow future optimisations (e.g. leveraging packs, streaming).

This version describes an MVP protocol suitable for Milestone 4. Later milestones
extend the server with optional capabilities (e.g. a remote HPC gateway) while
keeping the core sync protocol stable.

---

## 2. Server Responsibilities and Architecture

### 2.1 High-Level Responsibilities

The ChemVCS server is responsible for:

1. Hosting one or more ChemVCS repositories.
2. Exposing an HTTP API to:
   - Create and query repositories (administrative);
   - Upload and download objects (blobs, objects, snapshots);
   - Read and update refs;
   - Negotiate which objects are missing during push/pull.
3. Enforcing authentication and basic access control for repositories.
4. Maintaining integrity and consistency of on-disk repository data.

The server does **not**:

- Execute user code or scripts;
- Interpret domain-specific metadata beyond what is needed for consistency.

Note: ChemVCS may additionally expose repo-scoped HPC gateway endpoints that run
scheduler commands server-side for no-SSH workflows. Those endpoints are treated
as an extension and are documented separately.

### 2.2 Process and Deployment Model

The initial implementation is a single binary (e.g. `chemvcs-server`) that:

- Runs as a long-lived process;
- Listens on a configurable TCP port (HTTP);
- Reads configuration (repository root directory, authentication settings) from
  a config file and/or environment variables.

Typical deployment:

- A single server per group/lab/organisation;
- Fronted by a reverse proxy (e.g. Nginx) for TLS termination, logging, and
  further access control;
- File-based storage under a root such as `/var/lib/chemvcs/repos/`.

### 2.3 Internal Structure

The server reuses core concepts and components from the Go codebase:

- `model`, `objectstore`: to manage object encoding and storage;
- Repository abstraction to map a repo ID to a storage directory;
- Separate layers for HTTP routing, authentication, and request processing.

---

## 3. Repository Model on the Server

### 3.1 Repository Identification

Each repository is identified on the server by a **repository ID**, which is:

- A string used in URLs, e.g. `mygroup/project-a` or `user1/test-repo`;
- Mapped server-side to a filesystem directory, e.g.:

  ```text
  /var/lib/chemvcs/repos/mygroup/project-a/.chemvcs/
  ```

The mapping from repository ID to filesystem path is defined by server configuration
and internal logic, and must prevent directory traversal and similar attacks.

### 3.2 On-Disk Layout

Within each repository directory, the on-disk layout is equivalent to a local
ChemVCS repository’s `.chemvcs/` directory, as defined in the object model and
storage specification:

```text
/var/lib/chemvcs/repos/<repo-id>/
  ├─ objects/
  ├─ refs/
  ├─ HEAD
  └─ config
```

The server ensures that:

- All integrity constraints from the local repository spec are upheld;
- Object and ref writes follow the same ordering rules (write objects first,
  refs last).

### 3.3 Repository Creation and Administration

Repository creation and deletion are administrative operations. For the MVP:

- Repositories may be created:
  - Manually by the administrator (creating directories and initial config); or
  - Via an authenticated HTTP API endpoint (see §4).
- Permission models (owner, groups, public/private) can be introduced gradually.

---

## 4. HTTP API Overview

### 4.1 Base URL and Versioning

The server exposes an HTTP API under a base path, e.g.:

```text
http(s)://<host>/chemvcs/v1/
```

All endpoints described here assume a prefix of `/chemvcs/v1/`.

Versioning:

- `v1` is the initial API version;
- Backwards-compatible changes may be introduced within `v1`;
- Breaking changes will require `v2` or higher versions.

### 4.2 Authentication

For the MVP, the server supports token-based authentication:

- Clients include an `Authorization` header:

  ```http
  Authorization: Bearer <token>
  ```

- The server validates tokens against a configuration or a simple token store.
- Anonymous access (no token) may be allowed for read-only operations on
  selected repositories, depending on policy.

Authentication and token issuance are deployment-specific and not fully defined
here; this document assumes valid tokens can be obtained via separate procedures.

### 4.3 Content Types

The API uses standard content types:

- JSON: `application/json` for metadata requests and responses;
- Binary: `application/octet-stream` for object upload/download streams.

For simplicity, request and response bodies are small JSON objects, except for
bulk object transfer, which uses binary or streaming responses.

### 4.4 Error Responses

On errors, endpoints return an appropriate HTTP status code and a JSON body, e.g.:

```json
{
  "error": "not_found",
  "message": "Repository 'mygroup/project-a' not found"
}
```

Common `error` values:

- `not_found`
- `unauthorized`
- `forbidden`
- `bad_request`
- `conflict`
- `internal_error`

---

## 5. Core Endpoints

This section lists the main endpoints required for basic push/pull workflows.

### 5.1 Repository Info and Listing

#### 5.1.1 GET `/repos`

- Purpose: list repositories visible to the current user.
- Request: no body.
- Response `200 OK`:

  ```json
  {
    "repositories": [
      {
        "id": "mygroup/project-a",
        "description": "DFT study of system A"
      },
      {
        "id": "user1/test-repo",
        "description": "Personal test"
      }
    ]
  }
  ```

- Errors:
  - `401 Unauthorized` / `403 Forbidden` for unauthenticated or disallowed access.

#### 5.1.2 GET `/repos/{repoId}`

- Purpose: get basic information about a specific repository.
- Response `200 OK` example:

  ```json
  {
    "id": "mygroup/project-a",
    "description": "DFT study of system A",
    "default_branch": "main"
  }
  ```

- Errors:
  - `404 Not Found` if repository does not exist or is not visible;
  - `401` / `403` as appropriate.

### 5.2 Object Existence and Negotiation

Efficient push/pull requires determining which objects are already present on
the server or client.

#### 5.2.1 POST `/repos/{repoId}/objects/exists`

- Purpose: given a list of object hashes, determine which are already present.
- Request body:

  ```json
  {
    "ids": [
      "a1b2c3...",
      "d4e5f6..."
    ]
  }
  ```

- Response `200 OK`:

  ```json
  {
    "present": [
      "a1b2c3..."
    ],
    "missing": [
      "d4e5f6..."
    ]
  }
  ```

- Errors:
  - `404 Not Found` if repository does not exist;
  - `400 Bad Request` for invalid hash formats.

This endpoint is used by clients in push and pull operations to avoid redundant
transfers.

### 5.3 Object Upload

#### 5.3.1 POST `/repos/{repoId}/objects/upload`

- Purpose: upload one or more objects in a single request.
- Request:
  - Content-Type: `application/octet-stream` or a multipart format.
  - MVP approach: a simple framing format where each object is sent as:

    ```text
    <hash>\n
    <length>\n
    <raw-bytes>
    ```

    repeated for multiple objects.

  - Alternatively, a JSON envelope with base64-encoded content can be used, but
    is less efficient for large objects.

- Response `200 OK`:

  ```json
  {
    "stored": [
      "a1b2c3...",
      "d4e5f6..."
    ]
  }
  ```

- The server MUST verify that the hash matches the content before storing.

- Errors:
  - `400 Bad Request` for malformed frames or mismatched hashes;
  - `409 Conflict` if the server detects an existing object with a different
    content for the same hash (indicates corruption or collision);
  - `404 Not Found` if repository does not exist.

The exact framing and streaming details can be refined in implementation while
preserving the logical semantics: upload objects with known hashes, verify, store.

### 5.4 Object Download

#### 5.4.1 GET `/repos/{repoId}/objects/{hash}`

- Purpose: retrieve a single object or blob by hash.
- Response `200 OK`:

  - Headers:
    - `Content-Type: application/octet-stream`;
    - Optionally, custom header to indicate the logical type (blob/object/snapshot).
  - Body: raw bytes as stored.

- Errors:
  - `404 Not Found` if the object is missing;
  - `404 Not Found` if repository does not exist;
  - `400 Bad Request` for invalid hash format.

For bulk downloads, clients may either:

- Request known object hashes individually; or
- Use a future streaming/multipart endpoint (not defined in the MVP).

### 5.5 Refs

#### 5.5.1 GET `/repos/{repoId}/refs`

- Purpose: list refs (e.g. branches) in a repository.
- Response `200 OK` example:

  ```json
  {
    "refs": [
      {
        "name": "refs/heads/main",
        "target": "7f3c2b1a..."
      },
      {
        "name": "refs/heads/feature-x",
        "target": "ab12cd34..."
      }
    ]
  }
  ```

- Errors:
  - `404 Not Found` if repository does not exist.

#### 5.5.2 GET `/repos/{repoId}/refs/{refName}`

- Purpose: get a single ref’s target.
- Response `200 OK` example:

  ```json
  {
    "name": "refs/heads/main",
    "target": "7f3c2b1a..."
  }
  ```

- Errors:
  - `404 Not Found` if ref or repository does not exist.

#### 5.5.3 POST `/repos/{repoId}/refs/{refName}`

- Purpose: update a ref to point to a new snapshot, with optional optimistic
  concurrency control.
- Request body:

  ```json
  {
    "old_target": "ab12cd34...",  // optional, for concurrency control
    "new_target": "7f3c2b1a..."
  }
  ```

- Behaviour:
  - If `old_target` is provided and does not match the current ref on the server,
    the server MUST reject the update with `409 Conflict`.
  - If the `new_target` snapshot does not exist on the server, the server SHOULD
    reject the update with `400 Bad Request` or `409 Conflict`.

- Response `200 OK` example:

  ```json
  {
    "name": "refs/heads/main",
    "target": "7f3c2b1a..."
  }
  ```

- Errors:
  - `404 Not Found` if repository does not exist;
  - `409 Conflict` for concurrent updates or missing snapshots;
  - `400 Bad Request` for invalid hashes.

---

## 6. Push and Pull Workflows

This section describes how clients combine the endpoints to implement push and
pull operations.

### 6.1 Push (Client → Server)

Given a local branch `refs/heads/main` and a remote repository `origin`:

1. **Determine local objects for push**
   - The client computes the set of objects reachable from the local branch’s
     snapshot (root + closure over refs, including snapshots).

2. **Check server for existing objects**
   - The client sends the list of hashes (possibly batched) to
     `POST /repos/{repoId}/objects/exists`.
   - The server responds with `present` and `missing` lists.

3. **Upload missing objects**
   - The client uploads all `missing` objects using one or more calls to
     `POST /repos/{repoId}/objects/upload`.
   - The server validates and stores them.

4. **Update ref on server**
   - The client reads the current server ref (e.g. via
     `GET /repos/{repoId}/refs/heads/main`).
   - If the update should be fast-forward only, the client checks ancestry.
   - The client sends `POST /repos/{repoId}/refs/heads/main` with:
     - `old_target` = current server ref (if enforcing concurrency);
     - `new_target` = local snapshot hash.

5. **Handle conflicts**
   - If the server returns `409 Conflict`, the client reports that the remote
     has diverged and may require a pull/merge before pushing again.

### 6.2 Pull (Server → Client)

Given a local branch and a remote branch to track:

1. **Fetch remote ref**
   - Client calls `GET /repos/{repoId}/refs/heads/main` to get the remote
     snapshot hash.

2. **Determine remote objects for pull**
   - The client requests the closure of objects reachable from that snapshot.
     For MVP, this may involve:
       - Client calling `GET /repos/{repoId}/objects/{hash}` recursively; or
       - Using the `objects/exists` endpoint to compute what is missing locally.

3. **Download missing objects**
   - Client downloads all required objects (`GET /repos/{repoId}/objects/{hash}`)
     and stores them in its local object store.

4. **Update local ref**
   - Client updates local `refs/heads/main` to the remote snapshot hash (fast-forward
     logic handled locally).

Later optimisations may include server-side enumeration of reachable objects or
batch download endpoints.

---

## 7. Authentication, Authorisation, and Security

### 7.1 Authentication

The MVP uses bearer tokens in the `Authorization` header. The server must:

- Reject requests with missing or invalid tokens for protected operations;
- Optionally allow unauthenticated read-only access to selected repositories.

Token management (creation, revocation, expiry) is deployment-specific.

### 7.2 Authorisation

At minimum, the server enforces:

- Read vs write permissions per repository:
  - Read: allow `GET` on refs and objects.
  - Write: allow `POST` on objects and refs.

Repository-level ACLs (per-user, per-group) can be configured via server config
or an external identity system.

### 7.3 Transport Security

The recommended deployment is:

- Run `chemvcs-server` behind a reverse proxy that terminates TLS;
- Require HTTPS for all client–server communication;
- Use secure token storage and transmission practices.

### 7.4 Data Integrity and Validation

- The server must verify that object hashes match their content before storing.
- Ref updates must only point to existing snapshots.
- The server should implement basic sanity checks to prevent malicious or
  malformed requests from corrupting repositories.

---

## 8. Error Handling and Diagnostics

### 8.1 HTTP Status Codes

The server should use standard HTTP status codes:

- `200 OK` – request succeeded;
- `201 Created` – resource created (if applicable);
- `400 Bad Request` – invalid input, malformed hashes, etc.;
- `401 Unauthorized` – missing/invalid authentication;
- `403 Forbidden` – authenticated but not allowed;
- `404 Not Found` – repository or object not found;
- `409 Conflict` – ref concurrency conflict, object hash/content mismatch;
- `500 Internal Server Error` – unexpected server-side failure.

### 8.2 Error Payload

All error responses should include a JSON payload as described in §4.4, containing:

- `error` – short machine-oriented code;
- `message` – human-readable description.

### 8.3 Logging and Metrics

The server should log:

- Authentication failures;
- Repository-level operations (push/pull, ref updates);
- Integrity violations (hash mismatches, missing objects).

Optional metrics may include:

- Number of repositories;
- Number of objects per repository;
- Request rates per endpoint;
- Error rates and latency.

---

## 9. Deployment Considerations

### 9.1 Single-Server Deployment

The simplest deployment is:

- One `chemvcs-server` process;
- Repositories stored under `/var/lib/chemvcs/repos/`;
- Reverse proxy (e.g. Nginx) providing TLS and routing (`/chemvcs/v1/` prefix).

### 9.2 Backup and Recovery

Repository data is stored entirely on disk. Backup typically involves:

- Regularly backing up the repository root directory;
- Optionally snapshotting the filesystem (ZFS, LVM, etc.).

On restore:

- The server can operate directly on restored repositories, as long as the on-disk
  layout is intact.

### 9.3 Future Scalability

Potential future scaling strategies include:

- Sharding repositories across multiple servers;
- Using an object storage backend (e.g. S3-compatible) for `objects/`;
- Leveraging packfiles and caching to reduce I/O overhead.

These changes should not require altering the basic HTTP semantics defined here,
only implementation details and additional configuration.

---

This document defines the MVP remote protocol and server behaviour for ChemVCS.
It is the primary reference for implementing the `remote` client package in Go
and the `chemvcs-server` binary, and for integrating ChemVCS with existing
infrastructure in later phases.
