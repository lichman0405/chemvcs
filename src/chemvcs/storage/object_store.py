"""Content-addressable blob storage for ChemVCS.

This module implements a Git-like object store using SHA-256 hashing for
content addressing. Blobs are stored in .chemvcs/objects/ with automatic
deduplication and optional gzip compression for large files.
"""

import gzip
import hashlib
import os
import tempfile
from pathlib import Path
from typing import Optional

from chemvcs.constants import (
    GZIP_THRESHOLD,
    HASH_ALGORITHM,
    HASH_LENGTH,
    OBJECTS_DIR,
)


class BlobNotFoundError(Exception):
    """Raised when a blob cannot be found in the object store."""

    pass


class BlobCorruptedError(Exception):
    """Raised when a blob's hash doesn't match its content."""

    pass


class ObjectStore:
    """Content-addressable storage for file blobs.
    
    Stores file contents as blobs identified by SHA-256 hash. Provides
    automatic deduplication, atomic writes, and optional compression.
    
    Storage layout:
        .chemvcs/objects/<hash[:2]>/<hash[2:]>      # Raw blob
        .chemvcs/objects/<hash[:2]>/<hash[2:]>.gz   # Compressed blob
    
    Attributes:
        objects_dir: Path to the objects directory
        
    Example:
        >>> store = ObjectStore(Path(".chemvcs"))
        >>> blob_hash = store.write_blob(b"ENCUT = 520\\n")
        >>> content = store.read_blob(blob_hash)
        >>> assert content == b"ENCUT = 520\\n"
    """

    def __init__(self, chemvcs_dir: Path) -> None:
        """Initialize the object store.
        
        Args:
            chemvcs_dir: Path to .chemvcs directory
            
        Raises:
            ValueError: If chemvcs_dir doesn't exist
        """
        self.chemvcs_dir = Path(chemvcs_dir)
        self.objects_dir = self.chemvcs_dir / OBJECTS_DIR
        
        if not self.chemvcs_dir.exists():
            raise ValueError(f"ChemVCS directory not found: {chemvcs_dir}")

    def write_blob(
        self,
        content: bytes,
        compress: Optional[bool] = None,
    ) -> str:
        """Write a blob to the object store.
        
        If a blob with the same hash already exists, returns the hash without
        writing (deduplication). Uses atomic write (tmp file + rename) to
        prevent corruption.
        
        Args:
            content: Binary content to store
            compress: Force compression (True), no compression (False), or
                     auto-compress if >200MB (None, default)
        
        Returns:
            SHA-256 hash of the content (64 hex characters)
            
        Raises:
            OSError: If write fails (permissions, disk full, etc.)
            
        Example:
            >>> hash1 = store.write_blob(b"data")
            >>> hash2 = store.write_blob(b"data")
            >>> assert hash1 == hash2  # Deduplication
        """
        # Compute hash
        blob_hash = self._compute_hash(content)
        
        # Check if blob already exists (deduplication)
        if self.blob_exists(blob_hash):
            return blob_hash
        
        # Determine compression
        should_compress = compress
        if compress is None:
            should_compress = len(content) >= GZIP_THRESHOLD
        
        # Prepare content
        if should_compress:
            data_to_write = gzip.compress(content, compresslevel=6)
            use_gz_ext = True
        else:
            data_to_write = content
            use_gz_ext = False
        
        # Get target path
        blob_path = self._get_blob_path(blob_hash, compressed=use_gz_ext)
        blob_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Atomic write: tmp file -> rename
        tmp_fd, tmp_path = tempfile.mkstemp(
            dir=blob_path.parent,
            prefix=".tmp_",
            suffix=".blob",
        )
        try:
            with os.fdopen(tmp_fd, "wb") as f:
                f.write(data_to_write)
                f.flush()
                os.fsync(f.fileno())  # Ensure written to disk
            
            # Atomic rename (works on POSIX; on Windows may fail if target exists)
            # On Windows, we need to handle the case where target exists
            try:
                os.replace(tmp_path, blob_path)
            except OSError:
                # Race condition: another process created it first
                if blob_path.exists():
                    os.unlink(tmp_path)
                    return blob_hash
                raise
            
            return blob_hash
            
        except Exception:
            # Clean up temp file on error
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise

    def read_blob(self, blob_hash: str, verify_hash: bool = True) -> bytes:
        """Read a blob from the object store.
        
        Args:
            blob_hash: SHA-256 hash of the blob (64 hex characters)
            verify_hash: Whether to recompute and verify hash (default: True)
        
        Returns:
            Binary content of the blob
            
        Raises:
            BlobNotFoundError: If blob doesn't exist
            BlobCorruptedError: If hash verification fails
            ValueError: If blob_hash is invalid format
            
        Example:
            >>> content = store.read_blob("abc123...")
        """
        self._validate_hash(blob_hash)
        
        # Try compressed first, then uncompressed
        compressed_path = self._get_blob_path(blob_hash, compressed=True)
        uncompressed_path = self._get_blob_path(blob_hash, compressed=False)
        
        content: bytes
        is_compressed = False
        
        if compressed_path.exists():
            blob_path = compressed_path
            is_compressed = True
        elif uncompressed_path.exists():
            blob_path = uncompressed_path
        else:
            raise BlobNotFoundError(
                f"Blob not found: {blob_hash} (tried {uncompressed_path} and {compressed_path})"
            )
        
        # Read content
        with open(blob_path, "rb") as f:
            data = f.read()
        
        # Decompress if needed
        if is_compressed:
            content = gzip.decompress(data)
        else:
            content = data
        
        # Verify hash
        if verify_hash:
            actual_hash = self._compute_hash(content)
            if actual_hash != blob_hash:
                raise BlobCorruptedError(
                    f"Blob corrupted: expected {blob_hash}, got {actual_hash}"
                )
        
        return content

    def blob_exists(self, blob_hash: str) -> bool:
        """Check if a blob exists in the store.
        
        Args:
            blob_hash: SHA-256 hash of the blob
            
        Returns:
            True if blob exists (compressed or uncompressed)
            
        Example:
            >>> if store.blob_exists(hash):
            ...     content = store.read_blob(hash)
        """
        try:
            self._validate_hash(blob_hash)
        except ValueError:
            return False
        
        compressed_path = self._get_blob_path(blob_hash, compressed=True)
        uncompressed_path = self._get_blob_path(blob_hash, compressed=False)
        
        return compressed_path.exists() or uncompressed_path.exists()

    def get_blob_size(self, blob_hash: str) -> int:
        """Get the size of a blob on disk (after compression if applicable).
        
        Args:
            blob_hash: SHA-256 hash of the blob
            
        Returns:
            Size in bytes (compressed size if blob is compressed)
            
        Raises:
            BlobNotFoundError: If blob doesn't exist
        """
        compressed_path = self._get_blob_path(blob_hash, compressed=True)
        uncompressed_path = self._get_blob_path(blob_hash, compressed=False)
        
        if compressed_path.exists():
            return compressed_path.stat().st_size
        elif uncompressed_path.exists():
            return uncompressed_path.stat().st_size
        else:
            raise BlobNotFoundError(f"Blob not found: {blob_hash}")

    def _compute_hash(self, content: bytes) -> str:
        """Compute SHA-256 hash of content.
        
        Args:
            content: Binary data to hash
            
        Returns:
            Hex string of hash (64 characters for SHA-256)
        """
        hasher = hashlib.sha256()
        hasher.update(content)
        return hasher.hexdigest()

    def _get_blob_path(self, blob_hash: str, compressed: bool = False) -> Path:
        """Get the filesystem path for a blob.
        
        Uses Git-like sharding: objects/<hash[:2]>/<hash[2:]>[.gz]
        
        Args:
            blob_hash: SHA-256 hash (64 hex chars)
            compressed: Whether to add .gz extension
            
        Returns:
            Path to blob file
            
        Example:
            >>> path = store._get_blob_path("abc123...")
            >>> # .chemvcs/objects/ab/c123...
        """
        prefix = blob_hash[:2]
        suffix = blob_hash[2:]
        filename = f"{suffix}.gz" if compressed else suffix
        return self.objects_dir / prefix / filename

    def _validate_hash(self, blob_hash: str) -> None:
        """Validate that a hash string is properly formatted.
        
        Args:
            blob_hash: Hash string to validate
            
        Raises:
            ValueError: If hash is invalid format
        """
        if not isinstance(blob_hash, str):
            raise ValueError(f"Hash must be string, got {type(blob_hash)}")
        
        if len(blob_hash) != HASH_LENGTH:
            raise ValueError(
                f"Hash must be {HASH_LENGTH} characters, got {len(blob_hash)}"
            )
        
        # Check if valid hex
        try:
            int(blob_hash, 16)
        except ValueError as e:
            raise ValueError(f"Hash must be hexadecimal: {e}") from e
