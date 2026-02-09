"""Unit tests for ObjectStore."""

import gzip
import os
from pathlib import Path

import pytest

from chemvcs.storage.object_store import (
    BlobCorruptedError,
    BlobNotFoundError,
    ObjectStore,
)
from chemvcs.constants import GZIP_THRESHOLD


@pytest.fixture
def chemvcs_dir(tmp_path: Path) -> Path:
    """Create a temporary .chemvcs directory structure."""
    chemvcs = tmp_path / ".chemvcs"
    chemvcs.mkdir()
    (chemvcs / "objects").mkdir()
    (chemvcs / "commits").mkdir()
    return chemvcs


@pytest.fixture
def store(chemvcs_dir: Path) -> ObjectStore:
    """Create an ObjectStore instance."""
    return ObjectStore(chemvcs_dir)


class TestObjectStoreInit:
    """Test ObjectStore initialization."""

    def test_init_with_valid_dir(self, chemvcs_dir: Path) -> None:
        """Test initialization with valid .chemvcs directory."""
        store = ObjectStore(chemvcs_dir)
        assert store.chemvcs_dir == chemvcs_dir
        assert store.objects_dir == chemvcs_dir / "objects"

    def test_init_with_nonexistent_dir(self, tmp_path: Path) -> None:
        """Test initialization with non-existent directory fails."""
        nonexistent = tmp_path / "nonexistent"
        with pytest.raises(ValueError, match="not found"):
            ObjectStore(nonexistent)


class TestWriteBlob:
    """Test blob writing functionality."""

    def test_write_blob_basic(self, store: ObjectStore) -> None:
        """Test basic blob writing."""
        content = b"ENCUT = 520\nISMEAR = 0\n"
        blob_hash = store.write_blob(content)
        
        # Hash should be 64 hex characters (SHA-256)
        assert len(blob_hash) == 64
        assert all(c in "0123456789abcdef" for c in blob_hash)
        
        # Blob should exist
        assert store.blob_exists(blob_hash)

    def test_write_blob_deduplication(self, store: ObjectStore) -> None:
        """Test that identical content produces same hash and doesn't duplicate."""
        content = b"Test data for deduplication"
        
        hash1 = store.write_blob(content)
        hash2 = store.write_blob(content)
        
        # Hashes should be identical
        assert hash1 == hash2
        
        # Should only have one blob file on disk
        blob_path = store._get_blob_path(hash1)
        assert blob_path.exists()

    def test_write_blob_different_content(self, store: ObjectStore) -> None:
        """Test that different content produces different hashes."""
        content1 = b"First content"
        content2 = b"Second content"
        
        hash1 = store.write_blob(content1)
        hash2 = store.write_blob(content2)
        
        assert hash1 != hash2

    def test_write_blob_creates_sharded_directory(self, store: ObjectStore) -> None:
        """Test that blob storage uses 2-char sharding."""
        content = b"Test sharding"
        blob_hash = store.write_blob(content)
        
        # Check directory structure: objects/<first2chars>/<rest>
        shard_dir = store.objects_dir / blob_hash[:2]
        assert shard_dir.exists()
        assert shard_dir.is_dir()
        
        blob_file = shard_dir / blob_hash[2:]
        assert blob_file.exists()

    def test_write_blob_small_no_compression(self, store: ObjectStore) -> None:
        """Test that small files are not compressed by default."""
        content = b"Small file content"
        blob_hash = store.write_blob(content)
        
        # Should exist as uncompressed
        uncompressed_path = store._get_blob_path(blob_hash, compressed=False)
        assert uncompressed_path.exists()
        
        # Should NOT exist as compressed
        compressed_path = store._get_blob_path(blob_hash, compressed=True)
        assert not compressed_path.exists()

    def test_write_blob_large_auto_compression(self, store: ObjectStore) -> None:
        """Test that large files are automatically compressed."""
        # Create content larger than GZIP_THRESHOLD
        content = b"X" * (GZIP_THRESHOLD + 1000)
        blob_hash = store.write_blob(content)
        
        # Should exist as compressed
        compressed_path = store._get_blob_path(blob_hash, compressed=True)
        assert compressed_path.exists()
        
        # Should NOT exist as uncompressed
        uncompressed_path = store._get_blob_path(blob_hash, compressed=False)
        assert not uncompressed_path.exists()

    def test_write_blob_force_compression(self, store: ObjectStore) -> None:
        """Test forcing compression on small file."""
        content = b"Small but compressed"
        blob_hash = store.write_blob(content, compress=True)
        
        compressed_path = store._get_blob_path(blob_hash, compressed=True)
        assert compressed_path.exists()

    def test_write_blob_disable_compression(self, store: ObjectStore) -> None:
        """Test disabling compression on large file."""
        content = b"Y" * (GZIP_THRESHOLD + 1000)
        blob_hash = store.write_blob(content, compress=False)
        
        uncompressed_path = store._get_blob_path(blob_hash, compressed=False)
        assert uncompressed_path.exists()

    def test_write_blob_empty_content(self, store: ObjectStore) -> None:
        """Test writing empty blob."""
        content = b""
        blob_hash = store.write_blob(content)
        
        assert len(blob_hash) == 64
        assert store.blob_exists(blob_hash)
        
        # Should be able to read it back
        retrieved = store.read_blob(blob_hash)
        assert retrieved == b""


class TestReadBlob:
    """Test blob reading functionality."""

    def test_read_blob_basic(self, store: ObjectStore) -> None:
        """Test basic blob reading."""
        content = b"POSCAR content\nLi4 Co4 O8\n"
        blob_hash = store.write_blob(content)
        
        retrieved = store.read_blob(blob_hash)
        assert retrieved == content

    def test_read_blob_compressed(self, store: ObjectStore) -> None:
        """Test reading compressed blob."""
        content = b"Compressed test data" * 100
        blob_hash = store.write_blob(content, compress=True)
        
        retrieved = store.read_blob(blob_hash)
        assert retrieved == content

    def test_read_blob_with_hash_verification(self, store: ObjectStore) -> None:
        """Test that hash verification is performed by default."""
        content = b"Verify this hash"
        blob_hash = store.write_blob(content)
        
        # Reading should succeed with correct hash
        retrieved = store.read_blob(blob_hash, verify_hash=True)
        assert retrieved == content

    def test_read_blob_skip_verification(self, store: ObjectStore) -> None:
        """Test reading without hash verification."""
        content = b"No verification needed"
        blob_hash = store.write_blob(content)
        
        retrieved = store.read_blob(blob_hash, verify_hash=False)
        assert retrieved == content

    def test_read_blob_not_found(self, store: ObjectStore) -> None:
        """Test reading non-existent blob raises error."""
        fake_hash = "a" * 64
        
        with pytest.raises(BlobNotFoundError, match="Blob not found"):
            store.read_blob(fake_hash)

    def test_read_blob_corrupted(self, store: ObjectStore) -> None:
        """Test reading corrupted blob raises error."""
        content = b"Original content"
        blob_hash = store.write_blob(content)
        
        # Corrupt the blob by modifying file content
        blob_path = store._get_blob_path(blob_hash)
        with open(blob_path, "wb") as f:
            f.write(b"Corrupted content!")
        
        with pytest.raises(BlobCorruptedError, match="Blob corrupted"):
            store.read_blob(blob_hash, verify_hash=True)

    def test_read_blob_invalid_hash_format(self, store: ObjectStore) -> None:
        """Test reading with invalid hash format."""
        with pytest.raises(ValueError, match="must be 64 characters"):
            store.read_blob("tooshort")
        
        with pytest.raises(ValueError, match="hexadecimal"):
            store.read_blob("z" * 64)


class TestBlobExists:
    """Test blob existence checking."""

    def test_blob_exists_true(self, store: ObjectStore) -> None:
        """Test blob_exists returns True for existing blob."""
        content = b"Exists test"
        blob_hash = store.write_blob(content)
        
        assert store.blob_exists(blob_hash) is True

    def test_blob_exists_false(self, store: ObjectStore) -> None:
        """Test blob_exists returns False for non-existent blob."""
        fake_hash = "b" * 64
        assert store.blob_exists(fake_hash) is False

    def test_blob_exists_compressed(self, store: ObjectStore) -> None:
        """Test blob_exists returns True for compressed blob."""
        content = b"Compressed" * 1000
        blob_hash = store.write_blob(content, compress=True)
        
        assert store.blob_exists(blob_hash) is True

    def test_blob_exists_invalid_hash(self, store: ObjectStore) -> None:
        """Test blob_exists with invalid hash returns False."""
        assert store.blob_exists("invalid") is False
        assert store.blob_exists("") is False


class TestGetBlobSize:
    """Test blob size retrieval."""

    def test_get_blob_size_uncompressed(self, store: ObjectStore) -> None:
        """Test getting size of uncompressed blob."""
        content = b"Size test content"
        blob_hash = store.write_blob(content, compress=False)
        
        size = store.get_blob_size(blob_hash)
        assert size == len(content)

    def test_get_blob_size_compressed(self, store: ObjectStore) -> None:
        """Test getting size of compressed blob (returns compressed size)."""
        content = b"A" * 10000  # Highly compressible
        blob_hash = store.write_blob(content, compress=True)
        
        size = store.get_blob_size(blob_hash)
        # Compressed size should be much smaller than original
        assert size < len(content)
        
        # Verify it's actually the on-disk size
        blob_path = store._get_blob_path(blob_hash, compressed=True)
        assert size == blob_path.stat().st_size

    def test_get_blob_size_not_found(self, store: ObjectStore) -> None:
        """Test getting size of non-existent blob."""
        fake_hash = "c" * 64
        
        with pytest.raises(BlobNotFoundError):
            store.get_blob_size(fake_hash)


class TestComputeHash:
    """Test hash computation."""

    def test_compute_hash_deterministic(self, store: ObjectStore) -> None:
        """Test that same content produces same hash."""
        content = b"Deterministic hash test"
        
        hash1 = store._compute_hash(content)
        hash2 = store._compute_hash(content)
        
        assert hash1 == hash2

    def test_compute_hash_different_content(self, store: ObjectStore) -> None:
        """Test that different content produces different hash."""
        content1 = b"Content A"
        content2 = b"Content B"
        
        hash1 = store._compute_hash(content1)
        hash2 = store._compute_hash(content2)
        
        assert hash1 != hash2

    def test_compute_hash_length(self, store: ObjectStore) -> None:
        """Test that hash is correct length (SHA-256 = 64 hex chars)."""
        content = b"Hash length test"
        blob_hash = store._compute_hash(content)
        
        assert len(blob_hash) == 64


class TestGetBlobPath:
    """Test blob path generation."""

    def test_get_blob_path_uncompressed(self, store: ObjectStore) -> None:
        """Test path generation for uncompressed blob."""
        blob_hash = "abc123" + "d" * 58  # 64 chars total
        path = store._get_blob_path(blob_hash, compressed=False)
        
        assert path == store.objects_dir / "ab" / ("c123" + "d" * 58)

    def test_get_blob_path_compressed(self, store: ObjectStore) -> None:
        """Test path generation for compressed blob."""
        blob_hash = "xyz789" + "e" * 58
        path = store._get_blob_path(blob_hash, compressed=True)
        
        assert path == store.objects_dir / "xy" / ("z789" + "e" * 58 + ".gz")


class TestValidateHash:
    """Test hash validation."""

    def test_validate_hash_valid(self, store: ObjectStore) -> None:
        """Test validation of valid hash."""
        valid_hash = "a" * 64
        store._validate_hash(valid_hash)  # Should not raise

    def test_validate_hash_wrong_length(self, store: ObjectStore) -> None:
        """Test validation fails for wrong length."""
        with pytest.raises(ValueError, match="must be 64 characters"):
            store._validate_hash("abc")

    def test_validate_hash_not_hex(self, store: ObjectStore) -> None:
        """Test validation fails for non-hexadecimal."""
        with pytest.raises(ValueError, match="hexadecimal"):
            store._validate_hash("g" * 64)

    def test_validate_hash_not_string(self, store: ObjectStore) -> None:
        """Test validation fails for non-string."""
        with pytest.raises(ValueError, match="must be string"):
            store._validate_hash(123)  # type: ignore


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_concurrent_write_same_blob(self, store: ObjectStore) -> None:
        """Test that concurrent writes of same blob are handled gracefully."""
        content = b"Concurrent test"
        
        # Simulate race condition by writing same content twice
        hash1 = store.write_blob(content)
        hash2 = store.write_blob(content)
        
        # Should be same hash and no errors
        assert hash1 == hash2
        assert store.blob_exists(hash1)

    def test_binary_content_with_null_bytes(self, store: ObjectStore) -> None:
        """Test handling of binary content with null bytes."""
        content = b"\x00\x01\x02\xff\xfe\xfd"
        blob_hash = store.write_blob(content)
        
        retrieved = store.read_blob(blob_hash)
        assert retrieved == content

    def test_large_content_roundtrip(self, store: ObjectStore) -> None:
        """Test writing and reading large content (>200MB)."""
        # Create 250MB of data
        content = b"X" * (250 * 1024 * 1024)
        
        blob_hash = store.write_blob(content)
        retrieved = store.read_blob(blob_hash)
        
        assert retrieved == content
        
        # Should be compressed
        compressed_path = store._get_blob_path(blob_hash, compressed=True)
        assert compressed_path.exists()
        
        # Compressed size should be much smaller
        compressed_size = store.get_blob_size(blob_hash)
        assert compressed_size < len(content) / 10  # At least 10x compression
