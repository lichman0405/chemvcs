package model

import (
	"crypto/sha256"
	"encoding/hex"
)

// Blob represents raw binary data (e.g., file contents).
// Blobs are the leaves of the Merkle DAG and contain no metadata.
type Blob struct {
	data []byte
}

// NewBlob creates a new Blob with the given data.
func NewBlob(data []byte) *Blob {
	return &Blob{data: data}
}

// Hash computes the SHA-256 hash of the blob's data.
func (b *Blob) Hash() string {
	hash := sha256.Sum256(b.data)
	return hex.EncodeToString(hash[:])
}

// Data returns the blob's data.
func (b *Blob) Data() []byte {
	return b.data
}

// Size returns the size of the blob in bytes.
func (b *Blob) Size() int64 {
	return int64(len(b.data))
}

// ComputeBlobHash computes the SHA-256 hash of arbitrary data.
// This is a utility function for computing blob hashes without creating a Blob object.
func ComputeBlobHash(data []byte) string {
	hash := sha256.Sum256(data)
	return hex.EncodeToString(hash[:])
}
