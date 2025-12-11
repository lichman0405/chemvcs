package model

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"time"
)

// Snapshot represents a versioned state of the project at a point in time.
// Similar to a git commit, it captures the root object, parent snapshots,
// author, timestamp, and commit message.
type Snapshot struct {
	Version   int      `json:"version"`
	Root      string   `json:"root"`      // Hash of root Object
	Parents   []string `json:"parents"`   // Parent snapshot hashes
	Author    string   `json:"author"`    // Author name and email
	Timestamp string   `json:"timestamp"` // RFC3339 format
	Message   string   `json:"message"`   // Commit message
}

// NewSnapshot creates a new Snapshot with the given parameters.
func NewSnapshot(root, author, message string, parents []string) *Snapshot {
	if parents == nil {
		parents = []string{}
	}

	return &Snapshot{
		Version:   1,
		Root:      root,
		Parents:   parents,
		Author:    author,
		Timestamp: time.Now().UTC().Format(time.RFC3339),
		Message:   message,
	}
}

// Hash computes the SHA-256 hash of this Snapshot's canonical JSON representation.
func (s *Snapshot) Hash() (string, error) {
	canonical, err := s.MarshalCanonical()
	if err != nil {
		return "", fmt.Errorf("failed to marshal snapshot: %w", err)
	}

	hash := sha256.Sum256(canonical)
	return hex.EncodeToString(hash[:]), nil
}

// MarshalCanonical produces a canonical JSON representation of the Snapshot.
func (s *Snapshot) MarshalCanonical() ([]byte, error) {
	// Ensure deterministic output
	// Since Snapshot has no maps, standard JSON marshalling is sufficient
	// with keys in struct order
	return json.Marshal(s)
}

// UnmarshalSnapshot deserialises a Snapshot from JSON.
func UnmarshalSnapshot(data []byte) (*Snapshot, error) {
	var s Snapshot
	if err := json.Unmarshal(data, &s); err != nil {
		return nil, fmt.Errorf("failed to unmarshal snapshot: %w", err)
	}

	// Validate version
	if s.Version != 1 {
		return nil, fmt.Errorf("unsupported snapshot version: %d", s.Version)
	}

	return &s, nil
}

// ShortHash returns the first 8 characters of the snapshot hash for display.
func (s *Snapshot) ShortHash() string {
	hash, err := s.Hash()
	if err != nil {
		return "????????"
	}
	if len(hash) < 8 {
		return hash
	}
	return hash[:8]
}

// IsInitial returns true if this snapshot has no parents (initial commit).
func (s *Snapshot) IsInitial() bool {
	return len(s.Parents) == 0
}
