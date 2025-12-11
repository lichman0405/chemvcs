package model

// Reference represents a reference to another entity (Object or Blob).
// References are used to build the Merkle DAG structure.
type Reference struct {
	Kind string `json:"kind"` // "object" or "blob"
	ID   string `json:"id"`   // SHA-256 hash in hexadecimal
}

// NewObjectRef creates a reference to an Object.
func NewObjectRef(id string) Reference {
	return Reference{
		Kind: "object",
		ID:   id,
	}
}

// NewBlobRef creates a reference to a Blob.
func NewBlobRef(id string) Reference {
	return Reference{
		Kind: "blob",
		ID:   id,
	}
}

// IsObject returns true if this reference points to an Object.
func (r Reference) IsObject() bool {
	return r.Kind == "object"
}

// IsBlob returns true if this reference points to a Blob.
func (r Reference) IsBlob() bool {
	return r.Kind == "blob"
}
