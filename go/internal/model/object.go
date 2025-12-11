package model

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"sort"
)

// Object represents a typed node in the Merkle DAG with metadata and references.
type Object struct {
	Version int                    `json:"version"`
	Type    string                 `json:"type"`
	Meta    map[string]interface{} `json:"meta,omitempty"`
	Refs    []Reference            `json:"refs,omitempty"`
}

// NewObject creates a new Object with the specified type.
func NewObject(objType string) *Object {
	return &Object{
		Version: 1,
		Type:    objType,
		Meta:    make(map[string]interface{}),
		Refs:    []Reference{},
	}
}

// Hash computes the SHA-256 hash of this Object's canonical JSON representation.
// This ensures deterministic hashing regardless of construction order.
func (o *Object) Hash() (string, error) {
	canonical, err := o.MarshalCanonical()
	if err != nil {
		return "", fmt.Errorf("failed to marshal object: %w", err)
	}

	hash := sha256.Sum256(canonical)
	return hex.EncodeToString(hash[:]), nil
}

// MarshalCanonical produces a canonical JSON representation of the Object.
// Keys are sorted lexicographically to ensure deterministic output.
func (o *Object) MarshalCanonical() ([]byte, error) {
	// Create a temporary structure with ordered fields
	type orderedObject struct {
		Version int                    `json:"version"`
		Type    string                 `json:"type"`
		Meta    map[string]interface{} `json:"meta,omitempty"`
		Refs    []Reference            `json:"refs,omitempty"`
	}

	ordered := orderedObject{
		Version: o.Version,
		Type:    o.Type,
		Meta:    o.Meta,
		Refs:    o.Refs,
	}

	// Marshal to JSON with sorted keys
	data, err := json.Marshal(ordered)
	if err != nil {
		return nil, err
	}

	// Unmarshal and re-marshal to ensure key ordering
	// This is necessary because Go's json.Marshal doesn't guarantee map key order
	var intermediate interface{}
	if err := json.Unmarshal(data, &intermediate); err != nil {
		return nil, err
	}

	return marshalSorted(intermediate)
}

// marshalSorted recursively marshals data with sorted map keys.
func marshalSorted(data interface{}) ([]byte, error) {
	switch v := data.(type) {
	case map[string]interface{}:
		// Sort keys and build ordered structure
		keys := make([]string, 0, len(v))
		for key := range v {
			keys = append(keys, key)
		}
		sort.Strings(keys)

		// Build a slice of key-value pairs in sorted order
		type kvPair struct {
			Key   string
			Value interface{}
		}
		pairs := make([]kvPair, len(keys))
		for i, key := range keys {
			pairs[i] = kvPair{Key: key, Value: v[key]}
		}

		// Manually build JSON to ensure key order
		result := "{"
		for i, pair := range pairs {
			if i > 0 {
				result += ","
			}
			keyJSON, _ := json.Marshal(pair.Key)
			result += string(keyJSON) + ":"

			valueJSON, err := marshalSorted(pair.Value)
			if err != nil {
				return nil, err
			}
			result += string(valueJSON)
		}
		result += "}"
		return []byte(result), nil

	case []interface{}:
		// Handle arrays
		result := "["
		for i, item := range v {
			if i > 0 {
				result += ","
			}
			itemJSON, err := marshalSorted(item)
			if err != nil {
				return nil, err
			}
			result += string(itemJSON)
		}
		result += "]"
		return []byte(result), nil

	default:
		// Primitive types
		return json.Marshal(v)
	}
}

// AddRef adds a reference to this Object.
func (o *Object) AddRef(ref Reference) {
	o.Refs = append(o.Refs, ref)
}

// SetMeta sets a metadata value.
func (o *Object) SetMeta(key string, value interface{}) {
	if o.Meta == nil {
		o.Meta = make(map[string]interface{})
	}
	o.Meta[key] = value
}

// GetMeta retrieves a metadata value.
func (o *Object) GetMeta(key string) (interface{}, bool) {
	if o.Meta == nil {
		return nil, false
	}
	val, ok := o.Meta[key]
	return val, ok
}
