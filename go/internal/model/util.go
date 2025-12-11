package model

import (
	"encoding/json"
	"fmt"
)

// DecodeJSON is a generic helper for decoding JSON data into a target struct.
func DecodeJSON(data []byte, target interface{}) error {
	if err := json.Unmarshal(data, target); err != nil {
		return fmt.Errorf("failed to unmarshal JSON: %w", err)
	}
	return nil
}
