"""Core object representations matching Go's Object and Reference types."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class Reference:
    """Reference to another entity (Object or Blob)."""
    kind: str  # "object" or "blob"
    id: str    # SHA-256 hash

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Reference":
        """Create Reference from dictionary (parsed JSON)."""
        return cls(
            kind=data["kind"],
            id=data["id"]
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "kind": self.kind,
            "id": self.id
        }

    def is_object(self) -> bool:
        """Check if this is an object reference."""
        return self.kind == "object"

    def is_blob(self) -> bool:
        """Check if this is a blob reference."""
        return self.kind == "blob"


@dataclass
class CoreObject:
    """
    Core object representation matching Go's Object struct.
    
    This is a Python representation of the core VCS object,
    providing a bridge between Go's storage layer and Python's
    domain layer.
    """
    version: int
    type: str
    meta: Dict[str, Any] = field(default_factory=dict)
    refs: List[Reference] = field(default_factory=list)
    hash: Optional[str] = None  # Set when loaded from repository

    @classmethod
    def from_dict(cls, data: Dict[str, Any], hash: Optional[str] = None) -> "CoreObject":
        """
        Create CoreObject from dictionary (parsed JSON).
        
        Args:
            data: Dictionary containing object data
            hash: Optional hash of the object
        """
        refs = [Reference.from_dict(r) for r in data.get("refs", [])]
        return cls(
            version=data["version"],
            type=data["type"],
            meta=data.get("meta", {}),
            refs=refs,
            hash=hash
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            "version": self.version,
            "type": self.type,
        }
        if self.meta:
            result["meta"] = self.meta
        if self.refs:
            result["refs"] = [r.to_dict() for r in self.refs]
        return result

    def get_meta(self, key: str, default: Any = None) -> Any:
        """Get metadata value with optional default."""
        return self.meta.get(key, default)

    def set_meta(self, key: str, value: Any) -> None:
        """Set metadata value."""
        self.meta[key] = value

    def add_ref(self, ref: Reference) -> None:
        """Add a reference to this object."""
        self.refs.append(ref)

    def get_refs_by_kind(self, kind: str) -> List[Reference]:
        """Get all references of a specific kind."""
        return [r for r in self.refs if r.kind == kind]

    def get_object_refs(self) -> List[Reference]:
        """Get all object references."""
        return self.get_refs_by_kind("object")

    def get_blob_refs(self) -> List[Reference]:
        """Get all blob references."""
        return self.get_refs_by_kind("blob")
