"""Workflow domain object for representing computational workflows as DAGs."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

from ..core.objects import CoreObject, Reference
from ..util.errors import ValidationError


@dataclass
class WorkflowNode:
    """
    A node in the workflow DAG.
    
    Attributes:
        id: Unique identifier (typically hash of run or sub-workflow)
        kind: Type of node ("run" or "workflow")
        metadata: Additional metadata for this node
    """
    id: str
    kind: str
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate node data."""
        valid_kinds = {"run", "workflow"}
        if self.kind not in valid_kinds:
            raise ValidationError(
                f"Invalid node kind '{self.kind}'. Must be one of: {valid_kinds}"
            )


@dataclass
class Workflow:
    """
    Represents a computational workflow as a directed acyclic graph (DAG).
    
    A workflow consists of nodes (runs or sub-workflows) connected by
    dependencies (edges). It ensures that:
    - No cycles exist (DAG property)
    - All dependencies are satisfied before execution
    
    Attributes:
        name: Human-readable workflow name
        nodes: Dictionary of node_id -> WorkflowNode
        edges: List of (from_node_id, to_node_id) tuples representing dependencies
        metadata: Additional metadata
        id: Optional hash if loaded from repository
    """
    name: str
    nodes: Dict[str, WorkflowNode] = field(default_factory=dict)
    edges: List[Tuple[str, str]] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    id: Optional[str] = None

    def __post_init__(self):
        """Validate workflow data."""
        # Validate name
        if not self.name:
            raise ValidationError("Workflow name cannot be empty")
        
        # Validate edges reference existing nodes
        for from_node, to_node in self.edges:
            if from_node not in self.nodes:
                raise ValidationError(f"Edge references non-existent node: {from_node}")
            if to_node not in self.nodes:
                raise ValidationError(f"Edge references non-existent node: {to_node}")
        
        # Check for cycles
        if self._has_cycle():
            raise ValidationError("Workflow contains a cycle (must be a DAG)")

    def add_node(self, node_id: str, kind: str, metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Add a node to the workflow.
        
        Args:
            node_id: Unique identifier for the node
            kind: Node type ("run" or "workflow")
            metadata: Optional metadata dict
            
        Raises:
            ValidationError: If node_id already exists
        """
        if node_id in self.nodes:
            raise ValidationError(f"Node {node_id} already exists")
        
        self.nodes[node_id] = WorkflowNode(
            id=node_id,
            kind=kind,
            metadata=metadata or {}
        )

    def add_edge(self, from_node: str, to_node: str) -> None:
        """
        Add a dependency edge (from_node must complete before to_node).
        
        Args:
            from_node: ID of dependency node
            to_node: ID of dependent node
            
        Raises:
            ValidationError: If nodes don't exist or edge creates cycle
        """
        if from_node not in self.nodes:
            raise ValidationError(f"Node {from_node} does not exist")
        if to_node not in self.nodes:
            raise ValidationError(f"Node {to_node} does not exist")
        
        # Temporarily add edge to check for cycles
        self.edges.append((from_node, to_node))
        
        if self._has_cycle():
            self.edges.pop()  # Remove the edge
            raise ValidationError(
                f"Adding edge {from_node} -> {to_node} would create a cycle"
            )

    def get_dependencies(self, node_id: str) -> List[str]:
        """
        Get all nodes that must complete before the given node.
        
        Args:
            node_id: Node to get dependencies for
            
        Returns:
            List of node IDs that are dependencies
        """
        return [from_node for from_node, to_node in self.edges if to_node == node_id]

    def get_dependents(self, node_id: str) -> List[str]:
        """
        Get all nodes that depend on the given node.
        
        Args:
            node_id: Node to get dependents for
            
        Returns:
            List of node IDs that depend on this node
        """
        return [to_node for from_node, to_node in self.edges if from_node == node_id]

    def get_root_nodes(self) -> List[str]:
        """
        Get nodes with no dependencies (can be executed immediately).
        
        Returns:
            List of node IDs with no incoming edges
        """
        nodes_with_deps = {to_node for _, to_node in self.edges}
        return [node_id for node_id in self.nodes if node_id not in nodes_with_deps]

    def get_leaf_nodes(self) -> List[str]:
        """
        Get nodes with no dependents (final outputs).
        
        Returns:
            List of node IDs with no outgoing edges
        """
        nodes_with_dependents = {from_node for from_node, _ in self.edges}
        return [node_id for node_id in self.nodes if node_id not in nodes_with_dependents]

    def topological_sort(self) -> List[str]:
        """
        Get topologically sorted list of nodes (valid execution order).
        
        Returns:
            List of node IDs in topological order
            
        Raises:
            ValidationError: If workflow has cycles
        """
        # Build adjacency list and in-degree count
        adj_list = {node_id: [] for node_id in self.nodes}
        in_degree = {node_id: 0 for node_id in self.nodes}
        
        for from_node, to_node in self.edges:
            adj_list[from_node].append(to_node)
            in_degree[to_node] += 1
        
        # Start with nodes that have no dependencies
        queue = [node_id for node_id in self.nodes if in_degree[node_id] == 0]
        result = []
        
        while queue:
            # Process node
            current = queue.pop(0)
            result.append(current)
            
            # Reduce in-degree for dependent nodes
            for neighbor in adj_list[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)
        
        # If we didn't process all nodes, there's a cycle
        if len(result) != len(self.nodes):
            raise ValidationError("Workflow contains a cycle")
        
        return result

    def _has_cycle(self) -> bool:
        """Check if workflow has a cycle using DFS."""
        if not self.nodes:
            return False
        
        # Build adjacency list
        adj_list = {node_id: [] for node_id in self.nodes}
        for from_node, to_node in self.edges:
            adj_list[from_node].append(to_node)
        
        # Track visited and recursion stack
        visited = set()
        rec_stack = set()
        
        def has_cycle_dfs(node: str) -> bool:
            visited.add(node)
            rec_stack.add(node)
            
            for neighbor in adj_list[node]:
                if neighbor not in visited:
                    if has_cycle_dfs(neighbor):
                        return True
                elif neighbor in rec_stack:
                    return True
            
            rec_stack.remove(node)
            return False
        
        # Check all nodes
        for node_id in self.nodes:
            if node_id not in visited:
                if has_cycle_dfs(node_id):
                    return True
        
        return False

    def to_core_object(self) -> CoreObject:
        """
        Convert Workflow to CoreObject for storage.
        
        Returns:
            CoreObject with type="workflow"
        """
        # Serialize nodes
        nodes_data = {
            node_id: {
                "kind": node.kind,
                "metadata": node.metadata
            }
            for node_id, node in self.nodes.items()
        }
        
        # Serialize edges
        edges_data = [[from_node, to_node] for from_node, to_node in self.edges]
        
        meta = {
            "name": self.name,
            "nodes": nodes_data,
            "edges": edges_data,
        }
        
        # Merge additional metadata
        meta.update(self.metadata)
        
        # Create references to all node objects
        refs = [
            Reference(kind="object", id=node_id)
            for node_id in self.nodes.keys()
        ]
        
        obj = CoreObject(
            version=1,
            type="workflow",
            meta=meta,
            refs=refs
        )
        
        return obj

    @classmethod
    def from_core_object(cls, obj: CoreObject) -> "Workflow":
        """
        Create Workflow from CoreObject.
        
        Args:
            obj: CoreObject with type="workflow"
            
        Returns:
            Workflow instance
            
        Raises:
            ValidationError: If object type is not "workflow"
        """
        if obj.type != "workflow":
            raise ValidationError(
                f"Expected object type 'workflow', got '{obj.type}'"
            )
        
        # Extract required fields
        name = obj.get_meta("name")
        nodes_data = obj.get_meta("nodes", {})
        edges_data = obj.get_meta("edges", [])
        
        # Reconstruct nodes
        nodes = {}
        for node_id, node_info in nodes_data.items():
            nodes[node_id] = WorkflowNode(
                id=node_id,
                kind=node_info["kind"],
                metadata=node_info.get("metadata", {})
            )
        
        # Reconstruct edges
        edges = [tuple(edge) for edge in edges_data]
        
        # Extract additional metadata (exclude known fields)
        known_keys = {"name", "nodes", "edges"}
        metadata = {k: v for k, v in obj.meta.items() if k not in known_keys}
        
        return cls(
            name=name,
            nodes=nodes,
            edges=edges,
            metadata=metadata,
            id=obj.hash
        )

    def __repr__(self) -> str:
        id_str = f", id={self.id[:8]}" if self.id else ""
        return (
            f"Workflow('{self.name}', "
            f"nodes={len(self.nodes)}, "
            f"edges={len(self.edges)}{id_str})"
        )
