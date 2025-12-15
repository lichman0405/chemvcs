"""Tests for Workflow domain object."""

import pytest

from chemvcs_py.domain.workflow import Workflow, WorkflowNode
from chemvcs_py.core.objects import CoreObject
from chemvcs_py.util.errors import ValidationError


class TestWorkflowNode:
    """Test cases for WorkflowNode class."""
    
    def test_create_node(self):
        """Test creating a workflow node."""
        node = WorkflowNode(
            id="opt_1",
            kind="run",
            metadata={"description": "Optimization"}
        )
        
        assert node.id == "opt_1"
        assert node.kind == "run"
        assert node.metadata["description"] == "Optimization"
    
    def test_node_repr(self):
        """Test node string representation."""
        node = WorkflowNode(id="node_1", kind="run")
        repr_str = repr(node)
        
        assert "WorkflowNode" in repr_str
        assert "node_1" in repr_str
        assert "run" in repr_str


class TestWorkflow:
    """Test cases for Workflow class."""
    
    def test_create_empty_workflow(self):
        """Test creating an empty workflow."""
        workflow = Workflow(name="Test Workflow")
        
        assert workflow.name == "Test Workflow"
        assert len(workflow.nodes) == 0
        assert len(workflow.edges) == 0
    
    def test_add_single_node(self):
        """Test adding a single node."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run", {"description": "Optimize"})
        
        assert "opt" in workflow.nodes
        assert workflow.nodes["opt"].id == "opt"
        assert workflow.nodes["opt"].kind == "run"
        assert workflow.nodes["opt"].metadata["description"] == "Optimize"
    
    def test_add_multiple_nodes(self):
        """Test adding multiple nodes."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("sp", "run")
        workflow.add_node("freq", "run")
        
        assert len(workflow.nodes) == 3
        assert "opt" in workflow.nodes
        assert "sp" in workflow.nodes
        assert "freq" in workflow.nodes
    
    def test_add_duplicate_node_raises_error(self):
        """Test that adding duplicate node raises error."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        
        with pytest.raises(ValidationError, match="already exists"):
            workflow.add_node("opt", "run")
    
    def test_add_simple_edge(self):
        """Test adding an edge between two nodes."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "freq")
        
        assert ("opt", "freq") in workflow.edges
        assert len(workflow.edges) == 1
    
    def test_add_edge_nonexistent_from_node(self):
        """Test that adding edge with nonexistent from_node raises error."""
        workflow = Workflow(name="Test")
        workflow.add_node("freq", "run")
        
        with pytest.raises(ValidationError, match="does not exist"):
            workflow.add_edge("nonexistent", "freq")
    
    def test_add_edge_nonexistent_to_node(self):
        """Test that adding edge with nonexistent to_node raises error."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        
        with pytest.raises(ValidationError, match="does not exist"):
            workflow.add_edge("opt", "nonexistent")
    
    def test_cycle_detection_simple(self):
        """Test detection of simple cycle."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        
        workflow.add_edge("a", "b")
        
        # Try to create cycle: a -> b -> a
        with pytest.raises(ValidationError, match="cycle"):
            workflow.add_edge("b", "a")
    
    def test_cycle_detection_complex(self):
        """Test detection of cycle in complex DAG."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        workflow.add_node("c", "run")
        workflow.add_node("d", "run")
        
        workflow.add_edge("a", "b")
        workflow.add_edge("b", "c")
        workflow.add_edge("c", "d")
        
        # Try to create cycle: a -> b -> c -> d -> b (cycle)
        with pytest.raises(ValidationError, match="cycle"):
            workflow.add_edge("d", "b")
    
    def test_self_loop_detection(self):
        """Test detection of self-loop."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        
        with pytest.raises(ValidationError, match="cycle"):
            workflow.add_edge("a", "a")
    
    def test_get_dependencies(self):
        """Test getting node dependencies."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("sp", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "freq")
        workflow.add_edge("sp", "freq")
        
        deps = workflow.get_dependencies("freq")
        assert set(deps) == {"opt", "sp"}
    
    def test_get_dependencies_no_deps(self):
        """Test getting dependencies for node with no dependencies."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "freq")
        
        deps = workflow.get_dependencies("opt")
        assert deps == []
    
    def test_get_dependents(self):
        """Test getting node dependents."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("sp", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "sp")
        workflow.add_edge("opt", "freq")
        
        dependents = workflow.get_dependents("opt")
        assert set(dependents) == {"sp", "freq"}
    
    def test_get_dependents_no_dependents(self):
        """Test getting dependents for leaf node."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "freq")
        
        dependents = workflow.get_dependents("freq")
        assert dependents == []
    
    def test_get_root_nodes(self):
        """Test getting root nodes (no dependencies)."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt1", "run")
        workflow.add_node("opt2", "run")
        workflow.add_node("sp", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt1", "sp")
        workflow.add_edge("opt2", "freq")
        
        roots = workflow.get_root_nodes()
        assert set(roots) == {"opt1", "opt2"}
    
    def test_get_leaf_nodes(self):
        """Test getting leaf nodes (no dependents)."""
        workflow = Workflow(name="Test")
        workflow.add_node("opt", "run")
        workflow.add_node("sp", "run")
        workflow.add_node("freq", "run")
        
        workflow.add_edge("opt", "sp")
        workflow.add_edge("opt", "freq")
        
        leaves = workflow.get_leaf_nodes()
        assert set(leaves) == {"sp", "freq"}
    
    def test_topological_sort_simple(self):
        """Test topological sort on simple DAG."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        workflow.add_node("c", "run")
        
        workflow.add_edge("a", "b")
        workflow.add_edge("b", "c")
        
        order = workflow.topological_sort()
        
        # a must come before b, b must come before c
        assert order.index("a") < order.index("b")
        assert order.index("b") < order.index("c")
    
    def test_topological_sort_complex(self):
        """Test topological sort on complex DAG."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        workflow.add_node("c", "run")
        workflow.add_node("d", "run")
        workflow.add_node("e", "run")
        
        workflow.add_edge("a", "c")
        workflow.add_edge("b", "c")
        workflow.add_edge("c", "d")
        workflow.add_edge("c", "e")
        
        order = workflow.topological_sort()
        
        # Check dependencies are satisfied
        assert order.index("a") < order.index("c")
        assert order.index("b") < order.index("c")
        assert order.index("c") < order.index("d")
        assert order.index("c") < order.index("e")
    
    def test_topological_sort_empty(self):
        """Test topological sort on empty workflow."""
        workflow = Workflow(name="Test")
        
        order = workflow.topological_sort()
        assert order == []
    
    def test_topological_sort_single_node(self):
        """Test topological sort with single node."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        
        order = workflow.topological_sort()
        assert order == ["a"]
    
    def test_topological_sort_disconnected_components(self):
        """Test topological sort with disconnected components."""
        workflow = Workflow(name="Test")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        workflow.add_node("c", "run")
        workflow.add_node("d", "run")
        
        workflow.add_edge("a", "b")
        workflow.add_edge("c", "d")
        
        order = workflow.topological_sort()
        
        # Check each component is ordered correctly
        assert order.index("a") < order.index("b")
        assert order.index("c") < order.index("d")
        assert len(order) == 4
    
    def test_to_core_object(self):
        """Test conversion to CoreObject."""
        workflow = Workflow(name="Test Workflow")
        workflow.add_node("opt", "run", {"description": "Optimize"})
        workflow.add_node("freq", "run", {"description": "Frequency"})
        workflow.add_edge("opt", "freq")
        
        core_obj = workflow.to_core_object()
        
        assert isinstance(core_obj, CoreObject)
        assert core_obj.type == "workflow"
        assert core_obj.meta["name"] == "Test Workflow"
        assert "nodes" in core_obj.meta
        assert "edges" in core_obj.meta
        
        # Check nodes are serialized
        nodes = core_obj.meta["nodes"]
        assert "opt" in nodes
        assert nodes["opt"]["kind"] == "run"
        
        # Check edges are serialized
        edges = core_obj.meta["edges"]
        assert ["opt", "freq"] in edges
    
    def test_from_core_object(self):
        """Test reconstruction from CoreObject."""
        original = Workflow(name="Original Workflow")
        original.add_node("a", "run", {"desc": "Node A"})
        original.add_node("b", "run", {"desc": "Node B"})
        original.add_edge("a", "b")
        
        core_obj = original.to_core_object()
        reconstructed = Workflow.from_core_object(core_obj)
        
        assert reconstructed.name == original.name
        assert len(reconstructed.nodes) == len(original.nodes)
        assert len(reconstructed.edges) == len(original.edges)
        assert "a" in reconstructed.nodes
        assert "b" in reconstructed.nodes
        assert ("a", "b") in reconstructed.edges
    
    def test_round_trip_conversion(self):
        """Test that to/from CoreObject is lossless."""
        original = Workflow(
            name="Complex Workflow",
            metadata={"project": "test"}
        )
        original.add_node("opt1", "run", {"conformer": 0})
        original.add_node("opt2", "run", {"conformer": 1})
        original.add_node("analysis", "run", {"task": "compare"})
        
        original.add_edge("opt1", "analysis")
        original.add_edge("opt2", "analysis")
        
        core_obj = original.to_core_object()
        reconstructed = Workflow.from_core_object(core_obj)
        
        # Check all attributes match
        assert reconstructed.name == original.name
        assert len(reconstructed.nodes) == len(original.nodes)
        assert len(reconstructed.edges) == len(original.edges)
        
        # Check nodes
        for node_id in original.nodes:
            assert node_id in reconstructed.nodes
            assert reconstructed.nodes[node_id].kind == original.nodes[node_id].kind
        
        # Check edges
        assert set(reconstructed.edges) == set(original.edges)
        
        # Check metadata
        assert reconstructed.metadata == original.metadata
    
    def test_repr(self):
        """Test string representation."""
        workflow = Workflow(name="Test Workflow")
        workflow.add_node("a", "run")
        workflow.add_node("b", "run")
        workflow.add_edge("a", "b")
        
        repr_str = repr(workflow)
        assert "Workflow" in repr_str
        assert "Test Workflow" in repr_str
        assert "nodes=2" in repr_str
        assert "edges=1" in repr_str
