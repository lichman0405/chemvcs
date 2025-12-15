"""
Example demonstrating Run and Workflow usage with ChemVCS.

This example shows:
1. Creating computational run objects
2. Building multi-step workflows as DAGs
3. Tracking calculation status
4. Storing runs and workflows in repository
"""

import numpy as np
from chemvcs_py import Structure, Run, Workflow, WorkflowNode, Repo
from chemvcs_py.io import write_xyz

print("=" * 70)
print("ChemVCS Run and Workflow Example")
print("=" * 70)

# 1. Create a molecular structure
print("\n1. Creating H2O structure...")
water = Structure(
    formula="H2O",
    positions=np.array([
        [0.0, 0.0, 0.0],     # O
        [0.0, 0.757, 0.586],  # H
        [0.0, -0.757, 0.586]  # H
    ]),
    species=["O", "H", "H"],
    metadata={"method": "manual"}
)
print(f"   Created: {water}")

# 2. Create computational runs
print("\n2. Creating computational runs...")

# Geometry optimization run
opt_run = Run(
    structure_id="struct_001",  # Would be actual hash in real usage
    code="ORCA",
    code_version="5.0.3",
    parameters={
        "method": "DFT",
        "functional": "B3LYP",
        "basis": "def2-TZVP",
        "task": "opt"
    },
    status="planned",
    metadata={"description": "Geometry optimization"}
)
print(f"   Optimization: {opt_run}")

# Single point energy run
sp_run = Run(
    structure_id="struct_001",
    code="ORCA",
    code_version="5.0.3",
    parameters={
        "method": "DFT",
        "functional": "B3LYP",
        "basis": "def2-TZVP",
        "task": "sp"
    },
    status="planned",
    metadata={"description": "Single point energy"}
)
print(f"   Single point: {sp_run}")

# Frequency calculation run
freq_run = Run(
    structure_id="struct_002",  # Uses optimized structure
    code="ORCA",
    code_version="5.0.3",
    parameters={
        "method": "DFT",
        "functional": "B3LYP",
        "basis": "def2-TZVP",
        "task": "freq"
    },
    status="planned",
    metadata={"description": "Frequency calculation"}
)
print(f"   Frequency: {freq_run}")

# 3. Simulate run lifecycle
print("\n3. Simulating run lifecycle...")
print(f"   Initial status: {opt_run.status}")

opt_run.mark_submitted(job_id="job_12345")
print(f"   After submission: {opt_run.status} (job_id: {opt_run.metadata.get('job_id')})")

opt_run.mark_running()
print(f"   After start: {opt_run.status}")

opt_run.mark_finished()
opt_run.set_result("total_energy", -76.4321)
opt_run.set_result("optimized_structure_id", "struct_002")
print(f"   After completion: {opt_run.status}")
print(f"   Energy: {opt_run.get_energy()} Hartree")

# 4. Build a workflow DAG
print("\n4. Building computational workflow...")
workflow = Workflow(
    name="Water Calculation Workflow",
    metadata={"project": "H2O_study", "version": "1.0"}
)

# Add nodes
workflow.add_node("opt", "run", {"run_id": "run_001", "description": "Optimize geometry"})
workflow.add_node("sp", "run", {"run_id": "run_002", "description": "Single point energy"})
workflow.add_node("freq", "run", {"run_id": "run_003", "description": "Vibrational analysis"})

# Add dependencies: opt -> freq, opt -> sp
workflow.add_edge("opt", "freq")
workflow.add_edge("opt", "sp")

print(f"   Created workflow: {workflow.name}")
print(f"   Nodes: {list(workflow.nodes.keys())}")
print(f"   Edges: {workflow.edges}")

# 5. Query workflow structure
print("\n5. Analyzing workflow structure...")
root_nodes = workflow.get_root_nodes()
print(f"   Root nodes (no dependencies): {root_nodes}")

leaf_nodes = workflow.get_leaf_nodes()
print(f"   Leaf nodes (no dependents): {leaf_nodes}")

freq_deps = workflow.get_dependencies("freq")
print(f"   'freq' depends on: {freq_deps}")

opt_dependents = workflow.get_dependents("opt")
print(f"   'opt' is required by: {opt_dependents}")

# 6. Get execution order
print("\n6. Determining execution order...")
exec_order = workflow.topological_sort()
print(f"   Execution sequence: {' -> '.join(exec_order)}")

# 7. Convert to CoreObject for storage
print("\n7. Converting to CoreObjects...")
struct_obj = water.to_core_object()
opt_run_obj = opt_run.to_core_object()
workflow_obj = workflow.to_core_object()

print(f"   Structure object: type={struct_obj.type}, version={struct_obj.version}")
print(f"   Run object: type={opt_run_obj.type}, version={opt_run_obj.version}")
print(f"   Workflow object: type={workflow_obj.type}, version={workflow_obj.version}")

# 8. Demonstrate error handling
print("\n8. Testing workflow validation...")
try:
    # Try to create a cycle
    workflow.add_edge("freq", "opt")
    print("   ERROR: Cycle not detected!")
except Exception as e:
    print(f"   ✓ Cycle prevention works: {e}")

# 9. Complex workflow example
print("\n9. Creating complex multi-step workflow...")
complex_wf = Workflow(name="Multi-conformer Study")

# Add conformer search nodes
for i in range(3):
    complex_wf.add_node(f"conf_{i}", "run", {"conformer": i})
    complex_wf.add_node(f"opt_{i}", "run", {"conformer": i})
    complex_wf.add_node(f"energy_{i}", "run", {"conformer": i})
    
    # Dependencies: conf -> opt -> energy
    complex_wf.add_edge(f"conf_{i}", f"opt_{i}")
    complex_wf.add_edge(f"opt_{i}", f"energy_{i}")

# Add analysis node that depends on all energies
complex_wf.add_node("analysis", "run", {"task": "Boltzmann weighting"})
for i in range(3):
    complex_wf.add_edge(f"energy_{i}", "analysis")

print(f"   Complex workflow: {complex_wf.name}")
print(f"   Total nodes: {len(complex_wf.nodes)}")
print(f"   Total edges: {len(complex_wf.edges)}")
print(f"   Execution order has {len(complex_wf.topological_sort())} steps")

# 10. Try repository integration (if available)
print("\n10. Repository integration...")
try:
    repo = Repo()
    print(f"   ✓ Found repository at: {repo.repo_root}")
    print(f"   ✓ ChemVCS executable: {repo.chemvcs_exe}")
    
    # Could store objects here:
    # struct_hash = repo.store_object(struct_obj)
    # run_hash = repo.store_object(opt_run_obj)
    # workflow_hash = repo.store_object(workflow_obj)
    
except Exception as e:
    print(f"   ℹ No repository found (this is OK for demo): {e}")

print("\n" + "=" * 70)
print("Example completed successfully!")
print("=" * 70)
print("""
Summary:
- Created Structure, Run, and Workflow domain objects
- Demonstrated run lifecycle (planned → submitted → running → finished)
- Built workflow DAG with dependency management
- Validated workflow structure (cycle detection)
- Converted objects to CoreObject for storage
- Showed complex multi-conformer workflow pattern

Next steps:
- Initialize a ChemVCS repository: cd /path/to/project && chemvcs init
- Store objects using repo.store_object(core_obj)
- Query stored calculations: repo.list_objects(type_filter="run")
""")
