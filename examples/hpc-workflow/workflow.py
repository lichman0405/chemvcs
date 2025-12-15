"""
Complete example workflow demonstrating ChemVCS HPC integration.

This script shows how to use ChemVCS to manage a computational chemistry
workflow on an HPC cluster with SLURM. It demonstrates:

1. Repository initialization
2. Structure versioning
3. HPC job submission with provenance tracking
4. Job status monitoring
5. Result retrieval
6. Full workflow composition

Run this script in a directory containing:
- water.xyz: Input structure file
- vasp_relax.slurm: Geometry optimization submission script
- vasp_static.slurm: Single-point calculation submission script
"""

import time
from pathlib import Path
from chemvcs_py.api import Repo
from chemvcs_py.domain import Structure, Run, Workflow, RunStatus
from chemvcs_py.hpc import (
    SlurmAdapter,
    JobSubmitter,
    JobTracker,
    JobRetriever
)


def main():
    # Initialize repository
    print("=== Initializing ChemVCS repository ===")
    repo = Repo.init(".")
    print(f"Repository initialized at: {repo.repo_path}")
    
    # Step 1: Add water structure
    print("\n=== Step 1: Adding water structure ===")
    water_xyz = Path("water.xyz")
    if not water_xyz.exists():
        print("ERROR: water.xyz not found!")
        return
    
    structure = Structure(
        name="water",
        atoms=["O", "H", "H"],
        coordinates=[
            [0.0, 0.0, 0.119262],
            [0.0, 0.763239, -0.477047],
            [0.0, -0.763239, -0.477047]
        ],
        metadata={"source": "water.xyz", "description": "H2O molecule"}
    )
    
    struct_hash = repo.add_structure(structure)
    print(f"Structure added: {struct_hash[:8]}")
    
    # Step 2: Submit geometry optimization job
    print("\n=== Step 2: Submitting geometry optimization ===")
    
    # Initialize HPC components
    slurm = SlurmAdapter()
    submitter = JobSubmitter(repo, slurm)
    tracker = JobTracker(repo, slurm)
    retriever = JobRetriever(repo, slurm)
    
    # Create run for geometry optimization
    relax_run = Run(
        name="water-relax",
        status=RunStatus.PLANNED,
        structure_hash=struct_hash,
        command="mpirun -np 56 vasp_std",
        working_dir="/scratch/user/water-relax",
        metadata={
            "calculation_type": "geometry_optimization",
            "software": "VASP",
            "version": "6.3.0"
        }
    )
    
    # Add run to repository
    relax_hash = repo.add_run(relax_run)
    print(f"Relaxation run created: {relax_hash[:8]}")
    
    # Submit to SLURM
    script_path = Path("vasp_relax.slurm")
    if not script_path.exists():
        print("ERROR: vasp_relax.slurm not found!")
        return
    
    try:
        job_id = submitter.submit_run(
            run_hash=relax_hash,
            script_path=script_path,
            capture_env=True
        )
        print(f"Job submitted: {job_id}")
        print(f"Run status: {repo.get_run(relax_hash).status}")
    except Exception as e:
        print(f"Submission failed: {e}")
        return
    
    # Step 3: Monitor job status
    print("\n=== Step 3: Monitoring job status ===")
    print("(In real usage, you would poll periodically)")
    
    # Example: Check status once
    try:
        status = tracker.check_status(relax_hash)
        print(f"Current status: {status}")
        
        # Get detailed job info
        run = repo.get_run(relax_hash)
        if run.job_id:
            info = slurm.get_info(run.job_id)
            print(f"Job info: nodes={info.nodes}, state={info.state}")
    except Exception as e:
        print(f"Status check failed: {e}")
    
    # Step 4: Retrieve results (when job completes)
    print("\n=== Step 4: Retrieving results ===")
    print("(This would be done after job completion)")
    
    # Example retrieval (would fail if job not complete)
    try:
        result_files = retriever.retrieve_results(
            run_hash=relax_hash,
            output_patterns=["OUTCAR", "CONTCAR", "vasprun.xml"],
            destination=Path("./results")
        )
        print(f"Retrieved {len(result_files)} files")
    except Exception as e:
        print(f"Retrieval skipped: {e}")
    
    # Step 5: Create workflow for multi-step calculation
    print("\n=== Step 5: Creating multi-step workflow ===")
    
    # Create static calculation run
    static_run = Run(
        name="water-static",
        status=RunStatus.PLANNED,
        structure_hash=struct_hash,
        command="mpirun -np 28 vasp_std",
        working_dir="/scratch/user/water-static",
        metadata={
            "calculation_type": "single_point",
            "software": "VASP",
            "version": "6.3.0",
            "depends_on": relax_hash
        }
    )
    
    static_hash = repo.add_run(static_run)
    print(f"Static run created: {static_hash[:8]}")
    
    # Create workflow
    workflow = Workflow(
        name="water-opt-and-energy",
        run_hashes=[relax_hash, static_hash],
        metadata={
            "description": "Geometry optimization followed by single-point energy",
            "molecule": "water"
        }
    )
    
    workflow_hash = repo.add_workflow(workflow)
    print(f"Workflow created: {workflow_hash[:8]}")
    
    # Step 6: View provenance information
    print("\n=== Step 6: Provenance information ===")
    
    relax_run = repo.get_run(relax_hash)
    print(f"\nRelaxation run provenance:")
    print(f"  Job system: {relax_run.job_system}")
    print(f"  Queue: {relax_run.queue_name}")
    print(f"  Modules loaded: {relax_run.modules_loaded}")
    print(f"  Environment vars: {list(relax_run.environment_vars.keys())}")
    
    if relax_run.job_resources:
        print(f"  Resources:")
        for key, value in relax_run.job_resources.items():
            print(f"    {key}: {value}")
    
    print("\n=== Workflow complete ===")
    print(f"Repository contains:")
    print(f"  - 1 structure: {struct_hash[:8]}")
    print(f"  - 2 runs: {relax_hash[:8]}, {static_hash[:8]}")
    print(f"  - 1 workflow: {workflow_hash[:8]}")
    print(f"\nUse 'chemvcs log' to view full history")


if __name__ == "__main__":
    # Note: This example uses actual SLURM commands. To test without SLURM:
    # 1. Use MockAdapter from the user guide instead of SlurmAdapter
    # 2. Mock the job submission and status checks
    
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user")
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
