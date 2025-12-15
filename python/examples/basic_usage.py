"""
Example: Basic usage of ChemVCS Python layer

This script demonstrates how to:
1. Create a Structure from an XYZ file
2. Convert it to a CoreObject
3. Interact with a ChemVCS repository
"""

import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemvcs_py import Repo, Structure
from chemvcs_py.io import read_xyz, write_xyz
import numpy as np


def create_water_molecule():
    """Create a simple water molecule structure."""
    structure = Structure(
        formula="H2O",
        positions=np.array([
            [0.000, 0.000, 0.000],  # O
            [0.758, 0.587, 0.000],  # H
            [-0.758, 0.587, 0.000], # H
        ]),
        species=["O", "H", "H"],
        metadata={"title": "Water molecule"}
    )
    return structure


def main():
    print("ChemVCS Python Layer - Basic Example\n")
    print("=" * 50)
    
    # Create a water molecule
    print("\n1. Creating water molecule structure...")
    water = create_water_molecule()
    print(f"   {water}")
    print(f"   Number of atoms: {water.num_atoms}")
    print(f"   Formula: {water.formula}")
    print(f"   Is periodic: {water.is_periodic}")
    
    # Write to XYZ file
    print("\n2. Writing to XYZ file...")
    xyz_path = Path("water.xyz")
    write_xyz(water, xyz_path)
    print(f"   Written to {xyz_path}")
    
    # Read back from XYZ
    print("\n3. Reading from XYZ file...")
    water_loaded = read_xyz(xyz_path)
    print(f"   {water_loaded}")
    print(f"   Formula: {water_loaded.formula}")
    print(f"   Species: {water_loaded.species}")
    
    # Convert to CoreObject
    print("\n4. Converting to CoreObject...")
    obj = water.to_core_object()
    print(f"   Type: {obj.type}")
    print(f"   Version: {obj.version}")
    print(f"   Meta keys: {list(obj.meta.keys())}")
    
    # Try to connect to repository
    print("\n5. Attempting to connect to repository...")
    try:
        repo = Repo()
        print(f"   ✓ Found repository at: {repo.root}")
        print(f"   Current branch: {repo.get_current_branch()}")
        
        # List objects
        print("\n6. Listing objects in repository...")
        objects = repo.list_objects()
        if objects:
            print(f"   Found {len(objects)} object(s):")
            for obj_info in objects[:5]:  # Show first 5
                hash_short = obj_info['Hash'][:12]
                obj_type = obj_info['Type']
                print(f"     - {hash_short}  {obj_type}")
            if len(objects) > 5:
                print(f"     ... and {len(objects) - 5} more")
        else:
            print("   No objects found")
        
        # List structures specifically
        structures = repo.list_objects(type_filter="structure")
        print(f"\n7. Found {len(structures)} structure(s) in repository")
        
    except Exception as e:
        print(f"   ✗ Could not find repository: {e}")
        print("   (This is expected if not run from a ChemVCS repository)")
    
    # Demonstrate structure operations
    print("\n8. Structure operations...")
    center = water.get_center_of_mass()
    print(f"   Center of mass: [{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]")
    
    # Translate structure
    translated = water.translate([1.0, 0.0, 0.0])
    new_center = translated.get_center_of_mass()
    print(f"   After translation: [{new_center[0]:.3f}, {new_center[1]:.3f}, {new_center[2]:.3f}]")
    
    # Round-trip through CoreObject
    print("\n9. Round-trip conversion test...")
    obj = water.to_core_object()
    water_recovered = Structure.from_core_object(obj)
    print(f"   Original: {water.formula}, {water.num_atoms} atoms")
    print(f"   Recovered: {water_recovered.formula}, {water_recovered.num_atoms} atoms")
    print(f"   ✓ Positions match: {np.allclose(water.positions, water_recovered.positions)}")
    
    print("\n" + "=" * 50)
    print("Example complete!")
    
    # Clean up
    if xyz_path.exists():
        xyz_path.unlink()
        print(f"\nCleaned up {xyz_path}")


if __name__ == "__main__":
    main()
