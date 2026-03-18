"""LAMMPS data file parser.

Parses LAMMPS data files (``data.*``, ``*.lmp``, ``*.data``) that are
typically created by tools like moltemplate, msi2lmp, or VMD and read
by the LAMMPS ``read_data`` command.

Format reference: LAMMPS Manual §4.5.4 (22 Jul 2025 edition).

A data file has:
  * A title line (ignored by LAMMPS, used by us for identification)
  * A header section with keyword lines like ``100 atoms``,
    ``2 atom types``, ``0.0 10.0 xlo xhi``, etc.
  * Named sections: ``Atoms``, ``Masses``, ``Pair Coeffs``, ``Bonds``,
    ``Angles``, ``Dihedrals``, ``Impropers``, ``Velocities``, etc.

We intentionally do NOT load full atomic coordinates (potentially millions
of atoms), only the header metadata and topology counts.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class LammpsDataParser(BaseParser):
    """Parser for LAMMPS data files.

    Extracts header metadata suitable for version-tracking:

    * ``title``            – first line of the file
    * ``num_atoms``        – total number of atoms
    * ``num_atom_types``   – number of distinct atom types
    * ``num_bonds``        – number of bonds (0 for atomic systems)
    * ``num_bond_types``   – number of bond types
    * ``num_angles``       – number of angles
    * ``num_angle_types``  – number of angle types
    * ``num_dihedrals``    – number of dihedrals
    * ``num_dihedral_types``
    * ``num_impropers``    – number of impropers
    * ``num_improper_types``
    * ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, ``zhi``
    * ``xy``, ``xz``, ``yz``  – triclinic tilt factors (if present)
    * ``lx``, ``ly``, ``lz``  – box side lengths (derived)
    * ``volume``           – box volume (Å³, derived)
    * ``has_bonds``        – True if molecular topology is present
    * ``sections``         – list of section names present in the file
    * ``masses``           – dict of atom_type → mass (if Masses section present)
    """

    # ------------------------------------------------------------------ #
    # Significance classification                                          #
    # ------------------------------------------------------------------ #

    CRITICAL_FIELDS = frozenset({
        "num_atoms", "num_atom_types", "lx", "ly", "lz", "volume",
    })

    MAJOR_FIELDS = frozenset({
        "has_bonds", "num_bonds", "num_bond_types",
        "num_angles", "num_angle_types",
        "num_dihedrals", "num_dihedral_types",
        "num_impropers", "num_improper_types",
        "masses",
    })

    # ------------------------------------------------------------------ #
    # Header keyword patterns                                              #
    # ------------------------------------------------------------------ #

    # Integer count lines, e.g. "1000 atoms", "2 atom types", "500 bonds"
    _COUNT_RE = re.compile(
        r"^\s*(\d+)\s+"
        r"(atoms|atom\s+types|bonds|bond\s+types|angles|angle\s+types|"
        r"dihedrals|dihedral\s+types|impropers|improper\s+types"
        r")\s*(?:#.*)?$",
        re.IGNORECASE,
    )

    # Box boundary lines, e.g. "0.0 50.0 xlo xhi"
    _BOX_RE = re.compile(
        r"^\s*([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+"
        r"(xlo|ylo|zlo)\s+(xhi|yhi|zhi)\s*(?:#.*)?$",
        re.IGNORECASE,
    )

    # Triclinic tilt: "xy xz yz"
    _TILT_RE = re.compile(
        r"^\s*([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+"
        r"xy\s+xz\s+yz\s*(?:#.*)?$",
        re.IGNORECASE,
    )

    # Named section headers
    _SECTION_RE = re.compile(
        r"^(Atoms|Velocities|Masses|Pair\s+Coeffs|PairIJ\s+Coeffs|"
        r"Bond\s+Coeffs|Angle\s+Coeffs|Dihedral\s+Coeffs|Improper\s+Coeffs|"
        r"Bonds|Angles|Dihedrals|Impropers|"
        r"BondBond\s+Coeffs|BondAngle\s+Coeffs|"
        r"AngleAngle\s+Coeffs|AngleTorsion\s+Coeffs|"
        r"EndBondTorsion\s+Coeffs|MiddleBondTorsion\s+Coeffs|"
        r"BondBond13\s+Coeffs|AngleAngleTorsion\s+Coeffs"
        r")\s*(?:#.*)?$",
        re.IGNORECASE,
    )

    # ------------------------------------------------------------------ #
    # BaseParser interface                                                  #
    # ------------------------------------------------------------------ #

    def parse(self, content: str) -> Dict[str, Any]:
        """Parse LAMMPS data file content.

        Args:
            content: Full text of the data file.

        Returns:
            Dictionary of header metadata.

        Raises:
            ParserError: If the file cannot be parsed as a LAMMPS data file.
        """
        if not content or not content.strip():
            raise ParserError("LAMMPS data file is empty")

        lines = content.splitlines()

        # First non-empty line is the title
        title = ""
        for line in lines:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                title = stripped
                break

        data: Dict[str, Any] = {"title": title}

        # Parse header (lines before first named section)
        sections: List[str] = []
        in_header = True
        masses: Dict[int, float] = {}
        in_masses_section = False

        for line in lines[1:]:  # skip title
            stripped = line.strip()

            # Skip blank lines and comments
            if not stripped or stripped.startswith("#"):
                # Blank lines do NOT end a section — only the next section
                # header does. We simply skip blank/comment lines everywhere.
                continue

            # Detect section headers
            sec_match = self._SECTION_RE.match(stripped)
            if sec_match:
                in_header = False
                sec_name = " ".join(stripped.split()[:2])  # normalise
                sections.append(stripped.split("#")[0].strip())
                in_masses_section = sec_name.lower() == "masses"
                continue

            if in_header:
                # Try count line
                m = self._COUNT_RE.match(stripped)
                if m:
                    count = int(m.group(1))
                    key = m.group(2).lower().replace(" ", "_")
                    data[f"num_{key}"] = count
                    continue

                # Try box line
                m = self._BOX_RE.match(stripped)
                if m:
                    lo_val = float(m.group(1))
                    hi_val = float(m.group(2))
                    axis = m.group(3).lower()[:1]  # 'x', 'y', or 'z'
                    data[f"{axis}lo"] = lo_val
                    data[f"{axis}hi"] = hi_val
                    continue

                # Try tilt line
                m = self._TILT_RE.match(stripped)
                if m:
                    data["xy"] = float(m.group(1))
                    data["xz"] = float(m.group(2))
                    data["yz"] = float(m.group(3))
                    continue

            elif in_masses_section:
                # Masses section: "type_id mass  # optional comment"
                parts = stripped.split("#")[0].split()
                if len(parts) == 2:
                    try:
                        atom_type = int(parts[0])
                        mass = float(parts[1])
                        masses[atom_type] = mass
                    except ValueError:
                        pass

        # Derived quantities
        for axis in ("x", "y", "z"):
            lo = data.get(f"{axis}lo")
            hi = data.get(f"{axis}hi")
            if lo is not None and hi is not None:
                data[f"l{axis}"] = round(hi - lo, 10)

        lx = data.get("lx")
        ly = data.get("ly")
        lz = data.get("lz")
        if lx is not None and ly is not None and lz is not None:
            data["volume"] = round(lx * ly * lz, 6)

        data["has_bonds"] = data.get("num_bonds", 0) > 0
        data["sections"] = sections

        if masses:
            data["masses"] = {str(k): v for k, v in sorted(masses.items())}

        # Sanity check
        if "num_atoms" not in data:
            raise ParserError(
                "Could not find atom count in LAMMPS data file. "
                "Ensure the file starts with a title line and has a valid header."
            )

        return data

    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic diff between two parsed LAMMPS data files."""
        entries: List[DiffEntry] = []

        scalar_fields = [
            "num_atoms",
            "num_atom_types",
            "num_bonds",
            "num_bond_types",
            "num_angles",
            "num_angle_types",
            "num_dihedrals",
            "num_dihedral_types",
            "num_impropers",
            "num_improper_types",
            "lx", "ly", "lz",
            "volume",
            "has_bonds",
            "title",
        ]

        for field in scalar_fields:
            old_val = old_data.get(field)
            new_val = new_data.get(field)
            significance = self._classify(field)

            if old_val is None and new_val is not None:
                entries.append(DiffEntry(field, None, new_val, "added", significance))
            elif old_val is not None and new_val is None:
                entries.append(DiffEntry(field, old_val, None, "deleted", significance))
            elif old_val != new_val:
                # Suppress float noise for box dimensions
                if (
                    isinstance(old_val, float)
                    and isinstance(new_val, float)
                    and abs(old_val - new_val) < 1e-8
                ):
                    continue
                entries.append(
                    DiffEntry(field, old_val, new_val, "modified", significance)
                )

        # Diff masses dict
        old_masses: Dict[str, float] = old_data.get("masses", {})
        new_masses: Dict[str, float] = new_data.get("masses", {})
        if old_masses != new_masses:
            entries.append(
                DiffEntry("masses", old_masses, new_masses, "modified", "major")
            )

        return entries

    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate parsed LAMMPS data file."""
        errors: List[str] = []

        if data.get("num_atoms", 0) <= 0:
            errors.append("num_atoms must be positive")

        if data.get("num_atom_types", 0) <= 0:
            errors.append("num_atom_types must be positive")

        for axis in ("lx", "ly", "lz"):
            val = data.get(axis)
            if val is not None and val <= 0:
                errors.append(f"Box dimension {axis} must be positive (got {val})")

        return len(errors) == 0, errors

    # ------------------------------------------------------------------ #
    # Internal helpers                                                     #
    # ------------------------------------------------------------------ #

    def _classify(self, field: str) -> str:
        if field in self.CRITICAL_FIELDS:
            return "critical"
        if field in self.MAJOR_FIELDS:
            return "major"
        return "minor"
