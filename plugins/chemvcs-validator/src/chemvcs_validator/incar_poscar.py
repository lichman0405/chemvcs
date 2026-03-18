"""INCAR-POSCAR consistency validator.

Validates that INCAR parameters are consistent with the POSCAR structure.
Checks MAGMOM length, LDAU list lengths, ISPIN vs magnetism, NSW vs ISIF.
"""

import re
from pathlib import Path
from typing import Any, List, Optional, Tuple

from chemvcs.plugins.base import ValidationResult, ValidatorPlugin


def _parse_poscar_atoms(poscar_path: Path) -> Tuple[List[str], List[int]]:
    """Return (element_symbols, atom_counts) from a VASP 5+ POSCAR.

    Raises:
        ValueError: If the file cannot be parsed or is in VASP 4 format.
    """
    lines = poscar_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 7:
        raise ValueError("POSCAR has fewer than 7 lines")

    # Line 6 (index 5): element symbols (VASP 5+) or atom counts (VASP 4)
    # Line 7 (index 6): atom counts (VASP 5+) or S/D/C/K (VASP 4)
    tokens_5 = lines[5].split()
    tokens_6 = lines[6].split()

    # VASP 5+: line 5 has text symbols, line 6 has integers
    if all(t.lstrip("-").isdigit() for t in tokens_5):
        raise ValueError(
            "VASP 4.x POSCAR detected (no element symbols). "
            "Please use VASP 5+ format with element symbols on line 6."
        )

    try:
        counts = [int(t) for t in tokens_6]
    except ValueError:
        raise ValueError("Could not parse atom counts from POSCAR line 7")

    if len(tokens_5) != len(counts):
        raise ValueError(
            f"POSCAR element count ({len(tokens_5)}) does not match "
            f"atom-count entries ({len(counts)})"
        )

    return tokens_5, counts


def _expand_magmom(magmom_value: Any) -> Optional[List[float]]:
    """Expand a MAGMOM value (as parsed by pymatgen) into a flat list.

    pymatgen returns MAGMOM as a list that may contain floats or
    compressed entries like ``3*0.6``.  This function returns the fully
    expanded flat list, or None if expansion fails.
    """
    if magmom_value is None:
        return None

    # pymatgen already expands MAGMOM into a plain list of floats/ints
    # when it can.  Handle both string (raw) and list forms.
    if isinstance(magmom_value, (int, float)):
        return [float(magmom_value)]

    if isinstance(magmom_value, list):
        result: List[float] = []
        for item in magmom_value:
            if isinstance(item, (int, float)):
                result.append(float(item))
            else:
                # Try "n*val" compressed string form
                m = re.match(r"^(\d+)\*(.+)$", str(item).strip())
                if m:
                    n, val = int(m.group(1)), float(m.group(2))
                    result.extend([val] * n)
                else:
                    try:
                        result.append(float(item))
                    except (TypeError, ValueError):
                        return None  # Can't parse — skip check
        return result

    # Fallback: try parsing a raw string like "2*5 3*0"
    if isinstance(magmom_value, str):
        result = []
        for token in magmom_value.split():
            m = re.match(r"^(\d+)\*(.+)$", token)
            if m:
                n, val = int(m.group(1)), float(m.group(2))
                result.extend([val] * n)
            else:
                try:
                    result.append(float(token))
                except ValueError:
                    return None
        return result

    return None


class INCARPOSCARValidator(ValidatorPlugin):
    """Validates INCAR-POSCAR consistency.

    Checks performed
    ----------------
    1. **MAGMOM length** — must equal the total number of atoms in POSCAR.
    2. **LDAUU / LDAUJ / LDAUL length** — each list must equal the number of
       element species in POSCAR.
    3. **ISPIN = 2 without MAGMOM** — warns that initial moments are unset.
    4. **NSW = 0 with ISIF ≥ 3** — warns that cell/ionic relaxation is
       requested but no ionic steps are allowed.
    """

    @property
    def name(self) -> str:
        return "incar-poscar"

    @property
    def version(self) -> str:
        return "1.0.0"

    @property
    def description(self) -> str:
        return "Validates INCAR parameter consistency with POSCAR structure"

    @property
    def priority(self) -> int:
        return 20

    @property
    def enabled_by_default(self) -> bool:
        return True

    def can_validate(self, files: list[str]) -> bool:
        file_names = {Path(f).name for f in files}
        return "INCAR" in file_names and "POSCAR" in file_names

    def validate(
        self,
        workspace_root: Path,
        files: list[str],
        **kwargs: Any,
    ) -> ValidationResult:
        incar_path = workspace_root / "INCAR"
        poscar_path = workspace_root / "POSCAR"

        if not incar_path.exists() or not poscar_path.exists():
            return ValidationResult(passed=True, message="INCAR or POSCAR not found, skipping")

        errors: List[str] = []
        warnings: List[str] = []

        # ── Parse POSCAR ────────────────────────────────────────────────
        try:
            elements, counts = _parse_poscar_atoms(poscar_path)
        except ValueError as exc:
            return ValidationResult(
                passed=True,
                message="Could not parse POSCAR — skipping INCAR-POSCAR checks",
                warnings=[str(exc)],
            )

        n_atoms = sum(counts)
        n_species = len(elements)

        # ── Parse INCAR (via pymatgen) ───────────────────────────────────
        try:
            from pymatgen.io.vasp.inputs import Incar as PmgIncar
            incar = dict(PmgIncar.from_file(str(incar_path)))
        except Exception as exc:
            return ValidationResult(
                passed=True,
                message="Could not parse INCAR — skipping INCAR-POSCAR checks",
                warnings=[str(exc)],
            )

        # ── Check 1: MAGMOM length ───────────────────────────────────────
        if "MAGMOM" in incar:
            expanded = _expand_magmom(incar["MAGMOM"])
            if expanded is not None and len(expanded) != n_atoms:
                errors.append(
                    f"MAGMOM has {len(expanded)} value(s) but POSCAR has "
                    f"{n_atoms} atom(s) — they must match."
                )

        # ── Check 2: LDAUU / LDAUJ / LDAUL length ───────────────────────
        for tag in ("LDAUU", "LDAUJ", "LDAUL"):
            if tag in incar:
                val = incar[tag]
                if isinstance(val, list) and len(val) != n_species:
                    errors.append(
                        f"{tag} has {len(val)} value(s) but POSCAR has "
                        f"{n_species} element species ({', '.join(elements)}) — they must match."
                    )

        # ── Check 3: ISPIN=2 without MAGMOM ─────────────────────────────
        ispin = incar.get("ISPIN", 1)
        if ispin == 2 and "MAGMOM" not in incar:
            warnings.append(
                "ISPIN = 2 is set but MAGMOM is absent; "
                "VASP will use default initial moments (5 μB per atom for d/f elements)."
            )

        # ── Check 4: NSW=0 with ISIF≥3 ──────────────────────────────────
        nsw = incar.get("NSW", 0)
        isif = incar.get("ISIF", 2)
        if nsw == 0 and isinstance(isif, int) and isif >= 3:
            warnings.append(
                f"NSW = 0 but ISIF = {isif} (≥ 3) requests cell/volume relaxation; "
                "no ionic steps will run — check if this is intentional."
            )

        passed = len(errors) == 0
        if errors:
            message = f"INCAR-POSCAR consistency check failed ({len(errors)} error(s))"
        elif warnings:
            message = f"INCAR-POSCAR consistency check passed with {len(warnings)} warning(s)"
        else:
            message = (
                f"INCAR-POSCAR consistent: {n_atoms} atom(s), "
                f"{n_species} species ({', '.join(elements)})"
            )

        return ValidationResult(
            passed=passed,
            message=message,
            warnings=warnings,
            errors=errors,
            details={"n_atoms": n_atoms, "n_species": n_species, "elements": elements},
        )
