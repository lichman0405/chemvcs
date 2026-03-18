"""File format validator.

Validates basic VASP file format correctness by attempting to parse each
recognised file with the corresponding ChemVCS parser.  Parsing errors
are surfaced as validation errors; the raw file is never silently accepted.
"""

from pathlib import Path
from typing import Any, List

from chemvcs.parsers.base_parser import ParserError
from chemvcs.parsers.incar_parser import IncarParser
from chemvcs.parsers.kpoints_parser import KpointsParser
from chemvcs.plugins.base import ValidationResult, ValidatorPlugin

_incar_parser = IncarParser()
_kpoints_parser = KpointsParser()

# VASP files supported by this validator
_SUPPORTED = {"INCAR", "KPOINTS", "POSCAR"}


def _validate_poscar_format(path: Path) -> List[str]:
    """Return a list of format errors found in a POSCAR file.

    Checks (VASP 5+ format):
    - At least 8 non-empty lines
    - Line 2: scale factor is a non-zero float
    - Lines 3-5: 3 lattice-vector lines each with exactly 3 floats
    - Line 6: at least one element symbol (alpha text)
    - Line 7: non-negative integers, count matches line 6
    - Line 8: starts with S/s (selective dynamics) or C/c/D/d/K/k (coord mode)
    - Coordinate block: correct number of rows (sum of atom counts)
    """
    errors: List[str] = []
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return [f"Cannot read POSCAR: {exc}"]

    lines = [l for l in text.splitlines()]  # keep blank lines for line counting
    # Strip trailing blank lines but keep content lines
    content = [l for l in lines if l.strip()]
    if len(content) < 8:
        return [f"POSCAR has only {len(content)} non-empty lines; at least 8 required"]

    # Line 2: scale factor
    try:
        scale = float(content[1].split()[0])
        if scale == 0:
            errors.append("POSCAR scale factor (line 2) must not be zero")
    except (ValueError, IndexError):
        errors.append("POSCAR line 2 must be a numeric scale factor")

    # Lines 3-5: lattice vectors
    for i in range(2, 5):
        tokens = content[i].split()
        if len(tokens) < 3:
            errors.append(f"POSCAR lattice vector line {i+1} must have 3 components")
            continue
        try:
            [float(t) for t in tokens[:3]]
        except ValueError:
            errors.append(f"POSCAR lattice vector line {i+1} contains non-numeric values")

    # Line 6: element symbols
    elem_tokens = content[5].split()
    if not elem_tokens:
        errors.append("POSCAR line 6 (element symbols) is empty")
    elif all(t.lstrip("-").isdigit() for t in elem_tokens):
        errors.append(
            "POSCAR line 6 looks like atom counts (VASP 4 format). "
            "Use VASP 5+ format: add element symbols before atom counts."
        )

    # Line 7: atom counts
    count_tokens = content[6].split()
    atom_counts: List[int] = []
    for t in count_tokens:
        try:
            n = int(t)
            if n < 0:
                errors.append(f"POSCAR atom count '{t}' must be non-negative")
            atom_counts.append(n)
        except ValueError:
            errors.append(f"POSCAR atom count token '{t}' is not an integer")

    if elem_tokens and atom_counts and len(count_tokens) != len(elem_tokens):
        errors.append(
            f"POSCAR element count ({len(elem_tokens)}) does not match "
            f"atom-count entries ({len(count_tokens)})"
        )

    # Line 8: selective dynamics or coordinate mode
    coord_line_idx = 7
    if content[7][0].upper() == "S":
        coord_line_idx = 8  # selective dynamics present

    if coord_line_idx >= len(content):
        errors.append("POSCAR is missing coordinate mode line")
        return errors

    coord_mode = content[coord_line_idx].strip()[0].upper() if content[coord_line_idx].strip() else ""
    if coord_mode not in ("C", "D", "K"):
        errors.append(
            f"POSCAR coordinate mode line must start with C/D/K, got: '{content[coord_line_idx][:20]}'"
        )

    # Coordinate block: must have exactly sum(atom_counts) rows
    if atom_counts:
        expected = sum(atom_counts)
        coord_rows = content[coord_line_idx + 1 :]
        if len(coord_rows) < expected:
            errors.append(
                f"POSCAR has {len(coord_rows)} coordinate row(s) but "
                f"{expected} atom(s) declared"
            )

    return errors


class FileFormatValidator(ValidatorPlugin):
    """Validates basic VASP file format correctness.

    Validation strategy
    -------------------
    - **INCAR**: parsed by ``IncarParser`` (backed by pymatgen); any
      ``ParserError`` is reported as an error.
    - **KPOINTS**: parsed by ``KpointsParser``; ``ParserError`` is reported.
    - **POSCAR**: structural checks (scale, lattice vectors, element line,
      atom counts, coordinate block size) without an external library.

    This validator is **disabled by default** because it runs on every
    ``chemvcs add``.  Enable it via the plugin config if desired.
    """

    @property
    def name(self) -> str:
        return "file-format"

    @property
    def version(self) -> str:
        return "1.0.0"

    @property
    def description(self) -> str:
        return "Validates basic VASP file format (INCAR, KPOINTS, POSCAR)"

    @property
    def priority(self) -> int:
        return 5  # Run first

    @property
    def enabled_by_default(self) -> bool:
        return False  # Opt-in

    def can_validate(self, files: list[str]) -> bool:
        file_names = {Path(f).name for f in files}
        return bool(file_names & _SUPPORTED)

    def validate(
        self,
        workspace_root: Path,
        files: list[str],
        **kwargs: Any,
    ) -> ValidationResult:
        errors: List[str] = []
        checked: List[str] = []

        for rel in files:
            name = Path(rel).name
            if name not in _SUPPORTED:
                continue

            full = workspace_root / rel
            if not full.exists():
                continue

            checked.append(name)
            text = full.read_text(encoding="utf-8", errors="replace")

            if name == "INCAR":
                try:
                    _incar_parser.parse(text)
                except ParserError as exc:
                    errors.append(f"INCAR parse error: {exc}")

            elif name == "KPOINTS":
                try:
                    _kpoints_parser.parse(text)
                except ParserError as exc:
                    errors.append(f"KPOINTS parse error: {exc}")

            elif name == "POSCAR":
                poscar_errors = _validate_poscar_format(full)
                errors.extend(poscar_errors)

        if not checked:
            return ValidationResult(passed=True, message="No supported files to validate")

        passed = len(errors) == 0
        if passed:
            message = f"File format OK: {', '.join(sorted(checked))}"
        else:
            message = f"File format errors in: {', '.join(sorted(checked))}"

        return ValidationResult(passed=passed, message=message, errors=errors)
