"""INCAR file parser backed by pymatgen.

Uses pymatgen.io.vasp.inputs.Incar for robust, production-grade parsing
of VASP INCAR files.  All type inference, boolean handling, list values,
semicolon-separated entries, and comment stripping are delegated to pymatgen.
"""

from typing import Any, Dict, List, Tuple

from pymatgen.io.vasp.inputs import Incar as PmgIncar

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class IncarParser(BaseParser):
    """Parser for VASP INCAR files, powered by pymatgen.

    Pymatgen's ``Incar`` class provides:

    * Automatic type inference for every known VASP tag
    * Proper handling of booleans (``.TRUE.`` / ``.FALSE.`` / ``T`` / ``F``)
    * Multi-value parameters (``MAGMOM``, ``LDAUL``, etc.)
    * Semicolon-separated entries on a single line
    * Comment stripping (``#`` and ``!``)

    On top of pymatgen's parsed data we add **significance classification**
    (critical / major / minor) used by ``chemvcs diff``.
    """

    # ---------- significance classification ----------

    CRITICAL_TAGS = frozenset({
        "ENCUT", "PREC", "ALGO", "LREAL",
        "ISMEAR", "SIGMA", "EDIFF", "EDIFFG",
        "NSW", "IBRION", "POTIM", "ISIF",
        "NELM", "NELMIN",
    })

    MAJOR_TAGS = frozenset({
        "LWAVE", "LCHARG", "LORBIT", "LVTOT", "LVHAR",
        "ISPIN", "MAGMOM", "NCORE", "KPAR",
        "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ",
        "IVDW", "LUSE_VDW", "AGGAC",
    })

    # ---------- BaseParser interface ----------

    def parse(self, content: str) -> Dict[str, Any]:
        """Parse INCAR content into a typed parameter dictionary.

        Args:
            content: Raw INCAR file content.

        Returns:
            ``dict`` mapping uppercase tag names to their Python-typed values.

        Raises:
            ParserError: If pymatgen cannot parse *content*.
        """
        try:
            incar = PmgIncar.from_str(content)
            return dict(incar)
        except Exception as e:
            raise ParserError(f"Failed to parse INCAR: {e}") from e

    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic diff between two parsed INCAR dictionaries.

        Each changed tag is annotated with a *significance* level based on
        its impact on computational results.
        """
        diff_entries: List[DiffEntry] = []
        all_keys = set(old_data.keys()) | set(new_data.keys())

        for key in sorted(all_keys):
            old_value = old_data.get(key)
            new_value = new_data.get(key)

            significance = self._classify(key)

            if old_value is None:
                diff_entries.append(
                    DiffEntry(key, None, new_value, "added", significance)
                )
            elif new_value is None:
                diff_entries.append(
                    DiffEntry(key, old_value, None, "deleted", significance)
                )
            elif old_value != new_value:
                diff_entries.append(
                    DiffEntry(key, old_value, new_value, "modified", significance)
                )

        return diff_entries

    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate commonly-checked INCAR constraints."""
        errors: List[str] = []

        if "ENCUT" in data:
            if not isinstance(data["ENCUT"], (int, float)) or data["ENCUT"] <= 0:
                errors.append("ENCUT must be a positive number")

        if "ISMEAR" in data and isinstance(data["ISMEAR"], int):
            if data["ISMEAR"] < -5 or data["ISMEAR"] > 5:
                errors.append("ISMEAR should be in range -5 to 5")

        if "SIGMA" in data:
            if not isinstance(data["SIGMA"], (int, float)) or data["SIGMA"] <= 0:
                errors.append("SIGMA must be a positive number")

        if "NSW" in data:
            if not isinstance(data["NSW"], int) or data["NSW"] < 0:
                errors.append("NSW must be a non-negative integer")

        if "ISMEAR" in data and "SIGMA" in data:
            if (
                data["ISMEAR"] == 0
                and isinstance(data["SIGMA"], (int, float))
                and data["SIGMA"] > 0.2
            ):
                errors.append(
                    "Warning: ISMEAR=0 (Gaussian) with large SIGMA (>0.2) "
                    "may be inappropriate"
                )

        return len(errors) == 0, errors

    # ---------- helpers ----------

    def _classify(self, tag: str) -> str:
        """Return significance level for *tag*."""
        if tag in self.CRITICAL_TAGS:
            return "critical"
        if tag in self.MAJOR_TAGS:
            return "major"
        return "minor"
