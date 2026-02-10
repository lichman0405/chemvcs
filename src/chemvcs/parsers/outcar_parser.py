"""OUTCAR file parser backed by pymatgen.

Uses pymatgen.io.vasp.outputs.Outcar for parsing VASP output files.
For incomplete / truncated OUTCAR files (common during HPC runs) a
lightweight regex fallback extracts the most essential fields so that
chemvcs can still track partial results.
"""

import os
import re
import tempfile
from typing import Any, Dict, List, Tuple

from pymatgen.io.vasp.outputs import Outcar as PmgOutcar

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class OutcarParser(BaseParser):
    """Parser for VASP OUTCAR files, powered by pymatgen.

    Extracts key computed quantities that are meaningful for tracking
    the evolution of a VASP calculation across commits:

    * ``final_energy`` – total energy (eV)
    * ``efermi`` – Fermi energy (eV)
    * ``magnetization`` – per-atom magnetic moments
    * ``total_magnetization`` – total magnetization (μB)
    * ``charge`` – per-atom Bader-like charges
    * ``run_stats`` – timing / memory statistics
    * ``is_stopped`` – whether the run completed normally

    For incomplete OUTCAR files, a regex fallback extracts whatever
    essential fields are present (energy, Fermi level, timing).
    """

    # ---------- significance classification ----------

    CRITICAL_FIELDS = frozenset({
        "final_energy", "efermi", "is_stopped",
    })

    MAJOR_FIELDS = frozenset({
        "magnetization", "total_magnetization", "charge",
    })

    # ---------- BaseParser interface ----------

    def parse(self, content: str) -> Dict[str, Any]:
        """Parse OUTCAR content.

        Tries pymatgen first (requires a complete file).  Falls back to
        regex extraction for truncated / incomplete OUTCARs.

        Args:
            content: Raw OUTCAR file content.

        Returns:
            Dictionary of extracted quantities.

        Raises:
            ParserError: If *content* cannot be parsed at all.
        """
        # Try full pymatgen parse first
        try:
            return self._parse_with_pymatgen(content)
        except Exception:
            pass

        # Fallback: regex extraction for incomplete files
        try:
            data = self._parse_essential_fields(content)
            if data:
                return data
        except Exception:
            pass

        raise ParserError(
            "Failed to parse OUTCAR: file may be empty or not a valid VASP output"
        )

    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic diff between two parsed OUTCAR results."""
        entries: List[DiffEntry] = []

        fields = [
            "final_energy",
            "efermi",
            "is_stopped",
            "total_magnetization",
            "run_stats",
        ]

        for field in fields:
            old_val = old_data.get(field)
            new_val = new_data.get(field)
            significance = self._classify(field)

            if old_val is None and new_val is not None:
                entries.append(
                    DiffEntry(field, None, new_val, "added", significance)
                )
            elif old_val is not None and new_val is None:
                entries.append(
                    DiffEntry(field, old_val, None, "deleted", significance)
                )
            elif old_val != new_val:
                entries.append(
                    DiffEntry(field, old_val, new_val, "modified", significance)
                )

        # Filter out negligible energy diffs (float noise)
        entries = [
            e
            for e in entries
            if not (
                e.path == "final_energy"
                and e.old_value is not None
                and e.new_value is not None
                and abs(e.old_value - e.new_value) < 1e-8
            )
        ]

        return entries

    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate parsed OUTCAR data."""
        errors: List[str] = []

        if data.get("is_stopped"):
            errors.append("VASP run was stopped prematurely")

        if data.get("final_energy") is None:
            errors.append("No final energy found – run may be incomplete")

        return len(errors) == 0, errors

    # ---------- pymatgen-based parsing ----------

    def _parse_with_pymatgen(self, content: str) -> Dict[str, Any]:
        """Write *content* to a temp file and parse with pymatgen."""
        fd, tmppath = tempfile.mkstemp(suffix="_OUTCAR", text=True)
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(content)

            outcar = PmgOutcar(tmppath)
            return self._extract_pymatgen(outcar)
        finally:
            try:
                os.unlink(tmppath)
            except OSError:
                pass

    @staticmethod
    def _extract_pymatgen(outcar: PmgOutcar) -> Dict[str, Any]:
        """Extract our canonical fields from a pymatgen Outcar."""
        data: Dict[str, Any] = {}

        data["final_energy"] = getattr(outcar, "final_energy", None)
        data["efermi"] = getattr(outcar, "efermi", None)

        # Per-atom magnetization
        mag = getattr(outcar, "magnetization", None)
        data["magnetization"] = (
            [dict(m) for m in mag] if mag else []
        )

        data["total_magnetization"] = getattr(outcar, "total_mag", None)

        # Per-atom charge
        chg = getattr(outcar, "charge", None)
        data["charge"] = (
            [dict(c) for c in chg] if chg else []
        )

        data["run_stats"] = getattr(outcar, "run_stats", {}) or {}
        data["is_stopped"] = getattr(outcar, "is_stopped", False)

        return data

    # ---------- regex fallback for incomplete files ----------

    @staticmethod
    def _parse_essential_fields(content: str) -> Dict[str, Any]:
        """Extract essential fields via regex when pymatgen fails."""
        data: Dict[str, Any] = {}

        # Total energy (free energy TOTEN)
        m = re.search(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)\s*eV", content)
        if m:
            data["final_energy"] = float(m.group(1))

        # Fermi energy
        m = re.search(r"E-fermi\s*:\s*([-\d.]+)", content)
        if m:
            data["efermi"] = float(m.group(1))

        # Total magnetization
        m = re.search(
            r"number of electron\s+[\d.]+\s+magnetization\s+([-\d.]+)", content
        )
        if m:
            data["total_magnetization"] = float(m.group(1))

        # Run timing
        run_stats: Dict[str, Any] = {}
        m = re.search(r"Total CPU time used \(sec\):\s*([\d.]+)", content)
        if m:
            run_stats["Total CPU time used (sec)"] = float(m.group(1))

        m = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", content)
        if m:
            run_stats["Elapsed time (sec)"] = float(m.group(1))

        m = re.search(r"Maximum memory used \(kb\):\s*(\d+)", content)
        if m:
            run_stats["Maximum memory used (kb)"] = int(m.group(1))

        data["run_stats"] = run_stats

        # Defaults for missing fields
        data.setdefault("magnetization", [])
        data.setdefault("total_magnetization", None)
        data.setdefault("charge", [])
        data.setdefault("is_stopped", False)

        # If we found nothing useful at all, return empty to signal failure
        if data.get("final_energy") is None and data.get("efermi") is None:
            return {}

        return data

    # ---------- helpers ----------

    def _classify(self, field: str) -> str:
        if field in self.CRITICAL_FIELDS:
            return "critical"
        if field in self.MAJOR_FIELDS:
            return "major"
        return "minor"
