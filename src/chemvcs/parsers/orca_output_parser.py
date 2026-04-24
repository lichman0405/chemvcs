"""ORCA output file parser backed by cclib.

Uses ``cclib`` for primary parsing of ORCA ``.out`` files.  For
incomplete / truncated outputs (common on HPC clusters) a lightweight
regex fallback extracts the most essential fields so that ``chemvcs``
can still track partial results.

Architecture mirrors :class:`chemvcs.parsers.outcar_parser.OutcarParser`
which uses pymatgen as primary and regex as fallback.

Fields extracted
----------------
``final_energy``            – last total energy in Hartrees (float)
``termination_status``      – ``"normal"`` or ``"error"`` (str)
``scf_converged``           – whether the last SCF cycle converged (bool)
``scf_cycles``              – number of SCF iterations in last cycle (int)
``optimization_converged``  – whether geometry optimisation converged (bool)
``orca_version``            – version string from the ORCA header (str)
``charge``                  – molecular charge (int)
``multiplicity``            – spin multiplicity (int)
``dipole_moment``           – total dipole magnitude in Debye (float)
``wall_time``               – total wall-time string from footer (str)

Energy convention: Hartrees (ORCA native unit).
``cclib`` internally works in eV; conversion factor applied.
"""

from __future__ import annotations

import contextlib
import os
import re
import tempfile
from typing import Any

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError

# Hartree ↔ eV conversion factor (NIST 2018 CODATA)
_EV_TO_HARTREE = 1.0 / 27.211386245988


class OrcaOutputParser(BaseParser):
    """Parser for ORCA output (``.out``) files.

    Uses cclib as the primary parser with a regex fallback for robustness.
    """

    # ------------------------------------------------------------------ #
    # Significance classification                                         #
    # ------------------------------------------------------------------ #

    CRITICAL_FIELDS: frozenset[str] = frozenset(
        {
            "final_energy",
            "termination_status",
            "optimization_converged",
            "charge",
            "multiplicity",
        }
    )

    MAJOR_FIELDS: frozenset[str] = frozenset(
        {
            "dipole_moment",
            "scf_cycles",
            "scf_converged",
        }
    )

    # ------------------------------------------------------------------ #
    # BaseParser interface                                                  #
    # ------------------------------------------------------------------ #

    def parse(self, content: str) -> dict[str, Any]:
        """Parse ORCA output content.

        Tries cclib first (requires a reasonably complete file).
        Falls back to regex extraction for truncated / incomplete outputs.

        Args:
            content: Raw ORCA ``.out`` file content.

        Returns:
            Structured data dictionary.

        Raises:
            ParserError: If the content cannot be parsed at all.
        """
        if not content or not content.strip():
            raise ParserError("ORCA output file is empty")

        # Primary: cclib
        try:
            data = self._parse_with_cclib(content)
            if data:
                return data
        except Exception:
            pass

        # Fallback: regex
        try:
            data = self._parse_with_regex(content)
            # Only accept if at least one ORCA-specific field was found
            if data and (
                data.get("final_energy") is not None
                or data.get("termination_status") is not None
                or data.get("orca_version") is not None
            ):
                return data
        except Exception:
            pass

        raise ParserError(
            "Failed to parse ORCA output: file may be empty or not a valid ORCA output"
        )

    def diff(
        self,
        old_data: dict[str, Any],
        new_data: dict[str, Any],
    ) -> list[DiffEntry]:
        """Compute semantic diff between two parsed ORCA outputs."""
        entries: list[DiffEntry] = []

        fields = [
            "final_energy",
            "termination_status",
            "optimization_converged",
            "charge",
            "multiplicity",
            "scf_converged",
            "scf_cycles",
            "dipole_moment",
            "orca_version",
            "wall_time",
        ]

        for field in fields:
            old_val = old_data.get(field)
            new_val = new_data.get(field)
            significance = self._classify(field)

            if old_val is None and new_val is not None:
                entries.append(DiffEntry(field, None, new_val, "added", significance))
            elif old_val is not None and new_val is None:
                entries.append(DiffEntry(field, old_val, None, "deleted", significance))
            elif old_val != new_val:
                entries.append(DiffEntry(field, old_val, new_val, "modified", significance))

        # Suppress float noise for energy comparison (< 1e-10 Hartree)
        entries = [
            e
            for e in entries
            if not (
                e.path == "final_energy"
                and e.old_value is not None
                and e.new_value is not None
                and isinstance(e.old_value, float)
                and isinstance(e.new_value, float)
                and abs(e.old_value - e.new_value) < 1e-10
            )
        ]

        return entries

    def validate(self, data: dict[str, Any]) -> tuple[bool, list[str]]:
        """Validate parsed ORCA output data.

        Args:
            data: Parsed data dictionary from :meth:`parse`.

        Returns:
            ``(is_valid, list_of_error_messages)``
        """
        errors: list[str] = []

        if data.get("termination_status") != "normal":
            errors.append("ORCA run did not terminate normally")

        if data.get("final_energy") is None:
            errors.append("No final energy found — run may be incomplete")

        if data.get("optimization_converged") is False:
            errors.append("Geometry optimisation did not converge")

        return len(errors) == 0, errors

    # ------------------------------------------------------------------ #
    # cclib-based parsing                                                  #
    # ------------------------------------------------------------------ #

    def _parse_with_cclib(self, content: str) -> dict[str, Any] | None:
        """Write *content* to a temp file and parse with cclib.

        Returns ``None`` if cclib cannot parse the file.
        """
        try:
            import cclib  # noqa: PLC0415  (local import to keep optional)
        except ImportError as exc:
            raise ImportError(
                "cclib is required for ORCA output parsing: pip install cclib"
            ) from exc

        fd, tmppath = tempfile.mkstemp(suffix=".out", text=True)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(content)

            # Suppress cclib's own logging BEFORE calling ccopen/parse
            import logging

            logging.getLogger("cclib").setLevel(logging.ERROR)

            parser = cclib.io.ccopen(tmppath)
            if parser is None:
                return None

            ccdata = parser.parse()
        except Exception:
            return None
        finally:
            with contextlib.suppress(OSError):
                os.unlink(tmppath)

        return self._extract_from_cclib(ccdata, content)

    @staticmethod
    def _extract_from_cclib(ccdata: Any, raw_content: str) -> dict[str, Any]:
        """Extract canonical fields from a parsed cclib data object."""
        data: dict[str, Any] = {}

        # ---- Energy ------------------------------------------------- #
        # scfenergies is in eV; convert to Hartrees for ORCA convention
        scfenergies = getattr(ccdata, "scfenergies", None)
        mpenergies = getattr(ccdata, "mpenergies", None)
        ccenergies = getattr(ccdata, "ccenergies", None)

        final_ev: float | None = None
        if ccenergies is not None and len(ccenergies) > 0:
            final_ev = float(
                ccenergies[-1][-1] if hasattr(ccenergies[-1], "__len__") else ccenergies[-1]
            )
        elif mpenergies is not None and len(mpenergies) > 0:
            final_ev = float(
                mpenergies[-1][-1] if hasattr(mpenergies[-1], "__len__") else mpenergies[-1]
            )
        elif scfenergies is not None and len(scfenergies) > 0:
            final_ev = float(scfenergies[-1])

        data["final_energy"] = final_ev * _EV_TO_HARTREE if final_ev is not None else None

        # ---- Termination -------------------------------------------- #
        # cclib sets metadata["success"] = True when ORCA TERMINATED NORMALLY
        metadata = getattr(ccdata, "metadata", {}) or {}
        success = metadata.get("success", None)
        if success is True:
            data["termination_status"] = "normal"
        elif success is False:
            data["termination_status"] = "error"
        else:
            # Fallback: check raw content for the footer string
            if (
                "ORCA TERMINATED NORMALLY" in raw_content
                or "****ORCA TERMINATED NORMALLY****" in raw_content
            ):
                data["termination_status"] = "normal"
            else:
                data["termination_status"] = None

        # ---- SCF convergence ---------------------------------------- #
        scf_targets = getattr(ccdata, "scftargets", None)
        scf_values = getattr(ccdata, "scfvalues", None)
        if scf_targets is not None and scf_values is not None:
            data["scf_converged"] = True  # cclib only reports values when converged
        else:
            data["scf_converged"] = None

        # Number of SCF iterations in the last macro-cycle
        if scf_values is not None and len(scf_values) > 0:
            data["scf_cycles"] = len(scf_values[-1])
        else:
            data["scf_cycles"] = None

        # ---- Optimisation ------------------------------------------- #
        optdone = getattr(ccdata, "optdone", None)
        if optdone is not None:
            data["optimization_converged"] = (
                bool(optdone) if not hasattr(optdone, "__len__") else bool(len(optdone) > 0)
            )
        else:
            data["optimization_converged"] = None

        # ---- Molecular properties ----------------------------------- #
        data["charge"] = getattr(ccdata, "charge", None)
        data["multiplicity"] = getattr(ccdata, "mult", None)

        # Dipole moment magnitude (Debye) — cclib stores as 3-vector in Debye
        dipole = getattr(ccdata, "moments", None)
        if dipole is not None and len(dipole) > 1:
            import math

            dvec = dipole[1]
            data["dipole_moment"] = round(math.sqrt(sum(x**2 for x in dvec)), 6)
        else:
            data["dipole_moment"] = None

        # ---- Version / timing --------------------------------------- #
        data["orca_version"] = metadata.get("package_version", None)
        data["wall_time"] = _extract_wall_time(raw_content)

        return data

    # ------------------------------------------------------------------ #
    # Regex fallback for incomplete files                                  #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _parse_with_regex(content: str) -> dict[str, Any]:
        """Extract essential fields via regex when cclib fails or the
        output is truncated."""
        data: dict[str, Any] = {}

        # Final energy (Hartrees) — most recent occurrence wins
        for m in re.finditer(r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)", content):
            data["final_energy"] = float(m.group(1))

        # Termination status
        if re.search(r"ORCA TERMINATED NORMALLY", content):
            data["termination_status"] = "normal"
        elif re.search(r"ERROR\s*TERMINATION", content, re.IGNORECASE):
            data["termination_status"] = "error"
        else:
            data["termination_status"] = None

        # SCF convergence
        scf_m = re.search(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", content)
        if scf_m:
            data["scf_converged"] = True
            data["scf_cycles"] = int(scf_m.group(1))
        elif re.search(r"SCF NOT CONVERGED", content, re.IGNORECASE):
            data["scf_converged"] = False
            data["scf_cycles"] = None
        else:
            data["scf_converged"] = None
            data["scf_cycles"] = None

        # Geometry optimisation
        if re.search(r"THE OPTIMIZATION HAS CONVERGED", content):
            data["optimization_converged"] = True
        elif re.search(r"THE OPTIMIZATION HAS NOT CONVERGED", content):
            data["optimization_converged"] = False
        else:
            data["optimization_converged"] = None

        # ORCA version
        ver_m = re.search(
            r"Program Version\s+([\d.]+(?:-[\w]+)?)",
            content,
            re.IGNORECASE,
        )
        if ver_m:
            data["orca_version"] = ver_m.group(1)
        else:
            # Older format: "This ORCA version 4.x.y"
            ver_m2 = re.search(r"[Tt]his ORCA version\s+([\d.]+)", content)
            if ver_m2:
                data["orca_version"] = ver_m2.group(1)
            else:
                data["orca_version"] = None

        data["wall_time"] = _extract_wall_time(content)

        # Charge / multiplicity from header
        chrg_m = re.search(r"Total Charge\s+Charge\s+\.+\s+([+-]?\d+)", content)
        if chrg_m:
            data["charge"] = int(chrg_m.group(1))
        mult_m = re.search(r"Multiplicity\s+Mult\s+\.*\s+(\d+)", content)
        if mult_m:
            data["multiplicity"] = int(mult_m.group(1))

        data.setdefault("charge", None)
        data.setdefault("multiplicity", None)
        data.setdefault("dipole_moment", None)

        return data

    # ------------------------------------------------------------------ #
    # Significance helper                                                  #
    # ------------------------------------------------------------------ #

    def _classify(self, field: str) -> str:
        """Return significance level for *field*."""
        if field in self.CRITICAL_FIELDS:
            return "critical"
        if field in self.MAJOR_FIELDS:
            return "major"
        return "minor"


# ------------------------------------------------------------------ #
# Module-level helpers                                                  #
# ------------------------------------------------------------------ #


def _extract_wall_time(content: str) -> str | None:
    """Extract the TOTAL RUN TIME line from ORCA output."""
    m = re.search(
        r"TOTAL RUN TIME:\s*(\d+\s*days?\s*\d+\s*hours?\s*\d+\s*minutes?\s*[\d.]+\s*seconds?)",
        content,
        re.IGNORECASE,
    )
    if m:
        return m.group(1).strip()
    return None
