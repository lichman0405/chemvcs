"""ORCA input file parser (pure regex / static analysis).

Parses ORCA input scripts (``.inp``) by statically extracting:

* ``!`` keyword lines — method, basis set, run-type, convergence flags,
  dispersion, RI approximations, grids, relativistic options, etc.
* ``%block … end`` blocks — SCF, GEOM, PAL, CPCM, TDDFT, MDCI, etc.
* ``* xyz/int/intxyz/xyzfile charge mult … *`` — coordinate section
* Standalone ``%maxcore <MB>`` directives
* ``%pal nprocs <N> end`` — parallelism (both inline and block forms)

The parser is deliberately **case-insensitive** (ORCA itself ignores
case throughout).  All keyword strings are stored in upper-case.

Multi-step jobs separated by ``$new_job`` are partially supported: only
the first job block is parsed for MVP.

Format reference: ORCA 6.0 User Manual §4.
"""

from __future__ import annotations

import re
from typing import Any

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class OrcaInputParser(BaseParser):
    """Static parser for ORCA input files (``.inp``).

    Extracted data dictionary keys
    --------------------------------
    ``keywords``      – sorted list of upper-case keyword strings from ``!`` lines
    ``blocks``        – dict mapping lower-case block name → dict of key → value
    ``charge``        – molecular charge (int) or None
    ``multiplicity``  – spin multiplicity (int) or None
    ``coord_format``  – coordinate format string (``"xyz"``, ``"int"``, ``"xyzfile"``, …)
    ``atoms``         – list of raw atom lines (str) inside ``* … *`` block
    ``xyzfile``       – external XYZ file path (str) when ``coord_format`` is ``"xyzfile"``
    ``maxcore``       – memory per process in MB (int) or None
    ``nprocs``        – number of MPI processes (int) or None
    """

    # ------------------------------------------------------------------ #
    # Significance classification                                         #
    # ------------------------------------------------------------------ #

    #: Keywords that directly change the computed energy / property.
    CRITICAL_KEYWORDS: frozenset[str] = frozenset({
        # Wavefunction model — HF and post-HF
        "HF", "RHF", "UHF", "ROHF", "CASSCF",
        "MP2", "RI-MP2", "OO-MP2", "SCS-MP2",
        "DLPNO-MP2", "DLPNO-MP2-F12",
        "CCSD", "CCSD(T)", "CCSD-F12", "CCSD(T)-F12",
        "DLPNO-CCSD", "DLPNO-CCSD(T)", "DLPNO-CCSD(T1)",
        # DFT functionals
        "BLYP", "BP86", "PBE", "TPSS", "B97",
        "B3LYP", "PBE0", "TPSSh", "BHLYP", "X3LYP", "B1LYP",
        "B97-3C", "PBEh-3C", "R2SCAN-3C",
        "M06", "M06-2X", "M06-L", "M06-HF",
        "WB97X", "WB97X-D", "WB97X-D3", "WB97M-V", "WB97X-V",
        "CAM-B3LYP", "LC-BLYP", "LC-PBE", "LONG-RANGE-CORRECTED",
        "B2PLYP", "B2GP-PLYP", "DSDPBEP86", "PWPB95",
        "SCAN", "R2SCAN",
        # Run types — these change what property is computed
        "SP", "OPT", "COPT", "ZOPT", "GDIIS-OPT",
        "FREQ", "NUMFREQ", "NFREQ",
        "NMR",
        "SOCOUPLING",
        "MD",
        "IRC",
        "TDDFT",
        # Relativistic
        "DKH", "DKH2", "DKH3",
        "ZORA", "IORA",
        "X2C", "DKH-X2C",
    })

    #: Keywords that affect numerical accuracy / cost without changing the model.
    MAJOR_KEYWORDS: frozenset[str] = frozenset({
        # SCF convergence thresholds
        "SLOPPYSCF", "LOOSESCF", "NORMALSCF", "STRONGSCF",
        "TIGHTSCF", "VERYTIGHTSCF", "EXTREMESCF",
        # Quadrature grids
        "DEFGRID1", "DEFGRID2", "DEFGRID3",
        "GRID1", "GRID2", "GRID3", "GRID4", "GRID5",
        "FINALGRID1", "FINALGRID2", "FINALGRID3", "FINALGRID4", "FINALGRID5", "FINALGRID6",
        "GRIDX5",
        # RI approximations
        "RI", "RI-J", "RIJONX", "RIJK", "RIJCOSX", "NORIJCOSX", "NORI",
        # Dispersion corrections
        "D2", "D3", "D3BJ", "D3ZERO", "D4", "NODAMP", "NL",
        # Solvation
        "CPCM", "SMD",
        # Frozen core / ECP
        "FROZENCORE", "NOFROZENCORE", "SMALLCORE", "NOECP",
        # Memory / parallelism keywords (when specified on ! line)
        "LARGEMEM", "NORMALMEM", "SMALLMEM",
        # Auxiliary basis sets (affect RI accuracy)
        "AUTOAUX", "NOAUTOAUX",
    })

    # ------------------------------------------------------------------ #
    # BaseParser interface                                                  #
    # ------------------------------------------------------------------ #

    def parse(self, content: str) -> dict[str, Any]:
        """Parse ORCA ``.inp`` content into a structured dictionary.

        Args:
            content: Raw ORCA input file text.

        Returns:
            Structured data dictionary (see class docstring).

        Raises:
            ParserError: If the content is empty.
        """
        if not content or not content.strip():
            raise ParserError("ORCA input file is empty")

        # Multi-step jobs: keep only the first job block
        content = re.split(
            r"^\s*\$new_job\b",
            content,
            maxsplit=1,
            flags=re.IGNORECASE | re.MULTILINE,
        )[0]

        # Strip # comments; preserve line structure for block detection
        lines = self._strip_comments(content)
        joined = "\n".join(lines)

        data: dict[str, Any] = {
            "keywords": [],
            "blocks": {},
            "charge": None,
            "multiplicity": None,
            "coord_format": None,
            "atoms": [],
            "maxcore": None,
            "nprocs": None,
        }

        # ---- ! keyword lines ---------------------------------------- #
        keywords: set[str] = set()
        for m in re.finditer(r"^!(.+)$", joined, flags=re.MULTILINE):
            for token in m.group(1).split():
                keywords.add(token.upper())
        data["keywords"] = sorted(keywords)

        # ---- Standalone single-line directives ----------------------- #
        # %maxcore 4000
        for m in re.finditer(r"^%maxcore\s+(\d+)", joined, flags=re.MULTILINE | re.IGNORECASE):
            data["maxcore"] = int(m.group(1))
        # %pal nprocs N end   (inline form, single line)
        for m in re.finditer(
            r"^%pal\s+nprocs\s+(\d+)\s+end\b",
            joined,
            flags=re.MULTILINE | re.IGNORECASE,
        ):
            data["nprocs"] = int(m.group(1))

        # ---- %block … end blocks ------------------------------------ #
        for block_m in re.finditer(
            r"^%(\w+)(.*?)^end\b",
            joined,
            flags=re.MULTILINE | re.DOTALL | re.IGNORECASE,
        ):
            block_name = block_m.group(1).lower()
            block_body = block_m.group(2)

            kv: dict[str, Any] = {}
            for kv_m in re.finditer(r"^\s*(\w+)\s+([^\n]+)", block_body, flags=re.MULTILINE):
                k = kv_m.group(1).lower()
                v_raw = kv_m.group(2).strip()
                kv[k] = _try_numeric(v_raw)

            data["blocks"][block_name] = kv

            # Mirror %pal nprocs into top-level field
            if block_name == "pal" and "nprocs" in kv:
                data["nprocs"] = kv["nprocs"]

        # ---- Coordinate block --------------------------------------- #
        # * xyz charge mult
        #   atom lines …
        # *
        coord_m = re.search(
            r"^\*\s+(xyz|int|intxyz|gzmt)\s+(-?\d+)\s+(\d+)\s*\n(.*?)^\*\s*$",
            joined,
            flags=re.MULTILINE | re.DOTALL | re.IGNORECASE,
        )
        if coord_m:
            data["coord_format"] = coord_m.group(1).lower()
            data["charge"] = int(coord_m.group(2))
            data["multiplicity"] = int(coord_m.group(3))
            raw_atoms = [ln.strip() for ln in coord_m.group(4).splitlines() if ln.strip()]
            data["atoms"] = raw_atoms

        # * xyzfile charge mult filename   (no coordinate body)
        if data["charge"] is None:
            xyzfile_m = re.search(
                r"^\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)",
                joined,
                flags=re.MULTILINE | re.IGNORECASE,
            )
            if xyzfile_m:
                data["coord_format"] = "xyzfile"
                data["charge"] = int(xyzfile_m.group(1))
                data["multiplicity"] = int(xyzfile_m.group(2))
                data["xyzfile"] = xyzfile_m.group(3)

        return data

    def diff(
        self,
        old_data: dict[str, Any],
        new_data: dict[str, Any],
    ) -> list[DiffEntry]:
        """Compute semantic differences between two parsed ORCA inputs.

        Significance mapping:
        * **critical** — method, basis set, run type, charge, multiplicity,
          relativistic options, coordinate changes
        * **major**    — SCF convergence, grid, RI, dispersion, MaxCore,
          nprocs, CPCM settings
        * **minor**    — output verbosity, guess, symmetry, print options
        """
        entries: list[DiffEntry] = []

        # ---- Keyword-level diff ------------------------------------- #
        old_kw = set(old_data.get("keywords", []))
        new_kw = set(new_data.get("keywords", []))

        for kw in sorted(old_kw - new_kw):
            entries.append(
                DiffEntry(f"keyword.{kw}", kw, None, "deleted", self._classify_kw(kw))
            )
        for kw in sorted(new_kw - old_kw):
            entries.append(
                DiffEntry(f"keyword.{kw}", None, kw, "added", self._classify_kw(kw))
            )

        # ---- Scalar top-level fields -------------------------------- #
        _SCALAR_SIG: dict[str, str] = {
            "charge": "critical",
            "multiplicity": "critical",
            "coord_format": "critical",
            "maxcore": "major",
            "nprocs": "major",
        }
        for field, sig in _SCALAR_SIG.items():
            old_val = old_data.get(field)
            new_val = new_data.get(field)
            if old_val != new_val:
                change_type = (
                    "added" if old_val is None
                    else "deleted" if new_val is None
                    else "modified"
                )
                entries.append(DiffEntry(field, old_val, new_val, change_type, sig))

        # ---- Coordinate (atom) changes ------------------------------ #
        old_atoms = old_data.get("atoms", [])
        new_atoms = new_data.get("atoms", [])
        if old_atoms != new_atoms:
            entries.append(
                DiffEntry(
                    "coordinates",
                    f"{len(old_atoms)} atom(s)",
                    f"{len(new_atoms)} atom(s)",
                    "modified",
                    "critical",
                )
            )

        # ---- Block diffs -------------------------------------------- #
        old_blocks: dict[str, dict[str, Any]] = old_data.get("blocks", {})
        new_blocks: dict[str, dict[str, Any]] = new_data.get("blocks", {})
        all_block_names = set(old_blocks) | set(new_blocks)

        for block_name in sorted(all_block_names):
            old_block = old_blocks.get(block_name, {})
            new_block = new_blocks.get(block_name, {})
            all_keys = set(old_block) | set(new_block)

            for k in sorted(all_keys):
                old_val = old_block.get(k)
                new_val = new_block.get(k)
                if old_val != new_val:
                    sig = self._classify_block(block_name, k)
                    path = f"block.{block_name}.{k}"
                    change_type = (
                        "added" if old_val is None
                        else "deleted" if new_val is None
                        else "modified"
                    )
                    entries.append(DiffEntry(path, old_val, new_val, change_type, sig))

        return entries

    def validate(self, data: dict[str, Any]) -> tuple[bool, list[str]]:
        """Basic sanity checks on the parsed ORCA input.

        Args:
            data: Parsed data dictionary from :meth:`parse`.

        Returns:
            ``(is_valid, list_of_error_messages)``
        """
        errors: list[str] = []

        if not data.get("keywords"):
            errors.append(
                "No method/basis keywords found — expected at least one '!' line"
            )

        if data.get("charge") is None:
            errors.append(
                "No coordinate block found — expected '* xyz/int/xyzfile charge mult … *'"
            )

        if data.get("multiplicity") is not None and data["multiplicity"] < 1:
            errors.append("Spin multiplicity must be ≥ 1")

        return len(errors) == 0, errors

    # ------------------------------------------------------------------ #
    # Private helpers                                                      #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _strip_comments(content: str) -> list[str]:
        """Remove ``#``-prefixed comments from each line.

        Full-line ``#`` comments and inline ``# …`` suffixes are both removed.
        The line count is preserved so that multi-line regex patterns
        (``re.MULTILINE``) still match correctly.
        """
        lines = []
        for line in content.splitlines():
            stripped = re.sub(r"\s*#.*$", "", line)
            lines.append(stripped)
        return lines

    def _classify_kw(self, kw: str) -> str:
        """Return significance level for an ORCA keyword string."""
        ku = kw.upper()
        if ku in self.CRITICAL_KEYWORDS:
            return "critical"
        if ku in self.MAJOR_KEYWORDS:
            return "major"
        # Pattern-based fallback — recognize common basis set prefixes
        if re.match(
            r"(DEF2-|AUG-CC-P|CC-P|PC-\d|MINIX|STO-\d|6-3\d|MA-DEF2-|SARC-)"
            r"|(-SVP|-TZVP|-QZVP|-DZ|-TZ|-QZ)",
            ku,
        ):
            return "critical"
        # Pattern-based method recognition
        if re.match(r"(B\d|PBE\d?|M\d{2}[-\d]?|WB97|LC-|CAM-|RANGE-)", ku):
            return "critical"
        return "minor"

    @staticmethod
    def _classify_block(block_name: str, key: str) -> str:
        """Return significance for a ``%block`` key."""
        bn = block_name.lower()
        kl = key.lower()
        if bn == "pal" and kl == "nprocs":
            return "major"
        if bn == "scf" and kl in ("convergence", "maxiter", "thresh", "tcut"):
            return "major"
        if bn == "cpcm" and kl in ("epsilon", "refrac"):
            return "major"
        if bn == "tddft" and kl in ("nroots", "iroot"):
            return "major"
        if bn == "mdci" and kl in ("tcutpno", "tcutpairs", "locmet"):
            return "major"
        return "minor"


# ------------------------------------------------------------------ #
# Module-level helpers                                                  #
# ------------------------------------------------------------------ #

def _try_numeric(s: str) -> Any:
    """Attempt to parse *s* as int, then float; return str on failure."""
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return s
