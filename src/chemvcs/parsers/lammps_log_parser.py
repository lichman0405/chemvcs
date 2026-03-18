"""LAMMPS log file parser.

Parses ``log.lammps`` output files produced by LAMMPS MD simulations.
Supports both the traditional columnar thermo output and the YAML
format (available since LAMMPS 24Mar2022 with ``thermo_style yaml`` or
``thermo_modify line yaml``).

For incomplete log files (common during long HPC runs) a regex fallback
extracts the most recent thermodynamic values found in the file.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class LammpsLogParser(BaseParser):
    """Parser for LAMMPS ``log.lammps`` files.

    Extracts the following quantities for version-tracking purposes:

    * ``final_step``        – last timestep recorded
    * ``final_temp``        – final temperature (K)
    * ``final_press``       – final pressure (atm or bar, unit-set dependent)
    * ``final_pe``          – final potential energy
    * ``final_ke``          – final kinetic energy
    * ``final_etotal``      – final total energy
    * ``run_steps``         – total steps executed (sum of all ``run`` commands)
    * ``num_atoms``         – number of atoms (from header)
    * ``thermo_keywords``   – list of columns present in thermo output
    * ``completed``         – True if ``Total wall time:`` line was found
    * ``wall_time``         – wall-clock time string, e.g. ``"0:02:34"``
    """

    # ------------------------------------------------------------------ #
    # Significance classification                                          #
    # ------------------------------------------------------------------ #

    CRITICAL_FIELDS = frozenset({
        "final_etotal", "final_pe", "completed",
    })

    MAJOR_FIELDS = frozenset({
        "final_temp", "final_press", "run_steps", "num_atoms",
    })

    # ------------------------------------------------------------------ #
    # BaseParser interface                                                  #
    # ------------------------------------------------------------------ #

    def parse(self, content: str) -> Dict[str, Any]:
        """Parse LAMMPS log content.

        Tries YAML-format parsing first (modern LAMMPS), then falls back
        to columnar thermo parsing, then to a bare-minimum regex scan.

        Args:
            content: Full text of ``log.lammps``.

        Returns:
            Dictionary of extracted quantities.

        Raises:
            ParserError: If no thermodynamic data can be found at all.
        """
        if not content or not content.strip():
            raise ParserError("LAMMPS log file is empty")

        # 1. Try YAML thermo format (LAMMPS >= 24Mar2022)
        yaml_data = self._try_parse_yaml_thermo(content)
        if yaml_data:
            base = self._extract_base_fields(content)
            base.update(yaml_data)
            return base

        # 2. Try classic columnar thermo format
        col_data = self._try_parse_columnar_thermo(content)
        if col_data:
            base = self._extract_base_fields(content)
            base.update(col_data)
            return base

        # 3. Bare-minimum regex scan
        regex_data = self._try_regex_fallback(content)
        if regex_data:
            return regex_data

        raise ParserError(
            "Failed to parse LAMMPS log: no thermodynamic data found. "
            "File may be empty or not a valid LAMMPS log."
        )

    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic diff between two parsed LAMMPS log snapshots."""
        entries: List[DiffEntry] = []

        tracked_fields = [
            "final_etotal",
            "final_pe",
            "final_ke",
            "final_temp",
            "final_press",
            "run_steps",
            "num_atoms",
            "completed",
            "wall_time",
        ]

        for field in tracked_fields:
            old_val = old_data.get(field)
            new_val = new_data.get(field)
            significance = self._classify(field)

            if old_val is None and new_val is not None:
                entries.append(DiffEntry(field, None, new_val, "added", significance))
            elif old_val is not None and new_val is None:
                entries.append(DiffEntry(field, old_val, None, "deleted", significance))
            elif old_val != new_val:
                # Suppress float noise for energy fields
                if (
                    isinstance(old_val, float)
                    and isinstance(new_val, float)
                    and abs(old_val - new_val) < 1e-10
                ):
                    continue
                entries.append(DiffEntry(field, old_val, new_val, "modified", significance))

        return entries

    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate parsed LAMMPS log data."""
        errors: List[str] = []

        if not data.get("completed", False):
            errors.append("LAMMPS run did not complete normally (no 'Total wall time:' found)")

        if data.get("num_atoms") is not None and data["num_atoms"] <= 0:
            errors.append(f"Invalid atom count: {data['num_atoms']}")

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

    # -- YAML thermo format -------------------------------------------- #

    def _try_parse_yaml_thermo(self, content: str) -> Optional[Dict[str, Any]]:
        """Parse YAML-style thermo blocks (LAMMPS >= 24Mar2022).

        Expects blocks like::

            ---
            keywords: ['Step', 'Temp', 'TotEng', ...]
            data:
            - [100, 300.0, -1234.5, ...]
            - [200, 298.5, -1236.2, ...]
            ...

        Returns merged data from all blocks, using the *last* data row
        as the final state.
        """
        try:
            import yaml  # type: ignore[import-untyped]
        except ImportError:
            return None

        # Extract YAML documents delimited by '---' and either '...'
        # or the next '---'/EOF. LAMMPS commonly writes a single thermo
        # document followed by ordinary log lines such as wall time.
        yaml_blocks: List[str] = []
        current_block: List[str] = []
        in_block = False

        for line in content.splitlines():
            stripped = line.strip()

            if stripped == "---":
                if in_block and current_block:
                    yaml_blocks.append("\n".join(current_block))
                    current_block = []
                in_block = True
                continue

            if not in_block:
                continue

            if stripped == "...":
                yaml_blocks.append("\n".join(current_block))
                current_block = []
                in_block = False
                continue

            current_block.append(line)

        if in_block and current_block:
            yaml_blocks.append("\n".join(current_block))

        if not yaml_blocks:
            return None

        run_steps_total = 0
        last_row: Optional[Dict[str, float]] = None
        keywords: Optional[List[str]] = None

        for block in yaml_blocks:
            try:
                doc = yaml.safe_load("---\n" + block)
            except Exception:
                continue

            if not isinstance(doc, dict):
                continue
            if "keywords" not in doc or "data" not in doc:
                continue

            kws: List[str] = doc["keywords"]
            rows: List[List[Any]] = doc["data"]

            if not kws or not rows:
                continue

            keywords = kws
            last_row = dict(zip(kws, rows[-1]))

            # Accumulate run steps
            if "Step" in last_row and len(rows) > 1:
                first_step = rows[0][kws.index("Step")] if "Step" in kws else 0
                last_step = rows[-1][kws.index("Step")]
                run_steps_total += int(last_step) - int(first_step)

        if last_row is None:
            return None

        return self._rows_to_final(last_row, run_steps_total, keywords)

    # -- Columnar thermo format ----------------------------------------- #

    def _try_parse_columnar_thermo(self, content: str) -> Optional[Dict[str, Any]]:
        """Parse classic columnar thermo format.

        LAMMPS prints a header line like::

            Step Temp E_pair E_mol TotEng Press

        followed by numeric data rows. Multiple ``run``/``minimize``
        blocks are supported.
        """
        # Pattern: LAMMPS thermo header always starts with "Step" (case-sensitive)
        # followed by one or more column names (letters, digits, _, /, .).
        # This is intentionally broad to support user-defined thermo styles.
        header_re = re.compile(
            r"^(\s*Step(?:\s+[\w./]+)+)\s*$",
            re.MULTILINE,
        )

        run_steps_total = 0
        last_row: Optional[Dict[str, float]] = None
        keywords: Optional[List[str]] = None

        for match in header_re.finditer(content):
            header_line = match.group(0)
            kws = header_line.split()
            start = match.end()

            # Collect subsequent numeric rows
            rows: List[List[float]] = []
            for line in content[start:].splitlines():
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                try:
                    row = [float(p) for p in parts]
                    if len(row) == len(kws):
                        rows.append(row)
                    else:
                        break
                except ValueError:
                    break

            if not rows:
                continue

            keywords = kws
            last_row = dict(zip(kws, rows[-1]))

            if len(rows) > 1 and "Step" in kws:
                step_idx = kws.index("Step")
                run_steps_total += int(rows[-1][step_idx]) - int(rows[0][step_idx])

        if last_row is None:
            return None

        return self._rows_to_final(last_row, run_steps_total, keywords)

    # -- Regex fallback ------------------------------------------------- #

    def _try_regex_fallback(self, content: str) -> Optional[Dict[str, Any]]:
        """Last-resort regex scan for partial log files."""
        data = self._extract_base_fields(content)

        # Last occurrence of common energy patterns in per-step output
        toteng_matches = re.findall(r"(?:TotEng|E_total)\s+([-+]?\d+\.\d+)", content)
        if toteng_matches:
            data["final_etotal"] = float(toteng_matches[-1])

        temp_matches = re.findall(r"Temp\s+([\d.]+)", content)
        if temp_matches:
            data["final_temp"] = float(temp_matches[-1])

        if not toteng_matches and not temp_matches:
            return None

        return data

    # -- Shared field extraction ---------------------------------------- #

    def _extract_base_fields(self, content: str) -> Dict[str, Any]:
        """Extract fields common to all parse paths (num_atoms, wall_time...)."""
        data: Dict[str, Any] = {}

        # Number of atoms
        m = re.search(r"(\d+)\s+atoms$", content, re.MULTILINE)
        if m:
            data["num_atoms"] = int(m.group(1))

        # Wall time
        m = re.search(r"Total wall time:\s*(\S+)", content)
        if m:
            data["wall_time"] = m.group(1)
            data["completed"] = True
        else:
            data["completed"] = False

        return data

    @staticmethod
    def _rows_to_final(
        last_row: Dict[str, Any],
        run_steps_total: int,
        keywords: Optional[List[str]],
    ) -> Dict[str, Any]:
        """Convert last thermo row to our canonical dict."""
        data: Dict[str, Any] = {}

        # Map LAMMPS column names → our field names
        col_map = {
            "Step":    "final_step",
            "Temp":    "final_temp",
            "Press":   "final_press",
            "TotEng":  "final_etotal",
            "PotEng":  "final_pe",
            "E_pair":  "final_pe",     # fallback if PotEng absent
            "KinEng":  "final_ke",
        }

        for col, field in col_map.items():
            if col in last_row and field not in data:
                val = last_row[col]
                data[field] = float(val) if val is not None else None

        if run_steps_total > 0:
            data["run_steps"] = run_steps_total

        if keywords:
            data["thermo_keywords"] = keywords

        return data
