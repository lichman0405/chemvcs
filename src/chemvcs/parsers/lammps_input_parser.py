"""LAMMPS input script parser (static analysis only).

Parses LAMMPS input scripts (``in.*``, ``*.lammps``, ``lammps.in``) by
performing *static* extraction of key commands.  The LAMMPS scripting
language supports variables, loops, conditionals, and ``include``
directives; full interpretation would require an embedded LAMMPS engine.
This parser deliberately limits itself to what can be determined without
executing the script, providing useful semantics for version-tracking.

Format reference: LAMMPS Manual §5 (22 Jul 2025 edition).

Commands extracted (significance):
  CRITICAL
    * ``units``         – unit system (real, metal, lj, si, cgs, electron…)
    * ``pair_style``    – interatomic potential + cutoff
    * ``kspace_style``  – long-range solver (ewald, pppm…)
    * ``timestep``      – integration timestep
    * ``run``           – number of MD steps (last occurrence)
    * ``minimize``      – minimisation convergence (etol, ftol, N1, N2)

  MAJOR
    * ``boundary``      – boundary conditions (p/s/f for each of x/y/z)
    * ``atom_style``    – particle model (atomic, charge, full, molecular…)
    * ``fix_nvt``       – NVT thermostat settings (temp, Tdamp)
    * ``fix_npt``       – NPT barostat settings (temp, press, Pdamp)
    * ``fix_nve``       – NVE integrator (boolean presence)
    * ``dimension``     – system dimensionality (2 or 3)

  MINOR
    * ``thermo``        – output interval
    * ``dump_interval`` – dump output interval (first dump command)
    * ``restart``       – restart write interval (first restart command)
    * ``neigh_modify``  – neighbour-list settings
    * ``read_data``     – data file path(s) referenced
    * ``pair_coeff``    – first pair_coeff line (summary)
    * ``mass``          – mass definitions as dict {type: mass}
    * ``variable_names``– list of variable names defined

Notes
-----
* Lines beginning with ``#`` are comments; inline comments are stripped.
* Continuation lines ending with ``&`` are joined before parsing.
* Variable expansion (``${var}`` / ``$var``) is NOT performed — values
  containing unresolved variables are stored as raw strings.
"""

from __future__ import annotations

import contextlib
import shlex
from typing import Any

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class LammpsInputParser(BaseParser):
    """Static parser for LAMMPS input scripts.

    Extracts a structured summary suitable for semantic version-tracking.
    Does not execute the script or resolve LAMMPS variables.
    """

    # ------------------------------------------------------------------ #
    # Significance classification                                          #
    # ------------------------------------------------------------------ #

    CRITICAL_COMMANDS = frozenset(
        {
            "units",
            "pair_style",
            "kspace_style",
            "timestep",
            "run",
            "minimize",
        }
    )

    MAJOR_COMMANDS = frozenset(
        {
            "boundary",
            "atom_style",
            "dimension",
            "fix_nvt",
            "fix_npt",
            "fix_nve",
        }
    )

    # ------------------------------------------------------------------ #
    # BaseParser interface                                                  #
    # ------------------------------------------------------------------ #

    def parse(self, content: str) -> dict[str, Any]:
        """Parse LAMMPS input script and extract key settings.

        Args:
            content: Full text of the input script.

        Returns:
            Dictionary of extracted settings.

        Raises:
            ParserError: If content is completely unparseable.
        """
        if not content or not content.strip():
            raise ParserError("LAMMPS input script is empty")

        lines = self._preprocess(content)

        data: dict[str, Any] = {
            "read_data": [],
            "mass": {},
            "variable_names": [],
        }
        fix_list: list[dict[str, Any]] = []

        for line in lines:
            self._parse_line(line, data, fix_list)

        # Derive fix-specific summaries
        self._summarize_fixes(fix_list, data)

        return data

    def diff(
        self,
        old_data: dict[str, Any],
        new_data: dict[str, Any],
    ) -> list[DiffEntry]:
        """Compute semantic diff between two parsed LAMMPS input scripts."""
        entries: list[DiffEntry] = []

        scalar_fields = [
            "units",
            "pair_style",
            "kspace_style",
            "timestep",
            "run",
            "minimize_etol",
            "minimize_ftol",
            "minimize_maxiter",
            "minimize_maxeval",
            "boundary",
            "atom_style",
            "dimension",
            "fix_nvt",
            "fix_npt",
            "fix_nve",
            "thermo",
            "pair_coeff",
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
                # Suppress float noise for timestep
                if (
                    isinstance(old_val, float)
                    and isinstance(new_val, float)
                    and abs(old_val - new_val) < 1e-15
                ):
                    continue
                entries.append(DiffEntry(field, old_val, new_val, "modified", significance))

        # List-valued fields: read_data, variable_names
        for list_field, sig in [("read_data", "minor"), ("variable_names", "minor")]:
            old_list = set(old_data.get(list_field, []))
            new_list = set(new_data.get(list_field, []))
            if old_list != new_list:
                entries.append(
                    DiffEntry(
                        list_field,
                        sorted(old_list),
                        sorted(new_list),
                        "modified",
                        sig,
                    )
                )

        # Mass dict diff
        old_mass = old_data.get("mass", {})
        new_mass = new_data.get("mass", {})
        if old_mass != new_mass:
            entries.append(DiffEntry("mass", old_mass, new_mass, "modified", "minor"))

        return entries

    def validate(self, data: dict[str, Any]) -> tuple[bool, list[str]]:
        """Validate parsed LAMMPS input script."""
        errors: list[str] = []

        if "units" not in data:
            errors.append("'units' command not found – required for all simulations")

        if "pair_style" not in data:
            errors.append("'pair_style' command not found")

        ts = data.get("timestep")
        if ts is not None:
            try:
                if float(ts) <= 0:
                    errors.append(f"timestep must be positive (got {ts})")
            except (TypeError, ValueError):
                pass  # variable-containing timestep, skip

        return len(errors) == 0, errors

    # ------------------------------------------------------------------ #
    # Internal helpers                                                     #
    # ------------------------------------------------------------------ #

    def _classify(self, field: str) -> str:
        # Map field names back to command names for classification
        cmd_for_field = {
            "units": "units",
            "pair_style": "pair_style",
            "kspace_style": "kspace_style",
            "timestep": "timestep",
            "run": "run",
            "minimize_etol": "minimize",
            "minimize_ftol": "minimize",
            "minimize_maxiter": "minimize",
            "minimize_maxeval": "minimize",
            "boundary": "boundary",
            "atom_style": "atom_style",
            "dimension": "dimension",
            "fix_nvt": "fix_nvt",
            "fix_npt": "fix_npt",
            "fix_nve": "fix_nve",
        }
        cmd = cmd_for_field.get(field, field)
        if cmd in self.CRITICAL_COMMANDS:
            return "critical"
        if cmd in self.MAJOR_COMMANDS:
            return "major"
        return "minor"

    # -- Pre-processing ------------------------------------------------- #

    @staticmethod
    def _preprocess(content: str) -> list[str]:
        """Strip comments, join continuation lines, return clean lines."""
        result: list[str] = []
        pending = ""

        for raw_line in content.splitlines():
            # Strip inline comments (but not inside quotes – simplified approach)
            line = raw_line.strip()

            # Remove full-line comments
            if line.startswith("#"):
                continue

            # Remove inline comments
            comment_pos = -1
            in_single = False
            in_double = False
            for i, ch in enumerate(line):
                if ch == "'" and not in_double:
                    in_single = not in_single
                elif ch == '"' and not in_single:
                    in_double = not in_double
                elif ch == "#" and not in_single and not in_double:
                    comment_pos = i
                    break
            if comment_pos >= 0:
                line = line[:comment_pos].strip()

            if not line:
                continue

            # Continuation character
            if line.endswith("&"):
                pending += line[:-1].rstrip() + " "
                continue

            full_line = (pending + line).strip()
            pending = ""
            if full_line:
                result.append(full_line)

        if pending.strip():
            result.append(pending.strip())

        return result

    # -- Line dispatcher ------------------------------------------------ #

    def _parse_line(
        self,
        line: str,
        data: dict[str, Any],
        fix_list: list[dict[str, Any]],
    ) -> None:
        """Dispatch a single cleaned line to the appropriate handler."""
        try:
            tokens = self._tokenize(line)
        except ValueError:
            return

        if not tokens:
            return

        cmd = tokens[0].lower()

        if cmd == "units" and len(tokens) >= 2:
            data["units"] = tokens[1]

        elif cmd == "atom_style" and len(tokens) >= 2:
            data["atom_style"] = tokens[1]

        elif cmd == "dimension" and len(tokens) >= 2:
            try:
                data["dimension"] = int(tokens[1])
            except ValueError:
                data["dimension"] = tokens[1]

        elif cmd == "boundary" and len(tokens) >= 2:
            data["boundary"] = " ".join(tokens[1:4])  # e.g. "p p p"

        elif cmd == "pair_style" and len(tokens) >= 2:
            # Store full specification for semantic diff
            data["pair_style"] = " ".join(tokens[1:])

        elif cmd == "pair_coeff" and "pair_coeff" not in data:
            # Store only the first pair_coeff as a summary
            data["pair_coeff"] = " ".join(tokens[1:])

        elif cmd == "kspace_style" and len(tokens) >= 2:
            data["kspace_style"] = " ".join(tokens[1:])

        elif cmd == "timestep" and len(tokens) >= 2:
            raw = tokens[1]
            try:
                data["timestep"] = float(raw)
            except ValueError:
                data["timestep"] = raw  # keep raw if variable

        elif cmd == "run" and len(tokens) >= 2:
            # Keep the last run command (most representative)
            raw = tokens[1]
            try:
                data["run"] = int(raw)
            except ValueError:
                data["run"] = raw

        elif cmd == "minimize" and len(tokens) >= 5:
            try:
                data["minimize_etol"] = float(tokens[1])
                data["minimize_ftol"] = float(tokens[2])
                data["minimize_maxiter"] = int(tokens[3])
                data["minimize_maxeval"] = int(tokens[4])
            except (ValueError, IndexError):
                data["minimize_etol"] = tokens[1] if len(tokens) > 1 else None

        elif cmd == "thermo" and len(tokens) >= 2:
            try:
                data["thermo"] = int(tokens[1])
            except ValueError:
                data["thermo"] = tokens[1]

        elif cmd == "variable" and len(tokens) >= 2:
            vname = tokens[1]
            if vname not in data["variable_names"]:
                data["variable_names"].append(vname)

        elif cmd == "read_data" and len(tokens) >= 2:
            fpath = tokens[1]
            if fpath not in data["read_data"]:
                data["read_data"].append(fpath)

        elif cmd == "mass" and len(tokens) >= 3:
            try:
                atom_type = tokens[1]
                mass_val = float(tokens[2])
                data["mass"][atom_type] = mass_val
            except ValueError:
                pass

        elif cmd == "fix" and len(tokens) >= 4:
            fix_info: dict[str, Any] = {
                "fix_id": tokens[1],
                "group": tokens[2],
                "style": tokens[3].lower(),
                "args": tokens[4:],
            }
            fix_list.append(fix_info)

        elif cmd == "dump" and len(tokens) >= 5:
            # dump ID group style N file
            if "dump_interval" not in data:
                with contextlib.suppress(ValueError, IndexError):
                    data["dump_interval"] = int(tokens[4])

        elif cmd == "restart" and len(tokens) >= 2:
            if "restart_interval" not in data:
                with contextlib.suppress(ValueError):
                    data["restart_interval"] = int(tokens[1])

        elif cmd == "neigh_modify" and "neigh_modify" not in data:
            data["neigh_modify"] = " ".join(tokens[1:])

    # -- Fix summarization --------------------------------------------- #

    @staticmethod
    def _summarize_fixes(
        fix_list: list[dict[str, Any]],
        data: dict[str, Any],
    ) -> None:
        """Extract NVT/NPT/NVE ensemble settings from fix commands."""
        for fix in fix_list:
            style = fix["style"]
            args = fix["args"]

            if style in ("nvt", "nvt/omp", "nvt/gpu", "nvt/kk"):
                summary = _parse_fix_nvt_npt(style, args)
                data.setdefault("fix_nvt", summary)

            elif style in ("npt", "npt/omp", "npt/gpu", "npt/kk", "npt/sphere", "npt/asphere"):
                summary = _parse_fix_nvt_npt(style, args)
                data.setdefault("fix_npt", summary)

            elif style in ("nve", "nve/omp", "nve/gpu", "nve/kk", "nve/sphere", "nve/asphere"):
                data.setdefault("fix_nve", True)

            elif style in ("langevin",):
                # Langevin thermostat: fix ID group langevin Tstart Tstop damp seed
                if len(args) >= 3:
                    summary = {
                        "style": style,
                        "Tstart": _maybe_float(args[0]),
                        "Tstop": _maybe_float(args[1]),
                        "damp": _maybe_float(args[2]),
                    }
                    data.setdefault("fix_nvt", summary)

    # -- Tokenizer ----------------------------------------------------- #

    @staticmethod
    def _tokenize(line: str) -> list[str]:
        """Tokenize a LAMMPS line, respecting quoted strings."""
        try:
            return shlex.split(line)
        except ValueError:
            # shlex fails on unbalanced quotes (common with $-variables)
            return line.split()


# ------------------------------------------------------------------ #
# Module-level helpers                                               #
# ------------------------------------------------------------------ #


def _maybe_float(s: str) -> Any:
    """Try to convert to float, return raw string on failure."""
    try:
        return float(s)
    except (ValueError, TypeError):
        return s


def _parse_fix_nvt_npt(style: str, args: list[str]) -> dict[str, Any]:
    """Extract temp/press/damp from NVT or NPT fix args."""
    summary: dict[str, Any] = {"style": style}
    i = 0
    while i < len(args):
        kw = args[i].lower()
        if kw == "temp" and i + 3 < len(args):
            summary["Tstart"] = _maybe_float(args[i + 1])
            summary["Tstop"] = _maybe_float(args[i + 2])
            summary["Tdamp"] = _maybe_float(args[i + 3])
            i += 4
        elif kw in ("iso", "aniso", "tri", "x", "y", "z") and i + 3 < len(args):
            summary.setdefault("press_style", kw)
            summary["Pstart"] = _maybe_float(args[i + 1])
            summary["Pstop"] = _maybe_float(args[i + 2])
            summary["Pdamp"] = _maybe_float(args[i + 3])
            i += 4
        else:
            i += 1
    return summary
