"""KPOINTS file parser backed by pymatgen.

Uses pymatgen.io.vasp.inputs.Kpoints for robust parsing of all VASP
k-point formats: Automatic (Gamma / Monkhorst-Pack), Explicit, and
Line-mode for band-structure calculations.
"""

from typing import Any, Dict, List, Tuple

from pymatgen.io.vasp.inputs import Kpoints as PmgKpoints

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class KpointsParser(BaseParser):
    """Parser for VASP KPOINTS files, powered by pymatgen.

    Pymatgen's ``Kpoints`` class handles:

    * Automatic mesh (Gamma-centred and Monkhorst-Pack)
    * Explicit reciprocal / Cartesian k-point lists with weights
    * Line-mode for band-structure paths
    * All abbreviations (``G``, ``M``, ``Line``, etc.)

    The parsed result is normalised to a plain ``dict`` with a ``"type"``
    key (``"automatic"`` / ``"explicit"`` / ``"line"``), suitable for
    semantic diff.
    """

    # ---------- BaseParser interface ----------

    def parse(self, content: str) -> Dict[str, Any]:
        """Parse KPOINTS content using pymatgen.

        Args:
            content: Raw KPOINTS file content.

        Returns:
            Normalised dictionary describing the k-point scheme.

        Raises:
            ParserError: If pymatgen cannot parse *content*.
        """
        try:
            kpts = PmgKpoints.from_str(content)
            return self._to_dict(kpts)
        except Exception as e:
            raise ParserError(f"Failed to parse KPOINTS: {e}") from e

    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic diff between two parsed KPOINTS dictionaries."""
        diff_entries: List[DiffEntry] = []

        old_type = old_data.get("type")
        new_type = new_data.get("type")

        if old_type != new_type:
            diff_entries.append(
                DiffEntry("type", old_type, new_type, "modified", "critical")
            )
            return diff_entries

        if old_type == "automatic":
            diff_entries.extend(self._diff_automatic(old_data, new_data))
        elif old_type == "explicit":
            diff_entries.extend(self._diff_explicit(old_data, new_data))
        elif old_type == "line":
            diff_entries.extend(self._diff_line(old_data, new_data))

        return diff_entries

    # ---------- internal conversion ----------

    def _to_dict(self, kpts: PmgKpoints) -> Dict[str, Any]:
        """Convert a pymatgen ``Kpoints`` object to our canonical dict."""
        style = kpts.style.name  # "Gamma", "Monkhorst", "Reciprocal", "Cartesian", "Line_mode"

        if style in ("Gamma", "Monkhorst"):
            return self._auto_dict(kpts, style)
        if style == "Line_mode":
            return self._linemode_dict(kpts)
        # Reciprocal / Cartesian → explicit
        return self._explicit_dict(kpts, style)

    # -- automatic --

    def _auto_dict(self, kpts: PmgKpoints, style: str) -> Dict[str, Any]:
        grid = [int(g) for g in kpts.kpts[0]] if kpts.kpts else [1, 1, 1]
        shift = list(kpts.kpts_shift) if kpts.kpts_shift else [0.0, 0.0, 0.0]
        return {
            "type": "automatic",
            "comment": kpts.comment or "",
            "gamma_centered": style == "Gamma",
            "grid": grid,
            "shift": shift,
        }

    # -- explicit --

    def _explicit_dict(self, kpts: PmgKpoints, style: str) -> Dict[str, Any]:
        weights = kpts.kpts_weights or []
        kpoints = []
        for i, kpt in enumerate(kpts.kpts):
            w = weights[i] if i < len(weights) else 1.0
            kpoints.append({"coords": list(kpt), "weight": w})
        return {
            "type": "explicit",
            "comment": kpts.comment or "",
            "num_kpoints": kpts.num_kpts,
            "cartesian": style == "Cartesian",
            "kpoints": kpoints,
        }

    # -- line-mode --

    def _linemode_dict(self, kpts: PmgKpoints) -> Dict[str, Any]:
        segments: List[Dict[str, Any]] = []
        kpt_list = kpts.kpts or []
        labels = kpts.labels or []

        for i in range(0, len(kpt_list) - 1, 2):
            start_label = labels[i] if i < len(labels) else ""
            end_label = labels[i + 1] if (i + 1) < len(labels) else ""
            segments.append({
                "start": {
                    "coords": list(kpt_list[i]),
                    "label": start_label,
                },
                "end": {
                    "coords": list(kpt_list[i + 1]),
                    "label": end_label,
                },
            })

        return {
            "type": "line",
            "comment": kpts.comment or "",
            "num_points": kpts.num_kpts,
            "reciprocal": True,
            "segments": segments,
        }

    # ---------- diff helpers ----------

    def _diff_automatic(
        self, old: Dict[str, Any], new: Dict[str, Any]
    ) -> List[DiffEntry]:
        entries: List[DiffEntry] = []

        if old.get("gamma_centered") != new.get("gamma_centered"):
            entries.append(
                DiffEntry(
                    "gamma_centered",
                    old.get("gamma_centered"),
                    new.get("gamma_centered"),
                    "modified",
                    "major",
                )
            )

        if old.get("grid") != new.get("grid"):
            entries.append(
                DiffEntry(
                    "grid",
                    old.get("grid"),
                    new.get("grid"),
                    "modified",
                    "critical",
                )
            )

        if old.get("shift") != new.get("shift"):
            entries.append(
                DiffEntry(
                    "shift",
                    old.get("shift"),
                    new.get("shift"),
                    "modified",
                    "major",
                )
            )

        return entries

    def _diff_explicit(
        self, old: Dict[str, Any], new: Dict[str, Any]
    ) -> List[DiffEntry]:
        entries: List[DiffEntry] = []

        if old.get("num_kpoints") != new.get("num_kpoints"):
            entries.append(
                DiffEntry(
                    "num_kpoints",
                    old.get("num_kpoints"),
                    new.get("num_kpoints"),
                    "modified",
                    "critical",
                )
            )

        if old.get("cartesian") != new.get("cartesian"):
            entries.append(
                DiffEntry(
                    "cartesian",
                    old.get("cartesian"),
                    new.get("cartesian"),
                    "modified",
                    "major",
                )
            )

        return entries

    def _diff_line(
        self, old: Dict[str, Any], new: Dict[str, Any]
    ) -> List[DiffEntry]:
        entries: List[DiffEntry] = []

        if old.get("num_points") != new.get("num_points"):
            entries.append(
                DiffEntry(
                    "num_points",
                    old.get("num_points"),
                    new.get("num_points"),
                    "modified",
                    "minor",
                )
            )

        old_segs = len(old.get("segments", []))
        new_segs = len(new.get("segments", []))
        if old_segs != new_segs:
            entries.append(
                DiffEntry(
                    "num_segments",
                    old_segs,
                    new_segs,
                    "modified",
                    "major",
                )
            )

        return entries
