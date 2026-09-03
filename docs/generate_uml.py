#!/usr/bin/env python3
"""Regenerate CoreMS data-model UML class diagrams as SVG.

Requires:
  - pylint (provides ``pyreverse``), e.g. ``pip install -e ".[dev]"``
  - Graphviz ``dot`` on PATH

Usage (from repo root)::

    python3 docs/generate_uml.py
    # or
    make uml
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = REPO_ROOT / "docs" / "uml"


@dataclass(frozen=True)
class DiagramSpec:
    """One modality class diagram."""

    stem: str
    modules: Sequence[str]
    class_names: frozenset[str]


# Curated factory / parameter modules and entity classes (not full calc trees).
DIAGRAMS: tuple[DiagramSpec, ...] = (
    DiagramSpec(
        stem="Direct_Infusion_FTMS_Data_Model",
        modules=(
            "corems.mass_spectrum.factory.MassSpectrumClasses",
            "corems.ms_peak.factory.MSPeakClasses",
            "corems.molecular_formula.factory.MolecularFormulaFactory",
            "corems.encapsulation.factory.parameters",
        ),
        class_names=frozenset(
            {
                "MassSpecBase",
                "MassSpecProfile",
                "MassSpecfromFreq",
                "MassSpecCentroid",
                "MassSpecCentroidLowRes",
                "MSParameters",
                "_MSPeak",
                "ICRMassPeak",
                "TOFMassPeak",
                "OrbiMassPeak",
                "MolecularFormulaBase",
                "MolecularFormula",
                "MolecularFormulaIsotopologue",
            }
        ),
    ),
    DiagramSpec(
        stem="GC_MS_Data_Model",
        modules=(
            "corems.mass_spectra.factory.GC_Class",
            "corems.chroma_peak.factory.chroma_peak_classes",
            "corems.encapsulation.factory.parameters",
            "corems.mass_spectrum.factory.MassSpectrumClasses",
        ),
        class_names=frozenset(
            {
                "GCMSBase",
                "GCMSParameters",
                "ChromaPeakBase",
                "GCPeak",
                "GCPeakDeconvolved",
                "MSParameters",
                "MassSpecBase",
                "MassSpecCentroid",
                "MassSpecCentroidLowRes",
            }
        ),
    ),
    DiagramSpec(
        stem="LC_MS_Data_Model",
        modules=(
            "corems.mass_spectra.factory.lc_class",
            "corems.mass_spectra.factory.chromat_data",
            "corems.chroma_peak.factory.chroma_peak_classes",
            "corems.encapsulation.factory.parameters",
        ),
        class_names=frozenset(
            {
                "MassSpectraBase",
                "LCMSBase",
                "LCMSCollection",
                "ChromaPeakBase",
                "LCMSMassFeature",
                "LCMSParameters",
                "LCMSCollectionParameters",
                "TIC_Data",
                "EIC_Data",
                "MSParameters",
            }
        ),
    ),
)


def _simple_class_name(fqcn: str) -> str:
    return fqcn.rsplit(".", 1)[-1]


def _require_tools() -> None:
    if shutil.which("pyreverse") is None:
        sys.exit(
            "error: pyreverse not found on PATH. Install pylint "
            '(e.g. pip install -e ".[dev]").'
        )
    if shutil.which("dot") is None:
        sys.exit(
            "error: graphviz 'dot' not found on PATH. Install Graphviz "
            "(e.g. brew install graphviz)."
        )


def _run_pyreverse(
    project: str, modules: Sequence[str], work_dir: Path, env_pythonpath: str
) -> Path:
    cmd = [
        "pyreverse",
        "-o",
        "dot",
        "-p",
        project,
        "-d",
        str(work_dir),
        *modules,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = env_pythonpath
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout or "")
        sys.stderr.write(proc.stderr or "")
        sys.exit(f"error: pyreverse failed for project {project!r} (exit {proc.returncode})")
    classes_dot = work_dir / f"classes_{project}.dot"
    if not classes_dot.is_file():
        sys.exit(f"error: expected pyreverse output missing: {classes_dot}")
    return classes_dot


_NODE_RE = re.compile(
    r'^"(?P<id>[^"]+)"\s*\[(?P<attrs>.*)\]\s*;\s*$',
    re.MULTILINE,
)
_EDGE_RE = re.compile(
    r'^"(?P<src>[^"]+)"\s*->\s*"(?P<dst>[^"]+)"\s*\[(?P<attrs>.*)\]\s*;\s*$',
    re.MULTILINE,
)


def filter_classes_dot(dot_text: str, allowed: frozenset[str], digraph_name: str) -> str:
    """Keep only nodes whose simple class name is in ``allowed``, and edges between them."""
    keep_ids: set[str] = set()
    node_lines: list[str] = []
    for m in _NODE_RE.finditer(dot_text):
        node_id = m.group("id")
        if _simple_class_name(node_id) in allowed:
            keep_ids.add(node_id)
            node_lines.append(m.group(0))

    edge_lines: list[str] = []
    for m in _EDGE_RE.finditer(dot_text):
        src, dst = m.group("src"), m.group("dst")
        if src in keep_ids and dst in keep_ids:
            edge_lines.append(m.group(0))

    missing = sorted(allowed - {_simple_class_name(i) for i in keep_ids})
    if missing:
        sys.stderr.write(
            f"warning: classes not found in pyreverse graph for {digraph_name}: "
            f"{', '.join(missing)}\n"
        )
    if not node_lines:
        sys.exit(f"error: no allowed classes found in graph for {digraph_name}")

    body = "\n".join(node_lines + edge_lines)
    return f'digraph "{digraph_name}" {{\nrankdir=BT\ncharset="utf-8"\n{body}\n}}\n'


def _render_svg(dot_path: Path, svg_path: Path) -> None:
    proc = subprocess.run(
        ["dot", "-Tsvg", "-o", str(svg_path), str(dot_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout or "")
        sys.stderr.write(proc.stderr or "")
        sys.exit(f"error: dot failed for {dot_path} (exit {proc.returncode})")
    if not svg_path.is_file() or svg_path.stat().st_size == 0:
        sys.exit(f"error: empty or missing SVG: {svg_path}")


def generate_diagram(spec: DiagramSpec, out_dir: Path) -> Path:
    project = re.sub(r"[^A-Za-z0-9_]", "_", spec.stem)
    pythonpath = str(REPO_ROOT)
    with tempfile.TemporaryDirectory(prefix="corems_uml_") as tmp:
        work = Path(tmp)
        raw_dot = _run_pyreverse(project, spec.modules, work, pythonpath)
        filtered = filter_classes_dot(
            raw_dot.read_text(encoding="utf-8"),
            spec.class_names,
            digraph_name=spec.stem,
        )
        filtered_path = work / f"{spec.stem}.dot"
        filtered_path.write_text(filtered, encoding="utf-8")
        out_dir.mkdir(parents=True, exist_ok=True)
        svg_path = out_dir / f"{spec.stem}.svg"
        _render_svg(filtered_path, svg_path)
    return svg_path


def generate_all(out_dir: Path, diagrams: Iterable[DiagramSpec] = DIAGRAMS) -> list[Path]:
    _require_tools()
    written: list[Path] = []
    for spec in diagrams:
        path = generate_diagram(spec, out_dir)
        written.append(path)
        print(f"wrote {path.relative_to(REPO_ROOT)}")
    return written


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Directory for SVG outputs (default: {DEFAULT_OUT_DIR})",
    )
    args = parser.parse_args(argv)
    out_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    generate_all(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
