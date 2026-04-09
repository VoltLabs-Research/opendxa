#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


STRUCTURES = {
    "bct": {
        "title": "Body-Centered Tetragonal Lattice Reference",
        "coordination_number": 12,
        "a": 1.0,
        "c": 1.2,
    },
    "fct": {
        "title": "Face-Centered Tetragonal Lattice Reference",
        "coordination_number": 16,
        "a": 1.0,
        "c": 1.2,
    },
    "st": {
        "title": "Simple Tetragonal Lattice Reference",
        "coordination_number": 10,
        "a": 1.0,
        "c": 1.2,
    },
}


def format_number(value: float) -> str:
    if abs(value) < 1e-12:
        value = 0.0
    text = f"{value:.12f}".rstrip("0").rstrip(".")
    if "." not in text:
        text += ".0"
    return text


def run_atomsk(atomsk: Path, *args: str) -> None:
    subprocess.run([str(atomsk), *args], check=True)


def parse_cfg_positions(path: Path) -> tuple[list[tuple[float, float, float]], list[list[float]]]:
    lines = path.read_text().splitlines()
    h = [[0.0] * 3 for _ in range(3)]
    start = None

    for index, line in enumerate(lines):
        match = re.match(r"H0\((\d),(\d)\) =\s+([-+0-9.eE]+)", line)
        if match:
            row, column, value = match.groups()
            h[int(row) - 1][int(column) - 1] = float(value)
        if line.strip() == "entry_count = 3":
            start = index + 3
            break

    if start is None:
        raise RuntimeError(f"Unable to find atomic data in {path}")

    positions: list[tuple[float, float, float]] = []
    for line in lines[start:]:
        columns = line.split()
        if len(columns) < 3:
            continue
        rx, ry, rz = map(float, columns[:3])
        positions.append(
            (
                rx * h[0][0] + ry * h[1][0] + rz * h[2][0],
                rx * h[0][1] + ry * h[1][1] + rz * h[2][1],
                rx * h[0][2] + ry * h[1][2] + rz * h[2][2],
            )
        )

    return positions, h


def derive_neighbor_vectors(
    positions: list[tuple[float, float, float]],
    h: list[list[float]],
    coordination_number: int,
) -> list[tuple[float, float, float]]:
    center_guess = (h[0][0] / 2.0, h[1][1] / 2.0, h[2][2] / 2.0)
    center = min(
        positions,
        key=lambda pos: (
            (pos[0] - center_guess[0]) ** 2
            + (pos[1] - center_guess[1]) ** 2
            + (pos[2] - center_guess[2]) ** 2
        ),
    )

    candidates: list[tuple[float, tuple[float, float, float]]] = []
    for position in positions:
        if position == center:
            continue
        vector = (
            position[0] - center[0],
            position[1] - center[1],
            position[2] - center[2],
        )
        distance_sq = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
        candidates.append((distance_sq, vector))

    candidates.sort(
        key=lambda item: (
            round(item[0], 12),
            item[1][0],
            item[1][1],
            item[1][2],
        )
    )

    vectors: list[tuple[float, float, float]] = []
    for _, vector in candidates[:coordination_number]:
        vectors.append(
            tuple(0.0 if abs(component) < 1e-12 else round(component, 12) for component in vector)
        )

    return vectors


def build_yaml(name: str, title: str, coordination_number: int, vectors: list[tuple[float, float, float]]) -> str:
    lines = [
        "# -----------------------------------------------------------------------------",
        f"# {title}",
        "# -----------------------------------------------------------------------------",
        "#",
        "# For additional examples, detailed explanations, and the full specification,",
        "# please refer to the official documentation:",
        "#",
        "#   https://docs.voltcloud.dev",
        "#   https://github.com/voltlabs-research/",
        "#",
        "# -----------------------------------------------------------------------------",
        "",
        "# Unique identifier for the lattice.",
        f"name: {name}",
        "",
        "# Number of nearest neighbors per atom.",
        f"coordination_number: {coordination_number}",
        "",
        "# Relative vectors from the central atom to each of its neighbors.",
        "# These define the ideal local geometry.",
        "neighbor_vectors:",
    ]

    for x, y, z in vectors:
        lines.append(
            f"  - [{format_number(x)}, {format_number(y)}, {format_number(z)}]"
        )

    lines.append("")
    return "\n".join(lines)


def write_tetragonal_lattice(atomsk: Path, output_dir: Path, name: str, data: dict[str, object]) -> None:
    with tempfile.TemporaryDirectory(prefix=f"atomsk-{name}-") as temp_dir:
        temp = Path(temp_dir)
        unit_cfg = temp / f"{name}_unit.cfg"
        super_cfg = temp / f"{name}_555.cfg"

        run_atomsk(
            atomsk,
            "--create",
            name,
            format_number(float(data["a"])),
            format_number(float(data["c"])),
            "Fe",
            str(unit_cfg),
        )
        run_atomsk(
            atomsk,
            str(unit_cfg),
            "-duplicate",
            "5",
            "5",
            "5",
            str(super_cfg),
        )

        positions, h = parse_cfg_positions(super_cfg)
        vectors = derive_neighbor_vectors(positions, h, int(data["coordination_number"]))
        content = build_yaml(
            name,
            str(data["title"]),
            int(data["coordination_number"]),
            vectors,
        )
        (output_dir / f"{name}.yml").write_text(content)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate OpenDXA tetragonal lattice YAMLs from Atomsk."
    )
    parser.add_argument("atomsk", type=Path, help="Path to the atomsk executable")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "lattices",
        help="Directory where bct.yml, fct.yml and st.yml will be written",
    )
    args = parser.parse_args()

    atomsk = args.atomsk.resolve()
    if not atomsk.is_file():
        raise SystemExit(f"Atomsk executable not found: {atomsk}")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, data in STRUCTURES.items():
        write_tetragonal_lattice(atomsk, output_dir, name, data)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
