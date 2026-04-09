# Open-Source and Modular Implementation of the Dislocation Extraction Algorithm (DXA)

`OpenDXA` performs dislocation analysis from a crystal-state package.
It is decoupled from the upstream structure-identification algorithm, but it is not
agnostic to lattice geometry: the matrix phase must exist as a reference lattice
definition under `OpenDXA/lattices` or in a runtime lattice directory.

## Overview

OpenDXA consumes three files from the same snapshot:

1. An annotated LAMMPS dump
2. A `*_clusters.table`
3. A `*_cluster_transitions.table`

Those files are currently exported by:

- [CommonNeighborAnalysis](https://github.com/VoltLabs-Research/CommonNeighborAnalysis)
- [PolyhedralTemplateMatching](https://github.com/VoltLabs-Research/PolyhedralTemplateMatching)
- [PatternStructureMatching](https://github.com/VoltLabs-Research/PatternStructureMatching)

If the producer respects the contract documented below, OpenDXA can consume it.

## OpenDXA Reference Lattices

OpenDXA resolves the matrix phase by `topology_name`.
Reference lattices are loaded from YAML files under `OpenDXA/lattices`.

Each OpenDXA lattice YAML currently contains only:

- `name`
- `coordination_number`
- `neighbor_vectors`

### Lattice Search Path

OpenDXA resolves reference lattices in this order:

1. `--lattice-dir`, when provided
2. The compile-time source lattice directory
3. `share/volt/lattices` relative to the executable/package

This means you can test new OpenDXA lattice YAMLs at runtime without recompiling by
passing `--lattice-dir /path/to/lattices`.

## Required Contract

OpenDXA expects the annotated dump, clusters table, and transitions table to come
from the same frame and to reference the same cluster IDs and topology names.

### Annotated Dump

The annotated dump must be a valid single-frame LAMMPS dump with at least these
headers:

- `ITEM: TIMESTEP`
- `ITEM: NUMBER OF ATOMS`
- `ITEM: MAXIMUM NEIGHBOR DISTANCE`
- `ITEM: BOX BOUNDS`
- `ITEM: ATOMS ...`

`ITEM: MAXIMUM NEIGHBOR DISTANCE` is required by OpenDXA. It is used to derive
core geometric tolerances during reconstruction:

- Delaunay ghost layer size: `3.5 * maximumNeighborDistance`
- Interface mesh alpha: `5.0 * maximumNeighborDistance`

The `ITEM: ATOMS` section must include:

- atom coordinates
- `cluster_id`
- `neighbor_counts`
- `neighbor_indices_0` ... `neighbor_indices_17`
- `neighbor_lattice_x_0`, `neighbor_lattice_y_0`, `neighbor_lattice_z_0`, ... `neighbor_lattice_x_17`, `neighbor_lattice_y_17`, `neighbor_lattice_z_17`

Meaning of the reconstructed columns:

- `cluster_id`: cluster membership; `0` means no cluster / defect region
- `neighbor_counts`: number of valid reconstructed neighbors for that atom
- `neighbor_indices_k`: zero-based atom index of slot `k`, in dump order
- `neighbor_lattice_*_k`: ideal lattice-space vector assigned to slot `k`

The current reconstructed format reserves up to `18` neighbor slots per atom.

### `*_clusters.table`

The clusters table must contain:

- `cluster_id`
- `topology_name`

It may also contain the optional orientation matrix columns:

- `orientation_00` ... `orientation_22`

Rules:

- `cluster_id` must match the `cluster_id` written into the annotated dump
- `topology_name` must match an OpenDXA lattice `name`
- the matrix phase passed as `--reference-topology` must match the dominant matrix `topology_name`

Examples:

- FCC matrix: `--reference-topology fcc`
- BCC matrix: `--reference-topology bcc`
- FCT matrix: `--reference-topology fct`

### `*_cluster_transitions.table`

The transitions table must contain:

- `cluster1_id`
- `cluster2_id`
- `tm_00` ... `tm_22`
- `distance`

Where:

- `cluster1_id` and `cluster2_id` define a directed cluster-graph edge
- `tm_00` ... `tm_22` define the `3x3` transition matrix between clusters
- `distance` is the graph distance assigned to that transition

## OpenDXA CLI

Usage:

```bash
opendxa <annotated.dump> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<annotated.dump>` | Yes | Annotated dump exported by an upstream producer. | |
| `[output_base]` | No | Output basename. If omitted, OpenDXA derives it from the input dump path. | derived from input |
| `--clusters-table <path>` | Yes | Path to `*_clusters.table`. | |
| `--clusters-transitions <path>` | Yes | Path to `*_cluster_transitions.table`. | |
| `--cluster-transitions <path>` | No | Accepted alias for `--clusters-transitions`. | |
| `--reference-topology <name>` | Yes | Matrix-phase topology name resolved from OpenDXA lattice YAMLs. | |
| `--lattice-dir <path>` | No | Directory containing OpenDXA lattice YAMLs. | compiled/package lattice directory |
| `--maxTrialCircuitSize <int>` | No | Maximum Burgers circuit size. | `14` |
| `--circuitStretchability <int>` | No | Circuit stretchability factor. | `9` |
| `--lineSmoothingLevel <float>` | No | Smoothing applied to dislocation lines. | `1.0` |
| `--linePointInterval <float>` | No | Point spacing along exported lines. | `2.5` |
| `--help` | No | Print CLI help. | |

## Upstream Producers

### CommonNeighborAnalysis

Repository:
[VoltLabs-Research/CommonNeighborAnalysis](https://github.com/VoltLabs-Research/CommonNeighborAnalysis)

Supported input structures:

- `FCC`
- `BCC`
- `HCP`
- `CUBIC_DIAMOND`
- `HEX_DIAMOND`

### PolyhedralTemplateMatching

Repository:
[VoltLabs-Research/PolyhedralTemplateMatching](https://github.com/VoltLabs-Research/PolyhedralTemplateMatching)

Supported input structures:

- `SC`
- `FCC`
- `HCP`
- `BCC`
- `CUBIC_DIAMOND`
- `HEX_DIAMOND`

### PatternStructureMatching

Repository:
[VoltLabs-Research/PatternStructureMatching](https://github.com/VoltLabs-Research/PatternStructureMatching)

PatternStructureMatching uses dynamic lattice YAMLs instead of a fixed compiled list.
It still exports the same reconstructed-state contract that OpenDXA expects.

Supported CLI parameters:

- `--lattice-dir <path>`
- `--reference-lattice-dir <path>`
- `--patterns <csv>`
- `--dissolveSmallClusters`

#### Currently Supported PatternStructureMatching Lattices

| Lattice | Coordination Number |
| --- | ---: |
| `9R` | 12 |
| `L12` | 12 |
| `a7` | 6 |
| `bcc` | 14 |
| `bct` | 12 |
| `cubic_diamond` | 16 |
| `fcc` | 12 |
| `fct` | 16 |
| `hcp` | 12 |
| `hex_diamond` | 16 |
| `omega` | 14 |
| `sc` | 6 |
| `st` | 10 |

#### Adding a New PatternStructureMatching Lattice

To add a new lattice to PatternStructureMatching:

1. Create `plugins/PatternStructureMatching/lattices/<name>.yml`
2. Define:
   - `name`
   - `coordination_number`
   - `cell`
   - `coordinate_mode`
   - `basis`
3. Run PatternStructureMatching with:
   - `--lattice-dir /path/to/PatternStructureMatching/lattices`
   - `--patterns <name>`

Minimal example:

```yml
name: my_lattice
coordination_number: 8

cell:
  - [1.0, 0.0, 0.0]
  - [0.0, 1.0, 0.0]
  - [0.0, 0.0, 1.0]

coordinate_mode: fractional

basis:
  - species: 1
    position: [0.0, 0.0, 0.0]
```

#### Making a New PatternStructureMatching Lattice Usable by OpenDXA

If the new PatternStructureMatching lattice should also be consumable by OpenDXA:

1. Create `plugins/OpenDXA/lattices/<name>.yml`
2. Define:
   - `name`
   - `coordination_number`
   - `neighbor_vectors`
3. Run OpenDXA with `--lattice-dir /path/to/OpenDXA/lattices`
4. Ensure the upstream `topology_name` matches the same `<name>`

Minimal OpenDXA example:

```yml
name: my_lattice
coordination_number: 8

neighbor_vectors:
  - [1.0, 0.0, 0.0]
  - [-1.0, 0.0, 0.0]
  - [0.0, 1.0, 0.0]
  - [0.0, -1.0, 0.0]
  - [0.0, 0.0, 1.0]
  - [0.0, 0.0, -1.0]
  - [1.0, 1.0, 0.0]
  - [-1.0, -1.0, 0.0]
```

No OpenDXA code changes are required as long as:

- the producer emits a compatible reconstructed-state package
- `topology_name` matches the YAML `name`
- the matrix phase passed to `--reference-topology` exists in the OpenDXA lattice directory
