# Open-Source and Modular Implementation of the Dislocation Extraction Algorithm (DXA)

`OpenDXA` performs dislocation analysis from a crystal-state package.
It is decoupled from the upstream structure-identification algorithm, but it is not
agnostic to lattice geometry: the matrix phase must exist as a reference lattice
definition under `OpenDXA/lattices` or in a runtime lattice directory.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- OpenDXA
```

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
| `--reference-topology <name>` | Yes | Matrix-phase topology name resolved from OpenDXA lattice YAMLs. | |
| `--lattice-dir <path>` | No | Directory containing OpenDXA lattice YAMLs. | compiled/package lattice directory |
| `--max-trial-circuit-size <int>` | No | Maximum Burgers circuit size. | `14` |
| `--circuit-stretchability <int>` | No | Circuit stretchability factor. | `9` |
| `--line-smoothing-level <float>` | No | Smoothing applied to dislocation lines. | `1.0` |
| `--line-point-interval <float>` | No | Point spacing along exported lines. | `2.5` |
| `--ghost-layer-scale <float>` | No | Multiplier applied to the reconstructed maximum neighbor distance before building ghost atoms for the Delaunay tessellation. | `3.5` |
| `--interface-alpha-scale <float>` | No | Multiplier applied to the reconstructed maximum neighbor distance when running the interface alpha-shape filter. | `5.0` |
| `--inteface-alpha-scale <float>` | No | Accepted alias for `--interface-alpha-scale`. | |
| `--crystal-path-steps <int>` | No | Maximum crystal path depth used while assigning ideal edge vectors. | `4` |
| `--export-defect-mesh <bool>` | No | Enable or disable writing `*_defect_mesh.msgpack`. | `true` |
| `--export-interface-mesh <bool>` | No | Enable or disable writing `*_interface_mesh.msgpack`. | `false` |
| `--export-dislocations <bool>` | No | Enable or disable writing `*_dislocations.msgpack`. | `true` |
| `--export-circuit-information <bool>` | No | Include `circuit_information` inside the dislocations msgpack. | `true` |
| `--export-dislocation-network-stats <bool>` | No | Include `network_statistics` inside the dislocations msgpack. | `true` |
| `--export-junctions <bool>` | No | Include `junction_information` inside the dislocations msgpack. | `true` |
| `--clip-pbc-segments <bool>` | No | Clip exported dislocation polylines at periodic boundaries. If disabled, OpenDXA exports the raw traced lines. | `true` |
| `--cover-domain-with-finite-tets <bool>` | No | Add helper points so the Delaunay domain is fully covered by finite tetrahedra. | `false` |
| `--help` | No | Print CLI help. | |

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
- `ITEM: BOX BOUNDS`
- `ITEM: ATOMS ...`

OpenDXA now derives the maximum reconstructed neighbor distance at import time from
the annotated atom positions and neighbor graph. It no longer requires a dedicated
dump header for that value.

The `ITEM: ATOMS` section must include:

- atom coordinates
- `cluster_id`
- `neighbor_indices_0` ... `neighbor_indices_17`
- `neighbor_lattice_x_0`, `neighbor_lattice_y_0`, `neighbor_lattice_z_0`, ... `neighbor_lattice_x_17`, `neighbor_lattice_y_17`, `neighbor_lattice_z_17`

Meaning of the reconstructed columns:

- `cluster_id`: cluster membership; `0` means no cluster / defect region
- `neighbor_indices_k`: zero-based atom index of slot `k`, in dump order; unused slots must be `-1`
- `neighbor_lattice_*_k`: ideal lattice-space vector assigned to slot `k`

The current reconstructed format reserves up to `18` neighbor slots per atom.
OpenDXA derives the effective neighbor count for each atom from the first `-1`
sentinel in the `neighbor_indices_*` slots.

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

Where:

- `cluster1_id` and `cluster2_id` define an undirected direct cluster interface
- `tm_00` ... `tm_22` define the `3x3` transition matrix between clusters

Rules:

- only direct interfaces are exported
- reverse duplicate rows are not required
- self-transitions are not serialized
- OpenDXA reconstructs reverse edges and self-transitions internally at import time

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
