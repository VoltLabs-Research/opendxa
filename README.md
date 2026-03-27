# OpenDXA

OpenDXA is fully decoupled from the upstream structure-identification method. It expects a reconstructed-structure package generated from the same snapshot:

1. An annotated LAMMPS dump file
2. A `*_clusters.table` file
3. A `*_cluster_transitions.table` file

As long as those three files respect the contract below, OpenDXA does not need to know whether the upstream producer was CNA, PTM, or another compatible implementation.

## Required CLI Inputs

OpenDXA requires:

- `--clusters-table <path>`
- `--clusters-transitions <path>`
- `--reference-structure-label <int>`

`--reference-structure-label` must be the numeric `structure_type` label used by the matrix phase in `*_clusters.table`.

## Annotated Dump Requirements

The annotated dump must be a valid LAMMPS dump frame with atomic positions and these headers:

- `ITEM: TIMESTEP`
- `ITEM: NUMBER OF ATOMS`
- `ITEM: MAXIMUM NEIGHBOR DISTANCE`
- `ITEM: BOX BOUNDS`
- `ITEM: ATOMS ...`

`ITEM: MAXIMUM NEIGHBOR DISTANCE` is required by OpenDXA and must already be written by the upstream structural-analysis step.

This value is used as the characteristic neighborhood length scale of the reconstructed structure:

- Delaunay tessellation uses `ghostLayerSize = 3.5 * maximumNeighborDistance`
- Interface-mesh construction uses `alpha = 5.0 * maximumNeighborDistance`

Upstream producers currently generate it as follows:

- CNA: maximum local neighbor distance observed while classifying atoms and building ordered neighbor lists
- PTM: maximum minimum-image neighbor distance found by scanning the reconstructed PTM neighbor list

### Required Per-Atom Columns

The `ITEM: ATOMS` section must include the atom positions together with these reconstructed columns:

- `cluster_id`
- `neighbor_counts`
- `neighbor_indices_0` ... `neighbor_indices_15`
- `neighbor_lattice_x_0`, `neighbor_lattice_y_0`, `neighbor_lattice_z_0`, ... `neighbor_lattice_x_15`, `neighbor_lattice_y_15`, `neighbor_lattice_z_15`

The slot range ends at `15` because the reconstructed format reserves up to `16` neighbor slots per atom. `neighbor_counts` tells OpenDXA how many of those slots are valid for each atom.

These columns mean:

- `cluster_id`: cluster membership of the atom; `0` is the special no-cluster/defect label
- `neighbor_counts`: number of valid reconstructed neighbors for that atom
- `neighbor_indices_k`: zero-based atom index of neighbor slot `k` in dump order, not the LAMMPS `id`
- `neighbor_lattice_x_k`, `neighbor_lattice_y_k`, `neighbor_lattice_z_k`: components of the ideal neighbor lattice vector associated with slot `k`

OpenDXA consumes them in this order:

1. It rebuilds the reconstructed structure context from `cluster_id`, `neighbor_counts`, `neighbor_indices_*`, and `neighbor_lattice_*_*`
2. It links `cluster_id` to the imported cluster graph from `*_clusters.table` and `*_cluster_transitions.table`
3. It rebuilds the compact neighbor list used by neighbor queries and crystal-path traversal
4. It uses `neighbor_lattice_*_*` as ideal lattice-space bond vectors during crystal path finding and elastic mapping
5. It uses those bond vectors plus cluster transitions to construct the interface mesh and trace Burgers circuits

## `*_clusters.table` Requirements

The clusters table must contain:

- `cluster_id`
- `structure_type`

Where:

- `cluster_id`: unique cluster identifier referenced by the dump's per-atom `cluster_id`
- `structure_type`: structure label assigned to that cluster

`--reference-structure-label` must match the `structure_type` value used by the matrix phase. Common VOLT labels are:

- `1` = `SC`
- `2` = `FCC`
- `3` = `HCP`
- `4` = `BCC`
- `5` = `CUBIC_DIAMOND`
- `6` = `HEX_DIAMOND`

Examples:

- FCC matrix -> `--reference-structure-label 2`
- BCC matrix -> `--reference-structure-label 4`
- If the matrix clusters are written as `structure_type = 3`, use `--reference-structure-label 3`

The safe rule is simple: inspect the matrix-phase clusters in `*_clusters.table` and pass that same label to OpenDXA.

## `*_cluster_transitions.table` Requirements

The transitions table must contain:

- `cluster1_id`
- `cluster2_id`
- `tm_00`, `tm_01`, `tm_02`
- `tm_10`, `tm_11`, `tm_12`
- `tm_20`, `tm_21`, `tm_22`
- `distance`

Where:

- `cluster1_id` and `cluster2_id`: directed edge in the reconstructed cluster graph
- `tm_00` ... `tm_22`: entries of the `3 x 3` transition matrix between the two clusters
- `distance`: graph distance associated with that transition

## Installation

```bash
curl -sSL https://raw.githubusercontent.com/voltlabs-research/CoreToolkit/main/scripts/install.sh | bash
```
