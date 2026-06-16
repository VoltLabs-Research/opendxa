# Dislocation Analysis

Dislocation analysis (DXA) driven by the structure-identification algorithm selected by the user.

## Install

```bash
vpm install @voltlabs/opendxa
```

## CLI

```bash
opendxa <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump (annotated by an upstream structure-identification step). |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--clusters_table <path>` | yes | — | Cluster graph table from an upstream PTM/ACNA/PSM step. |
| `--clusters_transitions <path>` | yes | — | Cluster transitions table from the same upstream step. |
| `--neighbor_lattice <path>` | yes | — | Per-atom neighbor topology parquet from the same upstream step. |
| `--reference_topology <name>` | no | `fcc` | Matrix-phase lattice: `fcc`, `bcc`, `hcp`, `sc`, `cubic_diamond`, `hex_diamond`. |
| `--max_trial_circuit_size <int>` | no | `14` | Maximum Burgers circuit size. |
| `--circuit_stretchability <int>` | no | `9` | Circuit stretchability factor. |
| `--line_smoothing_level <float>` | no | `1.0` | Smoothing applied to dislocation lines. |
| `--line_point_interval <float>` | no | `2.5` | Point spacing along exported lines. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_dislocations.parquet` | Dislocations | LineExporter → glb |
| `{output_base}_defect_mesh.parquet` | Defect Mesh | MeshExporter → glb |
| `{output_base}_interface_mesh.parquet` | Interface Mesh | MeshExporter → glb |
| `{output_base}_atoms.parquet` | Structure Identification | AtomisticExporter → glb |
| `{output_base}_coherent_crystalline_regions.parquet` | Crystalline Coherent Regions | AtomisticExporter → glb |
| `{output_base}_dislocation_summary.parquet` | Burgers Segment Counts | ChartExporter → chart-png |
| `{output_base}_dislocation_summary.parquet` | Burgers Length Distribution | ChartExporter → chart-png |
| `{output_base}_dislocation_summary.parquet` | Network Statistics | — (listing-only) |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/open-dxa
