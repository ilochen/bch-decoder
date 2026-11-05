# Timing annotation archive

These SDFs were exported alongside the submitted design netlists. Each file annotates a specific design hierarchy and includes instance paths, interconnect delays, cell arcs, and timing checks. The synthesis and post-route results must be paired with their corresponding netlists. The submitted exports were copied byte for byte.

| Stage | File | CELL records | SHA-256 |
|---|---|---:|---|
| synthesis | [bch_syn.sdf](../implementation/synthesis/bch_syn.sdf) | 52,561 | `691a98f71dd7420bf855f40df6e5b394a825e23a529b43c1d11b45bc24e455c2` |
| postroute | [bch_apr.sdf](../implementation/postroute/bch_apr.sdf) | 51,982 | `5aaa40360d46d0cec035b3bca917ff14c18a1d13a1411ca0d96fb31a60bfb83f` |

The RTL regressions exercise the source design with functional timing. To perform a separate gate-level timing run, load the matching mapped netlist and the authorized simulation cell models in a licensed simulator, annotate the file on that design instance, and inspect the simulator’s annotation coverage and timing violations. Gate-level SDF simulation has not been run as part of this repository validation.
