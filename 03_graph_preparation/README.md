# Boundary-aware graph preparation

This step prepares graph inputs for the WDL after introgression intervals have been lifted to graph/sample-path coordinates.

## Why inject and then remove?

`odgi inject` is used to add introgression-path intervals to the graph. This is useful because injection splits graph nodes at the introgression boundaries, giving the graph a node partition that is aligned to tract boundaries.

After this boundary-aware node partition has been created, the injected introgression paths are removed before building the genotyping indexes. The resulting graph keeps the refined node granularity, but the injected archaic/introgression paths are not retained as haplotype or reference-like paths during downstream genotyping.

In practice:

- Use the injected graph as the WDL `injected_og_file`, so the workflow can recover introgressed paths and trace supporting topology.
- Use the cleaned graph/index outputs as the WDL graph/index inputs (`original_gbz_file`, `original_hapl_file`, `original_ri_file`, and `original_snarls_file`).

## Example

```bash
03_graph_preparation/inject_introgression_paths_odgi.sh \
  graph.og \
  lifted_introgression_intervals.bed \
  graph.with_introgression_paths.og \
  32

03_graph_preparation/remove_introgression_paths_and_build_indexes.sh \
  graph.with_introgression_paths.og \
  graph.boundary_aware_no_introgression_paths \
  graph_indexes \
  32
```

The second command writes a cleaned `.og` graph and the derived `.gbz`, `.ri`, `.hapl`, and `.snarls` files used by the WDL.
