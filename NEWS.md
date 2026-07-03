# treestructure 1.0.1

* The default target false discovery rate is now `fdr = 0.2` (was `0.1`). The stricter
  0.1 default was conservative on real trees and interacted awkwardly with
  `minCladeSize`; at `0.2` the number of designated clusters decreases monotonically as
  `minCladeSize` grows, and more genuine structure is recovered while the whole-tree
  error rate remains controlled.
* The node-support vignette was revised to the new default (`fdr = 0.2`).

# treestructure 1.0.0

## Calibrating to a false discovery rate

* `trestruct()` now calibrates the split threshold to a target false discovery rate
  by **default** (`fdr = 0.1`). At each scan the algorithm splits at the most extreme
  eligible candidate clade only if it clears a multiple-testing correction over the
  eligible candidates in that scan. The `fdr` is a property of the **whole tree**:
  under the global null of a single unstructured coalescent, the probability of
  designating any structure equals `fdr`; when real structure is present, `fdr`
  bounds the expected fraction of spurious splits among all splits. The analytic
  calibration requires no simulation and is deterministic. Supplying an explicit
  `level` (without `fdr`) selects the previous subjective-threshold behaviour, and
  `level = NULL` selects the CH index. Arguments were reordered
  (`fdr, level, minCladeSize, ...`) and `minCladeSize` now defaults to 10.

* The multiple-testing correction is selectable with `split`: `'bonferroni'` (the
  default) or `'bh'`, a Benjamini-Hochberg step-up that is less conservative and
  retains more power on large trees with abundant moderate structure, while still
  controlling the false discovery rate.

* The returned `TreeStructure` object now carries a **global-null test** in
  `$global.test` (the root-scan max|z|, the number of candidates, and a Bonferroni
  p-value for the presence of any structure at all), reported by `print()`.

* In `fdr` mode a **heterochronous-sampling diagnostic** is computed (`$hetero`) and
  printed. On serially sampled trees the coalescent statistic carries a small
  positive bias that can modestly inflate the realized FDR; this is reported as a
  caveat and discussed in the vignette. `fdr` composes with the existing
  `nodeSupportValues`/`nodeSupportThreshold` node-support filtering.

## Fast, simulation-free statistic

* The coalescent rank-sum null is characterized by exact analytic moments by
  default (`method = 'analytic'`) — a fast, deterministic alternative to Monte-Carlo
  simulation (`method = 'sim'`, the previous behaviour). `trestruct()` and
  `treestructure.test()` both gain the `method` argument.

* Fixed an error in the exact (Ei = 3) rank-sum transition used for
  monophyletic/monophyletic contrasts.

## Documentation and tests

* Added a `testthat` test suite covering the FDR calibration, the global-null test,
  the heterochronous diagnostic, and composition with node support.

* Vignettes updated to lead with the FDR calibration and document the precise
  meaning of the FDR, to show the isochronous FWER validation, and to present node
  support and the CH index as complementary ways of refining clusters.
