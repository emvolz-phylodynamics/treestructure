## Submission

treestructure 1.0.1 is a maintenance update.

## Changes in this version

* The default target false discovery rate is now `fdr = 0.2` (previously `0.1`). The
  stricter default was conservative on empirical trees and non-monotone in
  `minCladeSize`; at 0.2 the number of designated clusters decreases monotonically as
  `minCladeSize` grows, while the whole-tree error rate remains controlled.
* The vignettes were revised for the new default. The node-support vignette no longer
  accesses the network when it is built: the tree-construction code is shown for
  reference and the results are loaded from precomputed objects, so all vignettes build
  offline.
* Vignettes now use the `rmarkdown::html_vignette` output format; the unused `bookdown`
  suggestion was removed.

## Test environments

* win-builder: R-devel (x86_64-w64-mingw32, Windows Server 2022).
* Local: Ubuntu Linux (kernel 6.8), R 4.5.0 (release).

## R CMD check results

0 errors | 0 warnings | 0 notes.

`R CMD check --as-cran` is clean on win-builder (R-devel) and on the local Linux machine.
(The local run additionally warned that the system tool `qpdf` is not installed for the
PDF size-reduction check — a local tooling matter that does not arise on CRAN.)

## Reverse dependencies

There are no reverse dependencies on CRAN.
