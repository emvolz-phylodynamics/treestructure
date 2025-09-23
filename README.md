# treestructure

The R package treestructure contains function to detect population structure from 
the history of coalescent events recorded in phylogenetic trees. 
This method classifies each tip and internal node of a tree into disjoint sets 
characterized by similar coalescent patterns.

It also possible to use node support values (i.e. bootstrap values) to avoid 
designating population structure in badly supported clades.



## Installation

```r
# You will need to install the R package devtools 
# (https://github.com/r-lib/devtools)

install.packages("devtools")
devtools::install_github("emvolz-phylodynamics/treestructure")
```

Alternatively, you can download this repository from [here](https://github.com/emvolz-phylodynamics/treestructure/) and [follow instructions for your OS on installing packages from source](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

We also provide a command-line interface for `treestructure`. See `inst/tscl`.
If you would like to use the command-line interface you will need to install
the R package [getopt](https://github.com/trevorld/r-getopt).

You can then use the flag `-h` for how to use the `treestructure` command-line 
interface direct in a terminal.


## Tutorials 

* We recommend that you read the [Get started](http://emvolz-phylodynamics.github.io/treestructure/articles/treestructure.html) to 
understand the basic functions in `treestructure`.

* For an additional tutorial on how to include branch support to avoid designating 
population structure in badly supported clades, see this 
[example](http://emvolz-phylodynamics.github.io/treestructure/articles/supportValues.html) for Ebola.


## Author

treestructure has been developed by [Erik Volz](https://profiles.imperial.ac.uk/e.volz).


## Citation

Erik Volz, Wiuf Carsten, Yonatan Grad, Simon Frost, Ann Dennis, Xavier Didelot, 
(2020); [Identification of hidden population structure in time-scaled phylogenies](https://academic.oup.com/sysbio/article/69/5/884/5734655); 
Systematic Biology, 69: 884–896.
