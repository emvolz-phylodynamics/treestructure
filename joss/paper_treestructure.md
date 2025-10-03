---
title: 'treestructure: An R package to detect population structure in phylogenetic
  trees'
tags:
- R
- coalescence
- phylogenetic tree
date: "03 October 2025"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Fabricia F. Nascimento
  orcid: "0000-0001-9426-6835"
  affiliation: 1
- name: Vinicius B. Franceschi
  orcid: "0000-0002-0006-9337"
  affiliation: 1
- name: Erik M. Volz
  orcid: "0000-0001-6268-8937"
  corresponding: true
  affiliation: 1
bibliography: references.bib
affiliations:
- name: MRC Centre for Global Infectious Disease Analysis, School of Public Health,
    Imperial College London, London, UK
  index: 1
---

# Summary

How population structure can shape genetic diversity is a longstanding problem in 
population genetics. While the use of geographic locations, when available, can 
help answer some of these questions, it is still difficult to determine population 
structure when such metadata is not available or when the potential population 
structure is not easily observed. Methods developed to detect population 
structure have been developed, such as _CaveDive_ [@Helekal:2022] and 
_fastbaps_ [@Tonkin-Hill2019], and applied to the detection of outbreaks and 
variant surveillance [@Reimche:2023; @Binney:2025]. 

Here we present `treestructure`, an R package that has been previously 
described [@Volz:2020] and used in a variety of studies, such as detection of 
lineages with different demographic histories in SARS-CoV-2 [@Fountain-Jones:2020]; 
detection of fitness advantage on clades that showed similar demographic 
histories in _Neisseria gonorrhoea_ [@Joseph:2022]; and understanding a population 
of _Vibrio parahaemolyticus_ in Latin America [@Campbell:2024].

The `treestructure R package implements a statistical test based on the coalescent 
theory to detect unobserved population structure in a time-scaled phylogenetic tree. 
A time-scaled phylogenetic tree shows the evolutionary relationship between 
organisms in units of calendar time. We have now added new features to treestructure 
for the detection of population structure and for adding additional samples to a 
previous treestructure object. 


# Statement of need

`treestructure` is an R package developed to find clusters within a time-scaled 
phylogeny that are likely to show a distinct population structure, such as, 
demographic or epidemiological history [@Volz:2020]. The `treestructure` R package 
also groups clusters showing similar population structure into partitions. Here, 
we describe new functionalities added to the package, enhancing its practical 
utility and statistical robustness: 1) Methods to automatically choose clustering 
hyperparameters; 2) use of branch support values (_e.g._ bootstrap and posterior 
clade credibility) to filter out clusters with low phylogenetic confidence, 
and 3) adding new tips to a previous `treestructure` object, allowing clusters to 
be updated in an online fashion as new data becomes available.

For details of the main algorithm used in `treestructure` see @Volz:2020.

For details on installation, documentation and tutorial using the new features 
see the [package website](https://emvolz-phylodynamics.github.io/treestructure). 


## Clustering significance level 

Clustering methods require the specification of hyperparameters that specify how 
aggressively a method will partition data. The `treestructure` algorithm makes use 
of a _rank-sum significance level_, and clusters are defined when a coalescent-based 
statistical test detects a difference according to this level. Decreasing the 
significance level in the `treestructure` algorithm will increase the number of clusters detected. 
However, detecting more clusters will also increase the number of false positive detections. 

To determine the significance level, users can use additional metadata associated 
with each sample, and then select the significance level which gives a set of 
clusters that explains the most variance in the data of 
interest (_e.g._ use the cluster as a factor in an ANOVA). 

If metadata information is not available, users can use the new feature that 
implements the Caliński–Harabasz index or CH-index [@Calinski_and_Harabasz:1974] 
which is a metric based on within- and between-cluster variance in a given statistic 
to select a quasi-optimal significance level.  Within `treestructure` we use the 
node heights of the phylogeny itself, observed within each cluster, as the 
statistic that is used when computing the CH-index. Thus clusters are selected 
such that there is high between-cluster variance in phylogenetic node heights. 
If the user decided to use the CH-index, the option `level` in the _trestruct_ 
function should be set to NULL and a lower and upper bound for optimizing the 
significance level should be provided. A step-by-step tutorial on how to run such 
analysis can be found 
[here](https://emvolz-phylodynamics.github.io/treestructure/articles/supportValues.html).


## Branch support 

Whereas there is often a great deal of uncertainty about individual phylogenetic 
splits, a user may not want to cluster their data along branches which are poorly 
supported. We have also implemented the use of branch support (e.g. bootstrap 
and posterior probability) to refine clusters in `treestructure`. To use this 
functionality, the time-scaled tree should be annotated with node support values, 
and the user will need to define a node support threshold value between 0 and 100. 
Nodes with support value less than the threshold value will not be tested. 
This feature is very useful to filter out clusters that may not correspond to 
real phylogenetic splits.

\autoref{fig:strucfig} shows an example on how the use of node support can filter 
out clusters with low phylogenetic confidence.


![Down-sampled time-scaled phylogenetic tree for Ebola [publicly available](https://github.com/ebov/space-time/blob/master/Data/Makona_1610_cds_ig.GLM.MCC.tree) 
with 150 sequences. Clusters obtained by running `treestructure` **A.** without 
the use of node support and **B.** using node support threshold of 95. For both 
analyses we used a significance level of 0.01 and minimum clade size of 15 sequences. 
For an example analysing the complete Ebola dataset, a tutorial can be found [here](https://emvolz-phylodynamics.github.io/treestructure/articles/supportValues.html).
\label{fig:strucfig}](Figure/ebola_plots_joss_paper.png)


## Online inference by adding new samples to previous a treestructure object

Without the need to run multiple `treestructure` analyses, users can now update 
a `treestructure` object with new samples observed in a phylogenetic tree. The 
updated tree does not need to be time-scaled or binary, reducing the need for 
expensive computation.

This new feature is implemented in the _addtips_ function in the treestructure 
R package. The function _addtips_ will compare the new phylogenetic tree to the 
old `treestructure` object and it will merge the tips of the new tree into the 
`treestructure` object. Merging is carried out based on a phylogenetic criterion: 
New tips are added to the cluster which includes its most recent common ancestor 
in the new phylogeny. A step-by-step tutorial on how to use this feature can be 
found [here](https://emvolz-phylodynamics.github.io/treestructure/articles/updating_treestructure.html).


# Acknowledgements

We would like to acknowledge Oliver Stirrup for previous contribution to the package. 
EV acknowledges support from the UK Health Security Agency CARAA 104683ED 
“Development of phylogenetic analysis tools to track HIV transmissions for public 
health surveillance purposes” and CARAA 5126118 “Provision of expert advice and 
development support for emerging infections genomics and metagenomics analysis”.
We also thank support from the Wellcome Trust (Investigator Award 220885/Z/20/Z
“Population genomic analysis of HIV transmission fitness” awarded to EV).


# References


