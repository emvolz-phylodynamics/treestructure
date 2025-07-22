---
title: 'treestructure: An R package to detect population structure in phylogenetic
    trees'
tags:
  - R
  - coalescence
  - phylogenetic tree
authors:
  - name: Fabricia F. Nascimento
    orcid: "0000-0001-9426-6835"
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Vinicius B. Franceschi
    orcid: "0000-0002-0006-9337"
    affiliation: 1
  - name: Erik M. Volz
    orcid: "0000-0001-6268-8937"
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: MRC Centre for Global Infectious Disease Analysis, School of Public Health, 
     Imperial College London, London, UK
   index: 1

date: 19 May 2025
bibliography: references.bib
---

# Summary

How population structure can shape genetic diversity is a longstanding problem in population genetics. While the use of geographic locations, when available, can help answer some of these questions, it is still difficult to determine population structure when metadata (such as geographic location) is not available or when the potential population structure is not easily observed. 

The treestructure R package implements a statistical test based on the coalescent theory to detect unobserved population structure in a time-scaled phylogenetic tree. A time-scaled phylogenetic tree shows the evolutionary relationship between organisms in which time is in units of calendar time.
This R package has been previously described [@Volz:2020], but we have added new features for population structure detection and for adding additional samples to a previous treestructure object. 


# Statement of need

treestructure is an R package developed to find clusters within a time-scaled phylogeny that are likely to show a distinct population structure, such as, demographic or epidemiological history [@Volz:2020]. The treestructure R package also groups clusters showing similar population structure into partitions. Even though treestructure has been previously described [@Volz:2020], we have added new functionalities to the package that make it unique: choice of significance level, use of branch support to detect clusters and adding new tips to a previous treestructure object.

For a complete detail of the algorithm used in treestructure see @Volz:2020.

For details on installation, documentation and tutorial using the new features see the [package website](https://emvolz-phylodynamics.github.io/treestructure). 



## Significance level 

In practice, users will need to specify a significance level value. Increasing the significance level in the treestructure algorithm, will also increase the number of clusters detected. However, detecting more clusters will also increase the number of false positives. 

To determine the significance level, users can use additional metadata associated with each sample, and then select the significance level which gives a set of clusters that explains the most variance in the data of interest (e.g. use the cluster as a factor in an ANOVA). 

If metadata information is not available, users can use the new feature that implements the Caliński–Harabasz index or CH-index [@Calinski_and_Harabasz:1974] which is a metric based on within- and between-cluster variance in node heights to select a significance level within treestructure. If the user decided to use the CH-index, the option `level` in the _trestruct_ function should be set to NULL and a lower and upper bound for optimizing the significance level should be provided.


## Branch support 

We have also implemented the use of branch support (e.g. bootstrap and posterior probability) to detect clusters in treestructure. To use this functionality, the time-scale tree should have all node support values and the user will need to define a node support threshold value between 0 and 100. Nodes with support value less than the threshold value will not be tested. This feature is very useful to disregard low branch supports when assigning clusters using treestructure.


## Adding new samples to previous a treestructure object

Without the need to run multiple treestructure analyses, users can now update a treestructure object with new tips observed in a new phylogenetic tree. This new phylogenetic tree does not need to be time-scaled or binary.

This new feature is implemented in the _addtips_ function in the treestructure R package. The function _addtips_ will compare the new phylogenetic tree to the old treestructure object and it will merge the tips of the new tree into the treestructure object. Tips in the new tree that are not in the new treestructure object will be merged. Merging is carried out based on a phylogenetic criterion. The new tips are added to the cluster which shares the MRCA (most recent common ancestor). 


# Acknowledgements


# References
