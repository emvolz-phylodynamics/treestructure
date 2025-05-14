## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 11
)

## ----message=FALSE------------------------------------------------------------
library(treestructure)

## -----------------------------------------------------------------------------
( tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') ) ) 

## ----message=FALSE------------------------------------------------------------
s <- trestruct( tree ) 

## -----------------------------------------------------------------------------
print(s) 

## ----message=FALSE------------------------------------------------------------
plot(s)  + ggtree::geom_tiplab() 

## -----------------------------------------------------------------------------
plot( s, use_ggtree = FALSE )

## -----------------------------------------------------------------------------
structureData <- as.data.frame( s ) 
head( structureData )

## -----------------------------------------------------------------------------
with ( structureData, 
       ape::keep.tip(s$tree, taxon[ partition==1 ] ) 
       ) -> partition1
partition1
plot(partition1)

## ----message=FALSE------------------------------------------------------------
trestruct( tree, level = .05, minCladeSize = 5 )

