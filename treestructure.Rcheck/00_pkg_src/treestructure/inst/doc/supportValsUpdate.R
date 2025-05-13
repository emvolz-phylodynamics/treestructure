## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7,
  fig.height=11
)

## ----message=FALSE------------------------------------------------------------
library(treestructure)

## -----------------------------------------------------------------------------
mltr2_outl_rm <- readRDS( system.file('mltr2_outl_rm_sc2_feb2020.rds', 
                                      package='treestructure') )

ggtree::ggtree(mltr2_outl_rm)

## -----------------------------------------------------------------------------
timetr2_phylo <- readRDS( system.file('timetr2_phylo_sc2_feb2020.rds', 
                                      package='treestructure') )

ggtree::ggtree(timetr2_phylo)

## ----eval=FALSE---------------------------------------------------------------
#  trestruct_res_nobt <- trestruct(timetr2_phylo,
#                                  minCladeSize = 30,
#                                  nodeSupportValues = FALSE,
#                                  level = 0.01)

## -----------------------------------------------------------------------------
trestruct_res_nobt <- readRDS( system.file('trestruct_res_nobt.rds',
                                           package='treestructure') )

plot(trestruct_res_nobt, use_ggtree = T) + ggtree::geom_tippoint()

## -----------------------------------------------------------------------------
timetr2_boot <- as.integer(mltr2_outl_rm$node.label)
timetr2_boot[is.na(timetr2_boot)] <- 95
print(timetr2_boot)

## ----eval=FALSE---------------------------------------------------------------
#  trestruct_res <- trestruct(timetr2_phylo,
#                             minCladeSize = 30,
#                             nodeSupportValues = timetr2_boot,
#                             nodeSupportThreshold = 80,
#                             level = 0.01)

## -----------------------------------------------------------------------------
trestruct_res <- readRDS( system.file('trestruct_res.rds',
                                      package='treestructure') )

plot(trestruct_res, use_ggtree = T) + ggtree::geom_tippoint()

## -----------------------------------------------------------------------------
mltr_addtips <- readRDS( system.file('mltr_addtips_mar2020.rds', 
                                     package='treestructure') )
ggtree::ggtree(mltr_addtips)

## -----------------------------------------------------------------------------
trestruct_add_tips <- addtips(trst = trestruct_res, tre = mltr_addtips)
plot(trestruct_add_tips, use_ggtree = T) + ggtree::geom_tippoint()

