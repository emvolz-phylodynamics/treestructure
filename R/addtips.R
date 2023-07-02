#' Compares a new input tree to an old treestructure fit and merges tips into a new treestructure object. 
#' Tips in the new tree that are not in the new treestructure will be merged. 
#' 
#' Merging is carried out based on a phylogenetic criterion. The new tips are added to the cluster which shares its MRCA. 
#'
#' @author Erik Volz <erik.volz@gmail.com>
#'
#' @param trst Original treestructure fit that that will be updated. 
#' @param tre A new tree (ape::phylo) which may contain samples not in trst. This tree must be rooted, but does not need to be time-scaled or binary. 
#' @return A new treestructure fit. 
#' @examples 
#' set.seed(072023)
#' # simulate two trees and bind them to simulate structure
#' tr1 <- ape::rcoal( 50 )
#' tr2 <- ape::rcoal( 100 )
#' tr1$tip.label <- gsub(tr1$tip.label, patt = 't', rep = 's')
#' tr1$edge.length <- tr1$edge.length*.5 
#' tr1$root.edge <- 1 
#' tr2$root.edge <-1 
#' tr <- ape::bind.tree(tr1, tr2, position=.5 ) |> ape::multi2di()
#' # subsample the tree to simulating missing tips and estimate structure
#' ex <- sample( tr$tip.label, size = 30, replace=FALSE)
#' tr0 <- ape::drop.tip( tr, ex ) 
#' (s0 <- treestructure::trestruct( tr0 ))
#' # assign structure to the previously missing tips 
#' (s <- treestructure::addtips( s0, tr ))
#' @export
addtips <- function(trst, tre)
{
        notintree <- setdiff( trst$data$taxon , tre$tip.label)
        if ( length( notintree ) > 0 )
        {
                message( 'Some taxa in *trst* were not found in tre. These will be excluded from the output. 
                        ' )
        }
        clusterSets <- lapply( trst$clusterSets, function(x) intersect( x, tre$tip.label ))
        mrcas <- sapply( clusterSets, function(x) ape::getMRCA( tre, x )) 

        poe <- postorder( tre )
        desc <- lapply( 1:(Ntip(tre)+Nnode(tre)), function(x) c() )
        for ( ie in poe )
        {
                a <- tre$edge[ie,1]
                u <- tre$edge[ie,2]
                desc[[a]] <- c( desc[[a]], desc[[u]], u)
        }

        toadd <- setdiff( tre$tip.label ,  trst$data$taxon )
        toadd_nodes <- match( toadd, tre$tip.label )

        cluster_desc <- desc[mrcas]
        cluster_ranks <- sapply( cluster_desc, function(x) sum(mrcas %in% x))
        .assign <- function(u){
                whcl <- sapply( cluster_desc, function(x) u%in%x ) |> which() 
                v <- whcl[ which.min(cluster_ranks[whcl]) ]
                v
        }
        sapply( toadd_nodes, .assign) -> newclust
        stdf0 <- trst$data[ !duplicated(trst$data$cluster ), ]
        cl2pa <- setNames( stdf0$partition, stdf0$cluster )
        stdf <- rbind( 
                      trst$data 
                      , data.frame( 
                                   taxon = toadd 
                                   , cluster = newclust 
                                   , partition = cl2pa[ as.character(newclust ) ]
                                   , stringsAsFactors = FALSE 
                      )
                      )
        rownames(stdf) <- 1:nrow(stdf)
        stdf <- stdf[ stdf$taxon %in% tre$tip.label,  ]
        stdf <- stdf[ match(stdf$taxon, tre$tip.label), ]
        
        trst_out <- trst
        
        trst_out$data <- stdf 
        trst_out$clustering <- as.factor( stdf$cluster  ) |> setNames( stdf$taxon )
        trst_out$partition <- as.factor( stdf$partition ) |> setNames( stdf$taxon )
        trst_out$clusterSets <- lapply( split(stdf, stdf$cluster), function(x) x$taxon )
        trst_out$partitionSets <- lapply( split(stdf, stdf$partition), function(x) x$taxon )
        trst_out$D <- NULL 
        trst_out$tree <- tre 
        trst_out$cluster_mrca <- mrcas 
        trst_out
}


