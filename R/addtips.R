#' Compare and add tips into new treestructure object
#'
#' Compares a new input tree to an old treestructure fit and merges tips into a
#' new treestructure object.
#' Tips in the new tree that are not in the new treestructure will be merged.
#' Merging is carried out based on a phylogenetic criterion. The new tips are
#' added to the cluster which shares its MRCA (most recent common ancestor).
#'
#' @details
#' The MRCA of each cluster is located in the new tree, and every node then
#' inherits the cluster of its nearest ancestral MRCA in a single preorder pass.
#' Each new tip is assigned by direct lookup, so the cost is linear in the size
#' of \code{tre}. This makes it practical to place millions of tips, e.g. when
#' assigning lineages to a large maximum-likelihood tree from a fit on a
#' subsample.
#'
#' If several clusters share the same MRCA in the new tree, tips below it are
#' assigned to the cluster with the lowest index. New tips that do not descend
#' from the MRCA of any cluster are given \code{NA} for cluster and partition.
#'
#' @author Erik Volz, Vinicius Franceschi
#'
#' @param trst Original treestructure fit that that will be updated.
#' @param tre A new tree (ape::phylo) which may contain samples not in trst.
#'    This tree must be rooted, but does not need to be time-scaled or binary.
#' @param verbose If TRUE, report progress with \code{message()}.
#' @return A new treestructure fit.
#' @examples
#' set.seed(072023)
#' # simulate two trees and bind them to simulate structure
#' tr1 <- ape::rcoal( 50 )
#' tr2 <- ape::rcoal( 100 )
#' tr1$tip.label <- gsub(tr1$tip.label, patt = 't', rep = 's')
#' tr1$edge.length <- tr1$edge.length*.5
#' tr1$root.edge <- 1
#' tr2$root.edge <- 1
#' tr <- ape::bind.tree(tr1, tr2, position = .5 ) |> ape::multi2di()
#'
#' # subsample the tree to simulating missing tips and estimate structure
#' ex <- sample( tr$tip.label, size = 30, replace = FALSE)
#' tr0 <- ape::drop.tip( tr, ex )
#' (s0 <- treestructure::trestruct( tr0 ))
#'
#' # assign structure to the previously missing tips
#' (s <- treestructure::addtips( s0, tr ))
#' @export
addtips <- function(trst, tre, verbose = FALSE)
{
	notintree <- setdiff( trst$data$taxon , tre$tip.label)
	if ( length( notintree ) > 0 )
	{
		message( length( notintree ), ' taxa in *trst* were not found in tre. These will be excluded from the output.' )
	}
	clusterSets <- lapply( trst$clusterSets, function(x) intersect( x, tre$tip.label ))

	# MRCA of each cluster in the new tree; a single remaining tip is its own MRCA
	if ( verbose ) message( 'Computing MRCAs of ', length( clusterSets ), ' clusters' )
	mrcas <- vapply( clusterSets, function(x) {
		if ( length( x ) == 0 )
			return( NA_integer_ )
		if ( length( x ) == 1 )
			return( match( x, tre$tip.label ) )
		as.integer( ape::getMRCA( tre, x ) )
	}, integer(1) )

	# every node inherits the cluster of its nearest ancestral MRCA, in one preorder pass.
	# MRCAs are labelled in reverse so that the lowest cluster index wins when MRCAs coincide.
	if ( verbose ) message( 'Propagating cluster labels over ', nrow( tre$edge ), ' edges' )
	node_cluster <- integer( ape::Ntip(tre) + ape::Nnode(tre) ) # 0 = no ancestral MRCA
	ci <- rev( which( !is.na( mrcas ) ) )
	node_cluster[ mrcas[ci] ] <- ci
	preedges <- tre$edge[ rev( ape::postorder( tre ) ), , drop = FALSE ]
	for ( ie in seq_len( nrow( preedges ) ) )
	{
		u <- preedges[ie,2]
		if ( node_cluster[u] == 0 )
			node_cluster[u] <- node_cluster[ preedges[ie,1] ]
	}

	toadd <- setdiff( tre$tip.label ,  trst$data$taxon )
	toadd_nodes <- match( toadd, tre$tip.label )
	newclust <- node_cluster[ toadd_nodes ]
	newclust[ newclust == 0 ] <- NA
	if ( any( is.na( newclust ) ) )
		message( sum( is.na( newclust ) ), ' new tips do not descend from the MRCA of any cluster and are assigned NA.' )
	if ( verbose ) message( 'Assigned ', sum( !is.na( newclust ) ), ' of ', length( toadd ), ' new tips' )

	stdf0 <- trst$data[ !duplicated(trst$data$cluster ), ]
	cl2pa <- stats::setNames( stdf0$partition, stdf0$cluster )
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
	trst_out$clustering <- as.factor( stdf$cluster  ) |> stats::setNames( stdf$taxon )
	trst_out$partition <- as.factor( stdf$partition ) |> stats::setNames( stdf$taxon )
	trst_out$clusterSets <- lapply( split(stdf, stdf$cluster), function(x) x$taxon )
	trst_out$partitionSets <- lapply( split(stdf, stdf$partition), function(x) x$taxon )
	trst_out$D <- NULL
	trst_out$tree <- tre
	trst_out$cluster_mrca <- mrcas
	trst_out
}
