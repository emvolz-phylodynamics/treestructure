# Tests for addtips(): agreement with the original descendant-list algorithm, tie-breaking
# when cluster MRCAs coincide, and new tips outside every cluster MRCA.

# the original (quadratic) addtips assignment, kept as a reference: each new tip joins the
# ancestral cluster MRCA with the fewest other MRCAs below it, ie its nearest ancestral MRCA
.addtips_reference <- function(trst, tre){
	clusterSets <- lapply( trst$clusterSets, function(x) intersect( x, tre$tip.label ))
	mrcas <- sapply( clusterSets, function(x) ape::getMRCA( tre, x ))
	desc <- lapply( 1:(ape::Ntip(tre)+ape::Nnode(tre)), function(x) c() )
	for ( ie in ape::postorder( tre ) )
		desc[[ tre$edge[ie,1] ]] <- c( desc[[ tre$edge[ie,1] ]], desc[[ tre$edge[ie,2] ]], tre$edge[ie,2] )
	cluster_desc <- desc[mrcas]
	cluster_ranks <- sapply( cluster_desc, function(x) sum(mrcas %in% x))
	toadd <- setdiff( tre$tip.label, trst$data$taxon )
	newclust <- sapply( match( toadd, tre$tip.label ), function(u){
		whcl <- which( sapply( cluster_desc, function(x) u %in% x ) )
		whcl[ which.min( cluster_ranks[whcl] ) ]
	})
	stats::setNames( newclust, toadd )
}

.structured_tree <- function(n, seed){
	set.seed(seed)
	ks <- round( c(.4, .3, .2, .1) * n )
	trs <- lapply( seq_along(ks), function(i){
		t <- ape::rcoal( ks[i] )
		t$edge.length <- t$edge.length * i / 2
		t$tip.label <- paste0( 'p', i, '_', seq_len(ks[i]) )
		t$root.edge <- 1
		t
	})
	tr <- trs[[1]]
	for (i in 2:4) tr <- ape::bind.tree( tr, trs[[i]], position = .3 )
	ape::multi2di( tr )
}

test_that("addtips matches the original algorithm", {
	tr <- .structured_tree( 600, 11 )
	set.seed(1)
	s0 <- trestruct( ape::keep.tip( tr, sample( tr$tip.label, 300 ) ), verbosity = 0 )
	# restrict to tips below some cluster MRCA, where the original algorithm is defined
	m <- sapply( s0$clusterSets, function(x) ape::getMRCA( tr, x ) )
	tr <- ape::keep.tip( tr, unique( unlist( lapply( m, function(a) ape::extract.clade( tr, a )$tip.label ) ) ) )

	s <- addtips( s0, tr )
	ref <- .addtips_reference( s0, tr )
	expect_s3_class( s, 'TreeStructure' )
	expect_setequal( names( s$clustering ), tr$tip.label )
	expect_equal( as.integer( as.character( s$clustering[ names(ref) ] ) ), unname( ref ) )
})

test_that("tips under a shared MRCA join the lowest-index cluster", {
	set.seed(5)
	tr <- ape::rcoal( 40 )
	left <- ape::extract.clade( tr, tr$edge[ tr$edge[,1] == 41, 2 ][1] )$tip.label
	right <- setdiff( tr$tip.label, left )
	# both clusters straddle the root, so both have the root as MRCA
	d <- data.frame( taxon = c( left[1:3], right[1:3], left[4:6], right[4:6] )
		, cluster = rep( 1:2, each = 6 ), partition = rep( 1:2, each = 6 ) )
	trst <- list( data = d, clusterSets = split( d$taxon, d$cluster ) )
	s <- addtips( trst, tr )
	newtips <- setdiff( tr$tip.label, d$taxon )
	expect_true( all( as.character( s$clustering[ newtips ] ) == '1' ) )
	expect_equal( unname( .addtips_reference( trst, tr ) ), rep( 1L, length( newtips ) ) )
})

test_that("new tips outside every cluster MRCA get NA with a message", {
	set.seed(6)
	tr <- ape::rcoal( 40 )
	left <- ape::extract.clade( tr, tr$edge[ tr$edge[,1] == 41, 2 ][1] )$tip.label
	d <- data.frame( taxon = left[1:4], cluster = 1L, partition = 1L )
	trst <- list( data = d, clusterSets = split( d$taxon, d$cluster ) )
	expect_message( s <- addtips( trst, tr ), 'assigned NA' )
	outside <- setdiff( tr$tip.label, ape::extract.clade( tr, ape::getMRCA( tr, left[1:4] ) )$tip.label )
	expect_true( all( is.na( s$clustering[ outside ] ) ) )
	expect_false( anyNA( s$clustering[ left ] ) )
})
