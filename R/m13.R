#~ Treestructure
#~     Copyright (C) 2019  Erik Volz
#~     This program is free software: you can redistribute it and/or modify
#~     it under the terms of the GNU General Public License as published by
#~     the Free Software Foundation, either version 3 of the License, or
#~     (at your option) any later version.


.get_anc <- function(edge, u ){
	edge[ edge[,2]==u ,1]
}
.get_dgtrs <- function(edge, a){
	edge[ edge[,1] == a , 2 ]
}

.get_rootnode <- function(edge){
	setdiff( as.vector(edge), edge[,2] )
}

cocluster_accuracy <- function( x, y ){
# Statistic which measures overlap between equal length factors x and y
 n <- length(x)
 a <- matrix( FALSE, nrow =	n , ncol = n  )
 b <- matrix( FALSE, nrow =	n , ncol = n  )
 sx <- split( 1:n, x )
 sy <- split( 1:n, y )
 for ( z in sx )
   a[z,z] <- TRUE
 for ( z in sy )
   b[z,z] <- TRUE

  (sum( a == b ) - n ) / (n^2-n)
}




# DISSIMILARITY FUNCTIONS
# test statistic comparing two internals

# req ancestors
.ispara <- function(u, v, ancestors){
	(v %in% ancestors[[u]]) | (u %in% ancestors[[v]])
}


# requires tre
.ismonomono_cl <- function( tre, ancestors, cl1, cl2){
	t1 <- intersect( 1:(ape::Ntip(tre)), cl1 )
	t2 <- intersect( 1:(ape::Ntip(tre)), cl2 )
	u <- ape::getMRCA(tre,  t1 )
	v <- ape::getMRCA(tre,  t2 )
	if ( is.null( u ) | is.null(v) )
	  return(FALSE)
	!.ispara( u, v, ancestors)
}

# req nhs
.rank.sum <- function( nhs, uinternals, vinternals ) {
	x =rbind(
	 cbind(  nhs[ uinternals ], 1 )
	 , cbind(  nhs[ vinternals ], 0 )
	)
	x <- as.vector( x[ order(x[,1] ) , 2] )
	sum( x * (1:length(x)) )
}

# test stat null distribution
# E_i : event type (cf paper).
# 1: as wiuf and donelly; u is mono relative to v; mono/mono or mono/para
# 2: mono / mono
# 3: mono/mono OR mono/para OR para/mono
# requires nhs , tre , inodes
.uv.diss <- function(tre, nhs, inodes, n, ancestors, uset, vset, nsim , Ei = NA, returnabs=TRUE, detailedOut = FALSE){
	# 1 sample u
	# 0 co
	# -1 sample v

	if (is.na( Ei)){
		#nested = ifelse( .ismonomono( u, v), 1, 0 )
		Ei = ifelse( .ismonomono_cl( tre, ancestors, uset, vset ), 2, 3 )
	}

	uinternals <-  intersect( uset, inodes)
	vinternals <- setdiff(  intersect( vset, inodes), uinternals )
	utips <- intersect( uset, 1:n)
	vtips <- intersect( vset, 1:n)
	if ( (length( utips )==0) | (length(vtips) == 0) ){
		if ( detailedOut )
			return( list( statistic = 0 , null = c() , Z = x) )
		else
			return( 0 )
	}
	x <- rbind(
	 cbind(    nhs[ uinternals ], 0 )
	 , cbind(  nhs[ vinternals ], 0 )
	 , cbind(  nhs[ utips ], 1 )
	 , cbind(  nhs[ vtips ], -1 )
	)
	x <- as.vector( x[ order(x[,1] ) , 2] )
	nd <- Cuv_ranksum_nulldist(x, nsim, Ei  )

	rsuv <- .rank.sum(nhs, uinternals, vinternals)

	m_nd <- mean(nd)
	sd_nd <- stats::sd(nd)
	if ( returnabs )
		x = ( abs( (rsuv - m_nd) / sd_nd ) )
	else
		x = ( (rsuv - m_nd) / sd_nd )

	if ( detailedOut ){
		return( list ( statistic = unname( rsuv ), null = nd , Z = x )  )
	}
	return ( x )

	if (FALSE){ # alternative empirical measure:
		s <- sum(rsuv > nd )
		p <- min( s / nsim , (nsim - s) / nsim )

		p <- sum(rsuv > nd ) / nsim
		abs( qnorm( p )  )
	}
}


#


.tredat <- function ( tre ){
	n <- ape::Ntip( tre )
	nnode <- ape::Nnode(tre )
	rootnode <- .get_rootnode( tre$edge )
	inodes = iin <- (n+1):(n+nnode)

	D <- ape::node.depth.edgelength( tre )
	rh <- max( D[1:n] )
	sts <- D[1:n]
	shs <- rh - sts
	nhs = nodeheights <- rh - D

	poedges <- tre$edge[ ape::postorder( tre ), ]
	preedges <- tre$edge[ rev( ape::postorder( tre )), ]

	# PRECOMPUTE for each node
	# descendants; note also counts mrca node
	descendants <- lapply( 1:(n+nnode), function(u) u )
	for (ie in 1:nrow(poedges)){
		a <- poedges[ie,1]
		u <- poedges[ie,2]
		descendants[[a]] <- c( descendants[[a]], descendants[[u]] )
	}


	# descendant tips
	descendantTips <- lapply( 1:(n+nnode), function(u) c() )
	for (u in 1:n)
	  descendantTips[[u]] <- u
	for (ie in 1:nrow(poedges)){
		a <- poedges[ie,1]
		u <- poedges[ie,2]
		descendantTips[[a]] <- c( descendantTips[[a]], descendantTips[[u]] )
	}


	# descendant internal
	descInternal <- lapply( 1:(n+nnode), function(u) c() )
	for (u in (n+1):(n+nnode))
	  descInternal[[u]] <- u
	for (ie in 1:nrow(poedges)){
		a <- poedges[ie,1]
		u <- poedges[ie,2]
		descInternal[[a]] <- c( descInternal[[a]], descInternal[[u]] )
	}

	# ancestors
	ancestors <- lapply( 1:(n+nnode), function(u) c() )
	for (ie in 1:nrow(preedges)){
		a <- preedges[ie,1]
		u <- preedges[ie,2]
		ancestors[[u]] <- c( ancestors[[a]], a)
	}

	# number tips desc
	ndesc <- sapply( 1:(n+nnode), function(u) length( descendantTips[[u]] ) )

	# max dist to tip
	maxDistToTip <- rep( 0, n + nnode )
	for (ie in 1:nrow(poedges)){
		a <- poedges[ie,1]
		u <- poedges[ie,2]
		maxDistToTip[a] <- max( maxDistToTip[a] , maxDistToTip[u] + 1)
	}

	# vector heights of all descending
	descHeights <- lapply( 1:(n+nnode), function(u) nhs[ descendants[[u]] ] )
	descInternalHeights <- lapply( 1:(n+nnode), function(u) nhs[ descInternal[[u]] ] )

	# time range of each clade
	timerange <- t( sapply( 1:(n+nnode), function(u) range( descHeights[[u]]  ) ))

	# vector heights of tips descending
	descendantTipHeights <- lapply( 1:(n+nnode), function(u) nhs[ descendantTips[[u]] ] )

	prenodes <- unique( preedges[,1] )

	dgtrMat <- matrix( NA, nrow = ape::Nnode(tre)+ape::Ntip(tre), ncol =2)
	for (ie in 1:nrow(poedges))
	{
		if (is.na( dgtrMat[ poedges[ie,1] , 1] ) )
		  dgtrMat[ poedges[ie,1], 1] <- poedges[ie,2]
		else
		  dgtrMat[ poedges[ie,1], 2] <- poedges[ie,2]
	}

	# map node to branch length connecting to ancestor (useful for filtering multifurcations)
	node2edgelength <- rep(NA, n + nnode )
	node2edgelength[ tre$edge[,2] ] <-  tre$edge.length

	list (n = n , nnode = nnode , rootnode = rootnode, inodes = inodes
	  , D = D, rh = rh, sts = sts, shs = shs , nhs = nhs
	  , poedges = poedges, preedges = preedges, descendants = descendants
	  , descendantTips = descendantTips , descInternal = descInternal
	  , ancestors = ancestors
	  , ndesc = ndesc
	  , descHeights = descHeights
	  , descInternalHeights = descInternalHeights
	  , timerange = timerange
	  , descendantTipHeights = descendantTipHeights
	  , prenodes = prenodes
	  , dgtrMat = dgtrMat
	  , maxDistToTip = maxDistToTip
	  , node2edgelength = node2edgelength
	  , tre = tre
	)
}

#' Test treestructure hypothesis
#'
#' Test the hypothesis that two clades within a tree were generated by the same
#' coalescent process.
#'
#' @param tre An ape::phylo tree, must be binary and rooted
#' @param x A character vector of tip labels or numeric node numbers. If numeric,
#'    can include internal node numbers.
#' @param y as x, but must be disjoint with x
#' @param nsim Number of simulations (larger = slower and more accurate)
#'
#' @examples
#' tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') )
#'
#' (struc <- trestruct( tree ))
#'
#' #run the test
#'
#' results <- treestructure.test(tree, x = struc$clusterSets[[1]],
#'                               y = struc$clusterSets[[2]])
#'
#' print(results)
#' @export
treestructure.test <- function( tre, x, y, nsim = 1e4 )
{
	stopifnot( ape::is.rooted(tre))
	stopifnot( ape::is.binary(tre))
	stopifnot( length( intersect( x, y )) == 0 )
	if ( any(is.na(tre$tip.label)) | anyDuplicated(tre$tip.label)){
		message('Tree has NA or duplicated tip labels. Adding a unique id.')
		tre$tip.label <- paste0( paste0( 1:(ape::Ntip(tre)), '_' ), tre$tip.label )
	}
	tredat = .tredat( tre )
	#attach( tredat )

	if ( is.numeric( x ))
		uset = x
	else
		uset = match( x, tredat$tre$tip.label )
	if ( is.numeric(y))
		vset = y
	else
		vset = match( y, tredat$tre$tip.label )

	if ( any( is.na( uset ) ) | any (is.na( vset )) )
		stop ( 'Some tip labels in x or y could not be matched tre$tip.label. Check inputs.' )

	Uset = unique( c(uset,  do.call( c, tredat$ancestors[uset] ) ) )
	Uset <- setdiff( Uset, tredat$ancestors[[ ape::getMRCA( tredat$tre, uset ) ]] )
	Vset = unique( c( vset, do.call( c, tredat$ancestors[vset] )) )
	Vset = setdiff( Vset , tredat$ancestors[[ ape::getMRCA( tredat$tre, vset ) ]] )
	uvd = .uv.diss(tredat$tre, tredat$nhs, tredat$inodes, tredat$n, tredat$ancestors, Uset, Vset, nsim = nsim, returnabs=FALSE, detailedOut = TRUE)

	res = structure( list(
	  statistic = uvd$statistic
	  , p.value = with (uvd,  min( mean(statistic < null), mean( statistic> null) ) )
	  , estimate = with (uvd,  mean( statistic> null))
	  , std.err = stats::sd ( uvd$nd )
	  , conf.int = stats::quantile( uvd$null, c(0.025 , 0.975) )
	  , null.value = stats::median( uvd$null )
	  , alternative='Alternative hypothesis: Rank sum differs from coalescent distribution'
	  , method = 'Two tailed simulation quantiles'
	  , data.name = 'tre'
	  , data = tredat$tre
	  , nsim = nsim
	)
	, class = 'treestructure.htest'
	)

	#detach( tredat )
	res$null = uvd$null
	invisible(res)
}


#' @export
print.treestructure.htest <- function(x, ... ){
	stopifnot( inherits( x, 'treestructure.htest' ))
	#~ 		One Sample t-test

	#~ data:  rnorm(10)
	#~ t = -0.01619, df = 9, p-value = 0.9874
	#~ alternative hypothesis: true mean is not equal to 0
	#~ 95 percent confidence interval:
	#~  -0.6089170  0.6002629
	#~ sample estimates:
	#~    mean of x
	#~ -0.004327037

	cat( '     Treestructure rank-sum test\n'  )
	cat( '     \n'  )
	cat( 'data: '  )
	print( x$data )
	if ( x$p.value > 0 )
		cat( sprintf( 'Rank sum = %g, p-value = %g\n', x$statistic, x$p.value ))
	else
		cat( sprintf( 'Rank sum = %g, p-value < %g\n', x$statistic, 1/x$nsim ))
	cat( x$alternative ); cat('\n' )
	cat( '95 percent confidence interval:\n')
	cat( sprintf( '  %g  %g\n', x$conf.int[1], x$conf.int[2] ))
invisible(x)
}


#' Detect cryptic population structure in time trees
#'
#' Estimates a partition of a time-scaled tree by contrasting coalescent patterns.
#'
#' @param tre A tree of type ape::phylo. Must be rooted. If the tree has multifurcations,
#'    it will be converted to a binary tree before processing.
#' @param minCladeSize All clusters within partition must have at least this many
#'    tips.
#' @param minOverlap Threshold time overlap required to find splits in a clade.
#' @param nodeSupportValues Node support values such as produced by bootstrap or
#'    Bayesian credibility scores. Must be logical or vector with length equal
#'    to number of internal nodes in the tree. If nodeSupportValues = TRUE, then
#'    the function will get the information on node support from the tree.
#'    If numeric vector, these values should be between 0 and 100.
#' @param nodeSupportThreshold Threshold node support value between 0 and 100.
#'    Nodes with support lower than this threshold will not be tested.
#' @param nsim Number of simulations for computing null distribution of test
#'    statistics.
#' @param level Significance level for finding new split within a set of tips.
#'    Can also be NULL, in which case the optimal level is found according to the
#'    CH index (see details).
#' @param ncpu If > 1 will compute statistics in parallel using multiple CPUs.
#' @param verbosity If > 0 will print information about progress of the algorithm.
#' @param debugLevel If > 0 will produce additional data in return value.
#' @param levellb If optimising the `level` parameter, this is the lower bound
#'    for the search.
#' @param levelub If optimising the `level` parameter, this is the upper bound
#'    for the search.
#' @param res If optimising the `level` parameter, this is the number of values
#'    to test.
#' @return A TreeStructure object which includes cluster and partition assignment
#'    for each tip of the tree.
#'
#' @details
#' Estimates a partition of a time-scaled tree by contrasting coalescent patterns.
#' The algorithm is premised on a Kingman coalescent null hypothesis for the
#' ordering of node heights when contrasting two clades, and a test statistic is
#' formulated based on the rank sum of node times in the tree.
#' If node support values are available (as computed by bootstrap procedures),
#' the method can optionally exclude designation of structure on poorly supported
#' nodes. The method will not designate structure on nodes with zero branch length
#' relative to their immediate ancestor.
#' The significance level for detecting significant partitions of the tree can be
#' provided, or a range of values can be examined.
#' The \href{https://en.wikipedia.org/wiki/Calinski-Harabasz_index}{CH index}
#' based on within- and between-cluster variance in node heights can be used to
#' select a significance level if none is provided.
#'
#' @section References:
#' Volz EM, Carsten W, Grad YH, Frost SDW, Dennis AM, Didelot X.
#' \href{https://academic.oup.com/sysbio/article/69/5/884/5734655?login=false}{Identification of hidden population structure in time-scaled phylogenies.}
#' Systematic Biology 2020; 69(5):884-896.
#'
#' @author Erik M Volz
#'
#'
#' @examples
#' tree <- ape::rcoal(50)
#' struct <-  trestruct( tree )
#' print(struct)
#'
#' @export
trestruct <- function( tre, minCladeSize = 25, minOverlap = -Inf, nodeSupportValues = FALSE, nodeSupportThreshold = 95, nsim = 1e4, level = .01, ncpu = 1, verbosity = 1, debugLevel=0
	, levellb = 1e-3, levelub = 1e-1, res = 11)
{
	stopifnot( ape::is.rooted(tre))
	# stopifnot( ape::is.binary(tre))
	if ( minOverlap >= minCladeSize){
		stop('*minOverlap* should be < *minCladeSize*.')
	}
	if ( any(is.na(tre$tip.label)) | anyDuplicated(tre$tip.label)){
		if ( verbosity > 0 )
			cat('Tree has NA or duplicated tip labels. Adding a unique id.\n')
		else
			message('Tree has NA or duplicated tip labels. Adding a unique id.')
		tre$tip.label <- paste0( paste0( 1:(ape::Ntip(tre)), '_' ), tre$tip.label )
	}

	useNodeSupport <- FALSE
	if (!is.logical(nodeSupportValues) & is.vector(nodeSupportValues)) {
		# rewriting node.label here because if tree is multifurcating, support values would need to apply to most ancestral node in a multifurcation
		stopifnot( length(nodeSupportValues) == ape::Nnode( tre ) )
		stopifnot( is.numeric( nodeSupportValues ) )
		stopifnot( all( nodeSupportValues >= 0 ) )
		stopifnot( all( nodeSupportValues <= 100 ) )
		useNodeSupport <- TRUE
		tre$node.label <- as.character( nodeSupportValues )
		nodeSupportValues <- TRUE
	}
	if (!is.binary(tre)){# must be done after node labels are finalised
		tre <- multi2di( tre, random=FALSE) # maintain node order in case node labels present
	}
	if (is.logical(nodeSupportValues) & length(nodeSupportValues)==1){
		if ( nodeSupportValues ) {
			# stop('Not Implemented: Tree parsing node support values' )
			nodeSupportValues <- tre$node.label
			if ( is.null( nodeSupportValues ) )
				stop('*tre* must have node labels with numeric node support values.')
			nodeSupportValues <- as.numeric( nodeSupportValues )
			if ( all( is.na( nodeSupportValues)))
				stop('Failure to parse tree node labels as node support values. These should be provided as numbers between 0 and 100.')
			nodeSupportValues[ is.na( nodeSupportValues ) ] <- 0
			if ( max( nodeSupportValues ) <= 1 )
				nodeSupportValues <- round( 100 * nodeSupportValues )
			useNodeSupport <- TRUE
			nodeSupportValues <- c( rep(NA, ape::Ntip(tre)), nodeSupportValues )
		}
	}
	if ( is.logical(nodeSupportValues) & length(nodeSupportValues)!=1)
		stop('Failure to parse node support. Must be a single logical or a numeric vector with length equal to number of internal nodes in the tree. If numeric, these values should be between 0 and 100. If logical, node support should be included in labeled nodes of the tree.')

	tredat = .tredat( tre )
	stopifnot( is.null(level) | is.numeric(level) )
	if( !is.null( level ) & is.numeric(level))
		return( .trestruct( tre, minCladeSize, minOverlap , nodeSupportValues , nodeSupportThreshold , nsim , level[1] , ncpu , verbosity , debugLevel
		, useNodeSupport, tredat) )
	if (is.null(level))
	{
		levels <- seq( levellb, levelub, length.out=res )
		tss <- lapply( levels, function(l)
		{
			message( paste('Running treestructure with significance level', round(l,3) ))
			message('')
			.trestruct( tre, minCladeSize, minOverlap , nodeSupportValues , nodeSupportThreshold , nsim , l, ncpu , verbosity , debugLevel, useNodeSupport, tredat)
		})
		chs <- sapply( tss, .ch )
		chdf <- data.frame( level = levels, CH = chs, optimal= '' )
		chdf$optimal[ which.max( chs ) ] <- '***'
		ts <-  tss[[ which.max(chs) ]]
		ts$chdf <- chdf
		return(ts)
	}
}

.trestruct <- function( tre, minCladeSize = 25, minOverlap = -Inf, nodeSupportValues = FALSE, nodeSupportThreshold = 95, nsim = 1e3, level = .01, ncpu = 1, verbosity = 1, debugLevel=0
	, useNodeSupport, tredat)
{
	#attach( tredat )

	.ismonomono <- function( u, v){
		# NOTE if u is direct decendant of v or vice versa than the two tip sets are mono/mono related
		if ( utils::tail(tredat$ancestors[[u]],1) == v){
			return(TRUE)
		}
		if ( utils::tail(tredat$ancestors[[v]],1) == u){
			return(TRUE)
		}
		!((v %in% tredat$ancestors[[u]]) | (u %in% tredat$ancestors[[v]] ))
	}

	# DISSIMILARITY FUNCTIONS
	# do clades u and v overlap in time?
	.au.overlap <- function(a,u) {
		a_range <- range( tredat$nhs[ setdiff( tredat$descendants[[a]], tredat$descendants[[u]] ) ] )
		nu <- sum( (tredat$descInternalHeights[[u]]  > a_range[1]) & (tredat$descInternalHeights[[u]]  <= a_range[2]))
		( nu > minOverlap )
	}
	.overlap <- function( uset, vset){
		uhs <- tredat$nhs[ intersect( uset, tredat$inodes) ]
		vhs <- tredat$nhs[ intersect(vset, tredat$inodes) ]
		nu <- sum( (uhs > min(vhs)) & (uhs < max(vhs)))
		nv <- sum( (vhs > min(uhs)) & (vhs < max(uhs)))
		rv = ( min(nu,nv) > minOverlap )
		rv
	}
	# /DISSIMILARITY FUNCTIONS

	# FIND OUTLIERS
	debugdf = NULL
	zstar  <- stats::qnorm( 1-min(1,level)/2 )
	node2nodeset <- tredat$descendants
	shouldDig <- rep(FALSE, tredat$n + tredat$nnode )
	shouldDig[ tredat$rootnode ] <- TRUE

	# compute z score for u descendened from claderoot
	.calc.z <- function(u, v, Ei = 1, returnabs = TRUE){
		if ( u <= tredat$n )
		  return(0)
		if ( v <= tredat$n )
		  return(0)
		if ( tredat$ndesc[u] < minCladeSize )
		  return(0 )
		if ( tredat$ndesc[v] < minCladeSize )
		  return(0 )
		if ( u==v )
		  return(0)
		if ( is.na( tredat$node2edgelength[ u ] ) )
			return(0)
		if ( tredat$node2edgelength[ u ] == 0 )
			return(0)
		if ( useNodeSupport ){
			if (!is.na( nodeSupportValues[u] ) ){
				if ( nodeSupportValues[u] < nodeSupportThreshold ){
					return(0)
				}
			}
		}
		if (is.na( Ei)){
			#nested = ifelse( .ismonomono( u, v), 1, 0 )
			if ( all( node2nodeset[[u]] %in% node2nodeset[[v]] ))
			{
				nsu <- node2nodeset[[u]]
				nsv <- setdiff( node2nodeset[[v]] , node2nodeset[[u]] )
			} else if ( all( node2nodeset[[v]] %in% node2nodeset[[u]]) ){
				nsv <- node2nodeset[[v]]
				nsu <- setdiff( node2nodeset[[u]] , node2nodeset[[v]] )
			} else{
					nsu <-  setdiff( node2nodeset[[u]], node2nodeset[[v]] )
					nsv <- setdiff( node2nodeset[[v]], node2nodeset[[u]])
			}
			Ei = ifelse( .ismonomono_cl(tredat$tre, tredat$ancestors, nsu, nsv ), 2,1 )
		}
		if (Ei==1) {
			if (!.au.overlap(v, u) ){
				return(0)
			}
		} else {
			if (!.overlap( node2nodeset[[u]], node2nodeset[[v]]) ){
				return(0)
			}
		}
		if (Ei==1){ # u under v
			nsu <-   intersect(  node2nodeset[[u]], node2nodeset[[v]] )
			nsv <- setdiff( node2nodeset[[v]], node2nodeset[[u]])
			vtips <- intersect( nsv, 1:(ape::Ntip(tredat$tre)))
			utips <- intersect( nsu, 1:(ape::Ntip(tredat$tre)))
			if (length( vtips ) < minCladeSize)
			  return(0)
			if (length( utips ) < minCladeSize)
			  return(0)
			# v clade must include tips not in u
			if ( length( vtips ) == 0){
				return( 0 )
			}
			return( .uv.diss( tredat$tre, tredat$nhs, tredat$inodes, tredat$n, tredat$ancestors, nsu , nsv, nsim = nsim , Ei=Ei, returnabs = returnabs ) )
		} else{
			nsu <- setdiff( node2nodeset[[u]], node2nodeset[[v]] )
			nsv <- setdiff( node2nodeset[[v]], node2nodeset[[u]])
			.uv.diss( tredat$tre, tredat$nhs, tredat$inodes, tredat$n, tredat$ancestors, nsu, nsv, nsim = nsim , Ei=Ei , returnabs = returnabs)
		}
	}
	# find biggest outlier descend from a not counting rest of tree
	.find.biggest.outlier <- function(a, Ei = 1 ){
		zs <- if ( ncpu  > 1 ){
			unlist( parallel::mclapply(  node2nodeset[[a]], function(u) .calc.z( u, a, Ei = Ei )
			 , mc.cores = ncpu ) )
		} else{
			sapply( node2nodeset[[a]], function(u)  .calc.z( u, a, Ei = Ei ) )
		}
		wm <- which.max( zs )
		ustar <- node2nodeset[[a]][ wm ]
		if ( zs[wm] <= zstar )
		  return(NULL)
		ustar
	}
	# return ustar,  new exclude[a]
	.dig.clade <- function(a){
		u <- .find.biggest.outlier( a )
		if (is.null(u))
		  return(NULL)
		list( ustar = u
		  , newnodeset_a = setdiff( node2nodeset[[a]], tredat$descendants[[u]] )
		)
	}

	.init.nodeset <- function(u, digging){
		setdiff( tredat$descendants[[u]], do.call( c, lapply( digging, function(v) node2nodeset[[v]] ) ) )
	}

	if ( debugLevel > 0 ){
		# store z for each node compared to root
		debugdf = data.frame( z = sapply( node2nodeset[[tredat$rootnode]], function(u)  .calc.z( u, tredat$rootnode, Ei =1 , returnabs = FALSE) )
		 , cladeSize = sapply( node2nodeset[[tredat$rootnode]], function(u)  tredat$ndesc[u] )
		 , node = node2nodeset[[tredat$rootnode]]
		 )
	}

	node2cl <- c()
	clusterlist <- list()
	while( any(shouldDig) ){
		todig <- which( shouldDig )
		for (a in todig ){
			dc = .dig.clade( a )  # returns only u with highest z
			if ( is.null( dc )){
				shouldDig[a] <- FALSE
				if ( length( node2nodeset[[a]] ) > 0 ){
					#  test  if tips in nodeset before adding to clusterlist
					ans <- node2nodeset[[a]]
					atips <- setdiff( ans , 1:(ape::Ntip(tredat$tre)))
					if ( length( atips ) > 0 ){
						no <- length( clusterlist)
						clusterlist[[ no + 1 ]] <- node2nodeset[[a]]
						node2cl <- c( node2cl, a )
					}
				}
			}else{
				node2nodeset[[dc$ustar]] <- .init.nodeset( dc$ustar, setdiff(union( todig, node2cl ), a)  )
				if ( length(node2nodeset[[dc$ustar]]) > 0)
					shouldDig[ dc$ustar ] <- TRUE
				node2nodeset[[a]] <- dc$newnodeset_a
			}
		}
		if (verbosity > 0 ){
			cat(paste( 'Finding splits under nodes:',   paste(todig, collapse = ' ') , '\n') )
		}
	}
	# /OUTLIERS

	names( clusterlist ) <- node2cl
	no <- length(clusterlist)
	x <- setdiff( c( 1:tredat$n, tredat$inodes), do.call( c, clusterlist )) # remainder
	if (length(x) > 0){
		clusterlist[[no + 1]] <- x
	}



	# DISTANCE
	nc <- length( clusterlist )
	D <- matrix( NA, nrow = nc, ncol = nc )
	diag(D) <- 0
	# 1) fill in D where ancestry/size/overlap req's are met
	if ( nc  > 1){
		for (iu in 1:(nc-1)){
			for (iv in (iu+1):nc){
				# overlap
				if ( .overlap( clusterlist[[iu]] , clusterlist[[iv]] )) #
				{
				  D[iu,iv] <- .uv.diss ( tredat$tre, tredat$nhs, tredat$inodes, tredat$n, tredat$ancestors, clusterlist[[iu]] , clusterlist[[iv]], nsim = nsim,  Ei = 3)
				  D[iv, iu] <- D[iu, iv]
				}
			}
		}
	}
	# 2) impute any missing
	# impute and remaining missingness using weighted mean; weights based on tree distance
	.D <- D
	.D[ is.na(D) | is.infinite(D)] <- 0  #
	rownames(.D)  = colnames(.D) <- 1:nc
	D <- stats::as.dist(.D)

	# /DISTANCE



	# RETURN
	# for each tip:
	stats::setNames( rep(NA, ape::Ntip(tredat$tre)), tredat$tre$tip.label) -> clustervec
	for (k in 1:nc){
		cl <- clusterlist[[k]]
		itip <- intersect( 1:(ape::Ntip(tredat$tre)), cl)
		clustervec[ itip ] <- k
	}

	if ( nc  >  1){
		h <- ape::as.phylo( stats::hclust( D ) )
		#~ 	h$edge.length[ h$edge.length <= zstar ] <- 0
		partinds <- stats::cutree( stats::hclust( stats::as.dist( ape::cophenetic.phylo( h ) ) ), h = zstar)
		partition <- as.factor( stats::setNames( partinds[ clustervec ] , tredat$tre$tip.label ) )
		clustering <- as.factor( clustervec )
		clusters <- split( tredat$tre$tip.label, clustervec )
		partitionSets <- split( tredat$tre$tip.label, partition )

	} else{
		clustering <- stats::setNames(  as.factor( rep(1, tredat$n )), tredat$tre$tip.label )
		partition <- stats::setNames( as.factor( rep(1, tredat$n )), tredat$tre$tip.label )
		clusters <- list( tredat$tre$tip.label )
		partitionSets <- list( tredat$tre$tip.label )
		remainderClade <- NULL
	}

	rv = list(
	  clustering = clustering
	  , partition = partition
	  , clusterSets = clusters
	  , partitionSets = partitionSets
	  , D = D
	  , clusterList = clusterlist
	  , tree = tredat$tre
	  , level = level
	  , zstar = zstar
	  , cluster_mrca  = node2cl
	  , call = match.call()
	  , data =  data.frame( taxon = tre$tip.label
		   # , cluster = clustering
		   , cluster = clustering
		   , partition = partition
	     , row.names = 1:ape::Ntip(tredat$tre)
	     , stringsAsFactors=FALSE
		)
	  , debugdf = debugdf
	)
	class(rv) <- 'TreeStructure'
	#detach( tredat )
	rv
}


.ch <- function(trstr)
{
	cldf <- trstr$data # .computeclusters( tre, node2z, zth, rescale=FALSE )
	ints <- node.depth.edgelength( trstr$tree )[(ape::Ntip(trstr$tree)+1):(ape::Ntip(trstr$tree)+ape::Nnode(trstr$tree))] # internalnode times
	clnts <- lapply( split( cldf, cldf$cluster), function(d){
		tr1 <- ape::keep.tip( trstr$tree, d$taxon )
		ints1 <- node.depth.edgelength( tr1 )[(ape::Ntip(tr1)+1):(ape::Ntip(tr1)+ape::Nnode(tr1))] # internalnode times
	})
	cln <- sapply( clnts , length )
	clmeans <- sapply( clnts, mean )
	oz <- mean( ints )
	bcss <- sum( cln * (clmeans - oz )^2 )
	wcss <- sapply( clnts, function(x) mean( (x-mean(x))^2 ) ) |> sum()
	n <- sum( cln )
	k <- length( clnts )
	(bcss/(k-1)) / (wcss/(n-k))
}

.plot.TreeStructure.ggtree <- function(x, ... ){
	stopifnot(  'ggtree' %in% utils::installed.packages()[,1] )
	stopifnot( inherits( x, 'TreeStructure') )
	tre <- x$tree
	d <- x$data
	d$shape <-  rep('circle', ape::Ntip(tre))

	tre <- ggtree::groupOTU( tre, x$clusterSets)
	pl <- ggtree::`%<+%`( ggtree::ggtree( tre, ggplot2::aes_(color=~group), ... ) ,  d  )
	pl <- pl +  ggtree::geom_tippoint(ggplot2::aes_( color=~partition, shape=~shape, show.legend=TRUE), size = 2 )
	pl
}

#' Plot TreeStructure tree with cluster and partition variables
#' @param x  A TreeStructure object
#' @param use_ggtree Toggle ggtree or ape plotting behaviour
#' @param ... Additional arguments passed to ggtree or ape::plot.phylo
#' @export
#' @examples
#'
#' tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') )
#'
#' (struc <- trestruct( tree ))
#'
#' #plot treestructure object
#'
#' plot(struc)
plot.TreeStructure <- function(x, use_ggtree = TRUE , ... )
{
	stopifnot( inherits( x, 'TreeStructure') )
	if ( use_ggtree ){
		return( .plot.TreeStructure.ggtree (x , ... ) )
	} else{
		tr <- x$tree
		tr$tip.label <- as.character( x$partition  )
		ape::plot.phylo( tr , ... )
		ape::nodelabels( '', node = x$cluster_mrca , pch = 8, cex = 3, frame = 'none')
		ape::tiplabels( '', pch = 15, col   =  as.vector( x$partition ) , frame='none')
	}

	invisible(x)
}


#' @export
print.TreeStructure <- function(x, rows = 0, ...)
{
	stopifnot( inherits( x, 'TreeStructure') )

	npart <- nlevels( x$partition)
	nc <- nlevels( x$clustering )
	tre <- x$tree
	cat( 'Call: \r\n' )
	print( x$call )
	cat('\r\n')
	cat ( paste( 'Significance level:', x$level, '\n' ) )

	cat ( paste( 'Number of clusters:', nc , '\n') )
	cat ( paste( 'Number of partitions:', npart, '\n' ) )

	cat( 'Number of taxa in each cluster:\n' )
	print( table( x$clustering) )
	cat( 'Number of taxa in each partition:\n' )
	print( table( x$partition ))

	if ( rows > 0 ){
		cat ('Cluster and partition assignment: \n')
		print( x$data[1:min(nrow(x$data),rows), ] )
	}

	if (!is.null( x$chdf ))
	{
		cat('\r\n')
		cat( ' CH index ' )
		cat('\r\n')
		print( x$chdf )
	}

	cat( '...\r\n')
	cat( 'For complete data, use `as.data.frame(...)` \n' )

	invisible(x)
}

#' @export
as.data.frame.TreeStructure <- function(x, row.names=NULL, optional=FALSE, ...){
	stopifnot( inherits( x, 'TreeStructure') )
	if ( !is.null(row.names))
		rownames(x$data) <- row.names
	x$data
}

#' @export
coef.TreeStructure <- function(object, ... )
{
	stopifnot( inherits( object, 'TreeStructure') )
	object$partition
}
