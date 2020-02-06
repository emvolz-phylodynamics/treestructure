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

#' Detect cryptic population structure in time trees 
#'
#' @details 
#' Estimates a partition of a time-scaled tree by contrasting coalescent patterns. 
#' The algorithm is premised on a Kingman coalescent null hypothesis and a test statistic is formulated based on the rank sum of node times in the tree. 
#' 
#' @section References:
#' E.M. Volz, Wiuf, C., Grad, Y., Frost, S., Dennis, A., Didelot, X.D.  (2020) Identification of hidden population structure in time-scaled phylogenies.
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#' 
#' @param tre A tree of type ape::phylo. Must be rooted and binary. 
#' @param minCladeSize All clusters within parititon must have at least this many tips. 
#' @param minOverlap Threshold time overlap required to find splits in a clade  
#' @param nsim Number of simulations for computing null distribution of test statistics 
#' @param level Significance level for finding new split within a set of tips 
#' @param ncpu If >1 will compute statistics in parallel using multiple CPUs 
#' @param verbosity If > 0 will print information about progress of the algorithm 
#' @param debugLevel If > 0 will produce additional data in return value 
#' @return A TreeStructure object which includes cluster and partitition assignment for each tip of the tree. 
#' @examples 
#' tree <- ape::rcoal(50)
#' struct <-  trestruct( tree )
#' 
#' @export 
trestruct <- function( tre, minCladeSize = 25, minOverlap = -Inf, nsim = 1e3, level = .01, ncpu = 1, verbosity = 1, debugLevel=0 )
{
	stopifnot( ape::is.rooted(tre))
	stopifnot( ape::is.binary(tre))
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

.ismonomono <- function( u, v){
	# NOTE if u is direct decendant of v or vice versa than the two tip sets are mono/mono related 
	if ( utils::tail(ancestors[[u]],1) == v){
		return(TRUE)
	}
	if ( utils::tail(ancestors[[v]],1) == u){
		return(TRUE)
	}
	!((v %in% ancestors[[u]]) | (u %in% ancestors[[v]] ))
}

.ispara <- function(u, v){
	(v %in% ancestors[[u]]) | (u %in% ancestors[[v]])
}

.ismonomono_cl <- function( cl1, cl2){
	t1 <- intersect( 1:(ape::Ntip(tre)), cl1 )
	t2 <- intersect( 1:(ape::Ntip(tre)), cl2 )
	u <- ape::getMRCA(tre,  t1 )
	v <- ape::getMRCA(tre,  t2 )
	if ( is.null( u ) | is.null(v) ) 
	  return(FALSE)
	!.ispara( u, v)
}


# /PRECOMPUTE


# DISSIMILARITY FUNCTIONS
	# do clades u and v overlap in time? 
	.au.overlap <- function(a,u) {
		a_range <- range( nhs[ setdiff( descendants[[a]], descendants[[u]] ) ] )
		nu <- sum( (descInternalHeights[[u]]  > a_range[1]) & (descInternalHeights[[u]]  <= a_range[2]))
		( nu > minOverlap )
	}
	.overlap <- function( uset, vset){
		uhs <- nhs[ intersect( uset, inodes) ]
		vhs <- nhs[ intersect(vset, inodes) ]
		nu <- sum( (uhs > min(vhs)) & (uhs < max(vhs)))
		nv <- sum( (vhs > min(uhs)) & (vhs < max(uhs)))
		rv = ( min(nu,nv) > minOverlap )
		rv
	}
	
	# test statistic comparing two internals  
	.rank.sum <- function( uinternals, vinternals ) {
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
	.uv.diss <- function(uset, vset, Ei = NA, returnabs=TRUE){
		# 1 sample u
		# 0 co 
		# -1 sample v 
		
		if (is.na( Ei)){
			#nested = ifelse( .ismonomono( u, v), 1, 0 )
			Ei = ifelse( .ismonomono_cl( uset, vset ), 2, 3 )
		}
		
		uinternals <-  intersect( uset, inodes)
		vinternals <- setdiff(  intersect( vset, inodes), uinternals )
		utips <- intersect( uset, 1:n)
		vtips <- intersect( vset, 1:n)
		if ( (length( utips )==0) | (length(vtips) == 0) )
		 return( 0 )
		x <- rbind( 
		 cbind(    nhs[ uinternals ], 0 )
		 , cbind(  nhs[ vinternals ], 0 )
		 , cbind(  nhs[ utips ], 1 )
		 , cbind(  nhs[ vtips ], -1 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		nd <- Cuv_ranksum_nulldist(x, nsim, Ei  )
		
		rsuv <- .rank.sum( uinternals, vinternals)
		
		m_nd <- mean(nd)
		sd_nd <- stats::sd(nd)
		if ( returnabs )
			return( abs( (rsuv - m_nd) / sd_nd ) )
		else
			return( (rsuv - m_nd) / sd_nd )
		
#~ 		s <- sum(rsuv > nd )
#~ 		p <- min( s / nsim , (nsim - s) / nsim )
		
#~ 		p <- sum(rsuv > nd ) / nsim 
#~ 		abs( qnorm( p )  )
	}
	

# /DISSIMILARITY FUNCTIONS 


# FIND OUTLIERS
zstar  <- stats::qnorm( 1-min(1,level)/2 )
node2nodeset <- descendants 
shouldDig <- rep(FALSE, n + nnode )
shouldDig[ rootnode ] <- TRUE 

# compute z score for u descendened from claderoot 
.calc.z <- function(u, v, Ei = 1, returnabs = TRUE){ 
	if ( u <= n ) 
	  return(0)
	if ( v <= n ) 
	  return(0)
	if ( ndesc[u] < minCladeSize )
	  return(0 )
	if ( ndesc[v] < minCladeSize )
	  return(0 )
	if ( u==v )
	  return(0)
	
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
		Ei = ifelse( .ismonomono_cl(nsu, nsv ), 2,1 )
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
		vtips <- intersect( nsv, 1:(ape::Ntip(tre)))
		utips <- intersect( nsu, 1:(ape::Ntip(tre)))
		if (length( vtips ) < minCladeSize)
		  return(0)
		if (length( utips ) < minCladeSize)
		  return(0)
		# v clade must include tips not in u
		if ( length( vtips ) == 0){
			return( 0 )
		}
		return( .uv.diss( nsu , nsv, Ei=Ei, returnabs = returnabs ) ) 
	} else{
		nsu <- setdiff( node2nodeset[[u]], node2nodeset[[v]] )
		nsv <- setdiff( node2nodeset[[v]], node2nodeset[[u]])
		.uv.diss( nsu, nsv , Ei=Ei , returnabs = returnabs)
	}
}
# find biggest outlier descend from a not counting rest of tree 
.find.biggest.outlier <- function(a, Ei = 1 ){
	zs <- if ( ncpu  > 1 ){
		unlist( parallel::mclapply(  node2nodeset[[a]], function(u) .calc.z( u, a,Ei = Ei ) 
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
	  , newnodeset_a = setdiff( node2nodeset[[a]], descendants[[u]] )
	)
}

.init.nodeset <- function(u, digging){
	setdiff( descendants[[u]], do.call( c, lapply( digging, function(v) node2nodeset[[v]] ) ) )
}

debugdf = NULL 
if ( debugLevel > 0 ){
	# store z for each node compared to root 
	debugdf = data.frame( z = sapply( node2nodeset[[rootnode]], function(u)  .calc.z( u, rootnode, Ei =1 , returnabs = FALSE) )
	 , cladeSize = sapply( node2nodeset[[rootnode]], function(u)  ndesc[u] )
	 , node = node2nodeset[[rootnode]]
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
				atips <- setdiff( ans , 1:(ape::Ntip(tre)))
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
	x <- setdiff( c( 1:n, inodes), do.call( c, clusterlist )) # remainder 
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
					D[iu,iv] <- .uv.diss ( clusterlist[[iu]] , clusterlist[[iv]], Ei = 3)	
					#D[iu,iv] <- .uv.diss ( clusterlist[[iu]] , clusterlist[[iv]], Ei = 2)	
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
	stats::setNames( rep(NA, ape::Ntip(tre)), tre$tip.label) -> clustervec 
	for (k in 1:nc){
		cl <- clusterlist[[k]]
		itip <- intersect( 1:(ape::Ntip(tre)), cl)
		clustervec[ itip ] <- k 
	}
	
	if ( nc  >  1){
		h <- ape::as.phylo( stats::hclust( D ) ) 
		#~ 	h$edge.length[ h$edge.length <= zstar ] <- 0
		partinds <- stats::cutree( stats::hclust( stats::as.dist( ape::cophenetic.phylo( h ) ) ), h = zstar)
		partition <- as.factor( stats::setNames( partinds[ clustervec ] , tre$tip.label ) )
		clustering <- as.factor( clustervec )
		clusters <- split( tre$tip.label, clustervec )
		partitionSets <- split( tre$tip.label, partition )
		partitionNodes <-  lapply( 1:max(partinds), function(i) as.numeric(names(partinds))[ partinds==i ] )#
		names(partitionNodes) <- 1:max(partinds)
		
	} else{
		clustering <- stats::setNames(  as.factor( rep(1, n )), tre$tip.label )
		partition <- stats::setNames( as.factor( rep(1, n )), tre$tip.label )
		clusters <- list( tre$tip.label )
		partitionSets <- list( tre$tip.label )
		remainderClade <- NULL 
		partitionNodes <- NULL 
	}
	
	rv = list(
	  clustering = clustering 
	  , partition = partition 
	  , clusterSets  = clusters
	  , partitionSets = partitionSets
	  , partitionNodes =  partitionNodes
	  , D = D 
	  , clusterList = clusterlist 
	  , tree = tre
	  , level = level
	  , zstar = zstar 
	  , cluster_mrca  = node2cl 
	  , call = match.call() 
	  , data =  data.frame( taxon = tre$tip.label
		   , cluster = clustering
		   , partition = partition
	     , row.names = 1:ape::Ntip(tre)
	     , stringsAsFactors=FALSE
		)
	  , debugdf = debugdf 
	)
	class(rv) <- 'TreeStructure'
	rv
}


.plot.TreeStructure.ggtree <- function(x, ... ){
	stopifnot(  'ggtree' %in% utils::installed.packages()[,1] ) 
	stopifnot( inherits( x, 'TreeStructure') )
	tre <- x$tree 
	d <- x$data 
	d$shape <-  rep('circle', ape::Ntip(tre))
	
	tre <- ggtree::groupOTU( tre, x$clusterSets ) 
	pl <- ggtree::`%<+%`( ggtree::ggtree( tre, ggplot2::aes_(color=~group), ... ) ,  d  )
	pl <- pl +  ggtree::geom_tippoint(ggplot2::aes_( color=~partition, shape=~shape, show.legend=TRUE), size = 2 )	
	pl
}

#' Plot TreeStructure tree with cluster and partition variables 
#' @param x  A TreeStructure object 
#' @param use_ggtree Toggle ggtree or ape plotting behaviour 
#' @param ... Additional arguments passed to ggtree or ape::plot.phylo
#' @export 
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
	cat( 'Call: \n' )
	print( x$call )
	cat('\n')
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
	
	cat( '...\n') 
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
