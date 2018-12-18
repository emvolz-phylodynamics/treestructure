.get_anc <- function(edge, u ){
	edge[ edge[,2]==u ,1]
}
.get_dgtrs <- function(edge, a){
	edge[ edge[,1] == a , 2 ]
}

.get_rootnode <- function(edge){
	setdiff( as.vector(edge), edge[,2] )
}


#' Detect cryptic population structure in time trees 
#' 
#' @param tre A tree of type ape::phylo. Must be rooted and binary. 
#' @param minCladeSize All clusters within parititon must have at least this many tips. 
#' @param minOverlap Threshold time overlap required to find splits in a clade  
#' @param nsim Number of simulations for computing null distribution of test statistics 
#' @param level Significance level for finding new split within a set of tips 
#' @param ncpu If >1 will compute statistics in parallel using multiple CPUs 
#' @param verbosity If > 0 will print information about progress of the algorithm 
#' @return A TreeStructure object which includes cluster and partitition assignment for each tip of the tree. 
#' @examples 
#' library(ape)
#' tree <- rcoal(50)
#' struct <-  trestruct( tree )
#' @export 
trestruct <- function( tre, minCladeSize = 20, minOverlap = -Inf, nsim = 1e3, level = .01, ncpu = 1, verbosity = 1, ... )
{
	stopifnot( ape::is.rooted(tre))
	stopifnot( ape::is.binary(tre))
	if ( minOverlap >= minCladeSize){
		stop('*minOverlap* should be < *minCladeSize*.')
	}
	if ( any(is.na(tre$tip.label)) | anyDuplicated(tre$tip.label)){
		cat('Tree has NA or duplicated tip labels. Adding a unique id. \n')
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


wD <- ape::dist.nodes( tre )[-(1:n), -(1:n)]
wD <- pmax( wD, mean(upper.tri(wD)) / nnode  )
diag(wD) <- Inf

prenodes <- unique( preedges[,1] )

dgtrMat <- matrix( NA, nrow = ape::Nnode(tre)+ape::Ntip(tre), ncol =2)
for (ie in 1:nrow(poedges))
{
	if (is.na( dgtrMat[ poedges[ie,1] , 1] ) )
	  dgtrMat[ poedges[ie,1], 1] <- poedges[ie,2]
	else 
	  dgtrMat[ poedges[ie,1], 2] <- poedges[ie,2]
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
	# nested 0 if u is nested under v, 1 if u & v are both under root
	.uv.diss <- function(uset, vset, nested = 1){
		# 1 sample u
		# 0 co 
		# -1 sample v 
		
		uinternals <-  intersect( uset, inodes)
		vinternals <- setdiff(  intersect( vset, inodes), uinternals )
		x <- rbind( 
		 cbind(    nhs[ uinternals ], 0 )
		 , cbind(  nhs[ vinternals ], 0 )
		 , cbind(  nhs[ intersect( uset, 1:n) ], 1 )
		 , cbind(  nhs[ intersect( vset, 1:n) ], -1 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		nd <- Cuv_ranksum_nulldist(x, nsim, nested  )
		
		rsuv <- .rank.sum( uinternals, vinternals)
		
		m_nd <- mean(nd)
		sd_nd <- sd(nd)
		abs( (rsuv - m_nd) / sd_nd )
	}
	

# /DISSIMILARITY FUNCTIONS 


# FIND OUTLIERS
zstar  <- qnorm( 1-min(1,level)/2 )
node2nodeset <- descendants 
shouldDig <- rep(FALSE, n + nnode )
shouldDig[ rootnode ] <- TRUE 

# compute z score for u descendened from claderoot 
.calc.z <- function(u, v, nested = 0){
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
	if (nested==0) {
		if (!.au.overlap(v, u) ){
			return(0)
		}
	} else {
		if (!.overlap( node2nodeset[[u]], node2nodeset[[v]]) ){
			return(0)
		}
	}
	if (nested==0){ # u under v
		nsu <-  intersect( node2nodeset[[u]], node2nodeset[[v]] )
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
		.uv.diss( nsu , nsv, nested=nested )
	} else{
if ( length( intersect(   node2nodeset[[v]], node2nodeset[[u]]) ) > 0) stop('intersect error ')
		.uv.diss( node2nodeset[[u]] , node2nodeset[[v]], nested=nested )
	}
}
# find biggest outlier descend from a not counting rest of tree 
.find.biggest.outlier <- function(a ){
	zs <- if ( ncpu  > 1 ){
		unlist( parallel::mclapply(  node2nodeset[[a]], function(u)  .calc.z( u, a ) 
		 , mc.cores = ncpu ) )
	} else{
		sapply( node2nodeset[[a]], function(u)  .calc.z( u, a ) )
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

node2cl <- c()
clusterlist <- list()
while( any(shouldDig) ){
	todig <- which( shouldDig )
	for (a in todig ){
		dc = .dig.clade( a )  # returns only u with highest z  
		if ( is.null( dc )){
			shouldDig[a] <- FALSE
			if ( length( node2nodeset[[a]] ) > 0 ){
				no <- length( clusterlist)
				clusterlist[[ no + 1 ]] <- node2nodeset[[a]] 
				node2cl <- c( node2cl, a )
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
					D[iu,iv] <- .uv.diss ( clusterlist[[iu]] , clusterlist[[iv]], nested = 1)	
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
	D <- as.dist(.D)
# /DISTANCE 



# RETURN
	# for each tip: 
	setNames( rep(NA, ape::Ntip(tre)), tre$tip.label) -> clustervec 
	for (k in 1:nc){
		cl <- clusterlist[[k]]
		itip <- intersect( 1:(ape::Ntip(tre)), cl)
		clustervec[ itip ] <- k 
	}
	
	if ( nc  >  1){
		h <- ape::as.phylo( hclust( D ) ) 
		#~ 	h$edge.length[ h$edge.length <= zstar ] <- 0
		partinds <- cutree( hclust( as.dist( ape::cophenetic.phylo( h ) ) ), h = zstar)
		partition <- as.factor( setNames( partinds[ clustervec ] , tre$tip.label ) )
		clustering <- as.factor( clustervec )
		clusters <- split( tre$tip.label, clustervec )
		partitionSets <- split( tre$tip.label, partition )
		partitionNodes <-  lapply( 1:max(partinds), function(i) as.numeric(names(partinds))[ partinds==i ] )#
		names(partitionNodes) <- 1:max(partinds)
		
	} else{
		clustering <- as.factor( rep(1, n ))
		partition <- as.factor( rep(1, n ))
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
		)
	)
	class(rv) <- 'TreeStructure'
	rv
}


.plot.TreeStructure.ggtree <- function(x, ... ){
	stopifnot( require(ggtree ) )
	stopifnot( inherits( x, 'TreeStructure') )
	tre <- x$tree 
	d <- x$data 
	d$shape <-  rep('circle', ape::Ntip(tre))
	
	tre <- ggtree::groupOTU( tre, s$clusterSets ) 
	pl <- ggtree::ggtree( tre, aes(color=group), ... ) %<+% d + ggtree::geom_tippoint(aes( color=partition, shape=shape, show.legend=TRUE), size = 2 )
	#+ ggplot2::theme(legend.position="right")
	pl
}

#' Plot TreeStructure tree with cluster and partition variables 
#' @param x  A TreeStructure object 
#' @param use_ggtree Toggle ggtree or ape plotting behaviour 
#' @export 
plot.TreeStructure <- function(x, use_ggtree = TRUE , ... )
{
	stopifnot( inherits( x, 'TreeStructure') )
	if ( 'ggtree' %in% installed.packages() & use_ggtree ){
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
print.TreeStructure <- function(x, rows = 10, ...)
{
	stopifnot( inherits( x, 'TreeStructure') )
	
	npart <- nlevels( x$partition)
	nc <- nlevels( x$clustering )
	tre <- x$tree 
	cat( 'Call: \n' )
	print( x$call )
	cat('\n')
	
	cat ( paste( 'Number of clusters:', nc , '\n') )
	cat ( paste( 'Number of partitions:', npart, '\n' ) )
	cat ( paste( 'Significance level:', x$level, '\n' ) )
	cat ('Cluster and partition assignment: \n')
	
	print( x$data[1:min(nrow(x$data),rows), ] )
	cat( '...\n') 
	cat( 'For complete data, use `as.data.frame(...)` \n ' )
	
	invisible(x)
}

#' @export 
as.data.frame.TreeStructure <- function(x){
	stopifnot( inherits( x, 'TreeStructure') )
	x$data
}

#' @export 
coef.TreeStructure <- function(x, ... )
{
	stopifnot( inherits( x, 'TreeStructure') )
	x$partition
}
