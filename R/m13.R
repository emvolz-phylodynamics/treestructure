invisible('
- an option allow_nested which would prevent ident of nested clades 
- rather than classify outliers as internal nodes, classify them as sets of tips
	 this would allow nesting of mrcas, but tip sets would be disjoint 
') 

.get_anc <- function(edge, u ){
	edge[ edge[,2]==u ,1]
}
.get_dgtrs <- function(edge, a){
	edge[ edge[,1] == a , 2 ]
}

.get_rootnode <- function(edge){
	setdiff( as.vector(edge), edge[,2] )
}


#' Partition rooted tree into sets of clades
#' 
#' @param tre A tree of type ape::phylo. Must be rooted and binary. 
#' @param minCladeSize All clusters within parititon must have at least this many tips. 
#' @param minOverlap 
#' @param nsim Number of simulations for comparing clusters within partition
#' @param level Significance level for splitting partition to make new cluster 
#' @return A TreeStructure object which includes cluster and partitition assignment for each tip of the tree. 
#' @examples 
#' library(ape)
#' tree <- rcoal(100)
#' struct <-  trestruct( tree )
#' @export 
trestruct13 <- function( tre, minCladeSize = 20, minOverlap = -Inf, nsim = 1e3, level = .01,  ... )
{
	stopifnot( ape::is.rooted(tre))
	stopifnot( ape::is.binary(tre))
	if ( minOverlap >= minCladeSize){
		stop('*minOverlap* should be < *minCladeSize*.')
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
	#' @param nested 0 if u is nested under v, 1 if u & v are both under root
	.uv.diss <- function(uset,vset, nested = 1){
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

#' compute z score for u descendened from claderoot 
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
	if (nested==0){
		nsv <- setdiff( node2nodeset[[v]], node2nodeset[[u]])
		vtips <- intersect( nsv, 1:ape::Ntip(tre))
		if (length( vtips ) < minCladeSize)
		  return(0)
		# v clade must include tips not in u
		if ( length( vtips ) == 0){
			return( 0 )
		}
		.uv.diss( node2nodeset[[u]] , nsv, nested=nested )
	} else{
if ( length( intersect(   node2nodeset[[v]], node2nodeset[[u]]) ) > 0) stop('intersect error ')
		.uv.diss( node2nodeset[[u]] , node2nodeset[[v]], nested=nested )
	}
}
#' find biggest outlier descend from a not counting rest of tree 
.find.biggest.outlier <- function(a ){
	zs <- sapply( node2nodeset[[a]], function(u) 
	  .calc.z( u, a )
	)
	wm <- which.max( zs )
	ustar <- node2nodeset[[a]][ wm ]
	if ( zs[wm] <= zstar )
	  return(NULL)
	ustar 
}
#' return ustar,  new exclude[a] 
.dig.clade <- function(a){
	u <- .find.biggest.outlier( a )
	if (is.null(u)) 
	  return(NULL)
#~ print(a)
#~ print(u)
#~ browser()
	list( ustar = u
	  , newnodeset_a = setdiff( node2nodeset[[a]], descendants[[u]] )
	)
}





clusterlist <- list()
done <- FALSE
while( any(shouldDig) ){
	todig <- which( shouldDig )
	for (a in todig ){
		dc = .dig.clade( a )  # returns only u with highest z  
		
		if ( is.null( dc )){
			shouldDig[a] <- FALSE
			if ( length( node2nodeset[[a]] ) > 0 ){
				no <- length( clusterlist)
				clusterlist[[ no + 1 ]] <- node2nodeset[[a]] 
			}
		}else{
			shouldDig[ dc$ustar ] <- TRUE 
			node2nodeset[[a]] <- dc$newnodeset_a
		}
	}
cat('todig\n') 
print(todig)
}
# /OUTLIERS

	
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


#~ browser()


# RETURN
	# for each tip: 
	setNames( rep(NA, ape::Ntip(tre)), tre$tip.label) -> clustervec 
	for (k in 1:nc){
		cl <- clusterlist[[k]]
		itip <- intersect( 1:(ape::Ntip(tre)), cl)
		clustervec[ itip ] <- k 
	}
	
#~ browser()
# TODO nesting within clustervec! should be disjoint 
	
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
	)
	class(rv) <- 'TreeStructure'
	rv
}




summary.TreeStructure <- function(x, ... )
{
	invisible(x)
}

print.TreeStructure <- function(x, ...)
{
	invisible(x)
}

plot.TreeStructure <- function(x, ... )
{
	invisible(x)
}

coef.TreeStructure <- function(x, ... )
{
	invisible(x)
}

partition <- function(x, ... )
{
	invisible(x)
}

cluster <- function(x, ... )
{
	invisible(x)
}

distance <- function(x, ... )
{
	invisible(x)
}

