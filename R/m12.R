#~ library(ape)
#~ library(Rcpp)
#~ sourceCpp( 'uv_ranksum_nulldist0.cpp')

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
trestruct <- function( tre, minCladeSize = 20, minOverlap = -Inf, nsim = 1e3, level = .001,  ... )
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
	.uv.overlap <- function(u,v) {
		nu <- sum( (descInternalHeights[[u]]  > timerange[v,1]) & (descInternalHeights[[u]]  <= timerange[v,2]))
		nv <- sum( (descInternalHeights[[v]]  > timerange[u,1]) & (descInternalHeights[[v]]  <= timerange[u,2]))
		rv = ( min(nu,nv) > minOverlap )
		rv
	}
	.au.overlap <- function(a,u) {
		a_range <- range( nhs[ setdiff( descendants[[a]], descendants[[u]] ) ] )
		nu <- sum( (descInternalHeights[[u]]  > a_range[1]) & (descInternalHeights[[u]]  <= a_range[2]))
		( nu > minOverlap )
	}
	
	# test statistic comparing two internals  
	.rank.sum <- function( u, v ) {
		x <- rbind( 
		 cbind(  nhs[descInternal[[u]] ], 1 )
		 , cbind(  nhs[descInternal[[v]] ], 0 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		sum( x * (1:length(x)) )
	}
	
	# test stat null distribution 
	.uv.diss <- function(u,v){
		# 1 sample u
		# 0 co 
		# -1 sample v 
		x <- rbind( 
		 cbind(  nhs[descInternal[[u]] ], 0 )
		 , cbind(  nhs[descInternal[[v]] ], 0 )
		 , cbind(  nhs[descendantTips[[u]] ], 1 )
		 , cbind(  nhs[descendantTips[[v]] ], -1 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		nd <- Cuv_ranksum_nulldist(x, nsim, 1  )
		
		rsuv <- .rank.sum( u, v)
		#d <- sum( rsuv >= nd ) / nsim 
		#min( d, 1 - d)
		
		m_nd <- mean(nd)
		sd_nd <- sd(nd)
		abs( (rsuv - m_nd) / sd_nd )
	}
	
	# test statistic comparing internal with remaindrer 
	.rank.sum.u <- function( u , exclude = c() ){
		diroot <- setdiff( descInternal[[rootnode]], c(exclude, descInternal[[u]] ) )
		x <- rbind( 
		 cbind(  nhs[descInternal[[u]] ], 1 )
		 , cbind(  nhs[ diroot ], 0 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		sum( x * (1:length(x)) )
	}
	
	# null-with-remainder null distribution 
	# TODO 484/514 too many false positives? 
	.uv.diss.u <- function(u, exclude = c() ){
		# 1 sample u
		# 0 co 
		# -1 sample root 
		diroot <- setdiff( descInternal[[rootnode]], c(exclude, descInternal[[u]] ) )
		dtroot <- setdiff( descendantTips[[rootnode]], c( exclude, descendantTips[[u]] ) )
		if ( length( diroot) == 0 | length(dtroot)==0 )
		  return(0)
		x <- rbind( 
		 cbind(  nhs[descInternal[[u]] ], 0 )
		 , cbind(  nhs[ diroot ], 0 )
		 , cbind(  nhs[descendantTips[[u]] ], 1 )
		 , cbind(  nhs[ dtroot ], -1 )
		)
		x <- as.vector( x[ order(x[,1] ) , 2] )
		nd <- Cuv_ranksum_nulldist(x, nsim , 0 )
		
		rsuv <- .rank.sum.u( u, exclude = exclude )
		
		m_nd <- mean(nd)
		sd_nd <- sd(nd)
		
		abs( (rsuv - m_nd) / sd_nd )
	}
	

# /DISSIMILARITY FUNCTIONS 


# FIND OUTLIERS
zstar  <- qnorm( 1-min(1,level)/2 )
node2exclude <- lapply( 1:(n+nnode), function(u) setdiff( 1:(n+nnode),  descendants[[u]] ) ) # exclude all but descendant of u 
shouldDig <- rep(FALSE, n + nnode )
shouldDig[ rootnode ] <- TRUE 
#' compute z score for u descendened from claderoot 
.calc.z <- function(u, claderoot, z_exclude){
	if ( u <= n ) 
	  return(0)
	if ( u %in% z_exclude )
	  return(0)
	if ( ndesc[u] < minCladeSize )
	  return(0 )
	if ( (ndesc[claderoot] - ndesc[u] ) < minCladeSize )
	  return(0)
	if ( u==claderoot )
	  return(0)
	if (!.au.overlap(claderoot, u) ){
		return(0)
	}
	.uv.diss.u( u , exclude = z_exclude )
}
#' find biggest outlier descend from a not counting rest of tree 
.find.biggest.outlier <- function(a, z_exclude ){
	zs <- sapply( descendants[[a]], function(u) 
	  .calc.z( u, a, z_exclude )
	)
	wm <- which.max( zs )
	ustar <- descendants[[a]][ wm ]
	if ( zs[wm] <= zstar )
	  return(NULL)
#~ browser()
#~ .zs <- ( sapply( descendants[[a]], function(u) 
#~ 	  .calc.z( u, a, z_exclude )
#~ 	))
#~ o = descendants[[a]][ zs > zstar ]
	ustar 
}
#' return ustar,  new exclude[a] 
.dig.clade <- function(a){
	exclude <- node2exclude[[a]]
	u <- .find.biggest.outlier( a, exclude )
	if (is.null(u)) 
	  return(NULL)
	exclude <- unique( c( exclude, descendants[[u]] , tail(ancestors[[u]],1) )) #
	list( ustar = u, newexclude = exclude )
}


#' find biggest outlier(s) descended from a not counting rest of tree 
.find.biggest.outlier2 <- function(a, z_exclude ){
	zs <- sapply( descendants[[a]], function(u) 
	  .calc.z( u, a, z_exclude )
	)
	wm <- which( zs > zstar )
	if (length( wm)==0) 
	  return(NULL)
	ustars <- descendants[[a]][ wm ]
#~ e = new.env(); load('test.rda', envir=e)
#~ browser()
#~ do.call( save, c( as.list( c(ls(), ls(environment(.dig.clade), all.names=TRUE) ) ) , file = 'test.rda') ) 
	ustars
}
#' return nodes {u},  new exclude[a] 
.dig.clade2 <- function(a) {
	exclude <- node2exclude[[a]]
	us <- .find.biggest.outlier2( a, exclude )
	if (is.null(us)) 
	  return(NULL)
#~ browser()
	exclude <- c( exclude
	  , do.call( c
	    , lapply( us, function(u) c(descendants[[u]] , tail(ancestors[[u]],1) )) ))
	list( ustar = us, newexclude = exclude )
}




outliers <- c() 
done <- FALSE
while( any(shouldDig) ){
	todig <- which( shouldDig )
	for (a in todig ){
#~ 		dc = .dig.clade( a ) #TODO 
		dc = .dig.clade2( a )
		if ( is.null( dc )){
			shouldDig[a] <- FALSE
		}else{
			shouldDig[ dc$ustar ] <- TRUE 
			node2exclude[[a]] <- dc$newexclude 
			outliers <- c( outliers, dc$ustar )
		}
	}
cat('todig\n') 
print(todig)
}
# /OUTLIERS

	outliers <- unique( outliers ) #TODO should be unnecessary 
	clades <- lapply( outliers, function(u) ape::extract.clade( tre,u))
	names(clades) <- outliers 
		
	t2d <- do.call( c, lapply( clades, '[[', 'tip.label' ) )
	t2k <- setdiff( tre$tip.label, t2d  )
	#~ 	if ( length( t2k ) < minCladeSize ){
	#~ 		remainderClade = NULL 
	#~ 		remainderFit = ' ** Remainder of tree after removing tree splits is less than minimum clade size. ** ' 
	#~ 	} else 
	{
		remainderClade = ape::drop.tip ( tre, t2d )
		clades[[ as.character(rootnode) ]] <- remainderClade
	}
	
# DISTANCE
	nc = length( clades )
	D <- matrix( NA, nrow = nc, ncol = nc )
	diag(D) <- 0
	# 1) fill in D where ancestry/size/overlap req's are met 
	if ( nc  > 1){
		for (iu in 1:(nc-1)){
			u <- as.numeric( names(clades)[iu] )
			for (iv in (iu+1):nc){
				v <- as.numeric( names(clades)[iv] )
				# overlap
				if ( .uv.overlap(u,v)) # 
				{
					if ( u %in% descendants[[v]] ){
						D[iu,iv] <- .uv.diss.u( u , exclude = c( descendants[[u]], setdiff(1:(n+nnode), descendants[[v]]) ) )
					} else if(v %in% descendants[[u]]){
						D[iu,iv] <- .uv.diss.u( v , exclude = c( descendants[[v]], setdiff(1:(n+nnode), descendants[[u]]) ) )
					} else{
						D[iu, iv] <- .uv.diss(u,v) # TODO should treat remainder differently 
					}
					D[iv, iu] <- D[iu, iv] 
				}
			}
		} 
	}
	# 2) impute any missing 
	# impute and remaining missingness using weighted mean; weights based on tree distance 
	.D <- D 
	.D[ is.na(D) | is.infinite(D)] <- 0  #TODO shouldnt need to filter inf 
	rownames(.D)  = colnames(.D) <- names(clades)
	D <- as.dist(.D)
# /DISTANCE 


# RETURN
	clustering <- sapply( tre$tip.label, function(tip){
		x = c()
		k <- 0
		#nt <- -Inf 
		nt <- Inf 
		for ( clade in clades){
			k <- k + 1
			#if ( (tip %in% clade$tip.label) & (Ntip(clade) > nt)){ #TODO or maybe smaller level
			if ( (tip %in% clade$tip.label) & (ape::Ntip(clade) < nt)){ 
				x <- k 
				nt <- ape::Ntip(clade)
			}
		}
		x
	})
	if ( nc  >  1){
		h <- ape::as.phylo( hclust( D ) ) 
		#~ 	h$edge.length[ h$edge.length <= zstar ] <- 0
		partinds <- cutree( hclust( as.dist( ape::cophenetic.phylo( h ) ) ), h = zstar)
		partition <- as.factor( setNames( partinds[ clustering ] , tre$tip.label ) )
		clustering <- as.factor( clustering )
		clusters <- split( tre$tip.label, clustering )
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
#~ browser()
	
	rv = list(
	  clades = clades
	  , clustering = clustering 
	  , partition = partition 
	  , clusterSets  = clusters
	  , partitionSets = partitionSets
	  , partitionNodes =  partitionNodes
	  , D = D 
	  , outliers = outliers 
	  , remainderClade = remainderClade
	  , tree = tre
	)
	class(rv) <- 'TreeStructure'

#~ X11(); plot( tre ); nodelabels( outliers, outliers ); 
#~ X11(); plot( as.phylo( hclust( D )) )
#~ browser()
# structure( list(), class = 'TreeStructure' )
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

