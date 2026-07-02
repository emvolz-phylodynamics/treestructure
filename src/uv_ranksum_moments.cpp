#include <Rcpp.h>

using namespace Rcpp;

/*
 * Exact mean and standard deviation of the rank-sum null distribution, computed by a
 * deterministic recursion rather than by Monte-Carlo simulation (cf Cuv_ranksum_nulldist).
 * The pooled number of extant lineages A = Au + Av is deterministic along the event
 * sequence, so the null chain is one-dimensional in Au. We propagate the full
 * distribution of Au together with the first two moments of the accumulating rank sum,
 * so a single sweep returns E[R] and sd[R].
 * x=
 *   0 coalescent
 *   1 sample u
 *   -1 sample v
 * Ei : event type (as in Cuv_ranksum_nulldist). Ei==3 uses the EXACT transition of
 *      Equations 10-11 (Volz et al 2020), not the approximation which collapses to Ei==2.
 */

// R_{z,w} of Equation 11, defined for z,w >= 1
double Ruv_R( double z, double w ){
	double s = z + w ;
	return 1./(z*(z+1.)) + 1./(w*(w+1.)) - 1./(s*(s-1.)) ;
}

// probability that the next coalescent is within clade u given Au=z, Av=w
double Ruv_pu( double z, double w, int Ei ){
	if ( z <= 1. )                    // u needs >= 2 lineages to coalesce
		return 0. ;
	double p ;
	if ( Ei==1 ){
		p = (z + 1.) / (z + w) ;
	} else if ( Ei==2 ){
		p = (z - 1.) / (z + w - 2.) ;
	} else {                          // Ei==3, exact
		if ( w <= 1. )                // v cannot coalesce internally
			return 1. ;
		double num = (z - 1.) * Ruv_R( z - 1., w ) ;
		double den = num + (w - 1.) * Ruv_R( z, w - 1. ) ;
		p = num / den ;
	}
	if ( p < 0. ) p = 0. ;
	if ( p > 1. ) p = 1. ;
	return p ;
}

//[[Rcpp::export]]
NumericVector Cuv_ranksum_moments( NumericVector x, int Ei ){
	int U = 0 ;
	for (int i = 0; i < x.size(); i++)
		if ( x(i)==1 ) U++ ;

	// P[a] = Pr(Au=a) ; S[a] = E[R 1{Au=a}] ; Tm[a] = E[R^2 1{Au=a}]
	std::vector<double> P(U+1, 0.), S(U+1, 0.), Tm(U+1, 0.) ;
	P[0] = 1. ;
	double A = 0. ;
	int cocounter = 0 ;

	for (int i = 0; i < x.size(); i++){
		if ( x(i)==1 ){               // u sample: Au -> Au+1
			for (int a = U; a >= 1; a--){ P[a]=P[a-1]; S[a]=S[a-1]; Tm[a]=Tm[a-1]; }
			P[0]=0.; S[0]=0.; Tm[0]=0.;
			A += 1. ;
		} else if ( x(i)==-1 ){       // v sample
			A += 1. ;
		} else {                      // coalescent, pooled count = A
			cocounter++ ;
			std::vector<double> Pn(U+1, 0.), Sn(U+1, 0.), Tn(U+1, 0.) ;
			for (int a = 0; a <= U; a++){
				if ( P[a]==0. && S[a]==0. && Tm[a]==0. )
					continue ;
				double pui = Ruv_pu( (double)a, A - a, Ei ) ;
				double pvi = 1. - pui ;
				Pn[a] += P[a] * pvi ;      // v coalescent: Au unchanged
				Sn[a] += S[a] * pvi ;
				Tn[a] += Tm[a] * pvi ;
				if ( a >= 1 && pui > 0. ){ // u coalescent: Au -> Au-1, rank sum += cocounter
					double pm = P[a]*pui, sm = S[a]*pui, tm = Tm[a]*pui ;
					Pn[a-1] += pm ;
					Sn[a-1] += sm + cocounter * pm ;
					Tn[a-1] += tm + 2.*cocounter*sm + ((double)cocounter)*cocounter*pm ;
				}
			}
			P.swap(Pn); S.swap(Sn); Tm.swap(Tn);
			A -= 1. ;
		}
	}

	double ER = 0., ER2 = 0. ;
	for (int a = 0; a <= U; a++){ ER += S[a]; ER2 += Tm[a]; }
	double v = ER2 - ER*ER ;
	if ( v < 0. ) v = 0. ;
	return NumericVector::create( ER, sqrt(v) ) ;
}
