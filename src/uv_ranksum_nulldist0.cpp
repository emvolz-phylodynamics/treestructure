#include <Rcpp.h>

using namespace Rcpp;
//~ using namespace std;

/*
 * x=
 *   0 coalescent
 *   1 sample u
 *   -1 sample v
 */
//[[Rcpp::export]]
NumericVector Cuv_ranksum_nulldist( NumericVector x, int nsim, int Ei){
	NumericVector ranksums(nsim);
	int Aui, Avi;
	double pui, pvi;
	int cocounter;
	double Rikm1;
	double Rim1k;
	for (int k = 0; k < nsim;k++){
		ranksums(k) = 0;
		Aui = 0;
		Avi = 0;
		cocounter = 0;
		for (int i = 0; i < x.size(); i++){
			if (x(i)==1) {
				Aui++;
			} else if (x(i)==-1){
				Avi++;
			} else if (x(i)==0){
				cocounter++;
				if (Ei==2){
					pui = (double)(Aui - 1.) / (Aui + Avi - 2.);
				} else if (Ei==1) {
					pui = (double)(Aui + 1.) / (Aui + Avi );
				} else if (Ei==3){
					// exact E3 transition (Volz et al 2020, Eq 10-11); the previous
					// version computed Rim1k and Rikm1 identically as R_{z,w}, which
					// cancels and collapses this probability to the Ei==2 form.
					double z = Aui, w = Avi ;
					if ( z <= 1. ){
						pui = 0. ;
					} else if ( w <= 1. ){
						pui = 1. ;
					} else {
						Rim1k = 1./((z-1.)*z) + 1./(w*(w+1.)) - 1./((z+w-1.)*(z+w-2.)) ;
						Rikm1 = 1./(z*(z+1.)) + 1./((w-1.)*w) - 1./((z+w-1.)*(z+w-2.)) ;
						pui = (z-1.) * Rim1k / ( (z-1.)*Rim1k + (w-1.)*Rikm1 ) ;
					}
				}
				if ( (Aui > 1 ) && (Rf_runif(0,1) < pui) ){
					Aui--;
					ranksums(k) += cocounter;
				} else{
					Avi--;
				}
			}
		}
	}

	return ranksums;
}


