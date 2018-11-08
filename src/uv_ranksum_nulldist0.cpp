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
NumericVector Cuv_ranksum_nulldist( NumericVector x, int nsim, int monomono){
	NumericVector ranksums(nsim);
	int Aui, Avi; 
	double pui, pvi; 
	int cocounter; 
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
				if (monomono==1){
					pui = (double)(Aui - 1.) / (Aui + Avi - 2.); 
				} else{
					pui = (double)(Aui + 1.) / (Aui + Avi ); 
				}
//~ if (k==0)
  //~ cout << Aui << " " << Avi << " " << pui  << endl; 
				if (  Rf_runif(0,1) < pui ){
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


