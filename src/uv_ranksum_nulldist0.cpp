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
	int A;  
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
					A = Aui + Avi ;
					Rim1k = 1./(Aui*(Aui+1.)) + 1./(Avi*(Avi+1)) - 1./(A*(A-1.)) ;
					Rikm1 = 1./(Aui*(Aui+1.)) + 1./(Avi*(Avi+1.)) - 1. / (A*(A-1.) );
					pui = (Aui-1.) * Rim1k / ( (Aui-1.)*Rim1k + (Avi-1.)*Rikm1 ) ;
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


