// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; 
using namespace arma;


// [[Rcpp::export]]
Rcpp:: NumericVector Sample_vec(const  arma::vec & points, const  arma::vec & upper, const  arma::vec &lower, const  arma::vec &proba, const  arma::vec &val) {
    int N=points.size();
    Rcpp:: NumericVector J(N);
   

    for(int i = 0; i < N; i++) {
    	double s=0;
    	int j1=lower(i)-1;
    	
    	for(int j=j1; j<upper(i);j++){
    		s+=proba(j);
    	} 
    	double v=0.0;
    	int j=lower(i)-2;
    	while(s*points(i)>v){
    		j++;
    		v+=proba(j);
    	}
    	J(i)=val(j);
    	
    }
    return J;
}






















