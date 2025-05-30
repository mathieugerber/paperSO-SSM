// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; 
using namespace arma;


// [[Rcpp::export]]
Rcpp::List KF_fast(const arma::mat& mu, const arma::mat& Sigma, const arma::mat& theta_F, const arma::mat& theta_Q, const arma::vec& theta_R, const arma::mat& H, const arma::vec& y) {

  int N = Sigma.n_rows, s = mu.n_cols;
  int s2 = std::pow(s, 2);
  int count;

  arma::mat P(s, s);
  arma::mat F = arma::eye<arma::mat>(s, s);
  arma::mat I = arma::eye<arma::mat>(s, s);
  arma::mat resMu(N, s);
  arma::vec resW(N);
  arma::mat resP(N, s2);

  for(int i = 0; i < N; i++) {
	for(int j=0; j<s;j++){
  		F(j,j)=theta_F(i,j);
  	}
	arma::mat mu_pred = F * mu.row(i).t(); // E[X_{t+1}|Y_{1:t}]
    
	count = 0;
	for(int j = 0; j < s; j++) {
		for(int k = 0; k < s; k++) {
			P(j, k) = Sigma(i, count);
			count++;
		}
	}
    
	arma::mat P_pred = F * P * F.t();
	for(int j = 0; j < s; j++) {
		P_pred(j, j) += theta_Q(i, j);
	}
    
	double S =arma::as_scalar(H * P_pred * H.t());
	S += theta_R(i);
    
    	arma::mat K = P_pred * H.t() / S;
    
   	 double zhat = y(i) - arma::as_scalar(H * mu_pred);
    
    	resMu.row(i) = (mu_pred + K * zhat).t(); 
    
    	arma::mat P_new = (I - K * H) * P_pred;
    	count = 0;
    	for(int j = 0; j < s; j++) {
      		for(int k = 0; k < s; k++) {
        		resP(i, count) = P_new(j, k);
        		count++;
      		}
    	}
    	resW(i) = -0.5 * std::log(S) - 0.5 * std::pow(zhat, 2) / S;
  	}
  
  return Rcpp::List::create(
    Rcpp::Named("P") = resP,
    Rcpp::Named("MU") = resMu,
    Rcpp::Named("W") = resW
  );
}


















