
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
 NumericVector Stratified_Resampler(const NumericVector points,  const NumericVector W){
 
	int N = points.size();
  	NumericVector J(N);
  	int j = 0;
 	double s = 0.0;
  	double inv_N=1.0/(double)N;
  
	for(int i=0;i<N;i++)
	{
		s+= W[i];
		while(j < N && (j+points[j])*inv_N<=s)
		{
			J[j]=i+1;
			j++;
		}
	}
	return J;
}


// [[Rcpp::export]]
 NumericVector SSP_Resampler(const NumericVector points,  const NumericVector W){
 
	int  N=points.size();
	double bound;
	double cons=log(0.5);
	int test;
	NumericVector J(N), J2(N);
	
	for(int i=0;i<N;i++){
		J[i]=N*W[i];
	}
	
	int i=0;
	int n=0;
	int m=1;
	int count=0;
	while(m<N){
		double deltaN=ceil(J[n])-J[n];
		double deltaM= J[m]-floor(J[m]);
		double delta=deltaN;
		if(deltaM<deltaN){
			delta=deltaM;
		}
		double epsilonN= J[n]-floor(J[n]);
		double epsilonM=ceil(J[m])-J[m];
		double epsilon=epsilonN;
		if(epsilonM<epsilonN){
			epsilon=epsilonM;
		}
		if(delta==epsilon){
			bound=cons;
		}else{
			bound=log(epsilon)-log(delta+epsilon);
		}
		test=-1;
		if(log(points[i])<= bound){
			J[n]+= delta;
			J[m]+= -delta;
			if(delta==deltaN){
				count+=J[n];
				if(deltaM==deltaN){
					count+=J[m];
					n=m+1;
					m=m+2;
				}else{
					test=m;
					n=m;
					m=m+1;
				}
			}else{
			     count+=J[m];
			     test=n;
			     m=m+1;
			}
		}else{
			J[n]+= -epsilon;
			J[m]+= epsilon;
			if(epsilon==epsilonN){
				count+=J[n];
				if(epsilonM==epsilonN){
					count+=J[m];
					n=m+1;
					m=m+2;
				}else{
					test=m;
					n=m;
					m=m+1;
				}
			}else{
				count+=J[m];
				test=n;
				m=m+1;
			}
		}
		i++;
	}
	 
	if(test>=0){
		J[test]=N-count;
	}
	
	i=0;
	m=0;
	while(i<N && m<N){
		n=1;
		while( i<N && n<=J[m]){
			J2[i]=m+1;
			i++;
			n++;
		}
		m++;
		
	}
	return J2;
	
}



















