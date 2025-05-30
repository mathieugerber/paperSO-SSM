



##Adaptive slowly vanishing artificial dynamic

LG_A2<-function(y, splines, N, alpha, ch=1, sigma0, mu0, lb, ub, c_ess=0.7,  nu=100, t1=100, Delta=1)
{
	#pre-compute some quantities
	datalength<-length(y)
	ESS_bound<-c_ess*N 
	q<-ncol(splines)
	ns<-nrow(splines)
	d<-3*q+1 #total number of parameters
	
	
	#initializes some vectors and matrices
	ess.vec<-rep(0,datalength)
	X_est<-matrix(0,datalength,q)
	theta<-matrix(0,N,d)
	tilde_theta<-matrix(0,datalength,d)
	
	##sample from thata0 from mu_0(\dd theta) [here K_1(theta',dd theta)=delta_{\theta'}
	
	#beta parameters
	for(i in 1:q){
		theta[,i]<-runif(N,lb[1],ub[1])    
	}
	#rho parameters
	for(i in 1:q){
		theta[,q+i]<-runif(N,lb[2],ub[2])    
	}
	#sigma parameters
	for(i in 1:(q+1)){
		theta[,2*q+i]<-runif(N,lb[3],ub[3])    
	}
	
	
	#initial matrix of mean (mu_{1|1})
	mu<-matrix(mu0,N,q)
	
	#initial covariance matrix (P_{1|1})
	Sigma<-matrix(0,N,q^2)    #matrix(Sigma[1,],nrow=s,byrow=T)
	work<-diag(sigma0^2,q,q)
	for(i in 1:N){
		Sigma[i,]<-c(t(work))
	}
	
	#process first observation
	t<-1
	it<-1+24*( (t-1)/24-floor( (t-1)/24))
	mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
	H_mat<-matrix(splines[it,],ncol=q)
	
	res<-KF_fast(mu, Sigma, as.matrix(theta[,( (q+1):(2*q))]), as.matrix(theta[,( (2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
	mu<-as.matrix(res$MU)
	Sigma<-as.matrix(res$P)
	w<-c(res$W)
	
	w[is.nan(w)]<--Inf
	w1<- exp(w - max(w))
	W<- w1 / sum(w1)
	ESS<-1/sum(W^2)
	ess.vec[t]<-ESS
	
	#Estimate theta and E[X_t |Y{1:t}
	tilde_theta[t,]<-apply(W*theta,2,sum)
	theta_bar<-tilde_theta[t,]
	X_est[t,]<-apply(W*mu,2,sum) 
		 	
 	#iterate
	if(alpha>1){
		tp<--1
	}else{
		tp<-t1
	}
	for(t in 2:datalength){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		H_mat<-matrix(splines[it,],ncol=q)
		
		A<-1:N
		if(ESS<ESS_bound){
			A<-SSP_Resampler(runif(N), W)
			w<-rep(0,N)
		}
		ht<-ch*t^{-alpha}
      		
      		if(t==tp || ESS<ESS_bound){
		 	if(t==tp){
				for(i in 1:q){
					theta[,i]<-rtt(N,location= theta[A,i], scale=ht, df=nu, left=lb[1], right=ub[1])  
				}
				for(i in 1:q){
					theta[,q+i]<-rtt(N,location= theta[A,q+i], scale=ht, df=nu, left=lb[2], right=ub[2])  
				}
				for(i in 1:q){
					theta[,2*q+i]<-rtt(N,location= theta[A,2*q+i], scale=ht, df=nu, left=lb[3], right=ub[3])  
				}
				theta[,d]<-rtt(N,location= theta[A,d], scale=ht, df=nu, left=lb[3], right=ub[3])
				
				##update tp
				tp<-tp+Delta*ceiling(log(tp)^2)
				
			}else{
				for(i in 1:q){
					theta[,i]<-msm::rtnorm(N, theta[A,i], ht,lb[1],ub[1]) 
				}	
				for(i in 1:q){
					theta[,q+i]<-msm::rtnorm(N, theta[A,q+i], ht,lb[2],ub[2]) 
				}
				for(i in 1:q){
					theta[,2*q+i]<-msm::rtnorm(N, theta[A,2*q+i], ht, lb[3],ub[3]) 
				}
				theta[,d]<-msm::rtnorm(N, theta[A,d], ht, lb[3],ub[3]) 
			}
		}
		mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
		res<-KF_fast(mu[A,], Sigma[A,], as.matrix(theta[,((q+1):(2*q))]), as.matrix(theta[,((2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
		
		mu<-res$MU
		Sigma<-res$P
		w<-w+c(res$W)
	
		w[is.nan(w)]<--Inf
		w1<- exp(w - max(w))
		W<- w1 / sum(w1)
		ESS<-1/sum(W^2)
		ess.vec[t]<-ESS
	
		#Estimate theta and E[X_t |Y{1:t}
		tilde_theta[t,]<-apply(W*theta,2,sum)
		theta_bar<-( (t-1)*theta_bar+tilde_theta[t,])/t
		X_est[t,]<-apply(W*mu,2,sum) 

	}
	
   	return(list(THETA=tilde_theta, XF=X_est))
   
}


##Non-adaptive slowly vanishing artificial dynamic
LG_A1<-function(y, splines, N, alpha, ch=1, sigma0, mu0, lb, ub, c_ess=0.7,  nu=100, t1=100, Delta=1)
{
	#pre-compute some quantities
	datalength<-length(y)
	ESS_bound<-c_ess*N 
	q<-ncol(splines)
	ns<-nrow(splines)
	d<-3*q+1 #total number of parameters	
	
	#initializes some vectors and matrices
	ess.vec<-rep(0,datalength)
	X_est<-matrix(0,datalength,q)
	theta<-matrix(0,N,d)
	tilde_theta<-matrix(0,datalength,d)
	
	##sample from thata0 from mu_0(\dd theta) [here K_1(theta',dd theta)=delta_{\theta'}
	
	#beta parameters
	for(i in 1:q){
		theta[,i]<-runif(N,lb[1],ub[1])    
	}
	#rho parameters
	for(i in 1:q){
		theta[,q+i]<-runif(N,lb[2],ub[2])    
	}
	#sigma parameters
	for(i in 1:(q+1)){
		theta[,2*q+i]<-runif(N,lb[3],ub[3])    
	}
	
	
	#initial matrix of mean (mu_{1|1})
	mu<-matrix(mu0,N,q)
	
	#initial covariance matrix (P_{1|1})
	Sigma<-matrix(0,N,q^2)    #matrix(Sigma[1,],nrow=s,byrow=T)
	work<-diag(sigma0^2,q,q)
	for(i in 1:N){
		Sigma[i,]<-c(t(work))
	}
	
	#process first observation
	t<-1
	it<-1+24*( (t-1)/24-floor( (t-1)/24))
	mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
	H_mat<-matrix(splines[it,],ncol=q)
	
	res<-KF_fast(mu, Sigma, as.matrix(theta[,( (q+1):(2*q))]), as.matrix(theta[,( (2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
	mu<-as.matrix(res$MU)
	Sigma<-as.matrix(res$P)
	w<-c(res$W)
	
	w[is.nan(w)]<--Inf
	w1<- exp(w - max(w))
	W<- w1 / sum(w1)
	ESS<-1/sum(W^2)
	ess.vec[t]<-ESS
	
	#Estimate theta and E[X_t |Y{1:t}
	tilde_theta[t,]<-apply(W*theta,2,sum)
	theta_bar<-tilde_theta[t,]
	X_est[t,]<-apply(W*mu,2,sum) 
		 	
 	#iterate
	if(alpha>1){
		tp<--1
	}else{
		tp<-t1
	}
	for(t in 2:datalength){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		H_mat<-matrix(splines[it,],ncol=q)
		
		A<-1:N
		if(ESS<ESS_bound){
			A<-SSP_Resampler(runif(N), W)
			w<-rep(0,N)
		}
		ht<-ch*t^{-alpha}
      		
		 if(t==tp){
			for(i in 1:q){
				theta[,i]<-rtt(N,location= theta[A,i], scale=ht, df=nu, left=lb[1], right=ub[1])  
			}
			for(i in 1:q){
				theta[, q+i]<-rtt(N,location= theta[A,q+i], scale=ht, df=nu, left=lb[2], right=ub[2])  
			}
			for(i in 1:q){
				theta[,2*q+i]<-rtt(N,location= theta[A,2*q+i], scale=ht, df=nu, left=lb[3], right=ub[3])  
			}
			theta[,d]<-rtt(N,location= theta[A,d], scale=ht, df=nu, left=lb[3], right=ub[3])
			
			##update tp
			tp<-tp+Delta*ceiling(log(tp)^2)
				
		}else{
			for(i in 1:q){
				theta[,i]<-msm::rtnorm(N, theta[A,i], ht,lb[1],ub[1]) 
			}	
			for(i in 1:q){
				theta[,q+i]<-msm::rtnorm(N, theta[A,q+i], ht,lb[2],ub[2]) 
			}
			for(i in 1:q){
				theta[,2*q+i]<-msm::rtnorm(N, theta[A,2*q+i], ht, lb[3],ub[3]) 
			}
			theta[,d]<-msm::rtnorm(N, theta[A,d], ht, lb[3],ub[3]) 
		}
		mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
		res<-KF_fast(mu[A,], Sigma[A,], as.matrix(theta[,((q+1):(2*q))]), as.matrix(theta[,((2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
		mu<-res$MU
		Sigma<-res$P
		w<-w+c(res$W)
	
		w[is.nan(w)]<--Inf
		w1<- exp(w - max(w))
		W<- w1 / sum(w1)
		ESS<-1/sum(W^2)
		ess.vec[t]<-ESS
	
		#Estimate theta and E[X_t |Y{1:t}
		tilde_theta[t,]<-apply(W*theta,2,sum)
		theta_bar<-( (t-1)*theta_bar+tilde_theta[t,])/t
		X_est[t,]<-apply(W*mu,2,sum) 
	
	}

   	return(list(THETA=tilde_theta, XF=X_est))
   
}




##Adaptive slowly vanishing artificial dynamic with variance reset times

LG_A2_Real<-function(yy, splines, N, alpha, ch=1, sigma0, mu0, lb, ub, k.val=1, reset.time, reset.val=1, c_ess=0.7,  nu=100, t1=100, Delta=1)
{
	#pre-compute some quantities
	y0<-yy[1]
	y<-yy[-1]
	datalength<-length(y)
	ESS_bound<-c_ess*N 
	nk<-length(k.val)
	q<-ncol(splines)
	ns<-nrow(splines)
	d<-3*q+1 #total number of parameters
	
	
	#initializes some vectors and matrices
	pred<-matrix(0,datalength, length(k.val))
	theta<-matrix(0,N,d)
	tilde_theta<-matrix(0,datalength,d)
	
	##sample from thata0 from mu_0(\dd theta) [here K_1(theta',dd theta)=delta_{\theta'}
	
	#beta parameters
	for(i in 1:q){
		theta[,i]<-runif(N,lb[1],ub[1])    
	}
	#rho parameters
	for(i in 1:q){
		theta[,q+i]<-runif(N,lb[2],ub[2])    
	}
	#sigma parameters
	for(i in 1:(q+1)){
		theta[,2*q+i]<-runif(N,lb[3],ub[3])    
	}
	
	
	#initial matrix of mean (mu_{1|1})
	mu<-matrix(mu0,N,q)
	
	#initial covariance matrix (P_{1|1})
	Sigma<-matrix(0,N,q^2)    #matrix(Sigma[1,],nrow=s,byrow=T)
	work<-diag(sigma0^2,q,q)
	for(i in 1:N){
		Sigma[i,]<-c(t(work))
	}
	
	#process first observation
	t<-1
	it<-1+24*( (t-1)/24-floor( (t-1)/24))
	mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
	H_mat<-matrix(splines[it,],ncol=q)
	
	res<-KF_fast(mu, Sigma, as.matrix(theta[,( (q+1):(2*q))]), as.matrix(theta[,( (2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-y0-mean_ft))
	mu<-as.matrix(res$MU)
	Sigma<-as.matrix(res$P)
	w<-c(res$W)
	
	w[is.nan(w)]<--Inf
	w1<- exp(w - max(w))
	W<- w1 / sum(w1)
	ESS<-1/sum(W^2)
	
	#Estimate theta
	tilde_theta[t,]<-apply(W*theta,2,sum)
 	
	#Compute predictions
	Mu.t<-apply(W*mu,2,sum) 
	theta_use<-tilde_theta[t,]
	for(k in 1:nk){
		work<-  Mu.t*theta_use[((q+1):(2*q))]^{k.val[k]}
		it<-1+24*( (t+k.val[k]-1)/24-floor( (t+k.val[k]-1)/24))
		mean_ft<-sum(splines[it,]*theta_use[1:q])
	 	pred[t,k]<-y0+mean_ft+sum(splines[it,]*work)
	}
		 	
 	#iterate
	tp<-t1
	count<-2
	for(t in 2:datalength){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		if(it==1){
			y0<-y[t-1]
		}
		H_mat<-matrix(splines[it,],ncol=q)
		
		A<-1:N
		if(ESS<ESS_bound){
			A<-SSP_Resampler(runif(N), W)
			w<-rep(0,N)
		}
		if(t %in% reset.time){
			count<-reset.val
		}
		ht<-ch*count^{-alpha}
		count<-count+1
      		
      		if(t==tp || ESS<ESS_bound){
		 	if(t==tp){
				for(i in 1:q){
					theta[,i]<-rtt(N,location= theta[A,i], scale=ht, df=nu, left=lb[1], right=ub[1])  
				}
				
				for(i in 1:q){
					theta[,q+i]<-rtt(N,location= theta[A,q+i], scale=ht, df=nu, left=lb[2], right=ub[2])  
				}
				
				for(i in 1:q){
					theta[,2*q+i]<-rtt(N,location= theta[A,2*q+i], scale=ht, df=nu, left=lb[3], right=ub[3])  
				}
				
				theta[,d]<-rtt(N,location= theta[A,d], scale=ht, df=nu, left=lb[3], right=ub[3])
				
				##update tp
				tp<-tp+Delta*ceiling(log(tp)^2)
				
			}else{
				for(i in 1:q){
					theta[,i]<-msm::rtnorm(N, theta[A,i], ht,lb[1],ub[1]) 
				}
					
				for(i in 1:q){
					theta[,q+i]<-msm::rtnorm(N, theta[A,q+i], ht,lb[2],ub[2]) 
				}
				
				for(i in 1:q){
					theta[,2*q+i]<-msm::rtnorm(N, theta[A,2*q+i], ht, lb[3],ub[3]) 
				}
				
				theta[,d]<-msm::rtnorm(N, theta[A,d], ht, lb[3],ub[3]) 
			}
		}
		
		mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
		
		res<-KF_fast(mu[A,], Sigma[A,], as.matrix(theta[,((q+1):(2*q))]), as.matrix(theta[,((2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-y0-mean_ft))
	
		mu<-res$MU
		Sigma<-res$P
		w<-w+c(res$W)
	
		w[is.nan(w)]<--Inf
		w1<- exp(w - max(w))
		W<- w1 / sum(w1)
		ESS<-1/sum(W^2)

		#Estimate theta
		tilde_theta[t,]<-apply(W*theta,2,sum)
	
		#Compute predictions
		Mu.t<-apply(W*mu,2,sum) 
		theta_use<-tilde_theta[t,]
		
		tk<-24*floor((t-1)/24)+24
	 	work<- Mu.t*theta_use[( (q+1):(2*q))]^{(tk-t)}
	 	mean_ft<-sum(splines[24,]*theta_use[1:q])
		yk<-y0+mean_ft+sum(splines[24,]*work)
		
		for(k in 1:nk){
			work<-  Mu.t*theta_use[((q+1):(2*q))]^{k.val[k]}
			itk<-1+24*( (t+k.val[k]-1)/24-floor( (t+k.val[k]-1)/24))
			mean_ft<-sum(splines[itk,]*theta_use[1:q])
			if(itk>it){
	 			pred[t,k]<-y0+mean_ft+sum(splines[itk,]*work)
	 		}else{
	 			pred[t,k]<-yk+mean_ft+sum(splines[itk,]*work)
	 		}
		}
	
		

	}
   	return(list(THETA=tilde_theta, PRED=pred))
   
}


##Kalman filter fixed theta

KM_filter<-function(y, splines, sigma0, mu0, theta_val)
{
	#pre-compute some quantities
	datalength<-length(y)
	q<-ncol(splines)
	ns<-nrow(splines)
	d<-3*q+1 #total number of parameters
	
	
	#initializes some vectors and matrices
	X_est<-matrix(0,datalength,q)
	N<-2
	theta<-matrix(rep(theta_val,N),N,d, byrow=TRUE)
	
	#initial matrix of mean (mu_{1|1})
	mu<-matrix(mu0,N,q)
	
	#initial covariance matrix (P_{1|1})
	Sigma<-matrix(0,N,q^2)    #matrix(Sigma[1,],nrow=s,byrow=T)
	work<-diag(sigma0^2,q,q)
	for(i in 1:N){
		Sigma[i,]<-c(t(work))
	}
	
	#process first observation
	t<-1
	it<-1+24*( (t-1)/24-floor( (t-1)/24))
	
	mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
	H_mat<-matrix(splines[it,],ncol=q)
	
	res<-KF_fast(mu, Sigma, as.matrix(theta[,( (q+1):(2*q))]), as.matrix(theta[,( (2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
	mu<-as.matrix(res$MU)
	Sigma<-as.matrix(res$P)
	
	X_est[t,]<-mu[1,]
	
	for(t in 2:datalength){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		H_mat<-matrix(splines[it,],ncol=q)
		mean_ft<-apply(t(splines[it,]*t(theta[,1:q])),1,sum)
		res<-KF_fast(mu, Sigma, as.matrix(theta[,((q+1):(2*q))]), as.matrix(theta[,((2*q+1):(3*q))]^2), theta[,d]^2, H_mat, c(y[t]-mean_ft))
		mu<-res$MU
		Sigma<-res$P
		X_est[t,]<-mu[1,]

	}
   	return(list(XF=X_est))
   
}

