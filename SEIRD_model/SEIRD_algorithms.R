


A3_SEIRD<- function(y, N, nit, alpha, ch=1, lb, ub, chi_theta, c_lambda, c_kappa, q.val=0, c_ess=0.7, nu=100, t1=100, Delta=1){
 
	datalength<-nrow(y) 
	tilde_theta<-matrix(0, nit*datalength,7)
	cons<-log(c_lambda)
  	ESS_bound<-c_ess*N
  	index<-0
  	ess.vec<-rep(0,nit*datalength)
  
  	if(alpha>1){
   		 tp <- -1
  	}else{
    		 tp<- 1+t1*datalength
  	}
  
  	for(m in 1:nit){
    		index<-index+1
    		
    		if(m==1){
      			w <- rep(0,N)
      			theta <- matrix(0,ncol=7,nrow=N)
      			theta[,2] <- runif(N,min=lb[2],max=1-theta[,2])
      			theta[,1] <- runif(N,min=lb[1],max=1-theta[,2])
      			theta[,3] <- runif(N,min=lb[3],max=1-theta[,2])
      			for(j in 4:7){
       			 theta[,j] <- runif(N,min=lb[j],max=ub[j])
      			}
   		}else{
     			 
     			if(ESS<ESS_bound){
        			A<-SSP_Resampler(runif(N), W)
        			w<-rep(0,N)
      			}
      			
     			if(ESS<ESS_bound || index==tp){
        			if(index==tp){
        			
        				theta[,2] <- rtt(N,theta[A,2], ch*index^(-alpha),df=nu,left=lb[2],right=ub[2])
        				theta[,1] <- rtt(N,theta[A,1], ch*index^(-alpha),df=nu,left=lb[1],right=1-theta[,2])
          				theta[,3] <- rtt(N,theta[A,3], ch*index^(-alpha),df=nu,left=lb[3],right=1-theta[,2])
          				for(j in 4:7){
            					theta[,j] <- rtt(N,theta[A,j], ch*index^(-alpha),df=nu,left=lb[j],right=ub[j])
         				}
          				tp <- tp+Delta*datalength*ceiling((log(tp))^2)
       			 }else{
       			 	theta[,2] <- rtnorm(N,theta[A,2], ch*index^(-alpha),lower=lb[2],upper=ub[2])
       			 	theta[,1] <- rtnorm(N,theta[A,1], ch*index^(-alpha),lower=lb[1],upper=1-theta[,2])
          				theta[,3] <- rtnorm(N,theta[A,3], ch*index^(-alpha),lower=lb[3],upper=1-theta[,2])
          				for(j in  4:7){
            					theta[,j] <- rtnorm(N,theta[A,j], ch*index^(-alpha),lower=lb[j],upper=ub[j])
          				}
          				
        			}
        			
      			}
    		}
      		logLambda1<-log(theta[,6])
    		logLambda2<-log(theta[,7])

    		############ Sample X1
    		
    		E <-  runif(N, min=chi_theta$lower[2], max = chi_theta$upper[2])
    		I <-  runif(N, min=chi_theta$lower[3], max = chi_theta$upper[3])
    		R <- rep(0,N)
    		D <- rep(0,N)
    		S  <- 1-I-E
    		upper<-4*(theta[,1]+theta[,2])/S
    		upper[upper>chi_theta$upper[1]]<-chi_theta$upper[1]
    		beta <- runif(N, min=0, max = upper)
    		X <- cbind(S,E,I,R,D,beta)
    		
   		d_newcase <- dbeta(y[1,1],shape1=exp(cons+log(theta[,1])+log(X[,2])-logLambda1) ,shape2=exp(cons+log(1-theta[,1]*X[,2])-logLambda1),log=TRUE)
    		d_newdeath <- dbeta(y[1,2],shape1=exp(cons+log(theta[,3])+log(X[,3])-logLambda2),shape2=exp(cons+log(1-theta[,3]*X[,3])-logLambda2),log=TRUE)
    
    		work <- d_newcase+d_newdeath
    		work[is.nan(work)]<- -Inf
    		w<- w+work
    		w1<- exp(w - max(w))
    		W<- w1 / sum(w1)
    		ESS<-1/sum(W^2)
    		tilde_theta[index,]<-apply(W*theta,2,sum)
      	 	
 		
    		for(t in 2:datalength) {
      			index<-index+1
      			
      			if(ESS<ESS_bound){
        			A<-SSP_Resampler(runif(N), W)
        			w<-rep(0,N)
        			
     				X<-X[A,]
     				theta[,2] <- rtnorm(N,theta[A,2], ch*index^(-alpha),lower=lb[2],upper=ub[2])
     				theta[,1] <- rtnorm(N,theta[A,1], ch*index^(-alpha),lower=lb[1],upper=1-theta[,2])
        			theta[,3] <- rtnorm(N,theta[A,3], ch*index^(-alpha),lower=lb[3],upper=1-theta[,2])
        			for(j in  4:7){
          				theta[,j] <- rtnorm(N,theta[A,j], ch*index^(-alpha),lower=lb[j],upper=ub[j])
        			}
        			 
        			
      			}
      			logKappa<-log(theta[,5])
      			logLambda1<-log(theta[,6])
    			logLambda2<-log(theta[,7])
    			
    			
    			p1<-X[,1]-X[,6]*X[,1]*(X[,2]+q.val*X[,3])
    			p2<-X[,2]+X[,6]*X[,1]*(X[,2]+q.val*X[,3])-(theta[,1]+theta[,2])*X[,2]
    			p3<-X[,3]+theta[,1]*X[,2]-theta[,2]*X[,3]-theta[,3]*X[,3]
    			p4<-X[,4]+theta[,2]*X[,2]+theta[,2]*X[,3]
    			p5<-X[,5]+theta[,3]*X[,3]
    			
    			concentration_parameter <- c_kappa*cbind(exp(log(p1)-logKappa),
                                    			 exp(log(p2)-logKappa),
                                    			 exp(log(p3)-logKappa),
                                    			 exp(log(p4)-logKappa),
                                     			 exp(log(p5)-logKappa))
      
      			X[,1:5] <- LaplacesDemon::rdirichlet(N,alpha=concentration_parameter)
      			X[,6]<-exp(rtnorm(N,  log(X[,6]), theta[,4],lower=-Inf,upper=0))
      			
   			d_newcase <- dbeta(y[t,1],shape1=exp(cons+log(theta[,1])+log(X[,2])-logLambda1) ,shape2=exp(cons+log(1-theta[,1]*X[,2])-logLambda1),log=TRUE)
    			d_newdeath <- dbeta(y[t,2],shape1=exp(cons+log(theta[,3])+log(X[,3])-logLambda2),shape2=exp(cons+log(1-theta[,3]*X[,3])-logLambda2),log=TRUE)
    		
      			work <- d_newcase+d_newdeath
      			work[is.nan(work)]<- -Inf
     			w<- w+work
      			w1<- exp(w - max(w))
      			W<- w1 / sum(w1)
      			ESS<-1/sum(W^2)
     			tilde_theta[index,]<-apply(W*theta,2,sum) 
		}
    
  	}
  	return(list(THETA=tilde_theta))
}



#Particle filter
PF <- function(y, theta, N, chi_theta,  c_lambda, c_kappa, c_ess=0.7,q.val=0){
	datalength<-nrow(y)
	Xmat<-matrix(0,datalength,6)
	cons<-log(c_lambda)

	ltheta1<-log(theta[1])
	ltheta3<-log(theta[3])
	ltheta4<-log(theta[4])
	logKappa<-log(theta[5])
	logLambda1<-log(theta[6])
    	logLambda2<-log(theta[7])
    	
  	ESS_bound<-c_ess*N
  	ess.vec <- rep(0,datalength)
 	R_val <- rep(0,datalength)
  
    	E <-  runif(N, min=chi_theta$lower[2], max = chi_theta$upper[2])
    	I <-  runif(N, min=chi_theta$lower[3], max = chi_theta$upper[3])
    	R <- rep(0,N)
    	D <- rep(0,N)
    	S  <- 1-I-E
  
  	upper<-4*(theta[1]+theta[2])/S
    	upper[upper>chi_theta$upper[1]]<-chi_theta$upper[1]
    	beta <- runif(N, min=0, max = upper)
    
    	X <- cbind(S,E,I,R,D,beta)
    		
  	## compute weight
  	t<-1
   	d_newcase <- dbeta(y[t,1],shape1=exp(cons+ltheta1+log(X[,2])-logLambda1) ,shape2=exp(cons+log(1-theta[1]*X[,2])-logLambda1),log=TRUE)
    	d_newdeath <- dbeta(y[t,2],shape1=exp(cons+ltheta3+log(X[,3])-logLambda2),shape2=exp(cons+log(1-theta[3]*X[,3])-logLambda2),log=TRUE)
 
  	work <- d_newcase+d_newdeath
  	work[is.nan(work)]<- -Inf
  	w <- work
  	w1<- exp(w - max(w)) 
 	W<- w1 / sum(w1)
 	
  	Xmat[t,]<-apply(W*X,2,sum)
  	ESS<-1/sum(W^2)
  	ess.vec[t]<-ESS

  	## Filtered Rnumber before resampling
  	aaa <- q.val*theta[1]/(theta[2]+theta[3])
	bbb <- X[,1]*X[,6]/(theta[2] + theta[1])
  	R_val[t] <- sum(W*(1+aaa)*bbb)
  
  	for(t in 2:datalength){
		A<-1:N
   	 	if(ESS<ESS_bound){
      			A<-SSP_Resampler(runif(N), W)
      			w<-rep(0,N) 
      			X <- X[A,]
    		}
    
    		## Sample X_{t}
    		p1<-X[,1]-X[,6]*X[,1]*(X[,2]+q.val*X[,3])
    		p2<-X[,2]+X[,6]*X[,1]*(X[,2]+q.val*X[,3])-(theta[1]+theta[2])*X[,2]
    		p3<-X[,3]+theta[1]*X[,2]-theta[2]*X[,3]-theta[3]*X[,3]
    		p4<-X[,4]+theta[2]*X[,2]+theta[2]*X[,3]
    		p5<-X[,5]+theta[3]*X[,3]
    		
    		concentration_parameter <- c_kappa*cbind(exp(log(p1)-logKappa),
                                                exp(log(p2)-logKappa),
                                                exp(log(p3)-logKappa),
                                                exp(log(p4)-logKappa),
                                                exp(log(p5)-logKappa))
                                                
               X[,1:5] <- LaplacesDemon::rdirichlet(N,alpha=concentration_parameter)
      	        X[,6]<-exp(rtnorm(N,  log(X[,6]), theta[4],lower=-Inf,upper=0))
    		
    		## Compute weight
    		d_newcase <- dbeta(y[t,1],shape1=exp(cons+ltheta1+log(X[,2])-logLambda1) ,shape2=exp(cons+log(1-theta[1]*X[,2])-logLambda1),log=TRUE)
    	        d_newdeath <- dbeta(y[t,2],shape1=exp(cons+ltheta3+log(X[,3])-logLambda2),shape2=exp(cons+log(1-theta[3]*X[,3])-logLambda2),log=TRUE)
 
    		work <- d_newcase+d_newdeath
    		work[is.nan(work)]<- -Inf
    		w<- w+work 
    		w1<- exp(w - max(w)) 
    		W<- w1 / sum(w1)
 
    		ESS<-1/sum(W^2)
    		ess.vec[t]<-ESS
    
    		## Filtered Rnumber before resampling
    		aaa <- q.val*theta[1]/(theta[2]+theta[3])
    		bbb <- X[,1]*X[,6]/(theta[2] + theta[1])
    		R_val[t] <- sum(W*(1+aaa)*bbb)
    		Xmat[t,]<-apply(W*X,2,sum) 
    
  	}
  	return(list(R=R_val))
  
}

