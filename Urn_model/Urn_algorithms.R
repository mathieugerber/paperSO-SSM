
##theta=(j,k,r)
IF_Urn<-function(y, N, nit, alpha, alpha2, beta, ch=rep(1,4), ub, lb, c_ess=0.7,  t1=100, Delta=1)
{
	datalength<-length(y)
	tilde_theta<-matrix(0,nit*(datalength-1),3)
	ESS_bound<-c_ess*N 
	index<-0
	 
	
	theta<-matrix(0,N,3)
	size<-c(ub[1]-lb[1],ub[2]-lb[2],ub[3]-lb[3])
	Mat1<-(-size[1]):size[1]
	Mat2<-(-size[2]):size[2]
	Mat3<-(-size[3]):size[3]
	
	tp<-1+t1*(datalength-1)
	for(m in 1:nit){ 
	        A<-1:N
	        if(m==1){
	        	index<-1
	        	#sample theta0 from mu_0
			theta[,1]<-sample(c(lb[1]:ub[1]),N,replace=TRUE)
			theta[,2]<-sample(c(lb[2]:ub[2]),N,replace=TRUE)
			for(n in 1:N){
				upper<-min(c(ub[3],theta[n,1]+theta[n,2]))
				theta[n,3]<-sample(c(lb[3]:upper),1)	
			}
			w<-rep(0,N)
		}else{
			index<-index+1
			if(ESS<ESS_bound){
				A<-SSP_Resampler(runif(N), W)
				w<-rep(0,N)
			}
      			if(index==tp){ 
      				alpha_t<-log(index)^{-beta}
      				ht<-ch*index^{-alpha2*alpha_t}
      			}else{
      				ht<-ch*index^{-alpha}
      			}
      			
      			if(index==tp || ESS<ESS_bound){
      				for(j in 1:2){
					w1<-dbinom(0, size[j],prob=ht[j])
					w2<-dbinom(1:size[j], size[j],prob=ht[j])
					proba_vec<-c(w2[size[j]:1],2*w1,w2)
					u<-(2*size[j]+1)-(theta[A,j]-lb[j])
					l<-ub[j]-theta[A,j]+1
					theta[,j]<-theta[A,j]+Sample_vec(runif(N),u,l, proba_vec, (-size[j]):size[j])
			
				} 
				w1<-dbinom(0, size[3],prob=ht[3])
				w2<-dbinom(1:size[3], size[3],prob=ht[3])
			
				proba_vec<-c(w2[size[3]:1],2*w1,w2)
				upper<-theta[,1]+theta[,2]
				upper[upper>ub[3]]<-ub[3]
				u<-(2*size[3]+1)-(theta[A,3]-lb[3])+upper-ub[3]
				l<-ub[3]-theta[A,3]+1
				theta[,3]<-theta[A,3]+Sample_vec(runif(N),u,l, proba_vec, (-size[3]):size[3])
			}
			if(index==tp){
				tp<-tp+Delta*(datalength-1)*ceiling(log(tp)^{1-beta/2})
			}
		
		} 
	
		#Process first observations y2/y1
		t<-1
		work<- theta[,1]*theta[,2]
		if(y[t+1]==(y[t]-1)){
			W<-((theta[,1]-theta[,3]+y[t])*y[t] )/work
		}else if(y[t+1]== y[t]){
			W<-((theta[,3]-y[t])*y[t]+(theta[,1]-theta[,3]+y[t])*(theta[,2]-y[t]))/work
		}else{
			W<-(theta[,3]-y[t])*(theta[,2]-y[t])/work
		}
		W[W<=0]<- -1
		W[W>1]<- -1
		work<-rep(-Inf, N)
		work[W!=-1]<-log(W[W!=-1])
		
		work[is.nan(work)]<--Inf
		work[is.na(work)]<--Inf

		w<-w+work
		w1<- exp(w - max(w))
		W<- w1 / sum(w1)
		ESS<-1/sum(W^2)
	 	tilde_theta[index,]<- round(apply(W*theta,2,sum) )
	
		for(t in 2:(datalength-1)){
			A<-1:N
			index<-index+1
			if(ESS<ESS_bound){
				A<-SSP_Resampler(runif(N), W)
				w<-rep(0,N)
			}
      			ht<-ch*index^{-alpha}
      			
      			if(ESS<ESS_bound){
      				for(j in 1:2){
					w1<-dbinom(0, size[j],prob=ht[j])
					w2<-dbinom(1:size[j], size[j],prob=ht[j])
					proba_vec<-c(w2[size[j]:1],2*w1,w2)
					u<-(2*size[j]+1)-(theta[A,j]-lb[j])
					l<-ub[j]-theta[A,j]+1
					theta[,j]<-theta[A,j]+Sample_vec(runif(N),u,l, proba_vec, (-size[j]):size[j])
			
				}	 
				w1<-dbinom(0, size[3],prob=ht[3])
				w2<-dbinom(1:size[3], size[3],prob=ht[3])
			
				proba_vec<-c(w2[size[3]:1],2*w1,w2)
				upper<-theta[,1]+theta[,2]
				upper[upper>ub[3]]<-ub[3]
				u<-(2*size[3]+1)-(theta[A,3]-lb[3])+upper-ub[3]
				l<-ub[3]-theta[A,3]+1
				theta[,3]<-theta[A,3]+Sample_vec(runif(N),u,l, proba_vec, (-size[3]):size[3])	
			}
			#compute weights
			work<- theta[,1]*theta[,2]
			if(y[t+1]==(y[t]-1)){
				W<-((theta[,1]-theta[,3]+y[t])*y[t] )/work
			}else if(y[t+1]== y[t]){
				W<-((theta[,3]-y[t])*y[t]+(theta[,1]-theta[,3]+y[t])*(theta[,2]-y[t]))/work
			}else{
				W<-(theta[,3]-y[t])*(theta[,2]-y[t])/work
			}
			W[W<=0]<- -1
			W[W>1]<- -1
			work<-rep(-Inf, N)
			work[W!=-1]<-log(W[W!=-1])
		
			work[is.nan(work)]<--Inf
			work[is.na(work)]<--Inf

			w<-w+work
			w1<- exp(w - max(w))
			W<- w1 / sum(w1)
			ESS<-1/sum(W^2)
		
			tilde_theta[index,]<- round(apply(W*theta,2,sum) )
		}	
	}
	return(list(THETA=tilde_theta))
	
}	



 

	



	
