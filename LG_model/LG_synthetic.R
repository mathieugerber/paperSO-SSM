	 

sourceCpp("resampling.cpp")
sourceCpp("LG_model/KF.cpp")
source("LG_model/LG_algorithms.R")


################################################################################################################
### q=2 --- Algorithm A2-0.4
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1, dots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.4
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.4, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-04-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A2-04-q2_X.txt")


################################################################################################################
### q=2 --- Algorithm A2-0.5
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.5
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.5, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-05-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A2-05-q2_X.txt")



################################################################################################################
### q=2 --- Algorithm A2-0.6
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.6
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.6, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-06-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A2-06-q2_X.txt")





################################################################################################################
### q=2 --- Algorithm A2-1.1
################################################################################################################

##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-1.1
	Est<-LG_A2(y, nsMat, N=N.val, alpha=1.1, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-11-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A2-11-q2_X.txt")


################################################################################################################
### q=2 --- Algorithm A1-0.5
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A1-0.5
	Est<-LG_A1(y, nsMat, N=N.val, alpha=0.5, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A1-05-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A1-05-q2_X.txt")

################################################################################################################
### q=2 --- Algorithm A1-1.1
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-2			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A1-1.1
	Est<-LG_A1(y, nsMat, N=N.val, alpha=1.1, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A1-11-q2_theta.txt")
write.table(resX, "LG_model/LG_results/A1-11-q2_X.txt")



################################################################################################################
### q=4 --- Algorithm A2-0.4
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.4
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.4, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-04-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A2-04-q4_X.txt")


################################################################################################################
### q=4 --- Algorithm A2-0.5
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.5
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.5, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-05-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A2-05-q4_X.txt")



################################################################################################################
### q=4 --- Algorithm A2-0.6
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-0.6
	Est<-LG_A2(y, nsMat, N=N.val, alpha=0.6, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-06-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A2-06-q4_X.txt")





################################################################################################################
### q=4 --- Algorithm A2-1.1
################################################################################################################

##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A2-1.1
	Est<-LG_A2(y, nsMat, N=N.val, alpha=1.1, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A2-11-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A2-11-q4_X.txt")


################################################################################################################
### q=4 --- Algorithm A1-0.5
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A1-0.5
	Est<-LG_A1(y, nsMat, N=N.val, alpha=0.5, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A1-05-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A1-05-q4_X.txt")

################################################################################################################
### q=4 --- Algorithm A1-1.1
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#Number of spline functions (q>1)
##define theta_star
sigma.y<-0.5
sigX_star<-rep(1,q)  
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
 

##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,-1,0)


#create matrices to store the results
d<-3*q+1	 
resTheta<-matrix(0,M.val,d)
resX<-matrix(0,M.val,n.obs)
 

##Evaluate splines at points 0,1,adots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
 

##generate random seed
seed.vec<-.Random.seed

#parameter intial distribution for X
sigma0<-2
mu0<-0

##running time!
for(m in 1:M.val){		
	set.seed(seed.vec[m])
	##generate the data
	y<-rep(0,n.obs)
	x<-rnorm(q,mu0,sigma0)
	for(t in 1:n.obs){
		it<-1+24*( (t-1)/24-floor( (t-1)/24))
		mean_yt<-sum(nsMat[it,]*(theta_star[1:q]+x))
		y[t]<-rnorm(1,mean_yt,sd=theta_star[d])
		x<-theta_star[(q+1):(2*q)]*x+rnorm(q,0,sd=theta_star[(2*q+1):(3*q)])
	}
	##Kalman filter
	X_KF<-KM_filter(y,nsMat, sigma0=2, mu0=0, theta_val=theta_star)
	##Algorithm A1-1.1
	Est<-LG_A1(y, nsMat, N=N.val, alpha=1.1, sigma0=2, mu0=0, lb=lower, ub=upper)
	resTheta[m,]<-Est$THETA[n.obs,] 
	resX[m,]<-sqrt(apply((Est$XF-X_KF$XF)^2,1,sum))
	print(m)
}

write.table(resTheta, "LG_model/LG_results/A1-11-q4_theta.txt")
write.table(resX, "LG_model/LG_results/A1-11-q4_X.txt")


 














