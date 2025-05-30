sourceCpp("resampling.cpp")
sourceCpp("LG_model/KF.cpp")
source("LG_model/LG_algorithms.R")



################################################################################################################
### Load and prepare data 
### source: https://www.kaggle.com/datasets/selfishgene/historical-hourly-weather-data?select=temperature.csv_
################################################################################################################
data<-read.csv("Data/temperature.csv")

#select city of Portland (for which there is no missing values)
k<-3 			
data<-data[,c(1,k)] 
#Make the dataset set starts with the temperature on the 1st of January 2013 at 00:00am
data<-data[-(1:2196), ]  

#keep first n.obs observations.
data<-data[1:n.obs,2]	 

################################################################################################################
### q=3 Algorithm 2 (alpha=0.5)
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-3			#number of spline function (q>1)
kk<-1:8		#prediction horizon (in hours)


##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,0,0)


#create matrices to store the results
d<-3*q+1	 
datalength<-length(data)-1
resPred<-array(0,c(M.val,datalength,length(kk)))
 
##Evaluate splines at points 0,1, dots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
ncol(nsMat)

##times when the sequence of scale matrices is reset
change_points<-24*c(30*3, (30*3+30*6*(1:5)))

##running time!
for(m in 1:M.val){		
	Est<-LG_A2_Real(data, nsMat, N=N.val, alpha=0.5, sigma0=2, mu0=0,  lb=lower, ub=upper, k.val=kk, reset.time=change_points)
	resPred[m,,]<-Est$PRED
	print(m)
}

for(k in 1:length(kk)){
	 write.table(resPred[,,k], paste("LG_model/LG_results/Real_q3_resPred_",kk[k],".txt", sep=''))
}


################################################################################################################
### q=4 Algorithm 2 (alpha=0.5)
################################################################################################################
##Set the seed for the RNG
set.seed(9585)

q<-4			#number of spline function (q>1)
kk<-1:8		#prediction horizon (in hours)


##upper/lower bound for the parmeters (beta, rho, std)
upper<-c(10,1,4)
lower<-c(-10,0,0)


#create matrices to store the results
d<-3*q+1	 
datalength<-length(data)-1
resPred<-array(0,c(M.val,datalength,length(kk)))
theta<-{}

##Evaluate splines at points 0,1, dots,24
nn<-24/q
nsMat <- naturalSpline(0:24, knots=c( nn*(1:(q-1)) ), Boundary.knots=c(0,24), intercept = FALSE)
nsMat<-nsMat[-1,]
ncol(nsMat)

##times when the sequence of scale matrices is reset
change_points<-24*c(30*3, (30*3+30*6*(1:5)))
##running time!
for(m in 1:M.val){		
	Est<-LG_A2_Real(data, nsMat, N=N.val, alpha=0.5,  sigma0=2, mu0=0,  lb=lower, ub=upper, k.val=kk, reset.time=change_points)
	resPred[m,,]<-Est$PRED
	if(m==1) theta<-Est$THETA
	print(m)
}


write.table(theta, "LG_model/LG_results/Real_q4_traj.txt")
for(k in 1:length(kk)){
	 write.table(resPred[,,k], paste("LG_model/LG_results/Real_q4_resPred_",kk[k],".txt", sep=''))
}

 






