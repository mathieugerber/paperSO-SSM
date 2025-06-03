
sourceCpp("resampling.cpp")
source("SEIRD_model/SEIRD_algorithms.R")

################################################################################################################
### Load and prepare data 
data <- read.csv('Data/covid.csv')
##estimate size of UK population in 2020
UK_pop <- 67886004
##observations used to fit the model
obs<-cbind(data$new_cases_smoothed/UK_pop,data$new_deaths_smoothed/UK_pop)
##R number given in the OWID datase
R<-data$reproduction_rate
datalength<-length(R)


##upper/lower bounds for  theta=(eta,gamma,mu,sigmaBeta,kappa,lambda,lambda2)
upper <- c(1,1,1,1,1,1,1) 
lower <- c(0,0,0,0,0.01,0.01,0.01)

##upper/lower bounds for X_1
upperX=c(1,10^{-4},10^{-4})
lowerX=c(0,0,0)

#store this as a list
param<-list(lower=lowerX, upper=upperX)


################################################################
# Algorithm 3--alpha=0.5
################################################################

##Set the seed for the RNG
set.seed(9585)
 
#create matrices to store the results
datalength<-nrow(obs)
d<-7		 
burnIn<-(1:(b.val*datalength))   
res1<-matrix(0,M.val,d)
res2<-matrix(0,M.val,d)
resR1<-matrix(0,M.val,datalength)
resR2<-matrix(0,M.val,datalength)

##generate random seed
seed.vec<-.Random.seed
##running time!
for(m in  1:M.val){	
	set.seed(seed.vec[m])
	Est<-A3_SEIRD(obs,N=N.val, nit=n.it, alpha=0.5,  chi_theta=param,  lb=lower, ub =upper, c_lambda=10^{2}, c_kappa=10^{2})
	if(m==1){
		write.table(Est$THETA,"SEIRD_model/SEIRD_results/A3-05_traj.txt")
	}
	res1[m,]<-Est$THETA[nrow(Est$THETA),]       	      
	res2[m,]<-apply(Est$THETA[-burnIn,],2,mean)	      
	Est<- PF(obs,theta=res1[m,], N=20000, chi_theta=param,c_lambda=10^{2}, c_kappa=10^{2})
	resR1[m,]<-Est$R #estimated R number
	Est<- PF(obs,theta=res2[m,], N=20000, chi_theta=param,c_lambda=10^{2}, c_kappa=10^{2})
	resR2[m,]<-Est$R #estimated R number
	print(m)
}

write.table(res1,"SEIRD_model/SEIRD_results/A3-05_theta1.txt")
write.table(res2,"SEIRD_model/SEIRD_results/A3-05_theta2.txt")
write.table(resR1,"SEIRD_model/SEIRD_results/A3-05_R1.txt")
write.table(resR2,"SEIRD_model/SEIRD_results/A3-05_R2.txt")

################################################################
# Algorithm 3--alpha=1.1  
################################################################

set.seed(9585)
 
#create matrices to store the results
datalength<-nrow(obs)
d<-7		  
burnIn<-(1:(b.val*datalength))  
res1<-matrix(0,M.val,d)
res2<-matrix(0,M.val,d)
resR1<-matrix(0,M.val,datalength)
resR2<-matrix(0,M.val,datalength)

m.save<-min(M.val, 14)
##generate random seed
seed.vec<-.Random.seed
##running time!
for(m in  1:M.val){		
	set.seed(seed.vec[m])
	Est<-A3_SEIRD(obs,N=N.val, nit=n.it, alpha=1.1,  chi_theta=param,  lb=lower, ub =upper, c_lambda=10^{2}, c_kappa=10^{2})
	if(m==m.save){
		write.table(Est$THETA,"SEIRD_model/SEIRD_results/A3-11_traj.txt")
	}
	res1[m,]<-Est$THETA[nrow(Est$THETA),]       	##last iterate 
	res2[m,]<-apply(Est$THETA[-burnIn,],2,mean)  ## average over all values
	Est<- PF(obs,theta=res1[m,], N=20000, chi_theta=param,c_lambda=10^{2}, c_kappa=10^{2})
	resR1[m,]<-Est$R
	Est<- PF(obs,theta=res2[m,], N=20000, chi_theta=param,c_lambda=10^{2}, c_kappa=10^{2})
	resR2[m,]<-Est$R       
	print(m)
}


write.table(res1,"SEIRD_model/SEIRD_results/A3-11_theta1.txt")
write.table(res2,"SEIRD_model/SEIRD_results/A3-11_theta2.txt")
write.table(resR1,"SEIRD_model/SEIRD_results/A3-11_R1.txt")
write.table(resR2,"SEIRD_model/SEIRD_results/A3-11_R2.txt")






















