


sourceCpp("resampling.cpp")
sourceCpp("Urn_model/Sample_fast.cpp")
source("Urn_model/Urn_algorithms.R")
 

 
################################################################
# (T,c_star,alpha)=(200,10,0.5)
################################################################
##Set the seed for the RNG
set.seed(9585)

c_star<-10

j<-c_star		#number balls in urn 1--> need to be set to one for identifiability
k<-c_star		#number balls in urn 2
r<-c_star	 	#number  red balls


theta_star<-c(j,k,r)
datalength<-1000	

 
##upper/lower bound for  theta=(j,k,r)
lower<-rep(1,3)
upper<-rep(200,3)


y<-rep(0,datalength+1) 
t<-1
for(t in 2:(datalength+1)){
	P1<-(j-r+y[t-1])*y[t-1]/(j*k)
	P2<-((r-y[t-1])*y[t-1]+(j-r+y[t-1])*(k-y[t-1]))/(j*k)
	P3<-(r-y[t-1])*(k-y[t-1])/(j*k)
	y[t]<-y[t-1]+sample(c(-1,0,1),1, prob=c(P1,P2,P3))	
}

datalength<-200
y<-y[1:(datalength+1)]
lower[c(2,3)]<-rep(max(y),2)
 
#create matrices to store the results
d<-3	#number of parameters
times<-floor(seq(10,n.it,length.out=1000)*datalength)
res<-array(0,c(M.val,length(times),d))

##running time!
for(m in 1:M.val){		
	Est<-IF_Urn(y, N= N.val, nit=n.it , alpha=rep(0.5,3), alpha2=rep(0.5,3), ch=rep(1,3), beta=0.01, ub=upper, lb=lower)
	res[m,,]<-Est$THETA[times,] 
	print(m)
}

for(k in 1:length(times)){
	 write.table(res[,k,], paste("Urn_model/Urn_results/Csmall_Tsmall_05_",k,".txt", sep=''))
}


################################################################
# (T,c_star,alpha)=(1000,10,0.5)
################################################################
##Set the seed for the RNG
set.seed(9585)

c_star<-10
j<-c_star		#number balls in urn 1--> need to be set to one for identifiability
k<-c_star		#number balls in urn 2
r<-c_star	 	#number  red balls


theta_star<-c(j,k,r)
datalength<-1000	

 
##upper/lower bound for theta=(j,k,r)
lower<-rep(1,3)
upper<-rep(200,3)


y<-rep(0,datalength+1) 
t<-1
for(t in 2:(datalength+1)){
	P1<-(j-r+y[t-1])*y[t-1]/(j*k)
	P2<-((r-y[t-1])*y[t-1]+(j-r+y[t-1])*(k-y[t-1]))/(j*k)
	P3<-(r-y[t-1])*(k-y[t-1])/(j*k)
	y[t]<-y[t-1]+sample(c(-1,0,1),1, prob=c(P1,P2,P3))	
}



lower[c(2,3)]<-rep(max(y),2)
 

#create matrices to store the results
d<-3	#number of parameters
times<-floor(seq(10,n.it,length.out=1000)*datalength)
res<-array(0,c(M.val,length(times),d))

##running time!
for(m in 1:M.val){		
	Est<-IF_Urn(y, N= N.val, nit=n.it , alpha=rep(0.5,3), alpha2=rep(0.5,3), ch=rep(1,3), beta=0.01, ub=upper, lb=lower)
	res[m,,]<-Est$THETA[times,] 
	print(m)
}


for(k in 1:length(times)){
	 write.table(res[,k,], paste("Urn_model/Urn_results/Csmall_Tlarge_05_",k,".txt", sep=''))
}

################################################################
# (T,c_star,alpha)=(200,80,0.5)
################################################################
##Set the seed for the RNG
set.seed(9585)
set.seed(9585)

c_star<-80

j<-c_star		#number balls in urn 1--> need to be set to one for identifiability
k<-c_star		#number balls in urn 2
r<-c_star	 	#number  red balls


theta_star<-c(j,k,r)
datalength<-1000	

 
##upper/lower bound for theta=(j,k,r)
lower<-rep(1,3)
upper<-rep(200,3)


y<-rep(0,datalength+1) 
t<-1
for(t in 2:(datalength+1)){
	P1<-(j-r+y[t-1])*y[t-1]/(j*k)
	P2<-((r-y[t-1])*y[t-1]+(j-r+y[t-1])*(k-y[t-1]))/(j*k)
	P3<-(r-y[t-1])*(k-y[t-1])/(j*k)
	y[t]<-y[t-1]+sample(c(-1,0,1),1, prob=c(P1,P2,P3))	
}

datalength<-200
y<-y[1:(datalength+1)]
lower[c(2,3)]<-rep(max(y),2)
 
#create matrices to store the results
d<-3	#number of parameters
times<-floor(seq(10,n.it,length.out=1000)*datalength)
res<-array(0,c(M.val,length(times),d))

##running time!
for(m in 1:M.val){		
	Est<-IF_Urn(y, N= N.val, nit=n.it , alpha=rep(0.5,3), alpha2=rep(0.5,3), ch=rep(1,3), beta=0.01, ub=upper, lb=lower)
	res[m,,]<-Est$THETA[times,] 
	print(m)
}

for(k in 1:length(times)){
	 write.table(res[,k,], paste("Urn_model/Urn_results/Clarge_Tsmall_05_",k,".txt", sep=''))
}

################################################################
# (T,c_star,alpha)=(1000,80,0.5)
################################################################
##Set the seed for the RNG
set.seed(9585)

c_star<-80
j<-c_star		#number balls in urn 1--> need to be set to one for identifiability
k<-c_star		#number balls in urn 2
r<-c_star	 	#number  red balls


theta_star<-c(j,k,r)
datalength<-1000	

 
##upper/lower bound for the parmeters  
lower<-rep(1,3)
upper<-rep(200,3)


y<-rep(0,datalength+1) 
t<-1
for(t in 2:(datalength+1)){
	P1<-(j-r+y[t-1])*y[t-1]/(j*k)
	P2<-((r-y[t-1])*y[t-1]+(j-r+y[t-1])*(k-y[t-1]))/(j*k)
	P3<-(r-y[t-1])*(k-y[t-1])/(j*k)
	y[t]<-y[t-1]+sample(c(-1,0,1),1, prob=c(P1,P2,P3))	
}



lower[c(2,3)]<-rep(max(y),2)
 

#create matrices to store the results
d<-3	#number of parameters
times<-floor(seq(10,n.it,length.out=1000)*datalength)
res<-array(0,c(M.val,length(times),d))

##running time!
for(m in 1:M.val){		
	Est<-IF_Urn(y, N= N.val, nit=n.it , alpha=rep(0.5,3), alpha2=rep(0.5,3), ch=rep(1,3), beta=0.01, ub=upper, lb=lower)
	res[m,,]<-Est$THETA[times,] 
	print(m)
}


for(k in 1:length(times)){
	 write.table(res[,k,], paste("Urn_model/Urn_results/Clarge_Tlarge_05_",k,".txt", sep=''))
}



################################################################
# (T,c_star,alpha)=(200,10,0.4)
################################################################
##Set the seed for the RNG
set.seed(9585)

c_star<-10

j<-c_star		#number balls in urn 1--> need to be set to one for identifiability
k<-c_star		#number balls in urn 2
r<-c_star	 	#number  red balls


theta_star<-c(j,k,r)
datalength<-1000	

 
##upper/lower bound for the parmeters  
lower<-rep(1,3)
upper<-rep(200,3)


y<-rep(0,datalength+1) 
t<-1
for(t in 2:(datalength+1)){
	P1<-(j-r+y[t-1])*y[t-1]/(j*k)
	P2<-((r-y[t-1])*y[t-1]+(j-r+y[t-1])*(k-y[t-1]))/(j*k)
	P3<-(r-y[t-1])*(k-y[t-1])/(j*k)
	y[t]<-y[t-1]+sample(c(-1,0,1),1, prob=c(P1,P2,P3))	
}

datalength<-200
y<-y[1:(datalength+1)]
lower[c(2,3)]<-rep(max(y),2)
 
#create matrices to store the results
d<-3	#number of parameters
times<-floor(seq(10,n.it2,length.out=1000)*datalength)
res<-array(0,c(M.val,length(times),d))

##running time!
for(m in 1:M.val){		
	Est<-IF_Urn(y, N= N.val, nit=n.it2 , alpha=rep(0.4,3), alpha2=rep(0.4,3), ch=rep(1,3), beta=0.01, ub=upper, lb=lower)
	res[m,,]<-Est$THETA[times,] 
	print(m)
}

for(k in 1:length(times)){
	 write.table(res[,k,], paste("Urn_model/Urn_results/Csmall_Tsmall_04_",k,".txt", sep=''))
}





