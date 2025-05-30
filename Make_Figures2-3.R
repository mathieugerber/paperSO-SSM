################################################################################################################
##  Running this R script generates the plots in Figures 2-3 of the paper and save them in the folder "Figures"
##  WARNING: The data used in the paper are not open accessed, and in thus R script data from Kaggle are used
##           The data used below are not fro Bristol and cover another time period. 
## Remark: this R scirpt must be placed in you working directory.
################################################################################################################

##Load packages
library(gtools)
library(ggplot2)
library(scales)
library(reshape2)
library(ggpubr)
library(crch) 					
library(msm) 					
library(Rcpp) 	
library(splines2)
################################################################################################################
### Generate and store in the folder "LG_model/LG_results" the results needed to make Figures 2-3
################################################################################################################

#number of runs of the algorithms (M.val>1)
M.val<- 50 		#(M.val=50 in the paper)
#number of data points
n.obs<-20000 		#(1<n.obs<30000)
#number of particles
N.val<-10^4 		#(N.val=10^4  in the paper)

source("LG_model/LG_real.R")

################################################################################################################
### Load the observations (Temperature is in Kelvin)
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
y<-data[-1]
  

################################################################################################################
### Load predictions obtained for q=3
################################################################################################################
pred_mat<-array(0,c(M.val,length(y),8))
##Model based predictions with p=3######
for(k in 1:8){
 	pred_mat[,,k]<-t(apply(read.table(paste("LG_model/LG_results/Real_q3_resPred_",k,".txt", sep='')),1,as.numeric))
}

################################################################################################################
### Load predictions obtained for q=4
################################################################################################################
pred_mat2<-array(0,c(M.val,length(y),8))
for(k in 1:8){
 	pred_mat2[,,k]<-t(apply(read.table(paste("LG_model/LG_results/Real_q4_resPred_",k,".txt", sep='')),1,as.numeric))
}

################################################################################################################
### Compute prediction errors for q=3
################################################################################################################
v1<-matrix(0,8,length(y))
v2<-matrix(0,8,length(y))
np<-rep(0,8)

for(k in 1:8){
	pred<-pred_mat[,,k]
	m<-nrow(pred)
	obs<-y[-(1:k)] 	#contains values of Y_{1+k}, Y_{2+k}.... (observations predicted by the model) 
	np[k]<- min(length(obs),ncol(pred))		  
	xx<-1:np[k]
	error<-matrix(0,m,length(xx))
	for(i in 1:m){
		error[i,]<-cumsum(abs(pred[i,xx]-obs[xx]))/(1:length(xx))
	}
	#Compute min/max prediction errors
	v1[k,1:np[k]]<-apply(error,2,min)
	v2[k,1:np[k]]<-apply(error,2,max)
}
################################################################################################################
### Compute prediction errors for q=4
################################################################################################################
vv1<-matrix(0,8,length(y))
vv2<-matrix(0,8,length(y))
 
for(k in 1:8){
	pred2<-pred_mat2[,,k]
	m<-nrow(pred2)
	obs<-y[-(1:k)] 	            #contains values of Y_{1+k}, Y_{2+k}.... (observations predicted by the model). 
		
	xx<-1:np[k]
	error0<-matrix(0,m,length(xx))
	for(i in 1:m){
		error0[i,]<-cumsum(abs(pred2[i,xx]-obs[xx]))/(1:length(xx))
	}
	#Compute min/max prediction errors
	vv1[k,1:np[k]]<-apply(error0,2,min)
	vv2[k,1:np[k]]<-apply(error0,2,max)
}

################################################################################################################
### Compute prediction  errors when current temperature is used as prediction
################################################################################################################
error1<-matrix(0,8,length(y))
for(k in 1:8){
	y<-data[-1]		#contains values of Y_{1}, Y_{2}....
	naive1<-y     					 
	obs<-y[-(1:k)] 	#contains values of Y_{1+k}, Y_{2+k}.... 
	xx<-1:np[k]
	error1[k,1:np[k]]<-cumsum(abs(naive1[xx]-obs[xx]))/(1:length(xx))
}


################################################################################################################
### Compute predictions errors when temperature the day before is used as prediction
################################################################################################################
error2<-matrix(0,8,length(y))
np2<-rep(0,8)
for(k in 1:8){
	y<-data[-1]		#contains values of Y_{1}, Y_{2}....
	obs<-y[-(1:k)]		#contains values of Y_{1+k}, Y_{2+k}.... 
	obs2<-y[-(1:(k+24))]	#contains values of Y_{1+k+24}, Y_{2+k+24}....  
	naive2<-obs
	np2[k]<- min(length(obs2),length(naive2)) 
	xx<-1:np2[k]
	error2[k,1:np2[k]]<-cumsum(abs(naive2[xx]-obs2[xx]))/(1:length(xx))
}

################################################################################################################
### Make first row of plots in Figure 2
################################################################################################################

##number of points in the plot
ll<-min(10^3,ncol(v1), ncol(v2), ncol(vv1), ncol(vv2), ncol(error1), ncol(error2))

long.v1<-rep(0,4*ll)
long.v2<-rep(0,4*ll)
long.x<-rep(0,4*ll)
long.k<-rep(0,4*ll)
count<-1
for(k in 1:4){
	ss<-floor(seq(1,np[k], length.out=ll))
	for(i in 1:ll){
		long.v1[count]<-v1[k,ss[i]]
		long.v2[count]<-v2[k,ss[i]]
		long.x[count]<-ss[i] 
		long.k[count]<-paste("k=",k , sep='')
		count<-count+1
	}
}


long.vv1<-rep(0,4*ll)
long.vv2<-rep(0,4*ll)
count<-1;
for(k in 1:4){
	ss<-floor(seq(1,np[k], length.out=ll))
	for(i in 1:ll){
		long.vv1[count]<-vv1[k,ss[i]]
		long.vv2[count]<-vv2[k,ss[i]]
		count<-count+1
	}
}

long.e1<-rep(0,4*ll)
count<-1;
for(k in 1:4){
	ss<-floor(seq(1,np[k], length.out=ll))
	for(i in 1:ll){
		long.e1[count]=error1[k,ss[i]]
		count<-count+1
	}
}


long.e2<-rep(0,4*ll)
long.x2<-rep(0,4*ll)
long.k2<-rep(0,4*ll)
count<-1;
for(k in 1:4){
	x.val2<-24:(23+np2[k]) 
	ss<-floor(seq(1,np2[k], length.out=ll))
	for(i in 1:ll){
		long.e2[count]=error2[k,ss[i]]
		long.x2[count]=x.val2[ss[i]] 
		long.k2[count]=paste("k=",k , sep='')
		count<-count+1
	}
}


df<- data.frame(x=rep(long.x/1000,2), kval=rep(long.k,2), y1=c(long.v1,long.vv1), y2=c(long.v2,long.vv2), 
	Algorithm=factor(c(rep("Algo.2 q=3", length(long.x)),rep("Algo.2 q=4", length(long.x))), levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))
	

p1 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Algorithm, color = Algorithm), alpha=1)+
	xlab("time (in thousands)") +ylab("prediction error")+ylim(0,max(c(max(long.v2),max(long.vv2), max(long.e1), max(long.e2))))+theme_bw()+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("Algo.2 q=3"="green", "Algo.2 q=4"="black", "Previous hour"="red","Previous day"="blue"))+
	scale_color_manual(values=c("Algo.2 q=3"="green", "Algo.2 q=4"="black","Previous hour"="red","Previous day"="blue"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) 


additional_line <- data.frame(x = long.x/1000, yval = long.e1,  kval=long.k, Algorithm=factor(c(rep("Previous hour",length(long.x))),
				levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))
 

p1<-p1+geom_line(data = additional_line, aes(x =x, y = yval,   color = Algorithm),  linewidth = 1, linetype=1,inherit.aes = FALSE)
p1<-p1+geom_ribbon(data = additional_line, aes(x =x, ymin = yval, ymax=yval,  fill = Algorithm), inherit.aes = FALSE)

additional_line2 <- data.frame(x = long.x2/1000, yval = long.e2,  kval=long.k2, Algorithm=factor(c(rep("Previous day",length(long.x2))),
				levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))


p1<-p1+geom_line(data = additional_line2, aes(x =x, y = yval, color = Algorithm), linewidth = 1,inherit.aes = FALSE)
p1<-p1+geom_ribbon(data = additional_line2, aes(x =x, ymin = yval, ymax=yval, fill = Algorithm),inherit.aes = FALSE)
  


p1<-p1+facet_grid(.~kval)
p1<-p1+  theme(strip.text = element_text(size=25))
p1<-p1+xlab(" ")+theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
         axis.ticks.x=element_blank())

  
ggsave(p1,file="Figures/Figure2A.pdf", width=15, height=4)


################################################################################################################
### Make second row of plots in Figure 2
################################################################################################################

long.v1<-rep(0,4*ll)
long.v2<-rep(0,4*ll)
long.x<-rep(0,4*ll)
long.k<-rep(0,4*ll)
start<-1
count<-1
for(k in 5:8){
	ss<-floor(seq(start,np[k], length.out=ll))
	for(i in 1:ll){
		long.v1[count]<-v1[k,ss[i]]
		long.v2[count]<-v2[k,ss[i]]
		long.x[count]<-ss[i] 
		long.k[count]<-paste("k=",k , sep='')
		count<-count+1
	}
}


long.vv1<-rep(0,4*ll)
long.vv2<-rep(0,4*ll)
count<-1;
for(k in 5:8){
	ss<-floor(seq(start,np[k], length.out=ll))
	for(i in 1:ll){
		long.vv1[count]<-vv1[k,ss[i]]
		long.vv2[count]<-vv2[k,ss[i]]
		count<-count+1
	}
}


long.e1<-rep(0,4*ll)
count<-1;
for(k in 5:8){
	ss<-floor(seq(1,np[k], length.out=ll))
	for(i in 1:ll){
		long.e1[count]=error1[k,ss[i]]
		count<-count+1
	}
}

long.e2<-rep(0,4*ll)
long.x2<-rep(0,4*ll)
long.k2<-rep(0,4*ll)
count<-1;
for(k in 5:8){
	x.val2<-24:(23+np2[k]) 
	ss<-floor(seq(1,np2[k], length.out=ll))
	for(i in 1:ll){
		long.e2[count]=error2[k,ss[i]]
		long.x2[count]=x.val2[ss[i]] 
		long.k2[count]=paste("k=",k , sep='')
		count=count+1
	}
}


df<- data.frame(x=rep(long.x/1000,2), kval=rep(long.k,2), y1=c(long.v1,long.vv1), y2=c(long.v2,long.vv2), 
	Algorithm=factor(c(rep("Algo.2 q=3", length(long.x)),rep("Algo.2 q=4", length(long.x))), levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))
	

p2 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Algorithm, color = Algorithm), alpha=1)+
	xlab("time (in thousands)") +ylab("prediction error")+ylim(0,max(c(max(long.v2),max(long.vv2), max(long.e1), max(long.e2))))+theme_bw()+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("Algo.2 q=3"="green", "Algo.2 q=4"="black", "Previous hour"="red","Previous day"="blue"))+
	scale_color_manual(values=c("Algo.2 q=3"="green", "Algo.2 q=4"="black","Previous hour"="red","Previous day"="blue"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25))  


additional_line <- data.frame(x = long.x/1000, yval = long.e1,  kval=long.k, Algorithm=factor(c(rep("Previous hour",length(long.x))),
				levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))
 

p2<-p2+geom_line(data = additional_line, aes(x =x, y = yval,   color = Algorithm),  linewidth = 1, linetype=1,inherit.aes = FALSE)
p2<-p2+geom_ribbon(data = additional_line, aes(x =x, ymin = yval, ymax=yval,  fill = Algorithm), inherit.aes = FALSE)


additional_line2 <- data.frame(x = long.x2/1000, yval = long.e2,  kval=long.k2, Algorithm=factor(c(rep("Previous day",length(long.x2))),
				levels=c("Algo.2 q=3","Algo.2 q=4", "Previous hour", "Previous day")))


p2<-p2+geom_line(data = additional_line2, aes(x =x, y = yval, color = Algorithm), linewidth = 1,inherit.aes = FALSE)
p2<-p2+geom_ribbon(data = additional_line2, aes(x =x, ymin = yval, ymax=yval, fill = Algorithm),inherit.aes = FALSE)
  


p2<-p2+facet_grid(.~kval)
p2<-p2+  theme(strip.text = element_text(size=25))

ggsave(p2,file="Figures/Figure2B.pdf", width=15, height=4)



################################################################################################################
### Make Figure 3
################################################################################################################


theta<-t(apply(read.table("LG_model/LG_results/Real_q4_traj.txt"),1,as.numeric))
change_points<-24*c(30*3, (30*3+30*6*(1:5)))

T_end<-length(c(theta[,1]))

change_points<-change_points[change_points<=T_end]


df<-data.frame(res=c(c(theta[,1]),c(theta[,2]) ,c(theta[,3]),c(theta[,4])),
		 pp=c(rep("j=1",T_end), rep("j=2", T_end), rep("j=3", T_end),rep("j=4", T_end)),
		 axis=rep(1:T_end,4)/1000 )




p1<-ggplot(data=df,  aes(x=axis, y=res)) + geom_line(linewidth=1.5)+
	xlab("time (in thousands)") +ylab(expression(bold(hat(beta)["t,j"]^N))) +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="none")+
	theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), axis.title.y = element_text(vjust =-2))+theme(legend.text=element_text(size=25))+
	scale_linetype_manual(values=c(1, 1,1,1))+scale_color_manual(values=c("black", "red", "blue", "green"))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
	geom_vline(xintercept=change_points/1000, linetype="dashed")
	
p1<-p1+facet_grid(.~pp)
p1<-p1+  theme(strip.text = element_text(size=25))
 


ggsave(p1,file="Figures/Figure3.pdf", width=15.5, height=4.5)
























