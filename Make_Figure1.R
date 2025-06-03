################################################################################################################
### Running this R script generates the plots in Figure 1 of the paper and save them in the folder "Figures".
### Remark: this R scirpt must be placed in you working directory.
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
### Generate and store in the folder "LG_model/LG_results" the results needed to make Figure 1.
################################################################################################################

#number of runs of the algorithms (M.val>1)
M.val<-50   #(M.val=50 in the paper)
#number of data points (n.obs>1)
n.obs<-10^4 #(n.obs=10^4 in the paper)
#number of particles
N.val<-10^4 #(N.val=10^4 in the paper)

source("LG_model/LG_synthetic.R")

################################################################################################################
### Load the results for parameter inference
################################################################################################################

##Results for q=2
q<-2			 	 	 
##theta_star
sigma.y<-0.5
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
sigX_star<-rep(1,q)  
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
d<-length(theta_star)



SS1<-t(apply(read.table("LG_model/LG_results/A1-05-q2_theta.txt"),1,as.numeric))
error1<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error1[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A1-11-q2_theta.txt"),1,as.numeric))
error2<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error2[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-04-q2_theta.txt"),1,as.numeric))
error3<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error3[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-05-q2_theta.txt"),1,as.numeric))
error4<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error4[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-06-q2_theta.txt"),1,as.numeric))
error5<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error5[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-11-q2_theta.txt"),1,as.numeric))
error6<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error6[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}


##Results for q=4
q<-4					
##define theta_star
sigma.y<-0.5
set.seed(123)
beta_star<-runif(q,-2,2)
rho_star<-runif(q,-1,1)
sigX_star<-rep(1,q)  
theta_star<-c(beta_star,rho_star,sigX_star,sigma.y)
d<-length(theta_star)


SS1<-t(apply(read.table("LG_model/LG_results/A1-05-q4_theta.txt"),1,as.numeric))
error11<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error11[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A1-11-q4_theta.txt"),1,as.numeric))
error22<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error22[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-04-q4_theta.txt"),1,as.numeric))
error33<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error33[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-05-q4_theta.txt"),1,as.numeric))
error44<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error44[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-06-q4_theta.txt"),1,as.numeric))
error55<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error55[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}
SS1<-t(apply(read.table("LG_model/LG_results/A2-11-q4_theta.txt"),1,as.numeric))
error66<-rep(0,nrow(SS1))
for(i in 1:nrow(SS1)){
	error66[i]<- sqrt(sum(SS1[i,]-theta_star)^2)/d
}


T_end<-length(error1)
##Put all the results in a dataframe
all <- data.frame(
  Algorithm = factor(c(rep('A1-0.5', 2*T_end), rep('A1-1.1',2*T_end),
  			rep('A2-0.4', 2*T_end), rep('A2-0.5',2*T_end),
  			rep('A2-0.6', 2*T_end),
                       rep('A2-1.1', 2*T_end)), 
		levels = c('A1-0.5', 'A1-1.1', 'A2-0.4','A2-0.5','A2-0.6', 'A2-1.1' )),
  		q=c(rep(c(rep("q=2",T_end), rep("q=4",T_end)),6)),
		value = c( error1, error11, error2, error22, error3, error33, error4, error44, error5,error55, error6, error66)
)

################################################################################################################
### Load results for state inference
################################################################################################################

##Results for q=2
error<-read.table("LG_model/LG_results/A1-05-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error1<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A1-11-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error2<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-04-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error3<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-05-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error4<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-06-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error5<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-11-q2_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error6<-apply(error,2,mean)

##Results for q=4
error<-read.table("LG_model/LG_results/A1-05-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error11<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A1-11-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error22<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-04-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error33<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-05-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error44<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-06-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error55<-apply(error,2,mean)

error<-read.table("LG_model/LG_results/A2-11-q4_X.txt")
error<-apply(error,2,cumsum)
error<-t(t(error)/(1:ncol(error)))
error66<-apply(error,2,mean)

##Put all the results  dataframe
T_end<-length(error1)
all2 <- data.frame(
  Algorithm = factor(c(rep('A1-0.5', 2*T_end), rep('A1-1.1', 2*T_end), rep('A2-0.4', 2*T_end), rep('A2-0.5',2*T_end),
  			rep('A2-0.6',2*T_end), rep('A2-1.1',2*T_end)), 
		levels = c('A1-0.5', 'A1-1.1', 'A2-0.4', 'A2-0.5', 'A2-0.6','A2-1.1')),
 		q=c( rep(c(rep("q=2",T_end),rep("q=4",T_end)),6)),
  		value = c(error1,error11, error2, error22, error3, error33, error4, error44, error5, error55, error6, error66),
  		axis= c(rep(1:T_end,12)/1000)
)

################################################################################################################
### Make Figure 1
################################################################################################################

p1<-ggplot(data = all, aes(x=Algorithm, y=value, group=Algorithm)) + geom_boxplot(aes(fill=Algorithm))+theme_bw()+
    labs(y = "estimation error") +xlab(" ")+
    scale_y_log10(limits =  c(0.0001, 10^0.5), breaks = c(10^{-3}, 10^{-2}, 10^{-1}, 10^{0},10),labels = trans_format("log10", math_format(10^.x)))+
    scale_fill_manual(values=c("darkolivegreen","orange", "deeppink", "black","blue",  "red"))+
    theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25), axis.ticks.x=element_blank()) + 
    theme(legend.text=element_text(size=25), legend.title=element_text(size=25))+  scale_x_discrete(labels= c("","","", "","",""))

    
p1<-p1+facet_grid(.~q)
p1<-p1+ theme(strip.text = element_text(size=25))


p2<-ggplot(data=all2,  aes(x=axis, y=value, group=Algorithm,   colour=Algorithm, linetype=Algorithm)) + geom_line(linewidth=1.5)+theme_bw()+
	xlab("time t (in thousands)")+ylab("estimation error")+
	scale_y_log10(limits =  c(0.0001, 15), breaks = c(10^{-3}, 10^{-2}, 10^{-1},10^{0},10),labels = trans_format("log10", math_format(10^.x)))+#,labels = c("","","","",""))+
	scale_x_continuous(breaks = c(1,3,5,7,9))+
	scale_linetype_manual(values=rep(1,8))+scale_color_manual(values=c("darkolivegreen","orange", "deeppink", "black","blue",  "red"))+
        theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25))+
       theme(legend.text=element_text(size=25), legend.title=element_text(size=25))
       

p2<-p2+facet_grid(.~q)
p2<-p2+  theme(strip.text = element_text(size=25)) 

p1<- p1+ theme(legend.position = 'none')
p1<-p1+theme(plot.margin = margin(r=0))


p2<-p2+theme(plot.margin = margin(l=0))



ggsave(p1,file="Figures/Figure1A.pdf", width=6, height=5)
ggsave(p2, file="Figures/Figure1B.pdf", width=7, height=5)


