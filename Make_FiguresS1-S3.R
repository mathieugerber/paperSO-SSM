################################################################################################################
### Running this R script generates the plots in Figures S1-S3 in the supplementary material of the paper and 
##  save them in the folder "Figures"
##  Remark: this R scirpt must be placed in you working directory.
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

################################################################################################################
### Generate and store in the folder "Urn_model/Urn_results" the results needed to make Figures S1-S3
################################################################################################################

#number of runs of the algorithms (M.val>1)
M.val<-2	#(M.val=50 in the paper)

#number of particles
N.val<-1000	#(N.val=1000  in the paper)

#number of passes through the data for Figures S1-S2 (n.it>1)
n.it<-500	#(n.it=500 in the paper)

#number of passes through the data for Figures S3 (n.it>1)
n.it2<-1500	#(n.it2=1500 in the paper)

source("Urn_model/Urn_synthetic.R")

################################################################################################################
## Figure S1
################################################################################################################


times1<-floor(seq(min(10, n.it),n.it,length.out=1000)*200)
times2<-floor(seq(min(10, n.it),n.it,length.out=1000)*1000)


res1<-array(0,c(M.val,length(times1),3))
for(k in 1:length(times1)){
 	res1[,k,]<-t(apply(read.table(paste("Urn_model/Urn_results/Csmall_Tsmall_05_",k,".txt", sep='')),1,as.numeric))
}
 
res2<-array(0,c(M.val,length(times2),3))
for(k in 1:length(times2)){
 	res2[,k,]<-t(apply(read.table(paste("Urn_model/Urn_results/Csmall_Tlarge_05_",k,".txt", sep='')),1,as.numeric))
}


ub1<-matrix(0,length(times1),3)
lb1<-matrix(0,length(times1),3)
for(k in 1:length(times1)){
	for(j in 1:3) ub1[k,j]<-max(res1[,k,j]) 
	for(j in 1:3) lb1[k,j]<-min(res1[,k,j]) 
}
ss<-1:length(times)
x1<- times1/200


 

ub11<-matrix(0,length(times2),3)
lb11<-matrix(0,length(times2),3)
for(k in 1:length(times2)){
	for(j in 1:3) ub11[k,j]<-max(res2[,k,j]) 
	for(j in 1:3) lb11[k,j]<-min(res2[,k,j]) 
}
ss2<-1:length(times2)
x2<- times2/1000

   
df<- data.frame(x=c(rep(x1[ss], each = 3),rep(x2[ss2], each = 3)),
		Parameter=factor(c(rep(c("j", "k", "r"),length(ss)),rep(c("j", "k", "r"),length(ss2))),
				  levels=c("j", "k", "r")),
		y1=c( c(t(lb1[ss,])), c(t(lb11[ss2,]))), 
		y2=c( c(t(ub1[ss,])), c(t(ub11[ss2,]))), 
		T=c(rep("T=200",3*length(ss)), rep("T=1000",3*length(ss2)))
		)
		

df$T<-factor(df$T, levels=c('T=200', 'T=1000'))	 

p1 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Parameter, color =Parameter,  group =Parameter), alpha=0.5)+
	xlab("# passes through the data") +ylab("estimated value")+ theme_bw()+scale_x_continuous(breaks=c(100,200,300,400,500))+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("j"="blue", "k"="black", "r"="red"))+
	scale_color_manual(values=c("j"="blue", "k"="black", "r"="red"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) +geom_hline(yintercept = 10) 
	
	

	
p1<-p1+facet_grid(.~T)
p1<-p1+  theme(strip.text = element_text(size=25))    


ggsave(p1,file="Figures/FigureS1.pdf", width=15, height=5)	
 

################################################################################################################
## Figure S2
################################################################################################################

times1<-floor(seq(min(10, n.it),n.it,length.out=1000)*200)
times2<-floor(seq(min(10, n.it),n.it,length.out=1000)*1000)


res1<-array(0,c(M.val,length(times1),3))
for(k in 1:length(times1)){
 	res1[,k,]<-t(apply(read.table(paste("Urn_model/Urn_results/Clarge_Tsmall_05_",k,".txt", sep='')),1,as.numeric))
}
 
res2<-array(0,c(M.val,length(times2),3))
for(k in 1:length(times2)){
 	res2[,k,]<-t(apply(read.table(paste("Urn_model/Urn_results/Clarge_Tlarge_05_",k,".txt", sep='')),1,as.numeric))
}
 

ub1<-matrix(0,length(times1),3)
lb1<-matrix(0,length(times1),3)
for(k in 1:length(times1)){
	for(j in 1:3) ub1[k,j]<-max(res1[,k,j]) 
	for(j in 1:3) lb1[k,j]<-min(res1[,k,j]) 
}

ss<-1:length(times1)
x1<- times1/200 


 

ub11<-matrix(0,length(times2),3)
lb11<-matrix(0,length(times2),3)
for(k in 1:length(times2)){
	for(j in 1:3) ub11[k,j]<-max(res2[,k,j]) 
	for(j in 1:3) lb11[k,j]<-min(res2[,k,j]) 
}
ss2<-1:length(times2)
x2<- times2/1000


df<- data.frame(x=c(rep(x1[ss], each = 3),rep(x2[ss2], each = 3)),
		Parameter=factor(c(rep(c("j", "k", "r"),length(ss)),rep(c("j", "k", "r"),length(ss2))),
				  levels=c("j", "k", "r")),
		y1=c( c(t(lb1[ss,])), c(t(lb11[ss2,]))), 
		y2=c( c(t(ub1[ss,])), c(t(ub11[ss2,]))), 
		T=c(rep("T=200",3*length(ss)), rep("T=1000",3*length(ss2)))
		)
		

df$T<-factor(df$T, levels=c('T=200', 'T=1000'))	 

p1 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Parameter, color =Parameter,  group =Parameter), alpha=0.5)+
	xlab("# passes through the data") +ylab("estimated value")+ theme_bw()+scale_x_continuous(breaks=seq(n.it/5,n.it, length.out=5))+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("j"="blue", "k"="black", "r"="red"))+
	scale_color_manual(values=c("j"="blue", "k"="black", "r"="red"))+ geom_hline(yintercept = 80)+ 
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25))  
	
	

	
p1<-p1+facet_grid(.~T)
p1<-p1+  theme(strip.text = element_text(size=25))    


ggsave(p1,file="Figures/FigureS2.pdf", width=15, height=5)	
 

################################################################################################################
## Figure S3
#######################################################

times1<-floor(seq(min(10, n.it2),n.it2,length.out=1000)*200)
times2<-floor(seq(min(10, n.it2),n.it2,length.out=1000)*1000)


res1<-array(0,c(M.val,length(times1),3))
for(k in 1:length(times1)){
 	res1[,k,]<-t(apply(read.table(paste("Urn_model/Urn_results/Csmall_Tsmall_04_",k,".txt", sep='')),1,as.numeric))
}
 

ub1<-matrix(0,length(times1),3)
lb1<-matrix(0,length(times1),3)
for(k in 1:length(times1)){
	for(j in 1:3) ub1[k,j]<-max(res1[,k,j]) 
	for(j in 1:3) lb1[k,j]<-min(res1[,k,j]) 
}
ss<-1:length(times1)
x1<- times1/200


df<- data.frame(x=rep(x1[ss], each = 3),
		Parameter=factor(c(rep(c("j", "k", "r"),length(ss))),
				  levels=c("j", "k", "r")),
		y1=  c(t(lb1[ss,])),
		y2= c(t(ub1[ss,]))
		)
		
p1 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Parameter, color =Parameter,  group =Parameter), alpha=0.5)+
	xlab("# passes through the data") +ylab("estimated value")+ theme_bw()+scale_x_continuous(breaks=c(0,seq(n.it2/3,n.it2, length.out=3)))+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("j"="blue", "k"="black", "r"="red"))+
	scale_color_manual(values=c("j"="blue", "k"="black", "r"="red"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25))  + geom_hline(yintercept = 10)
	
ggsave(p1,file="Figures/FigureS3.pdf", width=15, height=5)	
 	
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

