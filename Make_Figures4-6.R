################################################################################################################
### Running this R script generates the plots in Figures 4-6 of the paper and save them in the folder "Figures"
##  Remark: this R scirpt must be placed in you working directory.
################################################################################################################

##Load packages
library(gtools)
library(ggplot2)
library(scales)
library(reshape2)
library(ggforce)
library(ggtext)
library(ggpubr)
library(crch) 					
library(msm) 
library(LaplacesDemon) 
library(Rcpp) 

################################################################################################################
### Generate and store in the folder "SEIRD_model/SEIRD_results" the results needed to make Figures 4-6
################################################################################################################

#number of runs of the algorithms (M.val>1)
M.val<-50	#(M.val=50 in the paper)

#number of particles
N.val<-5000	#(N.val=5000 in the paper)

#number of passes through the data (n.it>1)
n.it<-1000	#(n.it=1000 in the paper)

source("SEIRD_model/SEIRD_real.R")


################################################################################################################
### Load R number from OWID website 
### source: https://docs.owid.io/projects/etl/api/covid/#publications
################################################################################################################

data <- read.csv('Data/covid.csv')
data<-data[data$country=='United Kingdom',]
data<-data[-(1:63),]
data<-data[1:120,]
UK_pop <- 67886004
WTO_R<-data$reproduction_rate
datalength<-length(WTO_R)

################################################################################################################
## Figure 4
################################################################################################################

res1<-t(apply(read.table("SEIRD_model/SEIRD_results/A3-05_traj.txt"),1,as.numeric))
nn<- nrow(res1)
res2<-apply(res1,2,cumsum)/(1:nn)
res3<-t(apply(read.table("SEIRD_model/SEIRD_results/A3-11_traj.txt"),1,as.numeric))
npoints<-min(10000,nn)
xaxis<-rep(0,nn)
for(i in 1:n.it){
	xaxis[(120*(i-1)+1):(120*i)]<-(i-1)+1/(120:1)
}
select<-floor(seq(1,nn,length.out=npoints))
nn<-npoints
res1<-res1[select,]
res2<-res2[select,]
res3<-res3[select,]
xaxis<-xaxis[select]
 
all <- data.frame(
  Algorithm = factor(c(rep('A3-0.5', 7*nn), rep('AA3-0.5', 7*nn),rep('A3-1.1', 7*nn)), 
                   levels = c('A3-0.5', 'AA3-0.5', 'A3-1.1')),
  par=rep( c(rep("eta",nn), rep("gamma", nn), rep("mu", nn),rep("sigma[beta]", nn),
  	 rep("kappa", nn),rep("lambda[1]", nn),rep("lambda[2]", nn)),3),
  value = c(c(res1),c(res2),c(res3)),
  axis	 	= c(rep(xaxis,21)/100)
)

all$par<-factor(all$par, levels=c("eta", "gamma", "mu", "sigma[beta]", "kappa","lambda[1]","lambda[2]"))

p1<-ggplot(data=all,  aes(x=axis, y=value, group=Algorithm,   colour=Algorithm, linetype=Algorithm)) + geom_line(linewidth=1.5)+theme_bw()+
	xlab("number of passes through the data (in hundreds)")+ylab("estimated value")+
	scale_x_continuous(breaks = seq(n.it/500,n.it/100, length.out=5)
)+
	scale_y_log10(
       limits = c(10^-2.7, 1),
       breaks = c(1e-3, 1e-2, 1e-1, 1e0),
       labels = trans_format("log10", math_format(10^.x))
      ) +
	scale_linetype_manual(values=rep(1,8))+scale_color_manual(values=c("black", "red", "blue"))+
        theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25))+
       theme(legend.text=element_text(size=25), legend.title=element_text(size=25))
       
p1<-p1+facet_grid(.~par,labeller = label_parsed)
p1<-p1+  theme(strip.text = element_text(size=25))     



ggsave(p1,file="Figures/Figure4.pdf", width=15, height=5)

    
################################################################################################################
## Figure 5
################################################################################################################

res1<-t(apply(read.table("SEIRD_model/SEIRD_results/A3-05_theta1.txt"),1,as.numeric))
res2<-t(apply(read.table("SEIRD_model/SEIRD_results/A3-05_theta2.txt"),1,as.numeric))
res3<-t(apply(read.table("SEIRD_model/SEIRD_results/A3-11_theta1.txt"),1,as.numeric))


nn<-M.val
all <- data.frame(
  Algorithm = factor(c(rep('A3-0.5', 7*nn), rep('AA3-0.5', 7*nn),rep('A3-1.1', 7*nn)), 
                   levels = c('A3-0.5', 'AA3-0.5', 'A3-1.1')),
  par=rep( c(rep("eta",nn), rep("gamma", nn), rep("mu", nn),rep("sigma[beta]", nn),
  	 rep("kappa", nn),rep("lambda[1]", nn),rep("lambda[2]", nn)),3),
  value = c(c(res1),c(res2),c(res3)),
  axis	 	= c(rep(1,7*nn), rep(2,7*nn),rep(3,7*nn))
)

all$par<-factor(all$par, levels=c("eta", "gamma", "mu", "sigma[beta]", "kappa","lambda[1]","lambda[2]"))
all$axis<-factor(all$axis, levels=c('1', '2', '3'))

p1 <- ggplot(all, aes(x = axis, y = value, fill = Algorithm)) +
  scale_y_log10(
    limits = c(10^-4, 1),
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  geom_boxplot() +
  scale_fill_manual(values = c("black", "red", "blue", "green")) +
  scale_x_discrete(labels = c(
    '1' = "",  # black square
    '2' = "",    # red square
    '3' = ""    # blue square
  )) +
  facet_grid(. ~ par, labeller = label_parsed) +
  theme_bw() +
  labs(title = "", x = "", y = "estimated value") +
  theme(
    axis.text.y = element_text(size = 25, colour = "black"),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    strip.text = element_text(size = 25)
  )
 

circle_data <- data.frame(
  par = factor(c("eta", "gamma", "mu", "sigma[beta]", "kappa","lambda[1]","lambda[2]"), levels = levels(all$par)),  # match facet level
  x0 = rep(2,7),
  y0 = apply(res2,2,mean),
  r = rep(0.3,7)
)
p1 <- p1 + geom_circle(
  data = circle_data,
  aes(x0 = x0, y0 = y0, r = r),
  color = "red",
  inherit.aes = FALSE
)

circle_data <- data.frame(
  par = factor(c("eta", "mu", "kappa", "lambda[1]","lambda[2]"), levels = levels(all$par)),  # match facet level
  x0 = rep(3,5),
  y0 = apply(res3[,c(1,3,5:7)],2,mean),
  r = rep(0.3,5)
)
p1 <- p1 + geom_circle(
  data = circle_data,
  aes(x0 = x0, y0 = y0, r = r),
  color = "blue",
  inherit.aes = FALSE
)

ggsave(p1,file="Figures/Figure5.pdf", width=15.2, height=5)

################################################################################################################
## Figure 6
################################################################################################################

R1<-t(apply( read.table("SEIRD_model/SEIRD_results/A3-05_R1.txt"),1,as.numeric))
R11<-t(apply( read.table("SEIRD_model/SEIRD_results/A3-05_R2.txt"),1,as.numeric))


R2<-t(apply( read.table("SEIRD_model/SEIRD_results/A3-11_R1.txt"),1,as.numeric))

ub1<-rep(0,datalength)
lb1<-rep(0,datalength)
for(k in 1:datalength){
	ub1[k]<-max(R1[,k]) 
	lb1[k]<-min(R1[,k])
}



ub11<-rep(0,datalength)
lb11<-rep(0,datalength)
for(k in 1:datalength){
	ub11[k]<-max(R11[,k]) 
	lb11[k]<-min(R11[,k])
}



ub2<-rep(0,datalength)
lb2<-rep(0,datalength)
for(k in 1:datalength){
	ub2[k]<-max(R2[,k]) 
	lb2[k]<-min(R2[,k])
}


x1<-1:datalength

df<- data.frame(x=rep(x1,2), y1=c(lb1,lb11), y2=c(ub1, ub11),
		Estimate=factor(c(rep("A3-0.5", datalength),rep("AA3-0.5", datalength)),
				levels=c("A3-0.5", "AA3-0.5", "A3-1.1", "OWID"))
		)
		

p1 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Estimate, color = Estimate), alpha=0.5)+
	xlab("date") +ylab("R number")+ theme_bw()+
	scale_y_continuous(limits = c(0.5,3.03), breaks=c(0.5,1,1.5,2,2.5,3))+scale_x_continuous(breaks=c(15,45,75,105), labels=c("18/03", "17/04", "17/05", "16/06"))+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("A3-0.5"="red", "AA3-0.5"="black", "OWID"="green", "A3-1.1"="blue"))+
	scale_color_manual(values=c("A3-0.5"="red", "AA3-0.5"="black","OWID"="green","A3-1.1"="blue"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25))  


additional_line <- data.frame(x = x1, yval = WTO_R, Estimate=factor(c(rep("OWID",datalength)),
				levels=c("A3-0.5", "AA3-0.5", "A3-1.1", "OWID")))
 

p1<-p1+geom_line(data = additional_line, aes(x =x, y = yval,   color = Estimate),  linewidth = 1, linetype=1,inherit.aes = FALSE)
p1<-p1+geom_ribbon(data = additional_line, aes(x =x, ymin = yval, ymax=yval,  fill = Estimate), inherit.aes = FALSE)
p1<- p1+ theme(legend.position = 'none') 

 


x1<-1:datalength
		
df<- data.frame(x=rep(x1,3), y1=c(-ub1,-ub11,lb2), y2=c(-lb1, -lb11,ub2),
		Estimate=factor(c(rep("A3-0.5", datalength),rep("AA3-0.5", datalength),rep("A3-1.1", datalength)),
				levels=c("A3-0.5", "AA3-0.5", "A3-1.1", "OWID"))
		)
		

p2 <- ggplot(df, aes(x = x, ymin = y1, ymax = y2)) +geom_ribbon(aes(fill = Estimate, color = Estimate), alpha=0.5)+
	xlab("date") +ylab("R number")+ theme_bw()+
	scale_y_continuous(limits = c(0,9), breaks=c(1,3,5,7,9))+scale_x_continuous(breaks=c(15,45,75,105), labels=c("18/03", "17/04", "17/05", "16/06"))+
       theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25), legend.text=element_text(size=25))+
	scale_fill_manual(values=c("A3-0.5"="red", "AA3-0.5"="black", "OWID"="green", "A3-1.1"="blue"))+
	scale_color_manual(values=c("A3-0.5"="red", "AA3-0.5"="black","OWID"="green","A3-1.1"="blue"))+
	theme(legend.text=element_text(size=25), legend.title=element_text(size=25))  


additional_line <- data.frame(x = x1, yval = WTO_R, Estimate=factor(c(rep("OWID",datalength)),
				levels=c("A3-0.5", "AA3-0.5", "A3-1.1","AA3-1.1", "OWID")))
 



p2<-p2+geom_line(data = additional_line, aes(x =x, y = yval,   color = Estimate),  linewidth = 1, linetype=1,inherit.aes = FALSE)
p2<-p2+geom_ribbon(data = additional_line, aes(x =x, ymin = yval, ymax=yval,  fill = Estimate), inherit.aes = FALSE)
p2<-p2+ylab("")+theme(plot.margin = margin(l=0))

ggsave(p1,file="Figures/Figure6A.pdf", width=6.5, height=5)	
suppressWarnings(ggsave(p2,file="Figures/Figure6B.pdf", width=8.5, height=5))


 
 






