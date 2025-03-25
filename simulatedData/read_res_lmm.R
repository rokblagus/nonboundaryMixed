
library(ggplot2)
library(gridExtra)


######cnvergence params

##binom
lim_sing<-1e-5 

lim_sing_rho<-1e-5 

lim_sing_s<--log(sqrt(1e-5))


########################


###func to sum res

ent_loss<-function(true,est){
  q<-ncol(true)
  SS<-solve(true)%*%est
  sum(diag(SS))-log(det(SS))-q
}

get_loss<-function(Sigma,ests){
  Sigma_est<-matrix(c(ests[1],ests[2],ests[2],ests[3]),2,2)
  ent_loss(Sigma,Sigma_est)
}

min_eigen<-function(vr1,cvr,vr2){
  if (is.na(vr1)) NA else {
    Sigma<-matrix(c(vr1,cvr,cvr,vr2),2,2)
    min(eigen(Sigma)$values)
  }
}

mean_eigen<-function(vr1,cvr,vr2){
  if (is.na(vr1)) NA else {
    Sigma<-matrix(c(vr1,cvr,cvr,vr2),2,2)
    mean(eigen(Sigma)$values)
  }
}

get_bias<-function(est,true){
  mean(est)-true
}

get_rmse<-function(est,true){
  sqrt( mean( (est-true)**2 )  )
}
dfmsi<-NULL

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="lmm"
          
          
          for (meth in c("ml","reml","bml", "tml","to","o")){
            
            
            sd_int<-sqrt(vr_int)  
            multiplier_int<-mult
            
            sd_slope<-sd_int*multiplier_int
            
            b0 <- -0.1
            
            # regression coefficients of fixed effects
            betas <- c(0.17, -0.17, 0.17, 0.15, -0.1)
            
            
            rho=rh
            
            # list of covariance matrices for random effects
            corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2) 
            true_sig<-corr_mat
            
            sigma_eps=2
            
            alpha=0.05
            
            
            nms<-paste("lmm",N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
            
            nmsfile<-paste0("lmm/results/",meth,nms,".txt")
            
            dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)
            
            if (class(dd)!="try-error"){
              zz=zz+1
              
              nsim<-nrow(dd)
              
              ddm<-as.matrix(dd)
              
              is_sing_s<-function(v1,v2,tol){
                if (abs(log(sqrt(v1)))>tol|abs(log(sqrt(v2)))>tol) 1 else 0
              }
              
              sing_s<-loss<-min_e<-mean_e<-rep(NA,nsim)
              for (i in 1:nsim){
                loss[i]<-get_loss(true_sig,ddm[i,7:9])
                min_e[i]<- min_eigen(ddm[i,7],ddm[i,8],ddm[i,9])
                mean_e[i]<- mean_eigen(ddm[i,7],ddm[i,8],ddm[i,9])
                sing_s[i]<-is_sing_s(ddm[i,7],ddm[i,9],lim_sing_s)
              }
              if (meth=="ml"){  
                dfmse<-apply(ddm[,23:28],2,mean)
              if(is.null(dfmsi)) {dfmsi<-dfmse} else {dfmsi<-rbind(dfmsi,dfmse)}  
                }
              dfi<-data.frame(
                N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,
                method=meth,
                mloss=mean(loss,na.rm=TRUE),
                na_mloss=mean(is.na(loss)),
                singular=mean(min_e<lim_sing,na.rm=TRUE),
                na_singular=mean(is.na(min_e)),
                mean_eigen=mean(mean_e),
                na_mean_eigen=mean(is.na(mean_e)),
                covr=mean(ddm[,16],na.rm=TRUE),
                na_covr=mean(is.na(ddm[,16])),nsim=nsim,
                bound_rho=mean(1-abs(ddm[,10])<lim_sing_rho,na.rm=TRUE),
                na_bound_rho=mean(is.na(ddm[,10])),
                bound_s=mean(sing_s,na.rm=TRUE),
                na_bound_s=mean(is.na(sing_s)),
                bound_t=mean(ifelse(1-abs(ddm[,10])<lim_sing_rho|sing_s==1,1,0),na.rm=TRUE),
                na_bound_t=mean(is.na(ddm[,10])|is.na(sing_s)) ,na_fit=mean(is.na(ddm[,1])),
                bias_b5=get_bias(ddm[,6],-0.1),rmse_b5=get_rmse(ddm[,6],-0.1)
              )
              
              if (zz==1) df<-dfi else df<-rbind(df,dfi)
              
            }
            
            
          }
        }
      }}}}

#
split(df$na_fit,df$method) #  always zero, hence no issue

 

##get taus

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="lmm"
          
          
          
          
          sd_int<-sqrt(vr_int)  
          multiplier_int<-mult
          
          sd_slope<-sd_int*multiplier_int
          
          b0 <- -0.1
          
          # regression coefficients of fixed effects
          betas <- c(0.17, -0.17, 0.17, 0.15, -0.05)
          
          
          rho=rh
          
          # list of covariance matrices for random effects
          corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2) 
          true_sig<-corr_mat
          
          sigma_eps=2
          
          alpha=0.05
          
          
          nms<-paste("lmm",N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
          
          nmsfile<-paste0("lmm/results/taus",nms,".txt")
          
          dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)
          
          if (class(dd)!="try-error"){
            zz=zz+1
            
            nsim<-nrow(dd)
            
            ddm<-as.matrix(dd)
            
            taus<-apply(ddm,2,mean)
            
            dfit<-data.frame(
              N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,
              
              tauml= taus[1],
              
              taoracle=taus[1],nsim=nsim
            )
            
            if (zz==1) dft<-dfit else dft<-rbind(dft,dfit)
            
          }
          
          
        }
      }
    }}}


 

txt1<-"sigma[0]^2"
txt2<-"sigma[1]^2"

df$Nn<-paste0("N=",df$N,",n=",df$n)
df$vars<-paste0(txt1,"==",df$vr_int,"*','*~",txt2,"==",df$vr_int*df$mult)

dft$Nn<-paste0("N=",dft$N,",n=",dft$n)
dft$vars<-paste0(txt1,"==",dft$vr_int,"*','*~",txt2,"==",dft$vr_int*dft$mult)

dft$Nn2<-as.numeric(as.factor(dft$Nn))


dft$r2<-paste0("rho==",dft$rho)

 
df$Nn2<-as.numeric(as.factor(df$Nn))


df$r2<-paste0("rho==",df$rho)

 



###reduced df (only for PRIAL)

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="bin"
          
          
          for (meth in c("reml","bml", "tml","to", "o")){
            zz=zz+1
            dfi<-df[df$N==N&df$n==n&df$vr_int==vr_int&df$mult==mult&df$rho==rh&df$method%in%c(meth,"ml"),c(1:7,3+(22:25))]
            
            prial<-(dfi$mloss[dfi$method=="ml"]-dfi$mloss[dfi$method==meth])/dfi$mloss[dfi$method=="ml"]*100
            
            dfii<-dfi[dfi$method==meth,]
            dfii$prial<-prial
            
            if (zz==1) df2<-dfii else df2<-rbind(df2,dfii)
          }
        }
      }}}}



 
################final figs as for the paper

#rename methods as in paper

df$method<-as.character(df$method)
 
df$method[df$method=="ml"]<-"ML"
df$method[df$method=="reml"]<-"REML"
df$method[df$method=="bml"]<-"BM"

df$method[df$method=="tml"]<-"PML"
df$method[df$method=="to"]<-"POR"
df$method[df$method=="o"]<-"ORACLE"

 
df$method<-factor(df$method,levels=c("ML","REML","BM","PML","POR","ORACLE"))

df2$method<-as.character(df2$method)
 
df2$method[df2$method=="reml"]<-"REML"
df2$method[df2$method=="bml"]<-"BM"

df2$method[df2$method=="tml"]<-"PML"
df2$method[df2$method=="to"]<-"POR"

df2$method[df2$method=="o"]<-"ORACLE"

 
df2$method<-factor(df2$method,levels=c("REML","BM", "PML","POR","ORACLE"))



dfpallete<-c("black","red","cadetblue","blue","blueviolet","deeppink")
dfpallete2<-c("red","cadetblue","blue","blueviolet","deeppink")
dfpallete3<-c("red","blue","blueviolet","deeppink")

#supplement: boundary est all criteria, b5, loss(?)

ploss<-ggplot(df)+
  geom_point(aes(x=Nn2,y=log(mloss,base=10),colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=log(mloss,base=10),colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab(expression(log[10](loss)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("lmm/figs/supp_lmm_loss.pdf",height=6,width=12)
ploss
dev.off()


p1<-ggplot(df)+
  geom_point(aes(x=Nn2,y=bound_rho*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_rho*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("boundary correlation (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2<-ggplot(df)+
  geom_point(aes(x=Nn2,y=bound_s*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_s*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("boundary standard deviations (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p0<-ggplot(df)+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())



pp<-ggplot(df2[df2$method!="BM",])+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7)+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df2$Nn)))+
  ylab("PRIAL (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete3)+
  theme(legend.title=element_blank()) 

#need to break the axis!
df2$breaks<-rep(NA,nrow(df2))
df2$breaks[df2$prial>(-100)]<-"(-100*','*100)"
df2$breaks[df2$prial<=(-100)&df2$prial>(-1000)]<-"(-1000*','*-100)"
df2$breaks[df2$prial<=(-1000)]<-"(-inf*','*-1000)"
df2$breaks<-factor(df2$breaks,levels=c("(-100*','*100)","(-1000*','*-100)","(-inf*','*-1000)"))
library(tidyverse)
df21<-df2[df2$rho==0.5,]

df21$r2<-factor(df21$r2)
pps1<-ggplot(df21 )+
  geom_hline(data = df21 %>% filter(breaks == "(-100*','*100)"),
             aes(yintercept = 0))+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7)+
  facet_grid(r2+breaks~vars,scales="free_y", labeller = label_parsed,drop=TRUE )+
  theme_light()+
  theme(panel.spacing.y=unit(0, "lines"),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    #guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df21$Nn)))+
  ylab("PRIAL (%)")+
   
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())  


df2$breaks<-rep(NA,nrow(df2))
df2$breaks[df2$prial>(-100)]<-"(-100*','*100)"
df2$breaks[df2$prial<=(-100)&df2$prial>(-2000)]<-"(-2000*','*-100)"
df2$breaks[df2$prial<=(-2000)]<-"(-inf*','*-2000)"
df2$breaks<-factor(df2$breaks,levels=c("(-100*','*100)","(-2000*','*-100)","(-inf*','*-2000)"))

df22<-df2[df2$rho==0.8,]
df22$r2<-factor(df22$r2)
pps2<-ggplot(df22 )+
  geom_hline(data = df22 %>% filter(breaks == "(-100*','*100)"),
             aes(yintercept = 0))+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7)+
  facet_grid(r2+breaks~vars,scales="free_y", labeller = label_parsed,drop=TRUE )+
  theme_light()+
  theme(panel.spacing.y=unit(0, "lines"))+
  #guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df22$Nn)))+
  ylab("PRIAL (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank()) 

library(grid)
library(ggpubr)
fig<-ggarrange(pps1+ rremove("ylab"),pps2+ rremove("ylab"),ncol=1,nrow=2,common.legend = TRUE,legend="right")
fig2<-annotate_figure(fig, left = textGrob("PRIAL (%)", rot = 90, vjust = 1, gp = gpar(cex = 1)) )


pdf("lmm/figs/supp_lmm_boundary.pdf",height=10,width=10)
grid.arrange(p0,p1,p2,nrow=3)
dev.off()

##bias,rmse, cover for beta5

p1b<-ggplot(df)+
  geom_point(aes(x=Nn2,y=bias_b5,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bias_b5,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("bias")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df)+
  geom_point(aes(x=Nn2,y=rmse_b5,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=rmse_b5,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("RMSE")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
   
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b<-ggplot(df)+
  geom_point(aes(x=Nn2,y=covr*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=covr*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)


pdf("lmm/figs/supp_lmm_beta5.pdf",height=10,width=10)
grid.arrange(p1b,p2b,p3b,nrow=3)
dev.off()


#main document, boundary

p0<-ggplot(df)+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(df$Nn)))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("lmm/figs/main_lmm_boundary.pdf",height=6,width=12)
p0
dev.off()


#main document, prial

pdf("lmm/figs/supp_lmm_prial.pdf",height=10,width=12)
fig2
dev.off()

pdf("lmm/figs/main_lmm_prial.pdf",height=6,width=12)
pp
dev.off()
