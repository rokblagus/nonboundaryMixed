
 
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
    for (vr_int in c(0.001,0.01,0.1)){ #changed to include added results after rev avg 25
if (vr_int==0.001) multi=0.1 else multi=c(0.5,1.5)     
 for (mult in multi){
        for (rh in c(0.5,0.8)){
          mod="pois"
          
          
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
            
            
            nms<-paste("po",N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
            
            nmsfile<-paste0("pois/results/",meth,nms,".txt")
            
            dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)
            
            if (class(dd)!="try-error"){
              zz=zz+1
              
              nsim<-nrow(dd)
              
              ddm<-as.matrix(dd)
              if (meth=="ml"){  
                dfmse<-apply(ddm[,23:28],2,mean)
                if(is.null(dfmsi)) {dfmsi<-dfmse} else {dfmsi<-rbind(dfmsi,dfmse)}  
              }
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

split(df$na_fit,df$method) # always zero, hence no issue



 

txt1<-"sigma[0]^2"
txt2<-"sigma[1]^2"

df$Nn<-paste0("N=",df$N,",n=",df$n)
df$vars<-paste0(txt1,"==",df$vr_int,"*','*~",txt2,"==",df$vr_int*df$mult)


 
df$Nn2<-as.numeric(as.factor(df$Nn))


df$r2<-paste0("rho==",df$rho)

 


###reduced df (only for PRIAL)

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
     for (vr_int in c(0.001,0.01,0.1)){ #changed to include added results after rev avg 25
if (vr_int==0.001) multi=0.1 else multi=c(0.5,1.5)     
 for (mult in multi){
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



df$group<-ifelse(df$N==25,"1","2")
df2$group<-ifelse(df2$N==25,"1","2")

nms<-c("N=25\n n=moderate","N=25\n n=large","N=50\n n=moderate","N=50\n n=large")



#supplement: boundary est all criteria, PRIAL, all methods


p1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_rho*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_rho*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary correlation (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_s*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_s*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary standard deviations (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())



pp<-ggplot(df2,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7)+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("PRIAL (%)")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())

pdf("pois/figs/supp_pois_boundary.pdf",height=10,width=12)
grid.arrange(p0,p1,p2,nrow=3)
dev.off()

ploss<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](loss)))+
  #ylab("loss")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("pois/figs/supp_pois_loss.pdf",height=6,width=12)
ploss
dev.off()

##bias,rmse, cover for beta5

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b5,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bias_b5,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=rmse_b5,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=rmse_b5,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("RMSE")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=covr*100,colour =  method),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)


pdf("pois/figs/supp_pois_beta5.pdf",height=10,width=12)
grid.arrange(p1b,p2b,p3b,nrow=3)
dev.off()

 
#supp  document, prial

pdf("pois/figs/main_pois_prial.pdf",height=6,width=12)
pp
dev.off()

 


 


df_pois<-df
df_pois$model<-"PMM"

df_pois2<-df2
df_pois2$model<-"PMM"

save(df_pois,df_pois2,file="res_pois.Rdata")

