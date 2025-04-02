##creates tables for the lmm, binom, pois examples as R objects res_lin, res_binom, res_pois

###packages

library(glmmTMB)

library(lme4)

library(xtable)

library(nlme)

library(blme)

library(car)

library(merDeriv)  

 

###aux functions for creating PD

##function to create PD for any param for any link for q>1
#psi is always for the variance!
make_pseudo_data_rand_eigen_general_psi_v3_glmm<-function(psi,nu,const=1e8,param="precision",link_fun=function(x) 1/(1+exp(-x))){
  if (is.null(match.arg(param,c("precision","variance")))) stop("param needs to be one of: precision,variance")
  
  q<-ncol(psi)
  if (param=="precision") cc<-(nu-q-1)/q
  if (param=="variance") cc<-(nu+q+1)/q
  
  cc<-max(c(floor(cc),1))
   true<-psi/cc
  ee<-eigen(true,TRUE)
  ui<-list()
  for (j in 1:q){
    ui[[j]]<-sqrt(ee$values[j])*ee$vectors[,j]
  }
  
   pi<-list()
  
  for (j in 1:length(ui)){
    I<-diag(rep(1,length(ui[[j]])))
       pi[[j]]<-link_fun(I%*%matrix(ui[[j]],ncol=1))
  }
    Y<-unlist(pi)
  
   id<-rep(1:q,each=q)
  n<-rep(const,length(id))
  
  Zi<-matrix(0,ncol=q,nrow=q)
  
  for (j in 1:q){
    Zi[j,j]<-1
  }
  for (j in 1:q){
    if (j==1) Z=Zi else Z<-rbind(Z,Zi)
  }
  
  fact<-cc
  if (fact>1){
    Y<-rep(Y,fact)
    n<-rep(n,fact)
    id<-rep(1:(q*fact),each=q)
    for (j in 1:(q*fact)){
      if (j==1) Z=Zi else Z<-rbind(Z,Zi)
    }
  }
  
    data0<-list(Y=Y,grouping=id,nn=n,Z=Z)
  
   
  list(data=data0) 
}


##for Poisson, q=1

make_pseudo_data_rand_eigen_inter_alpha_beta<-function(alpha,beta,param="psi",const=1e8){
  if (is.null(match.arg(param,c("psi","sigma2","logsigma2")))) stop("param needs to be one of: psi,sigma2,logsigma2")
  
  if (param=="psi") N<-max(c(floor(2*(alpha-1)),1))
  if (param=="sigma2") N<-max(c(floor(2*(alpha+1)),1))
  if (param=="logsigma2") N<-max(c(floor(2*(alpha)),1))
  
  var.int<-beta*2/N
  fact<-N
  
  true=matrix(var.int,ncol=1,nrow=1)
  
  ee<-eigen(true,TRUE)
  u1<-sqrt(ee$values[1])*ee$vectors[,1]
  #u2<-sqrt(ee$values[2])*ee$vectors[,2]
  
  #matrix(u1,ncol=1)%*%matrix(u1,nrow=1)+matrix(u2,ncol=1)%*%matrix(u2,nrow=1)
  
  
  pi0=exp(u1[1])
  
  Y<-rep(c(pi0),fact) #the constant improves the convergence!
  n<-rep(rep(const,1),fact)
  id<-c(1:fact)
  Z<-matrix(rep(1,fact),ncol=1)
  data0<-list(Y=Y,grouping=id,nn=n,Z=Z)
  
  #fit0<-glmer(cbind(Y,nn-Y)~-1+(-1+Z|grouping),data=data0,family=binomial)
  #est.vcv<-VarCorr(fit0)$grouping[1:2,1:2]
  
  list(data=data0)#,fit=fit0,vcv.re=est.vcv)
}



#########################

#########################

####poisson example 

 

grouseticks$HEIGHT_C<-grouseticks$HEIGHT-mean(grouseticks$HEIGHT)

 
X<-model.matrix(TICKS~YEAR+HEIGHT_C,data=grouseticks)
 
Z1<-model.matrix(~1,data=grouseticks)
Z2<-model.matrix(~1,data=grouseticks)
Z3<-model.matrix(~YEAR,data=grouseticks)
Y<-grouseticks$TICKS
grouping1<-as.numeric(grouseticks$BROOD)
grouping2<-as.numeric(grouseticks$INDEX)
grouping3<-as.numeric(grouseticks$LOCATION)

xdf<-list(Y=Y,X=X,Z1=Z1,Z2=Z2,Z3=Z3,grouping1=grouping1,grouping2=grouping2,grouping3=grouping3)
 
fit_glmer<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2)+(-1+Z3|grouping3),data=xdf,family=poisson)
fit_bglmer<-bglmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2)+(-1+Z3|grouping3),data=xdf,family=poisson) #pazi ta p$$#"$ zamenja vrstni red!

 

fit_glmer_r <- glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2)+(-1+Z3|grouping3),data=xdf,family=poisson,REML=TRUE)



###tau
fiter_pois_tau<-function(tau,D_est,xdf){
  q<-ncol(D_est)
  ee<-eigen(D_est)
  ee$values[ee$values<1e-4]<-1e-4
  ee$values[ee$values>1e4]<-1e4
  lm<-mean(ee$values)
  li<-ee$values+tau*(lm-ee$values)
  psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  
  nu=2*q-1
  
  pd2<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) exp(x))
  
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)))
  Z1a<-rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)))
  Z2a<-rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd2$data$Z)))
  Z3a<-rbind(xdf$Z3,pd2$data$Z)
  
  Ya<-c(xdf$Y,pd2$data$Y)
  weightsa<-c(rep(1,length(xdf$Y)),pd2$data$nn)
  
  
  grouping1a<-c(xdf$grouping1,max(xdf$grouping1)+pd2$data$grouping)
  grouping2a<-c(xdf$grouping2,max(xdf$grouping2)+pd2$data$grouping)
  grouping3a<-c(xdf$grouping3,max(xdf$grouping3)+pd2$data$grouping)
  
  Ya2<-floor(Ya*weightsa)
  offset<-log(weightsa)
  
  
  xdfa<-list(Y=Ya2,ofset=offset,X=Xa,Z1=Z1a,Z2=Z2a,Z3=Z3a,grouping1=grouping1a,grouping2=grouping2a,grouping3=grouping3a)
  
  
  tmp2 <- glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2)+(-1+Z3|grouping3), family = poisson,
                  offset=ofset,
                  data=xdfa)
  tmp2
}

get_marLik_pois<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2)+(-1+Z3|grouping3),data=xdf,family=poisson,
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                             theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) #point estimates seem ok, but logLik is NA! They have a trick where they dont want to report loogLik if the model does not converge (which in our case defacto holds), but we can still accesss it via object$fit$objective which seems to give -loglik so it should be minimized
  
  
  -tmp3$fit$objective
  
}

tau_finder_pois<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_pois_tau(tau,D_est,xdf)
  abs(get_marLik_pois(fit_tau,xdf)-get_marLik_pois(fit_ml,xdf))-qchisq(1-alpha,1)/2
}


##ML

D<-VarCorr(fit_glmer)$cond$grouping3[1:3,1:3]


sek<-seq(from=0,to=1,by=0.1)
y<-rep(NA,length(sek))
zz=0
for (i in sek){
  zz=zz+1
  y[zz]<-tau_finder_pois(i,xdf,D,fit_glmer)
}
plot(sek,y,type="l")
abline(h=0)

opt_tau_ml<-uniroot(tau_finder_pois,c(0,1),xdf=xdf,D_est=D,
                 fit_ml=fit_glmer,alpha=0.05)

fit_tau_ml<-fiter_pois_tau(opt_tau_ml$root,D_est=D,xdf=xdf)


#summarize the res



mod<-fit_tau_ml 

#sd for 1st RE
vr<-mod$sdr$cov.fixed[5,5]
est<-mod$sdr$par.fixed[5]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[7:9,7:9]
est<-mod$sdr$par.fixed[7:9]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[10:12,10:12]
est<-mod$sdr$par.fixed[10:12]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:4,1:4]
est<-mod$sdr$par.fixed[1:4]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")")
)

resp2_ml<-resi  


mod<-fit_glmer

#sd for 1st RE
vr<-mod$sdr$cov.fixed[5,5]
est<-mod$sdr$par.fixed[5]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[7:9,7:9]
est<-mod$sdr$par.fixed[7:9]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[10:12,10:12]
est<-mod$sdr$par.fixed[10:12]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:4,1:4]
est<-mod$sdr$par.fixed[1:4]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")")
)

resm<-resi


mod<-fit_glmer_r

#sd for 1st RE
vr<-mod$sdr$cov.fixed[5-4,5-4]
est<-mod$sdr$par.fixed[5-4]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6-4,6-4]
est<-mod$sdr$par.fixed[6-4]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[(7:9)-4,(7:9)-4]
est<-mod$sdr$par.fixed[(7:9)-4]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[(10:12)-4,(10:12)-4]
est<-mod$sdr$par.fixed[(10:12)-4]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef


vr<- vcov(mod)$cond 
est<-fixef(mod)$cond

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")")
)

resm_reml<-resi


#bglmer
res.bglmer.cf<-fit_bglmer@optinfo$val
h<-fit_bglmer@optinfo$derivs$Hessian
h<-solve(h)
v<-forceSymmetric(h+t(h))

res.bglmer.se<-sqrt(diag(v))
res.bglmer<-paste(round(res.bglmer.cf,3),paste(" (",round(res.bglmer.se,3),")",sep=""),se="")

th<-fit_bglmer@optinfo$val[3:8]
h<-fit_bglmer@optinfo$derivs$Hessian
h<-solve(h)
v<-forceSymmetric(h+t(h))
v<-v[3:8,3:8]
names(th)<-colnames(v)<-rownames(v)<-paste0("theta",1:6)



s11<-deltaMethod(th,"sqrt(theta1*theta1)",v)
r12<-deltaMethod(th,"(theta1*theta2)/sqrt(theta1**2*(theta2**2+theta4**2)  )   ",v)
r13<-deltaMethod(th,"(theta1*theta3)/sqrt(theta1**2*(theta3**2+theta5**2+theta6**2)  )   ",v)
s22<-deltaMethod(th,"sqrt(theta2**2+theta4**2)",v)
s33<-deltaMethod(th,"sqrt(theta3**2+theta5**2+theta6**2)",v)

r23<-deltaMethod(th,"(theta2*theta3+theta4*theta5)/sqrt((theta2**2+theta4**2)*(theta3**2+theta5**2+theta6**2)  )   ",v)

res.s3.bml<-c(
  paste0(round(s11[1],3), " (",round(s11[2],3),")"   ),
  paste0(round(s22[1],3), " (",round(s22[2],3),")"   ),
  paste0(round(s33[1],3), " (",round(s33[2],3),")"   ),
  paste0(round(r12[1],3), " (",round(r12[2],3),")"   ),
  paste0(round(r13[1],3), " (",round(r13[2],3),")"   ),
  paste0(round(r23[1],3), " (",round(r23[2],3),")"   )
)

resb<-c(res.bglmer[c(9:12,1:2)],res.s3.bml)


res<-cbind(resm,resm_reml,resb[c(1:4,6,5,7:12)],resp2_ml)

colnames(res)<-c("ML","REML","BM",paste0("PML(",round(opt_tau_ml$root,2),")"))

nmr<-c("Intercept","Year$[$1996$]$","Year$[$1997$]$","Height","Brood","Index","Location",rep("",5))

res<-cbind(nmr,res)

rownm<-c("$\\beta_0$","$\\beta_1$","$\\beta_2$","$\\beta_3$","$\\sigma_1$","$\\sigma_2$",
         "$\\sigma_{3,1}$","$\\sigma_{3,2}$","$\\sigma_{3,3}$","$\\rho_{12}$","$\\rho_{13}$","$\\rho_{23}$")

rownames(res)<-rownm
res_pois<-res

print(xtable(res),sanitize.text.function=function(x){x})

 



######################

######################


####linear example: use aids data with log transformed Y


data(aids, package = "JM")
aids$CD4[aids$CD4==0]<-1e-8

aids$CD4<-log(aids$CD4)





 
Y<-aids$CD4
X<-model.matrix(CD4~gender+prevOI+AZT+obstime + drug + obstime:drug,data=aids)
Z<-model.matrix(CD4~1 + obstime,data=aids)
grouping<-as.numeric(as.factor(aids$patient))

xdf<-list(Y=Y,X=X,Z=Z,grouping=grouping)

fit_glmer<-glmmTMB(Y~X-1 +
                     (Z-1 | grouping), data = xdf,family = gaussian(link="identity"))


fit_bglmer<-blmer(Y~X-1 +
                    (Z-1 | grouping), data = xdf)


fit_glmer_r<-glmmTMB(Y~X-1 +
                     (Z-1 | grouping), data = xdf,family = gaussian(link="identity"),REML=TRUE)

 
##using tau
fiter_lin_tau<-function(tau,D_est,xdf){
  q<-ncol(D_est)
  ee<-eigen(D_est)
  ee$values[ee$values<1e-4]<-1e-4
  ee$values[ee$values>1e4]<-1e4
  lm<-mean(ee$values)
  li<-ee$values+tau*(lm-ee$values)
  psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  
  nu=2*q-1
  
  d2222<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi=psi,nu=nu,const=1e8,param="variance",link_fun = function(x) x )
  
  
  
  Ya<-c(xdf$Y,d2222$data$Y)
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(d2222$data$Z)))
  Za<-rbind(xdf$Z,d2222$data$Z)
  groupa<-c(xdf$grouping,max(xdf$grouping)+d2222$data$grouping)
   weights<-c(rep(1,length(xdf$Y)),d2222$data$nn)
  
  xdfa<-list(Y=Ya,X=Xa,Z=Za,grouping=groupa,weights=weights)
  
  tmp2 <- glmmTMB(Y~-1+X+(-1+Z|grouping), family = gaussian(link = "identity"),
                  dispformula = ~offset(-log(weights)),
                  data=xdfa)
  tmp2
}

get_marLik_lmm<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family=gaussian(link = "identity"),
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                           betad=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="betad")],
                           theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) #point estimates seem ok, but logLik is NA! They have a trick where they dont want to report loogLik if the model does not converge (which in our case defacto holds), but we can still accesss it via object$fit$objective which seems to give -loglik so it should be minimized
  
  
  -tmp3$fit$objective
  
}

tau_finder<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_lin_tau(tau,D_est,xdf)
  abs(get_marLik_lmm(fit_tau,xdf)-get_marLik_lmm(fit_ml,xdf))-qchisq(1-alpha,1)/2
}
 
#ML
D<-VarCorr(fit_glmer)$cond$group[1:2,1:2]
 
sek<-seq(from=0,to=0.5,by=0.05)
y<-rep(NA,length(sek))
zz=0
for (i in sek){
  zz=zz+1
  y[zz]<-tau_finder(i,xdf,VarCorr(fit_glmer)$cond$group[1:2,1:2],fit_glmer)
}
plot(sek,y,type="l")
abline(h=0)

opt_tau_ml<-uniroot(tau_finder,c(0,1),xdf=xdf,D_est=VarCorr(fit_glmer)$cond$group[1:2,1:2],
                 fit_ml=fit_glmer,alpha=0.05)

fit_tau_ml<-fiter_lin_tau(opt_tau_ml$root,D_est=VarCorr(fit_glmer)$cond$group[1:2,1:2],xdf=xdf)


##ml

mod<-fit_glmer

#for  RE
#diags
vr<-mod$sdr$cov.fixed[9:10,9:10]
est<-mod$sdr$par.fixed[9:10]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[11,11]
est<-mod$sdr$par.fixed[11]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#for sd: we need sqrt(exp(betad))#in the new release we need exp(betad)!

vr<-mod$sdr$cov.fixed[8,8]
est<-mod$sdr$par.fixed[8]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"sqrt(exp(theta))",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:7,1:7]
est<-mod$sdr$par.fixed[1:7]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

resm<-resi


##REML

mod<-fit_glmer_r

#for  RE
#diags
vr<-mod$sdr$cov.fixed[(9:10)-7,(9:10)-7]
est<-mod$sdr$par.fixed[(9:10)-7]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[11-7,11-7]
est<-mod$sdr$par.fixed[11-7]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#for sd: we need sqrt(exp(betad))#in the new release we need exp(betad)!

vr<-mod$sdr$cov.fixed[8-7,8-7]
est<-mod$sdr$par.fixed[8-7]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"sqrt(exp(theta))",vr)

#fixef
vr<-vcov(mod)$cond
est<-fixef(mod)$cond

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

resm_reml<-resi



##pml

mod<-fit_tau_ml 

#for  RE
#diags
vr<-mod$sdr$cov.fixed[9:10,9:10]
est<-mod$sdr$par.fixed[9:10]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[11,11]
est<-mod$sdr$par.fixed[11]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#for sd: we need sqrt(exp(betad))#in the new release we need exp(betad)

vr<-mod$sdr$cov.fixed[8,8]
est<-mod$sdr$par.fixed[8]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"sqrt(exp(theta))",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:7,1:7]
est<-mod$sdr$par.fixed[1:7]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

resp2_ml<-resi  

#bml

#bglmer
mod<-fit_bglmer
vv <- vcov(mod, full = TRUE)

#fixef

fixf<-paste(round(fixef(mod),3)," (", round(sqrt(diag(vv))[1:7],3),")",sep="")

#ranef

vr<-c(VarCorr(mod)$grouping)[-2]


vrv<-vv[8:10,8:10]

names(vr)<-colnames(vrv)<-rownames(vrv)<-c("sigma11","sigma12","sigma22")

s11<-deltaMethod(vr,"sqrt(sigma11)",vrv)
s22<-deltaMethod(vr,"sqrt(sigma22)",vrv)
r1<-deltaMethod(vr,"sigma12/sqrt(sigma11*sigma22)",vrv)

rnf<-c(paste0(round(s11$Estimate,3)," (",round(s11$SE,3),")"),
       paste0(round(s22$Estimate,3)," (",round(s22$SE,3),")"),
       paste0(round(r1$Estimate,3)," (",round(r1$SE,3),")")
)


#residual

st<-summary(mod)$sigma**2
vst<-vv[11,11]
names(st)<-names(vst)<-"sigma1"

s11<-deltaMethod(st,"sqrt(sigma1)",vst)

rs<-paste0(round(s11$Estimate,3)," (",round(s11$SE,3),")")

resi<-c(fixf,rnf,rs)

res_bglmer<-resi


res<-cbind(resm,resm_reml,res_bglmer, resp2_ml )
colnames(res)<-c("ML","REML","BM", paste0("PML(",round(opt_tau_ml$root,2),")") )

vrnm<-c("Intercept","Gender$[$male$]$","prevOI$[$AID$]$","AZT$[$failure$]$","Obstime","Drug$[$ddI$]$","Obstime$\\times$ Drug$[$ddI$]$",
        "Patient",rep("",3))

res<-cbind(vrnm,res)


rownm<-c("$\\beta_0$","$\\beta_1$","$\\beta_2$","$\\beta_3$","$\\beta_4$","$\\beta_5$","$\\beta_6$",
         "$\\sigma_1$","$\\sigma_2$",
         "$\\rho$","$\\sigma_\\epsilon$")

rownames(res)<-rownm
res_lin<-res

print(xtable(res),sanitize.text.function=function(x){x})




#########################

#########################

 
##binom, math data


dd<-read.csv2("MathEdataset.csv")

dd$Student.ID2<-paste0(dd$Student.ID,dd$Student.Country)

 
Y<-dd$Type.of.Answer
X<-model.matrix(Type.of.Answer~Question.Level+Topic,data=dd)
Z1<-model.matrix(~Question.Level,data=dd)
Z2<-model.matrix(~Question.Level,data=dd)
grouping1<-as.numeric(as.factor(dd$Student.Country))
grouping2<-as.numeric(as.factor(dd$Student.ID2))

xdf<-list(Y=Y,X=X,Z1=Z1,Z2=Z2,grouping1=grouping1,grouping2=grouping2)

fit_glmer<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family=binomial(link="logit"))
fit_bglmer<-bglmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family=binomial(link="logit")) #obrne REs!

 fit_glmer_r<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family=binomial(link="logit"),REML=TRUE)



##using tau
fiter_binom_tau<-function(tau,D_est,xdf){
  q<-ncol(D_est)
  ee<-eigen(D_est)
  ee$values[ee$values<1e-4]<-1e-4
  ee$values[ee$values>1e4]<-1e4
  lm<-mean(ee$values)
  li<-ee$values+tau*(lm-ee$values)
  psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  
  nu=2*q-1
  
  pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  
  
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
  Z1a<-rbind(xdf$Z1,pd1$data$Z)
  Z2a<-rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
  
  
  Ya<-c(xdf$Y,pd1$data$Y)
  weightsa<-c(rep(1,length(xdf$Y)),pd1$data$nn)
  
  
  grouping1a<-c(xdf$grouping1,max(xdf$grouping1)+pd1$data$grouping)
  grouping2a<-c(xdf$grouping2,max(xdf$grouping2)+pd1$data$grouping)
  
  
  
  
  xdfa<-list(Y=Ya,weights=weightsa,X=Xa,Z1=Z1a,Z2=Z2a,grouping1=grouping1a,grouping2=grouping2a)
  tmp2 <-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),weights = weights,data=xdfa,family=binomial(link="logit"))
  
  tmp2
}

get_marLik_binom<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family=binomial(link="logit"),
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                           
                           theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) #point estimates seem ok, but logLik is NA! They have a trick where they dont want to report loogLik if the model does not converge (which in our case defacto holds), but we can still accesss it via object$fit$objective which seems to give -loglik so it should be minimized
  
  
  -tmp3$fit$objective
  
}

tau_finder_binom<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_binom_tau(tau,D_est,xdf)
  abs(get_marLik_binom(fit_tau,xdf)-get_marLik_binom(fit_ml,xdf))-qchisq(1-alpha,1)/2
}


###ML
D1<-VarCorr(fit_glmer)$cond$grouping1[1:2,1:2]

ee<-eigen(D1)
ee$values[ee$values<1e-4]<-1e-4
ee$values[ee$values>1e4]<-1e4

l1<-mean(ee$values) 

psi_stein1<-diag(l1,2,2)*6




nu=2*2-1

pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi=psi_stein1,nu=nu,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))


Xa<-rbind(X,matrix(0,ncol=ncol(X),nrow=nrow(pd1$data$Z)))
Z1a<-rbind(Z1,pd1$data$Z)
Z2a<-rbind(Z2,matrix(0,ncol=ncol(Z2),nrow=nrow(pd1$data$Z)))


Ya<-c(Y,pd1$data$Y)
weightsa<-c(rep(1,length(Y)),pd1$data$nn)


grouping1a<-c(grouping1,max(grouping1)+pd1$data$grouping)
grouping2a<-c(grouping2,max(grouping2)+pd1$data$grouping)




xdfa<-list(Y=Ya,weights=weightsa,X=Xa,Z1=Z1a,Z2=Z2a,grouping1=grouping1a,grouping2=grouping2a)


fit_pen_ml <-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),weights = weights,data=xdfa,family=binomial(link="logit"))

sek<-seq(from=0,to=1,by=0.1)
y<-rep(NA,length(sek))
zz=0
for (i in sek){
  zz=zz+1
  y[zz]<-tau_finder_binom(i,xdf,D1,fit_glmer)
}
plot(sek,y,type="l")
abline(h=0)
#opt_tau can be set to 1 which gives the same res!
opt_tau_ml$root<-1


#####get res

mod<-fit_glmer

#for  RE
#diags
vr<-mod$sdr$cov.fixed[16:17,16:17]
est<-mod$sdr$par.fixed[16:17]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[18,18]
est<-mod$sdr$par.fixed[18]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re1<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


vr<-mod$sdr$cov.fixed[19:20,19:20]
est<-mod$sdr$par.fixed[19:20]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[21,21]
est<-mod$sdr$par.fixed[21]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re2<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


#fixef
vr<-mod$sdr$cov.fixed[1:15,1:15]
est<-mod$sdr$par.fixed[1:15]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,re1,re2)


resm<-resi

#reml

mod<-fit_glmer_r

#for  RE
#diags
vr<-mod$sdr$cov.fixed[(16:17)-15,(16:17)-15]
est<-mod$sdr$par.fixed[(16:17)-15]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[18-15,18-15]
est<-mod$sdr$par.fixed[18-15]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re1<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


vr<-mod$sdr$cov.fixed[(19:20)-15,(19:20)-15]
est<-mod$sdr$par.fixed[(19:20)-15]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[21-15,21-15]
est<-mod$sdr$par.fixed[21-15]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re2<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


#fixef
vr<-vcov(mod)$cond
est<-fixef(mod)$cond

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,re1,re2)


resm_reml<-resi


###
 

mod<-fit_pen_ml 

#for  RE
#diags
vr<-mod$sdr$cov.fixed[16:17,16:17]
est<-mod$sdr$par.fixed[16:17]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[18,18]
est<-mod$sdr$par.fixed[18]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re1<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


vr<-mod$sdr$cov.fixed[19:20,19:20]
est<-mod$sdr$par.fixed[19:20]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[21,21]
est<-mod$sdr$par.fixed[21]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

re2<-c(
  paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
  paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
  paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
)


#fixef
vr<-mod$sdr$cov.fixed[1:15,1:15]
est<-mod$sdr$par.fixed[1:15]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi<-c(fxf,re1,re2)


resp_ml<-resi  



#bglmer
res.bglmer.cf<-fit_bglmer@optinfo$val
h<-fit_bglmer@optinfo$derivs$Hessian
h<-solve(h)
v<-forceSymmetric(h+t(h))

res.bglmer.se<-sqrt(diag(v))
res.bglmer<-paste(round(res.bglmer.cf,3),paste(" (",round(res.bglmer.se,3),")",sep=""),se="")

fxf<-res.bglmer[-c(1:6)]

#re1: careful r1=r2!!
th<-fit_bglmer@optinfo$val[1:3]
h<-fit_bglmer@optinfo$derivs$Hessian
h<-solve(h)
v<-forceSymmetric(h+t(h))
v<-v[1:3,1:3]
names(th)<-colnames(v)<-rownames(v)<-paste0("theta",1:3)



s11<-deltaMethod(th,"sqrt(theta1*theta1)",v)
r12<-deltaMethod(th,"(theta1*theta2)/sqrt(theta1**2*(theta2**2+theta3**2)  )   ",v)
s22<-deltaMethod(th,"sqrt(theta2**2+theta3**2)",v)

res.s3.bml2<-c(
  paste0(round(s11[1],3), " (",round(s11[2],3),")"   ),
  paste0(round(s22[1],3), " (",round(s22[2],3),")"   ),
  paste0(round(r12[1],3), " (",round(r12[2],3),")"   )
)

#re2
th<-fit_bglmer@optinfo$val[4:6]
h<-fit_bglmer@optinfo$derivs$Hessian
h<-solve(h)
v<-forceSymmetric(h+t(h))
v<-v[4:6,4:6]
names(th)<-colnames(v)<-rownames(v)<-paste0("theta",1:3)



s11<-deltaMethod(th,"sqrt(theta1*theta1)",v)
r12<-deltaMethod(th,"(theta1*theta2)/sqrt(theta1**2*(theta2**2+theta3**2)  )   ",v)
s22<-deltaMethod(th,"sqrt(theta2**2+theta3**2)",v)

res.s3.bml1<-c(
  paste0(round(s11[1],3), " (",round(s11[2],3),")"   ),
  paste0(round(s22[1],3), " (",round(s22[2],3),")"   ),
  paste0(round(r12[1],3), " (",round(r12[2],3),")"   )
)

resb<-c(fxf,res.s3.bml1,res.s3.bml2)

res<-cbind(resm,resm_reml,resb,resp_ml )

colnames(res)<-c("ML","REML", "BM","PML(1)" )

 
vrnm<-c("Intercept","Question$[$basic$]$","Topic$[$complex no.$]$","$[$diff. eq.$]$","$[$differentiation$]$",
        "$[$fund. math.$]$","$[$graph th.$]$","$[$integration$]$","$[$lin. algebra$]$","$[$num. methods$]$",
        "$[$optimization$]$","$[$probability$]$","$[$fun. of single var.$]$","$[$set th.$]$","$[$statistics$]$",
        "Country",rep("",2),"StudentCountry",rep("",2))

res<-cbind(vrnm,res)


rownm<-c("$\\beta_0$","$\\beta_1$","$\\beta_2$","$\\beta_3$","$\\beta_4$","$\\beta_5$","$\\beta_6$",
         "$\\beta_7$","$\\beta_8$","$\\beta_9$","$\\beta_{10}$",
         "$\\beta_{11}$","$\\beta_{12}$","$\\beta_{13}$","$\\beta_{14}$",
         "$\\sigma_{1,1}$","$\\sigma_{1,2}$",
         "$\\rho_1$",
         "$\\sigma_{2,1}$","$\\sigma_{2,2}$",
         "$\\rho_2$")

rownames(res)<-rownm
res_binom<-res

print(xtable(res),sanitize.text.function=function(x){x})




