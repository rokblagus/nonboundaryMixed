##creates tables for the lmm, binom, pois examples as R objects res_lin, res_binom, res_pois

###packages

library(glmmTMB)

library(lme4)

library(xtable)

library(nlme)

library(blme)

library(car)

library(merDeriv)  

library(haven)

library(tidyverse)

library(ggplot2)

library(gridExtra)

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


#####Linear exampledata(ellenberg)

library(blmeco)
data(ellenberg)


ellenberg$gradient <- paste(ellenberg$Year, ellenberg$Soil)
table(ellenberg$Species, ellenberg$gradient)


ellenberg$water.z <- as.numeric(scale(ellenberg$Water))

ellenberg<-ellenberg[,names(ellenberg)%in%c("Yi.g","water.z","Species","gradient")]

ellenberg$Species<-as.numeric(as.factor(ellenberg$Species))
ellenberg$gradient<-as.numeric(as.factor(ellenberg$gradient))


ellenberg<-na.omit(ellenberg)
g1<-ggplot(ellenberg, aes(x = water.z, y = log(Yi.g))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "blue",fill="lightblue") +
  facet_wrap(~ Species) +
  labs(
    x = "water.s",
    y = expression(log(Yi.g)) 
  ) +
  theme_minimal()

 
pdf("building_1.pdf",height=6,width=8)
g1
dev.off()


mod <- lmer(log(Yi.g) ~ water.z + I(water.z^2) +
              (water.z + I(water.z^2)|Species) + (1|gradient),
            data=ellenberg)

mod <- glmmTMB(log(Yi.g) ~ water.z + I(water.z^2) +
              (water.z + I(water.z^2)|Species) + (1|gradient),
            data=ellenberg,family = gaussian(link="identity"))

 
Y<-log(ellenberg$Yi.g)
X<-model.matrix(log(Yi.g) ~ water.z + I(water.z^2),data=ellenberg)
Z1<-X
grouping1<- ellenberg$Species 
Z2<-matrix(1,ncol=1,nrow=nrow(X))
grouping2<-ellenberg$gradient


xdf<-list(Y=Y,X=X,Z1=Z1,grouping1=grouping1,Z2=Z2,grouping2=grouping2)

fit_glmer<-glmmTMB(Y~X-1 +
                     (Z1-1 | grouping1)+(Z2-1 | grouping2), data = xdf,family = gaussian(link="identity"))


fit_bglmer<-blmer(Y~X-1 +
                    (Z1-1 | grouping1)+(Z2-1 | grouping2), data = xdf)


fit_glmer_r<-glmmTMB(Y~X-1 +
                       (Z1-1 | grouping1)+(Z2-1 | grouping2), data = xdf,family = gaussian(link="identity"),REML=TRUE)



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
  
  pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) x  )
  
  
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
  Z1a<-rbind(xdf$Z1,pd1$data$Z)
  Z2a<-rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
  
  
  Ya<-c(xdf$Y,pd1$data$Y)
  weightsa<-c(rep(1,length(xdf$Y)),pd1$data$nn)
  
  
  grouping1a<-c(xdf$grouping1,max(xdf$grouping1)+pd1$data$grouping)
  grouping2a<-c(xdf$grouping2,max(xdf$grouping2)+pd1$data$grouping)
  
  
  
  
  xdfa<-list(Y=Ya,weights=weightsa,X=Xa,Z1=Z1a,Z2=Z2a,grouping1=grouping1a,grouping2=grouping2a)
   
  tmp2 <- glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                  dispformula = ~offset(-log(weights)),
                  data=xdfa)
  
  
  tmp2
}

get_marLik_lin<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
 
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family=gaussian(link = "identity"),
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                           betad=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="betad")],
                           theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) #point estimates seem ok, but logLik is NA! They have a trick where they dont want to report loogLik if the model does not converge (which in our case defacto holds), but we can still accesss it via object$fit$objective which seems to give -loglik so it should be minimized
  
  
  
  -tmp3$fit$objective
  
}

tau_finder_lin<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_lin_tau(tau,D_est,xdf)
  abs(get_marLik_lin(fit_tau,xdf)-get_marLik_lin(fit_ml,xdf))-qchisq(1-alpha,1)/2
}

D<-VarCorr(fit_glmer)$cond$grouping1[1:3,1:3]

sek<-seq(from=0,to=0.5,by=0.05)
y<-rep(NA,length(sek))
zz=0
for (i in sek){
  zz=zz+1
  y[zz]<-tau_finder_lin(i,xdf,VarCorr(fit_glmer)$cond$grouping1[1:3,1:3],fit_glmer)
}
plot(sek,y,type="l")
abline(h=0)

opt_tau_ml<-uniroot(tau_finder_lin,c(0,1),xdf=xdf,D_est=VarCorr(fit_glmer)$cond$grouping1[1:3,1:3],
                    fit_ml=fit_glmer,alpha=0.05)

fit_tau_ml<-fiter_lin_tau(opt_tau_ml$root,D_est=VarCorr(fit_glmer)$cond$grouping1[1:3,1:3],xdf=xdf)

plot_ml<-data.frame(fit=fitted(fit_glmer),res=residuals(fit_glmer))
plot_ml$method="ML"

plot_rml<-data.frame(fit=fitted(fit_glmer_r),res=residuals(fit_glmer_r))
plot_rml$method="REML"

plot_bml<-data.frame(fit=fitted(fit_bglmer),res=residuals(fit_bglmer))
plot_bml$method="BM"

plot_pml<-data.frame(fit=fitted(fit_tau_ml),res=residuals(fit_tau_ml))
plot_pml$method="PML"

plot_diag<-rbind(plot_ml,plot_bml,plot_rml,plot_pml[1:nrow(plot_ml),])

tuk<-ggplot(plot_diag, aes(x = fit, y = res)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  facet_wrap(~ method) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Fitted values",
    y = "Residuals" 
    
  ) +
  theme_minimal()

pdf("diag.pdf",height=6,width=6)
tuk
dev.off()

##plot coefs (at the internal scale!) with their Wald type CIs: you can make an argument that theta56 are unreasonable with ML/REML, further for theta456 you can get CIs; CIs for BM for some theta at this scale are very wide, resulting in very wide CIs 

mod<-fit_glmer

est<-mod$sdr$par.fixed
ses<-sqrt(diag(mod$sdr$cov.fixed))

res<-cbind(est,ses,est-qnorm(0.975)*ses,est+qnorm(0.975)*ses)
colnames(res)<-c("Estimate","SE","CI_low","CI_up")

df <- as.data.frame(res)
df$Coefficient <- rownames(res)  # use row names as labels
df$Coefficient[df$Coefficient=="beta"]<-paste0(df$Coefficient[df$Coefficient=="beta"],0:2)
df$Coefficient[df$Coefficient=="theta"]<-paste0(df$Coefficient[df$Coefficient=="theta"],1:7)


# Make Coefficient a factor to control order in plot
df$Coefficient <- factor(df$Coefficient) 
#df$group<-factor(ifelse(df$Estimate<(-100),1,2))


df$Method="ML"

df_ml<-df

##PML

mod<-fit_tau_ml

est<-mod$sdr$par.fixed
ses<-sqrt(diag(mod$sdr$cov.fixed))

res<-cbind(est,ses,est-qnorm(0.975)*ses,est+qnorm(0.975)*ses)
colnames(res)<-c("Estimate","SE","CI_low","CI_up")

df <- as.data.frame(res)
df$Coefficient <- rownames(res)  # use row names as labels
df$Coefficient[df$Coefficient=="beta"]<-paste0(df$Coefficient[df$Coefficient=="beta"],0:2)
df$Coefficient[df$Coefficient=="theta"]<-paste0(df$Coefficient[df$Coefficient=="theta"],1:7)


# Make Coefficient a factor to control order in plot
df$Coefficient <- factor(df$Coefficient) 
#df$group<-factor(ifelse(df$Coefficient%in%c("theta5","theta6"),1,2))



df$Method="PML"
df_pml<-df

##REML

mod<-fit_glmer_r

est<-c(fixef(mod)$cond,mod$sdr$par.fixed)
ses<-c(sqrt(diag(vcov(mod)$cond)),sqrt(diag(mod$sdr$cov.fixed)))

res<-cbind(est,ses,est-qnorm(0.975)*ses,est+qnorm(0.975)*ses)
colnames(res)<-c("Estimate","SE","CI_low","CI_up")

df <- as.data.frame(res)
df$Coefficient <- df_ml$Coefficient

# Make Coefficient a factor to control order in plot
#df$group<-factor(ifelse(df$Coefficient%in%c("theta5","theta6"),1,2))

df$Method="REML"

df_reml<-df


##BM

#estr<-getME(fit_bglmer, "theta") 

#hess<-fit_bglmer@optinfo$derivs$Hessian

#we need to put thetas on the same scale as used in TMB! the best way is via D and its SE, as we use in our Tables!

mod<-fit_bglmer
vv <- vcov(mod, full = TRUE)



#ranef, first RE

vr<-c(VarCorr(mod)$grouping1)[-c(4,7,8)]


vrv<-vv[-c(1:3,10:11),-c(1:3,10:11)]

names(vr)<-colnames(vrv)<-rownames(vrv)<-c("sigma11","sigma12","sigma13","sigma22","sigma23","sigma33")

#s11<-deltaMethod(vr,"log(sqrt(sigma11))",vrv)
#s22<-deltaMethod(vr,"log(sqrt(sigma22))",vrv)
#s33<-deltaMethod(vr,"log(sqrt(sigma33))",vrv)

#avoid deltaMethod so that we have the same func for all!

#theta1:
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  #theta1 <- L[3,1] / L[3,3]
  #theta2 <- L[3,2] / L[3,3]
  log(s1)
}
grad <- grad(func = theta_expr, x = vr)
var_theta1 <- t(grad) %*% vrv %*% grad
se_theta1 <- sqrt(var_theta1)

theta1<-c(theta_expr(vr),as.numeric(se_theta1))

#theta2:
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  #theta1 <- L[3,1] / L[3,3]
  #theta2 <- L[3,2] / L[3,3]
  log(s2)
}
grad <- grad(func = theta_expr, x = vr)
var_theta2 <- t(grad) %*% vrv %*% grad
se_theta2 <- sqrt(var_theta2)

theta2<-c(theta_expr(vr),as.numeric(se_theta2))

#theta3:
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  #theta1 <- L[3,1] / L[3,3]
  #theta2 <- L[3,2] / L[3,3]
  log(s3)
}
grad <- grad(func = theta_expr, x = vr)
var_theta3 <- t(grad) %*% vrv %*% grad
se_theta3 <- sqrt(var_theta3)

theta3<-c(theta_expr(vr),as.numeric(se_theta3))


#theta4:
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  theta1 <- L[3,1] / L[3,3]
  theta2 <- L[3,2] / L[3,3]
  theta0
}
grad <- grad(func = theta_expr, x = vr)
var_theta4 <- t(grad) %*% vrv %*% grad
se_theta4 <- sqrt(var_theta4)

theta4<-c(theta_expr(vr),as.numeric(se_theta4))

#theta5
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  theta1 <- L[3,1] / L[3,3]
  theta2 <- L[3,2] / L[3,3]
  theta1
}
grad <- grad(func = theta_expr, x = vr)
var_theta5 <- t(grad) %*% vrv %*% grad
se_theta5 <- sqrt(var_theta5)

theta5<-c(theta_expr(vr),as.numeric(se_theta5))



#theta6
theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  sigma12 <- vr["sigma12"]
  sigma13 <- vr["sigma13"]
  sigma22 <- vr["sigma22"]
  sigma23 <- vr["sigma23"]
  sigma33 <- vr["sigma33"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  s2 <- sqrt(sigma22)
  s3 <- sqrt(sigma33)
  
  # Correlations
  r12 <- sigma12 / (s1 * s2)
  r13 <- sigma13 / (s1 * s3)
  r23 <- sigma23 / (s2 * s3)
  
  # Correlation matrix
  R <- matrix(c(1,   r12, r13,
                r12, 1,   r23,
                r13, r23, 1), nrow=3, byrow=TRUE)
  
  # Cholesky factor (lower triangular)
  L <- t(chol(R))
  
  # Compute the theta parameters
  theta0 <- L[2,1] / L[2,2]
  theta1 <- L[3,1] / L[3,3]
  theta2 <- L[3,2] / L[3,3]
  theta2
}
grad <- grad(func = theta_expr, x = vr)
var_theta6 <- t(grad) %*% vrv %*% grad
se_theta6 <- sqrt(var_theta6)

theta6<-c(theta_expr(vr),as.numeric(se_theta6))


#theta7
vr<-c(VarCorr(mod)$grouping2) 


vrv<-vv[10,10,drop=FALSE]

names(vr)<-colnames(vrv)<-rownames(vrv)<-c("sigma11")


theta_expr <- function(vr){
  sigma11 <- vr["sigma11"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  log(s1)
}
grad <- grad(func = theta_expr, x = vr)
var_theta7 <- t(grad) %*% vrv %*% grad
se_theta7 <- sqrt(var_theta7)

theta7<-c(theta_expr(vr),as.numeric(se_theta7))

#residual

st<-summary(mod)$sigma**2
vst<-vv[11,11,drop=FALSE]
names(st)<-colnames(vst)<-rownames(vst)<-"sigma1"

s11<-deltaMethod(st,"log(sqrt(sigma1))",vst)

theta_expr <- function(vr){
  sigma11 <- vr["sigma1"]
  
  # Compute standard deviations
  s1 <- sqrt(sigma11)
  log(s1)
}
grad <- grad(func = theta_expr, x = st)
var_disp <- t(grad) %*% vst %*% grad
se_disp <- sqrt(var_disp)

disp<-c(theta_expr(st),as.numeric(se_disp))



est<-c(fixef(fit_bglmer),disp[1],
       c(theta1[1],theta2[1],theta3[1],theta4[1],theta5[1],theta6[1],theta7[1]) 
)
ses<-c(sqrt(diag(vcov(fit_bglmer))),disp[2],
       c(theta1[2],theta2[2],theta3[2],theta4[2],theta5[2],theta6[2],theta7[2])
)

res<-cbind(est,ses,est-qnorm(0.975)*ses,est+qnorm(0.975)*ses)
colnames(res)<-c("Estimate","SE","CI_low","CI_up")

df <- as.data.frame(res)
df$Coefficient <- df_ml$Coefficient 
#df$group<-factor(ifelse(df$Coefficient%in%c("theta5","theta6"),1,2))
df$Method="BM"

df_bm<-df

df_all <- bind_rows(df_ml, df_reml, df_pml, df_bm)
df_all$Coefficient<-as.character(df_all$Coefficient)
#df_all$Coefficient[df_all$Coefficient%in%paste0("theta",1:7)]<-paste0("  ",as.character(df_all$Coefficient[df_all$Coefficient%in%paste0("theta",1:7)]))
df_all$Coefficient<-factor(df_all$Coefficient)
df_all$Method<-factor(df_all$Method,levels=c("ML","REML","PML","BM"))



signed_log <- function(x,a=1/10) {
  #sign(x) * log10(abs(x) + 1)
  #asinh(x/a)
  #sign(x) * abs(x)^a
  x
}


df_all$Coef2<-as.character(df_all$Coefficient)
df_all$Coef2[df_all$Coef2=="betadisp"]<-"dispersion"




plot_coefs<-ggplot(df_all, 
                   aes(y = Method, x = signed_log(Estimate), color = Method)) +
  geom_pointrange(aes(xmin = signed_log(CI_low), xmax = signed_log(CI_up))) +
  scale_color_manual(values = dfpallete) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  facet_grid(Coef2 ~ ., scales = "free_y", switch = "y") + 
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    panel.spacing.y = unit(0.5, "cm"),
    panel.grid.major.y = element_blank(),    
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),    
    panel.grid.minor.x = element_blank()      
  ) +
  labs(x = "Estimate with 95% CI", y = "Coefficients", color = NULL) +
  coord_cartesian(xlim = c(-5, 5))






####get predictions with 95% confidence intervals 
#maybe we would need to avoid boot, since we did not evaluate this in sims, wald could maybe be justified?
#note: this are not prediction intervals (ie our uncertanty about the individual value), but CI for the mean predicted response (not for individual observation)

#make a new df for which predictions are to be calculated

#prediction intervals, for each method to show that (if) they are of different widths



tau<-opt_tau_ml$root
D_est<-D
q<-ncol(D_est)
ee<-eigen(D_est)
ee$values[ee$values<1e-4]<-1e-4
ee$values[ee$values>1e4]<-1e4
lm<-mean(ee$values)
li<-ee$values+tau*(lm-ee$values)
psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q

nu=2*q-1

pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) x  )


Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
Z1a<-rbind(xdf$Z1,pd1$data$Z)
Z2a<-rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))


Ya<-c(xdf$Y,pd1$data$Y)
weightsa<-c(rep(1,length(xdf$Y)),pd1$data$nn)


grouping1a<-c(xdf$grouping1,max(xdf$grouping1)+pd1$data$grouping)
grouping2a<-c(xdf$grouping2,max(xdf$grouping2)+pd1$data$grouping)




xdfa<-list(Y=Ya,weights=weightsa,X=Xa,Z1=Z1a,Z2=Z2a,grouping1=grouping1a,grouping2=grouping2a)


xdfadf<-data.frame(Y=xdfa$Y,X=xdfa$X,Z1=xdfa$Z1,Z2=xdfa$Z2,gr1=xdfa$grouping1,gr2=xdfa$grouping2,w=xdfa$weights)




newtime=seq(from=-5,to=5,by=0.01)

newx<-data.frame(X..Intercept.=1,X.water.z=newtime,
                 X.I.water.z.2.=newtime**2,
                 Z1..Intercept.=1,
                 Z1.water.z=newtime,
                 Z1.I.water.z.2.=newtime**2,
                 gr1=4,
                 Z2=1,
                 gr2=1,
                 w=1,
                 x=newtime
)

#PML

tmp2 <- glmmTMB(Y~-1+X..Intercept.+X.water.z+X.I.water.z.2.+
                  (-1+Z1..Intercept.+Z1.water.z+Z1.I.water.z.2.|gr1)+
                  (-1+Z2|gr2), family = gaussian(link = "identity"),
                dispformula = ~offset(-log(w)),
                data=xdfadf)

preds<-predict(tmp2,newdata=newx, re.form = NULL,se.fit=TRUE)

results <- cbind(
  
  Predicted = preds$fit,
  Lower = preds$fit-qnorm(0.975)*preds$se.fit,
  Upper = preds$fit+qnorm(0.975)*preds$se.fit
)

predint_pml<-cbind(newx,results)  


#ML

tmp2 <- glmmTMB(Y~-1+X..Intercept.+X.water.z+X.I.water.z.2.+
                  (-1+Z1..Intercept.+Z1.water.z+Z1.I.water.z.2.|gr1)+
                  (-1+Z2|gr2), family = gaussian(link = "identity"),
                
                data=xdfadf[1:nrow(xdf$X),])

preds<-predict(tmp2,newdata=newx, re.form = NULL,se.fit=TRUE)

results <- cbind(
  
  Predicted = preds$fit,
  Lower = preds$fit-qnorm(0.975)*preds$se.fit,
  Upper = preds$fit+qnorm(0.975)*preds$se.fit
)

predint_ml<-cbind(newx,results)  


#REML

tmp2 <- glmmTMB(Y~-1+X..Intercept.+X.water.z+X.I.water.z.2.+
                  (-1+Z1..Intercept.+Z1.water.z+Z1.I.water.z.2.|gr1)+
                  (-1+Z2|gr2), family = gaussian(link = "identity"),
                REML=TRUE,
                data=xdfadf[1:nrow(xdf$X),])

preds<-predict(tmp2,newdata=newx, re.form = NULL,se.fit=TRUE)

results <- cbind(
  
  Predicted = preds$fit,
  Lower = preds$fit-qnorm(0.975)*preds$se.fit,
  Upper = preds$fit+qnorm(0.975)*preds$se.fit
)

predint_reml<-cbind(newx,results)  

#BM
tmp2 <- blmer(Y~-1+X..Intercept.+X.water.z+X.I.water.z.2.+
                (-1+Z1..Intercept.+Z1.water.z+Z1.I.water.z.2.|gr1)+
                (-1+Z2|gr2) ,
              
              data=xdfadf[1:nrow(xdf$X),])

#newx$gr1<-factor(newx$gr1,levels=levels(factor(xdf$grouping1)))
#newx$gr2<-factor(newx$gr2,levels=levels(factor(xdf$grouping2)))

#preds<-predict(tmp2,newdata=newx, re.form = NULL,se.fit=TRUE)
#note that this gives error due to a bug: https://github.com/lme4/lme4/issues/815
#we need all levels of clusters to be present!

zz=0
for (gr1 in unique(xdf$grouping1)){
  for (gr2 in unique(xdf$grouping2)){
    zz=zz+1
    new_i<-newx[,c(1:6,8,10:11)]
    new_i$gr1<-gr1
    new_i$gr2<-gr2
    if (zz==1) new_blmer<-new_i else new_blmer<-rbind(new_blmer,new_i)
  }
}
preds<-predict(tmp2,newdata=new_blmer, re.form = NULL,se.fit=TRUE)


results <- cbind(
  
  Predicted = preds$fit,
  Lower = preds$fit-qnorm(0.975)*preds$se.fit,
  Upper = preds$fit+qnorm(0.975)*preds$se.fit
)

predint_bm<-cbind(newx,results[which(new_blmer$gr1==4&new_blmer$gr2==1),])




predint_bm$Method<-"BM"
predint_ml$Method<-"ML"
predint_reml$Method<-"REML"
predint_pml$Method<-"PML"

res_a<-rbind(predint_ml,predint_reml,predint_pml,predint_bm)

res_a$Method<-factor(res_a$Method,levels=c("ML","REML","PML","BM"))

dfpallete<-c("black","red","blue","cadetblue")



plot_pred<-ggplot(res_a, aes(x = x, y = Predicted, color = Method, fill = Method)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
  geom_line(size = 1,alpha=0.7) +
  theme_minimal() +
  scale_color_manual(values = dfpallete) +
  scale_fill_manual(values = dfpallete) +
  labs(x = "Groundwater depth (standardized)", y = "Predicted log-transformed yield (with 95% confidence interval for the mean response)")+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))




###combine coefs and predictions
pdf("fig_elenberg_coefs.pdf",height=8,width=15)
grid.arrange(plot_coefs,plot_pred,nrow=1)
dev.off()




##Model diagnostics

#C_thetai=(hatthtea_i-hattheta)^t cov(hattheta)^-1 (hatthtea_i-hattheta), where theta q-vector (we work at the scale of internal parametrization)
#to do this, we need hattheta_i, ie the estimate when i is ommited

#PML

vr<-fit_tau_ml$sdr$cov.fixed[5:10,5:10] #cov
est<-fit_tau_ml$sdr$par.fixed[5:10] #hattheta


esti<-list()
for (i in 1:6){
  
  omit_idx<-which(xdfa$grouping1==i)
  xdfs<- lapply(xdfa, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 dispformula = ~offset(-log(weights)),
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[5:10]
}
cooki<-list()
for (i in 1:6) cooki[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cooki))

# ML

vr<-fit_glmer$sdr$cov.fixed[5:10,5:10] #cov
est<-fit_glmer$sdr$par.fixed[5:10] #hattheta


esti<-list()
for (i in 1:6){
  
  omit_idx<-which(xdf$grouping1==i)
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[5:10]
}
cookio<-list()
for (i in 1:6) cookio[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookio))

#REML
vr<-fit_glmer_r$sdr$cov.fixed[2:7,2:7] #cov
est<-fit_glmer_r$sdr$par.fixed[2:7] #hattheta


esti<-list()
for (i in 1:6){
  
  omit_idx<-which(xdf$grouping1==i)
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 REML=TRUE,
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[2:7]
}
cookior<-list()
for (i in 1:6) cookior[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookior))

#BML

est<-getME(fit_bglmer, "theta")[1:6]

hess<-fit_bglmer@optinfo$derivs$Hessian



esti<-list()
for (i in 1:6){
  
  omit_idx<-which(xdf$grouping1==i)
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-blmer(Y~X-1 +
                 (Z1-1 | grouping1)+(Z2-1 | grouping2), data = xdfs)
  esti[[i]]<-getME(fit_a, "theta")[1:6]
}

cookiobm<-list()
for (i in 1:6) cookiobm[[i]]<-(esti[[i]]-est)%*%hess[1:6,1:6]%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookiobm))

cooks<-cbind(1:6,unlist(cookio),unlist(cookior),unlist(cooki),unlist(cookiobm))
colnames(cooks)<-c("Species","ML","REML","PML","BM")

#what about omitting each row (slow? okish)

#PML
vr<-fit_tau_ml$sdr$cov.fixed[5:10,5:10] #cov
est<-fit_tau_ml$sdr$par.fixed[5:10] #hattheta


esti<-list()
for (i in 1:nrow(xdf$X)){
  
  omit_idx<-i
  xdfs<- lapply(xdfa, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 dispformula = ~offset(-log(weights)),
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[5:10]
}
cooki_i<-list()
for (i in 1:nrow(xdf$X)) cooki_i[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cooki_i),type="b")



# ML

vr<-fit_glmer$sdr$cov.fixed[5:10,5:10] #cov
est<-fit_glmer$sdr$par.fixed[5:10] #hattheta


esti<-list()
for (i in 1:nrow(xdf$X)){
  
  omit_idx<-i
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[5:10]
}
cookio_i<-list()
for (i in 1:nrow(xdf$X)) cookio_i[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookio_i),type="b")

#REML
vr<-fit_glmer_r$sdr$cov.fixed[2:7,2:7] #cov
est<-fit_glmer_r$sdr$par.fixed[2:7] #hattheta


esti<-list()
for (i in 1:nrow(xdf$X)){
  
  omit_idx<-i
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), family = gaussian(link = "identity"),
                 REML=TRUE,
                 data=xdfs)
  esti[[i]]<-fit_a$sdr$par.fixed[2:7]
}
cookior_i<-list()
for (i in 1:nrow(xdf$X)) cookior_i[[i]]<-(esti[[i]]-est)%*%solve(vr)%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookior_i),type="b")

#BML

est<-getME(fit_bglmer, "theta")[1:6]

hess<-fit_bglmer@optinfo$derivs$Hessian



esti<-list()
for (i in 1:nrow(xdf$X)){
  
  omit_idx<-i
  xdfs<- lapply(xdf, function(x) {
    if (is.vector(x) && !is.matrix(x)) {
      x[-omit_idx]
    } else if (is.matrix(x)) {
      x[-omit_idx, , drop = FALSE]
    } else {
      x
    }
  })
  fit_a<-blmer(Y~X-1 +
                 (Z1-1 | grouping1)+(Z2-1 | grouping2), data = xdfs)
  esti[[i]]<-getME(fit_a, "theta")[1:6]
}

cookiobm_i<-list()
for (i in 1:nrow(xdf$X)) cookiobm_i[[i]]<-(esti[[i]]-est)%*%hess[1:6,1:6]%*%matrix(esti[[i]]-est,ncol=1)
plot(unlist(cookiobm_i))


cooks_i<-cbind(1:nrow(xdf$X),unlist(cookio_i),unlist(cookior_i),unlist(cooki_i),unlist(cookiobm_i))
colnames(cooks_i)<-c("Species","ML","REML","PML","BM")

library(tidyverse)

# Convert to data frame
df_cooks <- as.data.frame(cooks_i)

# Make sure Species is treated as a factor or integer
df_cooks$Species <- as.integer(df_cooks$Species)

# Reshape to long format: key = Method, value = Cook
df_long <- df_cooks %>%
  pivot_longer(cols = -Species, names_to = "Method", values_to = "Cook")

df_long$Method <- factor(df_long$Method, levels = c("ML", "REML", "PML", "BM"))


dfpallete<-c("black","red","blue","cadetblue")
#cutoff <- qchisq(0.95, df = 6)

p1<-ggplot(df_long, aes(x = Species, y = Cook, color = Method)) +
  geom_point() +
  geom_line() +
  # geom_hline(yintercept = cutoff, linetype = "dashed", color = "darkgrey") +
  facet_wrap(~ Method, scales = "free_y",nrow=1) +  
  scale_color_manual(values = dfpallete) +
  theme_minimal() +
  labs(x = "Observation", y = "Cook's Distance") +
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(from=0,by=20,to=200))

#by cluster
df_cookss <- as.data.frame(cooks)

# Make sure Species is treated as a factor or integer
df_cookss$Species <- as.integer(df_cookss$Species)

# Reshape to long format: key = Method, value = Cook
df_longs <- df_cookss %>%
  pivot_longer(cols = -Species, names_to = "Method", values_to = "Cook")

df_longs$Method <- factor(df_longs$Method, levels = c("ML", "REML", "PML", "BM"))

p2<-ggplot(df_longs, aes(x = Species, y = Cook, color = Method)) +
  geom_point() +
  geom_line() +
  # geom_hline(yintercept = cutoff, linetype = "dashed", color = "darkgrey") +
  facet_wrap(~ Method, scales = "free_y",nrow=1) +  
  scale_color_manual(values = dfpallete) +
  theme_minimal() +
  labs(x = "Species", y = "Cook's Distance") +
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 1:6)

library(gridExtra)

pdf("fig_elenberg_cook.pdf",width = 15,height=10)
grid.arrange(p1, p2, nrow = 2)
dev.off()




#####sum res, Table


#####sum res

#PML

mod<-fit_tau_ml 



#sd for 2nd RE
vr<-mod$sdr$cov.fixed[11,11]
est<-mod$sdr$par.fixed[11]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[5:7,5:7]
est<-mod$sdr$par.fixed[5:7]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[8:10,8:10]
est<-mod$sdr$par.fixed[8:10]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:3,1:3]
est<-mod$sdr$par.fixed[1:3]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")


#disp

vr<-mod$sdr$cov.fixed[4,4]
est<-mod$sdr$par.fixed[4]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"(exp(theta))",vr)



resi<-c(fxf,
        
        
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")"),
        
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

resp2_ml<-resi


#ML

mod<-fit_glmer 



#sd for 2nd RE
vr<-mod$sdr$cov.fixed[11,11]
est<-mod$sdr$par.fixed[11]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[5:7,5:7]
est<-mod$sdr$par.fixed[5:7]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[8:10,8:10]
est<-mod$sdr$par.fixed[8:10]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:3,1:3]
est<-mod$sdr$par.fixed[1:3]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")


#disp

vr<-mod$sdr$cov.fixed[4,4]
est<-mod$sdr$par.fixed[4]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"(exp(theta))",vr)



resi<-c(fxf,
        
        
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")"),
        
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

res_ml<-resi



###REML
 

mod<-fit_glmer_r 



#sd for 2nd RE
vr<-mod$sdr$cov.fixed[11-3,11-3]
est<-mod$sdr$par.fixed[11-3]
names(est)<-"theta"

s2<-deltaMethod(est,"exp(theta)",vr)

#for 3rd RE
#diags
vr<-mod$sdr$cov.fixed[(5:7)-3,(5:7)-3]
est<-mod$sdr$par.fixed[(5:7)-3]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)
s3.3<-deltaMethod(est,"exp(theta3)",vr)

#off diags

vr<-mod$sdr$cov.fixed[(8:10)-3,(8:10)-3]
est<-mod$sdr$par.fixed[(8:10)-3]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",0:2)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)
s3.13<-deltaMethod(est,"theta1/sqrt(theta1^2+theta2^2+1)",vr)
s3.23<-deltaMethod(est,"(theta0*theta1+theta2)/sqrt(theta0^2+1)/sqrt(theta1^2+theta2^2+1)",vr)

#fixef
est<-fixef(mod)$cond
vr<-vcov(mod)$cond

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")


#disp

vr<-mod$sdr$cov.fixed[4-3,4-3]
est<-mod$sdr$par.fixed[4-3]
names(est)<-names(vr)<-"theta"

s.r<-deltaMethod(est,"(exp(theta))",vr)



resi<-c(fxf,
        
        
        paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
        paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
        paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")"),
        
        paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")"),
        paste0(round(s3.13$Estimate,3)," (",round(s3.13$SE,3),")"),
        paste0(round(s3.23$Estimate,3)," (",round(s3.23$SE,3),")"),
        
        paste0(round(s2$Estimate,3)," (",round(s2$SE,3),")"),
        
        paste0(round(s.r$Estimate,3)," (",round(s.r$SE,3),")")
)

res_reml<-resi


##BM



mod<-fit_bglmer
vv <- vcov(mod, full = TRUE)

#fixef

fixf<-paste(round(fixef(mod),3)," (", round(sqrt(diag(vv))[1:3],3),")",sep="")

#ranef

vr<-c(VarCorr(mod)$grouping1)[-c(4,7,8)]


vrv<-vv[-c(1:3,10:11),-c(1:3,10:11)]

names(vr)<-colnames(vrv)<-rownames(vrv)<-c("sigma11","sigma12","sigma13","sigma22","sigma23","sigma33")

s11<-deltaMethod(vr,"sqrt(sigma11)",vrv)
s22<-deltaMethod(vr,"sqrt(sigma22)",vrv)
s33<-deltaMethod(vr,"sqrt(sigma33)",vrv)
r12<-deltaMethod(vr,"sigma12/sqrt(sigma11*sigma22)",vrv)
r13<-deltaMethod(vr,"sigma13/sqrt(sigma11*sigma33)",vrv)
r23<-deltaMethod(vr,"sigma23/sqrt(sigma33*sigma22)",vrv)

rnf<-c(paste0(round(s11$Estimate,3)," (",round(s11$SE,3),")"),
       paste0(round(s22$Estimate,3)," (",round(s22$SE,3),")"),
       paste0(round(s33$Estimate,3)," (",round(s33$SE,3),")"),
       paste0(round(r12$Estimate,3)," (",round(r12$SE,3),")"),
       paste0(round(r13$Estimate,3)," (",round(r13$SE,3),")"),
       paste0(round(r23$Estimate,3)," (",round(r23$SE,3),")")
)

#R2

vr<-c(VarCorr(mod)$grouping2) 


vrv<-vv[10,10,drop=FALSE]

names(vr)<-colnames(vrv)<-rownames(vrv)<-c("sigma11")

s11<-deltaMethod(vr,"sqrt(sigma11)",vrv)

rnf2<-c(paste0(round(s11$Estimate,3)," (",round(s11$SE,3),")"))

#residual

st<-summary(mod)$sigma**2
vst<-vv[11,11]
names(st)<-names(vst)<-"sigma1"

s11<-deltaMethod(st,"sqrt(sigma1)",vst)

rs<-paste0(round(s11$Estimate,3)," (",round(s11$SE,3),")")

resi<-c(fixf,rnf,rnf2,rs)

res_bglmer<-resi




res<-cbind(res_ml,res_reml,res_bglmer, resp2_ml )
colnames(res)<-c("ML","REML","BM", "PML" )

vrnm<-c("Intercept","water.z","I(water.z^2)","species (Intercept)","(water.z)","(I(water.z^2))",rep("",3),"gradient (Intercept)","")

res<-cbind(vrnm,res)


rownm<-c("$\\beta_0$","$\\beta_1$","$\\beta_2$",
         "$\\sigma_{1,1}$","$\\sigma_{1,2}$","$\\sigma_{1,3}$",
         "$\\rho_{12}$","$\\rho_{13}$","$\\rho_{23}$","$\\sigma_{2}$","$\\sigma_\\epsilon$")

rownames(res)<-rownm
res_lin<-res

print(xtable(res),sanitize.text.function=function(x){x})



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

colnames(res)<-c("ML","REML","BM","PML")

nmr<-c("Intercept","Year$[$1996$]$","Year$[$1997$]$","Height","Brood","Index","Location",rep("",5))

res<-cbind(nmr,res)

rownm<-c("$\\beta_0$","$\\beta_1$","$\\beta_2$","$\\beta_3$","$\\sigma_1$","$\\sigma_2$",
         "$\\sigma_{3,1}$","$\\sigma_{3,2}$","$\\sigma_{3,3}$","$\\rho_{12}$","$\\rho_{13}$","$\\rho_{23}$")

rownames(res)<-rownm
res_pois<-res

print(xtable(res),sanitize.text.function=function(x){x})

 




#########################

#########################

 
##binom, math data


dd<-read.csv2("MathEdataset.csv")

#dd$Student.ID2<-paste0(dd$Student.ID,dd$Student.Country)

 
Y<-dd$Type.of.Answer
X<-model.matrix(Type.of.Answer~Question.Level+Topic,data=dd)
Z1<-model.matrix(~Question.Level,data=dd)
Z2<-model.matrix(~Question.Level,data=dd)
grouping1<-as.numeric(as.factor(dd$Student.Country))
#grouping2<-as.numeric(as.factor(dd$Student.ID2))
grouping2<-as.numeric(as.factor(dd$Student.ID))



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

colnames(res)<-c("ML","REML", "BM","PML" )

 
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




