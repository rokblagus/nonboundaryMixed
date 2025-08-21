################################################

##packages

library(mvtnorm) #for simulating u
library(nlme) #to fit non-linear model
library(pedigreemm) #to fit models with depended clusters due to a known pedigree
library(glmmTMB) 
library(spaMM)  
 


################################################


##aux functions

expit<-function(x) 1/(1+exp(-x))
"%^%" <- function(x, n)   with(eigen(x), vectors %*% (values^n * t(vectors)))

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

#univariate logistic

make_pseudo_data_rand_eigen_inter_alpha_beta<-function(alpha,beta,param="sigma2",const=1e8){
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
  
  
  pi0=1/(1+exp(-u1[1]))
  
  Y<-rep(c(pi0),fact) #the constant improves the convergence!
  n<-rep(rep(const,1),fact)
  id<-c(1:fact)
  Z<-matrix(rep(1,fact),ncol=1)
  data0<-list(Y=Y,grouping=id,nn=n,Z=Z)
  
  #fit0<-glmer(cbind(Y,nn-Y)~-1+(-1+Z|grouping),data=data0,family=binomial)
  #est.vcv<-VarCorr(fit0)$grouping[1:2,1:2]
  
  list(data=data0)#,fit=fit0,vcv.re=est.vcv)
}



################################################


###Extension to non-linear MM as discussed by Vones and Carter (Biometrics, 1992)

## Simulate data

set.seed(10)
n_subj <- 20
times <- seq(0, 10, length.out = 8)
dat <- expand.grid(subject = factor(1:n_subj), t = times)
 
alpha1 <- 5
alpha2 <- 0.4
   
#true covariance          
vr_int=0.1**2
vr_slope1=0.01**2
vr_slope2=0.05**2
rho<-0.8
D<-diag(c(vr_int,vr_slope1,vr_slope2))
D[1,2]<-D[2,1]<-rho*sqrt(vr_int*vr_slope1)
D[1,3]<-D[3,1]<-rho*sqrt(vr_int*vr_slope2)
D[2,3]<-D[3,2]<-rho*sqrt(vr_slope2*vr_slope1)
      
    
      
bn<-rmvnorm(n_subj,mu=rep(0,3),sigma=D)
b_i <- bn[,2]
b_i2 <- bn[,3]
b_0 <- bn[,1]
      
epsilon <- rnorm(nrow(dat), mean = 0, sd = 0.5)
        
dat$y <- alpha1 * exp(-alpha2 * dat$t) + b_0[as.numeric(dat$subject)]+
          b_i[as.numeric(dat$subject)] * dat$t +b_i2[as.numeric(dat$subject)] * dat$t^2+
              epsilon
        
         
##fit the model to original data

fit_model_orig <- nlme(
               y ~ intercept + alpha1 * exp(-alpha2 * t),
               data = dat,
               fixed =  alpha1 + alpha2 ~ 1,
               random = intercept ~ t+I(t^2) | subject,   
               start = c(5,0.4),control=nlmeControl(returnObject = TRUE,msMaxIter=500)
           )
         
summary(fit_model_orig)  #note boundary estimate of D




##fit the same model, but using the revised dataset which is compatible with our approach
 
dat_tr<-dat
dat_tr$w<-rep(1,nrow(dat_tr))
dat_tr$subject<-as.numeric(dat_tr$subject)
dat_tr$zt<-dat_tr$t
dat_tr$zt2<-dat_tr$t^2
dat_tr$zt0<-rep(1,nrow(dat_tr))
dat_tr$t2<-rep(1,nrow(dat_tr)) #this is used to make sure that pds do not affect estimation of beta1 and beta2

fit_model_rev <- nlme(
  y ~ intercept + alpha1*t2 * exp(-alpha2*t ),
  data = dat_tr,
  fixed =  alpha1 + alpha2 ~ 1,
  random = intercept ~ zt0+zt+zt2-1 | subject,  
  start = c(5,0.4),weights = varFixed(~I(1/w)),control=nlmeControl(returnObject = TRUE)
)           

#note that results are exactly the same as when using dat

##add pseudo-clusters
           
#set prior parameters

nu0<-2*3-1
           
psi0<-D*3*3 #oracle

pd<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi=psi0,nu=nu0,const=1e8,param="variance",link_fun=function(x) x)
           
           
#dataset for PD only

dat_pd<-data.frame(subject=max(dat_tr$subject)+pd$data$grouping,t=rep(0,length(pd$data$Y)),t2=rep(0,length(pd$data$Y)),y=pd$data$Y,w=pd$data$nn,zt=pd$data$Z[,2],
                             zt2=pd$data$Z[,3],zt0=pd$data$Z[,1])

#augmented dataset

dat_aug<-rbind(dat_tr,dat_pd)


#fit the model to augmented dataset

fit_model_aug <- nlme(
             y ~ intercept + alpha1*t2 * exp(-alpha2*t ),
             data = dat_aug,
             fixed =  alpha1 + alpha2 ~ 1,
             random = intercept ~ zt0+zt+zt2-1 | subject,   
             start = c(5,0.4),weights = varFixed(~I(1/w)),control=nlmeControl(returnObject = TRUE,msMaxIter=500)
           )
           
summary(fit_model_aug)           
           
     

###Extension to mixed models with dependent clusters as discussed by Jiang et al. (JASA, 2021)


 
#  Create a pedigree for 20 individuals
 

# Founders: I1–I8 (no parents)
founders <- data.frame(
  label = paste0("I", 1:8),
  dam = NA,
  sire = NA
)

# Offspring: I9–I20 (randomly assign valid parents from founders)
set.seed(123)
offspring <- data.frame(
  label = paste0("I", 9:20),
  dam = sample(founders$label[1:4], 12, replace = TRUE),
  sire = sample(founders$label[5:8], 12, replace = TRUE)
)

# Combine to full pedigree
pedigree_df <- rbind(founders, offspring)

# Create pedigree object
ped <- with(pedigree_df, pedigree(sire = sire, dam = dam, label = label))

 
#  Simulate random effects from pedigree-based covariance
 

A <- as.matrix(getA(ped))   
n <- nrow(A)
set.seed(42)
u <- t(chol(A)) %*% rnorm(n, sd = 1)  # Random effects ~ N(0, A)

#simulate other data

X <- model.matrix(~1, data = data.frame(id = pedigree_df$label))  # Intercept only
Z <- diag(n)
eta <- X %*% 0.5 + Z %*% u  # Linear predictor
prob <- 1 / (1 + exp(-eta))  # Inverse logit
y <- rbinom(n, size = 1, prob = prob)

df <- data.frame(y = y, id = pedigree_df$label,M=rep(1,length(pedigree_df$label)))


#   Fit the GLMM with pedigree random effects


fit <- pedigreemm(y ~ 1 + (1 | id),
                  family = binomial,
                  data = df,
                  weights=M,
                  pedigree = list(id = ped))

 
#restructure the dataset to be compatible with the PD approach
xdf<-df
xdf$x<-rep(1,nrow(xdf))
xdf$z<-rep(1,nrow(xdf))

fit <- pedigreemm(y ~ x-1 + (z-1 | id),
                  family = binomial,
                  data = xdf,
                  weights=M,
                  pedigree = list(id = ped))


##add PD (note logistic model and univariate RE)

alpha0=0.5
beta0=1*(alpha0+1) #oracle

pd<-make_pseudo_data_rand_eigen_inter_alpha_beta(alpha=alpha0,beta=beta0,param="sigma2",const=1e8)

df_pd<-data.frame(y=pd$data$Y,id=paste0("I",pd$data$grouping+20),M=pd$data$nn,x=rep(0,nrow(pd$data$Z)),z=pd$data$Z)

df_aug<-rbind(xdf,df_pd)


#augment pedigree: new clusters are "parents" with no "ofsprings"
pedigree_pd<-data.frame(label=paste0("I",pd$data$grouping+20),dam=NA,sire=NA)

pedigree_df_aug <- rbind(pedigree_df,pedigree_pd)

 ped_aug <- with(pedigree_df_aug, pedigree(sire = sire, dam = dam, label = label))

#fit the model to augmented data
 
 fit_aug <- pedigreemm(y ~ x-1 + (z-1 | id),
                   family = binomial,
                   data = df_aug,
                   weights=M,
                   pedigree = list(id = ped_aug))
 
 summary(fit)
 summary(fit_aug)

###Extension to mixed models with dependent clusters with spatial covariance across a closed area
 
 

 set.seed(100)
 n <- 50
 coords <- cbind(runif(n ), runif(n ))
 dist_mat <- as.matrix(dist(coords))
 
 # Parameters
 beta0<-0
 nu <- 1.5
 rho <- 0.8
 sigma2 <- 0.2
 
 # Get correlation matrix
 R <- MaternCorr(dist_mat, nu = nu, rho = rho)
 
 # Simulate random effects
 u <- t(chol(R )) %*% rnorm(n,sd=sqrt(sigma2))
 
# Simulate binary response
 eta <- beta0+ u
 prob <- plogis(eta)
 y <- rbinom(n, size = 1, prob = prob)
 
 df <- data.frame(y = y, x = coords[, 1], ycoord = coords[, 2],w=rep(1,length(y)),xx=rep(1,length(y)),id=1:length(y))
 
  
 df$pos <- numFactor( df$x ,  df$y )
 df$ID <- factor(rep(1, nrow(df)))
 
# fit the model
 m_tmb <- glmmTMB(cbind( y*w , (1-y)*w)  ~ -1+xx + mat(pos + 0 | ID), df,family = binomial()) 
 

 ##add PD
 
 alpha0=0.5
 beta0=sigma2*(alpha0+1) #oracle
 
 pd<-make_pseudo_data_rand_eigen_inter_alpha_beta(alpha=alpha0,beta=beta0,param="sigma2",const=1e8)
 n_pseudo<-nrow(pd$data$Z)
 x_pseudo <- seq(from=1e3,to=1e3+100, length.out=n_pseudo) 
 y_pseudo <- seq(from=1e3,to=1e3+100, length.out = n_pseudo)
 
  
 df_pd<-data.frame(y=pd$data$Y,w=pd$data$nn,xx=rep(0,nrow(pd$data$Z)),x=x_pseudo,ycoord=y_pseudo,id=pd$data$grouping+nrow(df))
 df_pd$pos <- numFactor( df_pd$x ,  df_pd$y )
 df_pd$ID <- factor(rep(1, nrow(df_pd)))
 
 
 df_aug<-rbind(df,df_pd)
 
 # fit the model
 m_tmb_a <- glmmTMB(cbind(floor(y*w),floor((1-y)*w)) ~ -1+xx + mat(pos + 0 | ID), df_aug,family = binomial()) 
 
 #increase confidence in our prior belief
 
 alpha0=10
 beta0=sigma2*(alpha0+1) #oracle
 
 pd<-make_pseudo_data_rand_eigen_inter_alpha_beta(alpha=alpha0,beta=beta0,param="sigma2",const=1e8)
 n_pseudo<-nrow(pd$data$Z)
 x_pseudo <- seq(from=1e3,to=1e3+100, length.out=n_pseudo) 
 y_pseudo <- seq(from=1e3,to=1e3+100, length.out = n_pseudo)
 
 
 df_pd<-data.frame(y=pd$data$Y,w=pd$data$nn,xx=rep(0,nrow(pd$data$Z)),x=x_pseudo,ycoord=y_pseudo,id=pd$data$grouping+nrow(df))
 df_pd$pos <- numFactor( df_pd$x ,  df_pd$y )
 df_pd$ID <- factor(rep(1, nrow(df_pd)))
 
 
 df_aug<-rbind(df,df_pd)
 
 # fit the model
 m_tmb_a2<- glmmTMB(cbind(floor(y*w),floor((1-y)*w)) ~ -1+xx + mat(pos + 0 | ID), df_aug,family = binomial()) 
 
 
 