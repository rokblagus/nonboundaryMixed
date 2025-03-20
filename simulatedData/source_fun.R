####packages

library(glmmTMB)
library(mvtnorm)
library(blme)


 
####


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

###LMM functions

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



fiter_binom_tau<-function(tau,D_est,xdf){
  q<-ncol(D_est)
  ee<-eigen(D_est)
  ee$values[ee$values<1e-4]<-1e-4
  ee$values[ee$values>1e4]<-1e4
  lm<-mean(ee$values)
  
  li<-ee$values+tau*(lm-ee$values)
  psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  
  nu=2*q-1
  
  d2222<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi=psi,nu=nu,const=1e8,param="variance",link_fun = function(x) 1/(1+exp(-x)) )
  
  
  
  Ya<-c(xdf$Y,d2222$data$Y)
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(d2222$data$Z)))
  Za<-rbind(xdf$Z,d2222$data$Z)
  groupa<-c(xdf$grouping,max(xdf$grouping)+d2222$data$grouping)
  weights<-c(rep(1,length(xdf$Y)),d2222$data$nn)
  
  xdfa<-list(Y=Ya,X=Xa,Z=Za,grouping=groupa,weights=weights)
  
  tmp2 <- glmmTMB(Y~-1+X+(-1+Z|grouping), family = binomial(link = "logit"),
                  weights = weights,
                  data=xdfa)
  tmp2
}

get_marLik_binom<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family=binomial(link = "logit"),
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                           theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) #point estimates seem ok, but logLik is NA! They have a trick where they dont want to report loogLik if the model does not converge (which in our case defacto holds), but we can still accesss it via object$fit$objective which seems to give -loglik so it should be minimized
  
  
  -tmp3$fit$objective
  
}

tau_finder_binom<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_binom_tau(tau,D_est,xdf)
  abs(get_marLik_binom(fit_tau,xdf)-get_marLik_binom(fit_ml,xdf))-qchisq(1-alpha,1)/2
}

######Poisson functions


fiter_pois_tau<-function(tau,D_est,xdf){
  q<-ncol(D_est)
  ee<-eigen(D_est)
  ee$values[ee$values<1e-4]<-1e-4
  ee$values[ee$values>1e4]<-1e4
  lm<-mean(ee$values)
  
  li<-ee$values+tau*(lm-ee$values)
  psi<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  
  nu=2*q-1
  
  d2222<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi=psi,nu=nu,const=1e8,param="variance",link_fun = function(x) exp(x) )
  
  
  
  Ya<-c(xdf$Y,d2222$data$Y)
  Xa<-rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(d2222$data$Z)))
  Za<-rbind(xdf$Z,d2222$data$Z)
  groupa<-c(xdf$grouping,max(xdf$grouping)+d2222$data$grouping)
  weights<-c(rep(1,length(xdf$Y)),d2222$data$nn)
  
  
  Ya2<-floor(Ya*weights)
  offset<-log(weights)
  
  xdfa<-list(Y=Ya2,ofset=offset,X=Xa,Z=Za,grouping=groupa)
  
  
  tmp2 <- glmmTMB(Y~-1+X+(-1+Z|grouping), family = poisson, offset=ofset, data=xdfa)
  
   
  tmp2
}



 
get_marLik_pois<-function(fited_model,xdf){
  
  tmp2<-fited_model
  
  tmp3<-glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family=poisson(link = "log"),
                
                start=list(beta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="beta")],
                           theta=tmp2$sdr$par.fixed[which(names(tmp2$sdr$par.fixed)=="theta")]),
                control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")) 
  
  
  -tmp3$fit$objective
  
}


 
tau_finder_pois<-function(tau,xdf,D_est,fit_ml,alpha=0.05){
  fit_tau<-fiter_pois_tau(tau,D_est,xdf)
  abs(get_marLik_pois(fit_tau,xdf)-get_marLik_pois(fit_ml,xdf))-qchisq(1-alpha,1)/2
}





######## Supporting functions for data generation in simulation study #######

# Generate grouping variables and determine the number of observations in every cluster
generate_grouping <-
  function(
    # number of clusters for each level of random effects
    num_cluster,
    # a scalar or vector of length(num_cluster);
    # a fixed number of observations within clusters or the expected value for Poisson variable
    num_obs,
    # should the number of observations be sampled from Poisson distribution?
    # In that case, the num_obs is used as lambda in zero-truncated poisson
    pois = FALSE) {
    N <- num_cluster
    
    if (!(length(num_obs) == 1 | length(num_obs) == num_cluster)) {
      stop("Parameter num_obs should be of length 1 or num_cluster")
    }
    
    if (pois) {
      flag <- TRUE
      while (flag) {
        n <- rpois(N, num_obs)
        flag <- sum(n >= 1) != N
      }
    } else{
      n <- num_obs
    }
    
    gr <- cbind(cl = 1:N, n_j = n)
    return(gr)
  }


# Generate fixed effects model matrix according to Geroldinger et. al. (2022)

generate_fixed_effects <-
  function(
    # matrix of grouping factors and number of observations
    grouping,
    # integer scalar; how many columns in X will be considered cluster level
    cluster_level = NULL,
    # correlation matrix of covariates in X; should be identity in case of independent variables
    corr,
    # list of functions to transform to transform the variables
    trans) {
    # if (ncol(corr) < cluster_level)
    #   stop(
    #     "Number of variables according to the covariate matrix and cluster_level argument does not match."
    #   )
    
    chol_cov_mat <- chol(corr)
    N <- max(grouping[, 1])
    n <- grouping[, 2]
    
    # generate cluster level variables (they have constant values for observations within the same cluster)
    if (cluster_level > 0)
      z_cl <-
      matrix(rnorm(N * cluster_level), ncol = cluster_level)
    # generate within level variables (the values for observations differ also within cluster)
    if ((ncol(corr) - cluster_level) > 0)
      z_with <-
      matrix(rnorm(sum(n) * (ncol(corr) - cluster_level)), ncol = (ncol(corr) -
                                                                     cluster_level)) # different within cluster
    
    # cluster level variables are always at the beginning
    z_x <- cbind(z_cl[rep(1:N, n), ], z_with)
    # introduce the covariance between covariates
    Z_x <- z_x %*% chol_cov_mat
    
    # transform the standard normal variables according to defined transformations
    x <- lapply(1:ncol(Z_x), function(x)
      trans[[x]](Z_x[, x]))
    X <- cbind(1, do.call(cbind, x))
    
    return(X)
  }


### dimension of the random effects covariance matrix
d = 2

## define covariance matrix between standard normal distributions between fixed effects variables according to Geroldinger et. al. (2022)
q = 5
cov_mat <- matrix(0, ncol = q, nrow = q)
cov_mat[upper.tri(cov_mat)] <- c(0.5,
                                 0.5, 0,
                                 0, 0.5, 0,
                                 0, 0, -0.3, 0.5)
cov_mat <- cov_mat + t(cov_mat)
diag(cov_mat) <- 1

## define transformations of standard normal distributions for fixed effects variables according to Geroldinger et. al. (2022)
trans_x <- list(
  x1 = function(x) {
    as.integer(x < (0))
  },
  x2 = function(x) {
    as.integer(x < (0))
  },
  x3 = function(x) {
    as.integer(x < 0)
  },
  x4 = function(x) {
    as.integer(x >= 0.5) +
      as.integer(x >= 1.5)
  },
  x5 = function(x) {
    a <- 10 * x + 55
    q <- quantile(a, c(0.25, 0.75))
    iqr <- abs(diff(q))
    ifelse(a < q[1] - 3 * iqr, q[1] - 3 * iqr,
           ifelse(a > q[2] + 3 * iqr,  a > q[2] + 3 * iqr, a))
  }
)



get_param<-function(fit,b0,betas){
  
  if (class(fit)=="try-error"){rep(NA,6+3+1+6+6+6)} else {
  fixef<-fixef(fit)$cond
  vrcv<-VarCorr(fit)$cond$group[1:2,1:2]
  
  vars<-c(vrcv[1,],vrcv[2,2])
  rho<-cov2cor(vrcv)[1,2]
  
  cib0<-confint(fit)[1,1:2]
  
  if (sum(is.na(cib0))==0){
    hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
    w_b0<-cib0[2]-cib0[1]
  } else {
    hit_b0<- w_b0<-NA
  }
  
  hits<-wdths<-rep(NA,5)
  for (ii in 1:5){
    cib1<-confint(fit)[ii+1,1:2]
    if (sum(is.na(cib1))==0){
      hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
      wdths[ii]<-cib1[2]-cib1[1]
    }
    
  }
  
  c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,sqrt(diag(vcov(fit)$cond)))
  }
}




get_param_bglmer<-function(fit,b0,betas){
  
  if (class(fit)=="try-error"){rep(NA,6+3+1+6+6+6)} else {
    fixef<-fixef(fit)
    vrcv<-VarCorr(fit)$group[1:2,1:2]
    
    vars<-c(vrcv[1,],vrcv[2,2])
    rho<-cov2cor(vrcv)[1,2]
    
    ses<-sqrt(diag(vcov(fit)))
    
    make_wald_ci<-function(est,se,alpha=0.05){
      cr<-qnorm(1-alpha/2)
      cbind(est-cr*se,  est+cr*se )
    }
    cis<-make_wald_ci(fixef,ses)
    
    cib0<-cis[1,1:2]
    
    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }
    
    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-cis[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }
      
    }
    
    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,ses)
  }
}

 


funi_bin<-function(ii,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha) {    
  
  gr <- generate_grouping(num_cluster,num_subj , pois = TRUE)
  
  flag=FALSE
  while(flag==FALSE){
  X <- generate_fixed_effects(
    grouping = gr,
    cluster_level = 2,
    corr = cov_mat,
    trans = trans_x
  )
  
  X[,6]<-(X[,6]-mean(X[,6]))/sd(X[,6])
  if (length(unique(X[,2]))==1|length(unique(X[,3]))==1) flag=FALSE else flag=TRUE
    
  }
  
  
  n<-num_cluster
  
  for (i in 1:n){
    idi<-rep(gr[i,1],gr[i,2])
    if(i==1) id<-idi else id<-c(id,idi)
  }
  for (i in 1:n){
    ni<-gr[i,2]
    
    b<-rmvnorm(1, mean = rep(0, 2), sigma = corr_mat)
    
    
    lpfi<-X[which(id==i),]%*%matrix(c(b0,betas),ncol=1)
    lpri<-rep(b[1],ni)+X[which(id==i),6]*rep(b[2],ni)
    
    lp<-lpfi+lpri
    
    
    y<-rbinom(ni,size=1,prob=1/(1+exp(-lp)))
    
    
    dfi<-data.frame(y=y,x1=X[which(id==i),2],
                    x2=X[which(id==i),3],
                    x3=X[which(id==i),4],
                    x4=X[which(id==i),5],
                    x5=X[which(id==i),6],id=id[which(id==i)])
    if (i==1) df<-dfi else df<-rbind(df,dfi)
  }
  
  nms<-paste("bin",num_cluster,num_subj,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
  
  #saveRDS(df, file=paste0("data/",ii,nms,".Rda"))
  
  X<-model.matrix(~x1+x2+x3+x4+x5,data=df)
  Z<-model.matrix(~x5,data=df)
  
  xdf<-list(Y=df$y,X=X,Z=Z,grouping=df$id)  
  
  
  
  
  fit_bin <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family=binomial(link="logit")),silent=TRUE)
  
  
  
  #REML
  fit_bin_r <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family="binomial",REML=TRUE),silent=TRUE)
  
   
  
  #ML
  
  if (class(fit_bin)=="try-error") D<-diag(1,2,2) else D<-VarCorr(fit_bin)$cond$group[1:2,1:2]  
  
   
  opt_tau<-try(uniroot(tau_finder_binom,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1 }
  
  opt_ml_tau<-opt_tau$root
  
  fit_tau_ml<-try(fiter_binom_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
   
  #oracle, tau
  
  D<-corr_mat
  
  
  
   
  opt_tau<-try(uniroot(tau_finder_binom,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1} 
  opt_oracle_tau<-opt_tau$root
  
  fit_tau_oracle<-try(fiter_binom_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
   
  #oracle
  
  fit_oracle<-try(fiter_binom_tau(0,D_est=D,xdf=xdf),silent=TRUE)
  
  
  #####
  
  #bglmer
  
  fit_bglmer<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,family=binomial(link="logit")),silent=TRUE)
  
   
  
 
  
  #####
 
  fit_ml<-fit_bin
  fit_reml<-fit_bin_r
  
  resml<-get_param(fit_ml,b0=b0,betas=betas)
  resreml<-get_param(fit_reml,b0=b0,betas=betas)
   restauml<-get_param(fit_tau_ml,b0=b0,betas=betas)
  restauoracle<-get_param(fit_tau_oracle,b0=b0,betas=betas)
  
  
  resoracle<-get_param(fit_oracle,b0=b0,betas=betas)
  
  resbglmer<-get_param_bglmer(fit_bglmer,b0=b0,betas=betas)
  
  
  
    
  
   
  
  write(resml,paste0("results/ml",nms,".txt"),ncolumns=length(resml),append = TRUE,sep="\t")
  write(resreml,paste0("results/reml",nms,".txt"),ncolumns=length(resreml),append = TRUE,sep="\t")
  write(resbglmer,paste0("results/bml",nms,".txt"),ncolumns=length(resbglmer),append = TRUE,sep="\t")
  
    write(restauml,paste0("results/tml",nms,".txt"),ncolumns=length(restauml),append = TRUE,sep="\t")
  write(restauoracle,paste0("results/to",nms,".txt"),ncolumns=length(restauoracle),append = TRUE,sep="\t")
  
  
   
  write(resoracle,paste0("results/o",nms,".txt"),ncolumns=length(resoracle),append = TRUE,sep="\t")
  
  taus<-c(opt_ml_tau,opt_oracle_tau)
  
  write(taus,paste0("results/taus",nms,".txt"),ncolumns=length(taus),append = TRUE,sep="\t")
  
  print(ii)
}




funi_pois<-function(ii,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha) {    
  
  gr <- generate_grouping(num_cluster,num_subj , pois = TRUE)
  
  flag=FALSE
  while(flag==FALSE){
    X <- generate_fixed_effects(
      grouping = gr,
      cluster_level = 2,
      corr = cov_mat,
      trans = trans_x
    )
    
    X[,6]<-(X[,6]-mean(X[,6]))/sd(X[,6])
    if (length(unique(X[,2]))==1|length(unique(X[,3]))==1) flag=FALSE else flag=TRUE
    
  }
  
  
  n<-num_cluster
  
  for (i in 1:n){
    idi<-rep(gr[i,1],gr[i,2])
    if(i==1) id<-idi else id<-c(id,idi)
  }
  for (i in 1:n){
    ni<-gr[i,2]
    
    b<-rmvnorm(1, mean = rep(0, 2), sigma = corr_mat)
    
    
    lpfi<-X[which(id==i),]%*%matrix(c(b0,betas),ncol=1)
    lpri<-rep(b[1],ni)+X[which(id==i),6]*rep(b[2],ni)
    
    lp<-lpfi+lpri
    
    
    y<-rpois(ni,lambda=exp(lp))
    
    
    dfi<-data.frame(y=y,x1=X[which(id==i),2],
                    x2=X[which(id==i),3],
                    x3=X[which(id==i),4],
                    x4=X[which(id==i),5],
                    x5=X[which(id==i),6],id=id[which(id==i)])
    if (i==1) df<-dfi else df<-rbind(df,dfi)
  }
  
  nms<-paste("po",num_cluster,num_subj,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
  
  #saveRDS(df, file=paste0("data/",ii,nms,".Rda"))
  
  X<-model.matrix(~x1+x2+x3+x4+x5,data=df)
  Z<-model.matrix(~x5,data=df)
  
  xdf<-list(Y=df$y,X=X,Z=Z,grouping=df$id)  
  
  
  
  
  fit_bin <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family=poisson(link="log")),silent=TRUE)
  
  
  
  #REML
  fit_bin_r <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family=poisson(link="log"),REML=TRUE),silent=TRUE)
  
    
  
  #ML
  
  if (class(fit_bin)=="try-error") D<-diag(1,2,2) else D<-VarCorr(fit_bin)$cond$group[1:2,1:2]  
  
  
  opt_tau<-try(uniroot(tau_finder_pois,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1 }
  
  opt_ml_tau<-opt_tau$root
  
  fit_tau_ml<-try(fiter_pois_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
   
  #oracle, tau
  
  D<-corr_mat
  
  
  
   
  opt_tau<-try(uniroot(tau_finder_pois,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1} 
  opt_oracle_tau<-opt_tau$root
  
  fit_tau_oracle<-try(fiter_pois_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
   
  #oracle
  
  fit_oracle<-try(fiter_pois_tau(0,D_est=D,xdf=xdf),silent=TRUE)
  
  
  #####
  
  #bglmer
  
  fit_bglmer<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,family=poisson(link="log")),silent=TRUE)
  
  
  
  
  #####
  
  fit_ml<-fit_bin
  fit_reml<-fit_bin_r
  
  resml<-get_param(fit_ml,b0=b0,betas=betas)
  resreml<-get_param(fit_reml,b0=b0,betas=betas)
   restauml<-get_param(fit_tau_ml,b0=b0,betas=betas)
  restauoracle<-get_param(fit_tau_oracle,b0=b0,betas=betas)
  
   resoracle<-get_param(fit_oracle,b0=b0,betas=betas)
  
  resbglmer<-get_param_bglmer(fit_bglmer,b0=b0,betas=betas)
  
  write(resml,paste0("results/ml",nms,".txt"),ncolumns=length(resml),append = TRUE,sep="\t")
  write(resreml,paste0("results/reml",nms,".txt"),ncolumns=length(resreml),append = TRUE,sep="\t")
  write(resbglmer,paste0("results/bml",nms,".txt"),ncolumns=length(resbglmer),append = TRUE,sep="\t")
  
   write(restauml,paste0("results/tml",nms,".txt"),ncolumns=length(restauml),append = TRUE,sep="\t")
  write(restauoracle,paste0("results/to",nms,".txt"),ncolumns=length(restauoracle),append = TRUE,sep="\t")
  
  
   
  write(resoracle,paste0("results/o",nms,".txt"),ncolumns=length(resoracle),append = TRUE,sep="\t")
  
  taus<-c(opt_ml_tau,opt_oracle_tau)
  
  write(taus,paste0("results/taus",nms,".txt"),ncolumns=length(taus),append = TRUE,sep="\t")
  
  print(ii)
}


##lmm

get_param_lmm<-function(fit,b0,betas){
  
  if (class(fit)=="try-error") {rep(NA,6+3+1+6+6+6+1)} else {
  fixef<-fixef(fit)$cond
  vrcv<-VarCorr(fit)$cond$group[1:2,1:2]
  
  vars<-c(vrcv[1,],vrcv[2,2])
  rho<-cov2cor(vrcv)[1,2]
  
  cib0<-confint(fit)[1,1:2]
  
  if (sum(is.na(cib0))==0){
    hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
    w_b0<-cib0[2]-cib0[1]
  } else {
    hit_b0<- w_b0<-NA
  }
  
  hits<-wdths<-rep(NA,5)
  for (ii in 1:5){
    cib1<-confint(fit)[ii+1,1:2]
    if (sum(is.na(cib1))==0){
      hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
      wdths[ii]<-cib1[2]-cib1[1]
    }
    
  }
  
  sigma(fit)
  c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,sqrt(diag(vcov(fit)$cond)),sigma(fit))
  }
}



get_param_bglmer_lmm<-function(fit,b0,betas){
  
  if (class(fit)=="try-error"){rep(NA,6+3+1+6+6+6+1)} else {
    fixef<-fixef(fit)
    vrcv<-VarCorr(fit)$group[1:2,1:2]
    
    vars<-c(vrcv[1,],vrcv[2,2])
    rho<-cov2cor(vrcv)[1,2]
    
    ses<-sqrt(diag(vcov(fit)))
    
    make_wald_ci<-function(est,se,alpha=0.05){
      cr<-qnorm(1-alpha/2)
      cbind(est-cr*se,  est+cr*se )
    }
    cis<-make_wald_ci(fixef,ses)
    
    cib0<-cis[1,1:2]
    
    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }
    
    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-cis[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }
      
    }
    
    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,ses,sigma(fit))
  }
}


funi_lmm<-function(ii,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha,sigma_eps) {    
  
  gr <- generate_grouping(num_cluster,num_subj , pois = TRUE)
  
  flag=FALSE
  while(flag==FALSE){
    X <- generate_fixed_effects(
      grouping = gr,
      cluster_level = 2,
      corr = cov_mat,
      trans = trans_x
    )
    
    X[,6]<-(X[,6]-mean(X[,6]))/sd(X[,6])
    if (length(unique(X[,2]))==1|length(unique(X[,3]))==1) flag=FALSE else flag=TRUE
    
  }
  
  
  
  
  n<-num_cluster
  
  for (i in 1:n){
    idi<-rep(gr[i,1],gr[i,2])
    if(i==1) id<-idi else id<-c(id,idi)
  }
  for (i in 1:n){
    ni<-gr[i,2]
    
    b<-rmvnorm(1, mean = rep(0, 2), sigma = corr_mat)
    
    
    lpfi<-X[which(id==i),]%*%matrix(c(b0,betas),ncol=1)
    lpri<-rep(b[1],ni)+X[which(id==i),6]*rep(b[2],ni)
    
    lp<-lpfi+lpri
    
    
    y<-lp+rnorm(ni,sd=sigma_eps)
    
    
    dfi<-data.frame(y=y,x1=X[which(id==i),2],
                    x2=X[which(id==i),3],
                    x3=X[which(id==i),4],
                    x4=X[which(id==i),5],
                    x5=X[which(id==i),6],id=id[which(id==i)])
    if (i==1) df<-dfi else df<-rbind(df,dfi)
  }
  
  nms<-paste("lmm",num_cluster,num_subj,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
  
  #saveRDS(df, file=paste0("data/",ii,nms,".Rda"))
  
  
  X<-model.matrix(~x1+x2+x3+x4+x5,data=df)
  Z<-model.matrix(~x5,data=df)
  
  xdf<-list(Y=df$y,X=X,Z=Z,grouping=df$id)  
  
  
  
  
  fit_bin <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family=gaussian(link="identity")),silent=TRUE)
  
  
  
  #REML
  fit_bin_r <- try(glmmTMB(Y~-1+X+(-1+Z|grouping), data=xdf, family=gaussian(link="identity"),REML=TRUE),silent=TRUE)
  
   
  
  #ML
  
  if (class(fit_bin)=="try-error") D<-diag(1,2,2) else D<-VarCorr(fit_bin)$cond$group[1:2,1:2] 
  
    
  opt_tau<-try(uniroot(tau_finder,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1 }
  
  opt_ml_tau<-opt_tau$root
  
  fit_tau_ml<-try(fiter_lin_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
    
  #oracle, tau
  
  D<-corr_mat
  
  
  
    
  opt_tau<-try(uniroot(tau_finder,c(0,1),xdf=xdf,D_est=D,fit_ml=fit_bin,alpha=alpha),silent=TRUE)
  if (class(opt_tau)=="try-error") {opt_tau<-list();opt_tau$root<-1} 
  opt_oracle_tau<-opt_tau$root
  
  fit_tau_oracle<-try(fiter_lin_tau(opt_tau$root,D_est=D,xdf=xdf),silent=TRUE)
  
 #oracle
  
  fit_oracle<-try(fiter_lin_tau(0,D_est=D,xdf=xdf),silent=TRUE)
  
  
  #bglmer
  
  fit_bglmer<-try(blmer(Y~-1+X + (-1+Z | grouping), data = xdf),silent=TRUE)
  
  #####
  
  
  fit_ml<-fit_bin
  fit_reml<-fit_bin_r
  
  resml<-get_param_lmm(fit_ml,b0=b0,betas=betas)
  resreml<-get_param_lmm(fit_reml,b0=b0,betas=betas)
  restauml<-get_param_lmm(fit_tau_ml,b0=b0,betas=betas)
  restauoracle<-get_param_lmm(fit_tau_oracle,b0=b0,betas=betas)
  
   
  resoracle<-get_param_lmm(fit_oracle,b0=b0,betas=betas)
  
  resbglmer<-get_param_bglmer_lmm(fit_bglmer,b0=b0,betas=betas)
  
  write(resml,paste0("results/ml",nms,".txt"),ncolumns=length(resml),append = TRUE,sep="\t")
  write(resreml,paste0("results/reml",nms,".txt"),ncolumns=length(resreml),append = TRUE,sep="\t")
  write(resbglmer,paste0("results/bml",nms,".txt"),ncolumns=length(resbglmer),append = TRUE,sep="\t")
  
   write(restauml,paste0("results/tml",nms,".txt"),ncolumns=length(restauml),append = TRUE,sep="\t")
  write(restauoracle,paste0("results/to",nms,".txt"),ncolumns=length(restauoracle),append = TRUE,sep="\t")
  
  
     
  write(resoracle,paste0("results/o",nms,".txt"),ncolumns=length(resoracle),append = TRUE,sep="\t")
  
  taus<-c(opt_ml_tau,opt_oracle_tau)
  
  write(taus,paste0("results/taus",nms,".txt"),ncolumns=length(taus),append = TRUE,sep="\t")
  
  print(ii)
  
}

 

####################
####################




##sim params, fixed!
sd_slope<-sd_int*multiplier_int

b0 <- -0.1

# regression coefficients of fixed effects
betas <- c(0.17, -0.17, 0.17, 0.15, -0.1)




# list of covariance matrices for random effects
corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2) 

sigma_eps=2

alpha=0.05


B=1000

if (model=="lmm"){
  
  m=0
  
  while(m<B){
    
    m=m+1
    
    dd<-try( funi_lmm(m,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha,sigma_eps)  ,silent=TRUE)
    
    if (class(dd)=="try-error") m<-m-1
    
  }
  
}
  
if (model=="bin"){
  
  m=0
  
  while(m<B){
    
    m=m+1
    
    dd<-try( funi_bin(m,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha)  ,silent=TRUE)
    
    if (class(dd)=="try-error") m<-m-1
    
  }
  
}

if (model=="pois"){
  
  m=0
  
  while(m<B){
    
    m=m+1
    
    dd<-try( funi_pois(m,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas,alpha)  ,silent=TRUE)
    
    if (class(dd)=="try-error") m<-m-1
    
  }
  
}