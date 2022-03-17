## Supplementary information 2
# Fit to the simulated data. 
# Adapted from Appendix S1 from Hefley et al. 2017, Ecology Letters.
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.12763&file=ele12763-sup-0001-AppendixS1.pdf
# Alterations for this paper denoted with triple hash: ###

library(raster) 
library(fields) 
library(RColorBrewer) 
library(rasterVis) 
library(truncnorm) 
library(mvtnorm) 
library(gridExtra) 
library(fBasics) 
library(coda)
library("matrixcalc")
library("epsiwal")

##########################################################################

fit.model.mcmc_XtypesOfstartingPoints <- function(n.iter, print.iter, thin.u, save.u, 
                                                  alpha.start, beta.start, gamma.start, theta.start, phi.start,
                                                  alpha.prior.var, beta.prior.var, gamma.prior.var, theta.prior, phi.prior,  
                                                  alpha.tune, beta.tune, gamma.tune, theta.tune, phi.tune, 
                                                  y, X, K, t.stop, t.keep, dt, spatial.covariates, diffusion.coef, growth.coef, d, 
                                                  us.fact){
  
  ### Number of points of treatment access
  Nd                    <- nrow(d)
  
  # Link function 
  link                  <- function(nu){ptruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)}
  
  # Broad-scale grid cells 
  us.cells              <- aggregate(spatial.covariates$cell,fact=us.fact,FUN=mean) 
  us.cells[]            <- 1:length(us.cells[]) 
  us.res                <- res(us.cells)
  
  # Fine-scale grid cells 
  ds.cells              <- spatial.covariates$cell
  
  # First-order neighbourhood matrix 
  NN                    <- neighbourhood(us.cells)
  
  ### Distance matrix for initial state dist.d1.indiv <- matrix(NA, q, N.d1)
  D                     <- matrix(NA, q, Nd)
  for (kH in 1:Nd){
    dist.d              <- raster(, nrows = q^0.5, ncols = q^0.5, 
                                  xmn = 0, xmx = 1, 
                                  ymn = 0, ymx = 1, crs = NA) 
    D[, kH]             <- rdist(data.frame(SpatialPoints(dist.d)), matrix(d[kH,], 1, 2)) 
  }
  
  # Spatial covariates z(s) and w(s) 
  Z                     <- as.matrix(spatial.covariates[[which(diffusion.coef == 1)]]) 
  W                     <- as.matrix(spatial.covariates[[which(growth.coef    == 1)]])
  
  ### Variables that will be saved - including resistance hotspot magnitudes and dispersals
  theta                 <- matrix(, n.iter+1, Nd)
  phi                   <- matrix(, n.iter+1, Nd)
  alpha                 <- matrix(, n.iter+1, sum(diffusion.coef)) 
  beta                  <- matrix(, n.iter+1, dim(X)[2]) 
  gamma                 <- matrix(, n.iter+1, sum(growth.coef)) 
  accept                <- matrix(, n.iter+1, 4)
  
  if(save.u==TRUE){ 
    u.save              <- matrix(, n.iter/thin.u+1, max(spatial.covariates$cell[]) * length(t.keep)); 
    keep.u              <- seq(thin.u, n.iter, by=thin.u)}
  
  theta[1,]             <- theta.start 
  phi[1,]               <- phi.start 
  alpha[1,]             <- alpha.start 
  beta[1,]              <- beta.start 
  gamma[1,]             <- gamma.start 
  
  ### Names for variables - including resistance hotspot magnitudes and dispersals
  colnames_theta        <- c()
  colnames_phi          <- c()
  for(i in 1:Nd){
    colnames_theta[i]   <- c(paste0("theta", i))
    colnames_phi[i]     <- c(paste0("phi", i))
  }
  colnames(theta)       <- colnames_theta
  colnames(phi)         <- colnames_phi
  colnames(alpha)       <- names(spatial.covariates)[which(diffusion.coef==1)] 
  colnames(beta)        <- colnames(X) 
  colnames(gamma)       <- names(spatial.covariates)[which(growth.coef==1)] 
  colnames(accept)      <- c("accept.rate.thetas.phis",
                             "accept.rate.alpha", "accept.rate.beta","accept.rate.gamma")
  # Begin MCMC loop 
  for(i in 1:n.iter){
    # Sample u(s,t) on fine-scale grid 
    nu                  <- exp(X %*% beta[i,]) 
    mu                  <- exp(Z %*% alpha[i,]) 
    lambda              <- W %*% gamma[i,] 
    ds.cells$mu         <- mu 
    ds.cells$lambda     <- lambda/mu 
    mu.bar              <- aggregate(ds.cells$mu, fact = us.fact, 
                           fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))}) 
    lambda.bar          <- mu.bar*aggregate(ds.cells$lambda,fact=us.fact, 
                                   fun=function(x,na.rm){(mean(x,na.rm=TRUE))}) 
    H                   <- propagator.plain(NN, mu.bar[], lambda.bar[], us.res[1], us.res[2], dt=dt) 
    ### u0 based on distance from resistance hotspots
    ds.cells$u0         <- (exp(-D[,1]^2 / phi[i,1]^2) / sum(exp(-D[,1]^2 / phi[i,1]^2))) * theta[i,1] 
    for (kH in 2:Nd){
      ds.cells$u0[]     <- ds.cells$u0[] + (exp(-D[,kH]^2 / phi[i,kH]^2) / sum(exp(-D[,kH]^2 / phi[i,kH]^2))) * theta[i,kH] 
    }
    us.cells$c0         <- extract(ds.cells$mu * ds.cells$u0, SpatialPoints(us.cells)) 
    c.all               <- calc.c(H, us.cells$c0, t.stop, t.keep) 
    u.all               <- disaggregate(c.all, us.fact) / ds.cells$mu 
    u                   <- vec(u.all[]) 
    if(length(which(keep.u==i))==1 & save.u==TRUE){u.save[which(keep.u==i)+1,] <- u} 
    u                   <- u[K]
    
    ### Sample phi and theta for Nd hotspots
    theta.star          <- c()
    phi.star            <- c()
    for (kH in 1:Nd){
      theta.star[kH]    <- rnorm(1, theta[i,kH], theta.tune)  
      phi.star[kH]      <- rnorm(1, phi[i,kH],   phi.tune)  
    }
    if(min(phi.star, theta.star) > 0){ 
      ### u0 based on distance from resistance hotspots
      ds.cells$u0.star     <- (exp(-D[,1]^2 / phi.star[1]^2) / sum(exp(-D[,1]^2 / phi.star[1]^2))) * theta.star[1] 
      for (kH in 2:Nd){
        ds.cells$u0.star[] <- ds.cells$u0.star[] + (exp(-D[,kH]^2/phi.star[kH]^2)/sum(exp(-D[,kH]^2/phi.star[kH]^2))) * theta.star[kH]
      }
      us.cells$c0.star     <- extract(ds.cells$mu * ds.cells$u0.star, SpatialPoints(us.cells)) 
      c.all.star           <- calc.c(H, us.cells$c0.star, t.stop, t.keep) 
      u.all.star           <- disaggregate(c.all.star,us.fact) / ds.cells$mu 
      u.star               <- vec(u.all.star[])[K]
      mh1.sum              <- 0
      mh2.sum              <- 0
      for (kH in 1:Nd){
        mh1.sum            <- log(dtruncnorm(theta.star[kH], mean=0, sd=theta.prior^0.5, a=0, b=Inf)) 
        mh2.sum            <- log(dtruncnorm(theta[i,kH],    mean=0, sd=theta.prior^0.5, a=0, b=Inf))
      }
      mh1                  <- sum(dbinom(y, 1, link(nu*u.star), log=TRUE)) + mh1.sum 
      mh2                  <- sum(dbinom(y, 1, link(nu*u), log=TRUE))      + mh2.sum 
      mh                   <- exp(mh1 - mh2)
    } else {mh = 0} 
    if(mh > runif(1)){
      phi[i+1,]            <- phi.star; 
      theta[i+1,]          <- theta.star; 
      us.cells$c0          <- us.cells$c0.star; 
      accept[i+1,1]        <- 1; 
      u                    <- u.star
    } else { 
      phi[i+1,]           <- phi[i,]; 
      theta[i+1,]         <- theta[i,]; 
      accept[i+1,1]       <- 0}
    
    #Sample alpha
    alpha.star          <- rmvnorm(1,alpha[i,],diag(alpha.tune,dim(Z)[2])) 
    mu.star             <- exp(Z %*% t(alpha.star)) 
    lambda              <- W %*% gamma[i,] 
    ds.cells$mu         <- mu.star 
    ds.cells$lambda     <- lambda/mu.star 
    mu.bar              <- aggregate(ds.cells$mu, fact=us.fact, 
                           fun = function(x, na.rm){(1/mean(1/x, na.rm=TRUE))}) 
    lambda.bar          <- mu.bar*aggregate(ds.cells$lambda,fact=us.fact, 
                                    fun=function(x,na.rm){(mean(x,na.rm=TRUE))}) 
    H.star              <- propagator.plain(NN,mu.bar[],lambda.bar[],us.res[1],us.res[2],dt=dt) 
    c.all.star          <- calc.c(H.star,us.cells$c0,t.stop,t.keep) 
    u.all.star          <- disaggregate(c.all.star,us.fact) / ds.cells$mu 
    u.star              <- vec(u.all.star[])[K] 
    mh1                 <- sum(dbinom(y, 1, link(nu*u.star),log=TRUE)) + 
                           sum(dnorm(alpha.star, 0, alpha.prior.var^0.5, log=TRUE)) 
    mh2                 <- sum(dbinom(y, 1, link(nu*u),log=TRUE)) + 
                           sum(dnorm(alpha[i,], 0, alpha.prior.var^0.5, log=TRUE)) 
    mh                  <- exp(mh1 - mh2) 
    if(mh > runif(1)){
      alpha[i+1,]       <- alpha.star; 
      accept[i+1,2]     <- 1; 
      u                 <- u.star
      } else { 
      alpha[i+1,]       <- alpha[i,]; 
      accept[i+1,2]     <- 0}
    
    # Sample beta 
    beta.star           <- rmvnorm(1, beta[i,], diag(beta.tune,dim(X)[2])) 
    nu.star             <- exp(X %*% t(beta.star)) 
    mh1                 <- sum(dbinom(y, 1, link(nu.star*u), log=TRUE)) +
                           sum(dnorm(beta.star,0,beta.prior.var^0.5, log=TRUE)) 
    mh2                 <- sum(dbinom(y, 1, link(nu*u), log=TRUE)) + 
                          sum(dnorm(beta[i,], 0, beta.prior.var^0.5, log=TRUE)) 
    mh                  <- exp(mh1 - mh2) 
    if(mh > runif(1)){
      beta[i+1,]        <- beta.star; 
      accept[i+1,3]     <- 1;
      nu                <- nu.star
    } else{ 
      beta[i+1,]        <- beta[i,]; 
      accept[i+1,3]     <- 0}
    
    # Sample gamma 
    gamma.star          <- rmvnorm(1, gamma[i,], diag(gamma.tune, dim(W)[2])) 
    lambda.star         <- W %*% t(gamma.star) 
    mu                  <- exp(Z %*% alpha[i+1,]) 
    ds.cells$mu         <- mu 
    ds.cells$lambda     <- lambda.star/mu 
    mu.bar              <- aggregate(ds.cells$mu,fact=us.fact, 
                           fun = function(x, na.rm){(1/mean(1/x, na.rm=TRUE))}) 
    lambda.bar          <- mu.bar*aggregate(ds.cells$lambda,fact = us.fact, 
                                   fun = function(x, na.rm){(mean(x, na.rm=TRUE))}) 
    H.star              <- propagator.plain(NN, mu.bar[], lambda.bar[], us.res[1], us.res[2], dt=dt) 
    c.all.star          <- calc.c(H.star, us.cells$c0, t.stop, t.keep) 
    u.all.star          <- disaggregate(c.all.star, us.fact) / ds.cells$mu 
    u.star              <- vec(u.all.star[])[K] 
    mh1                 <- sum(dbinom(y, 1, link(nu*u.star), log=TRUE)) + 
                           sum(dnorm(gamma.star, 0, gamma.prior.var^0.5, log=TRUE)) 
    mh2                 <- sum(dbinom(y, 1, link(nu*u),      log=TRUE)) + 
                           sum(dnorm(gamma[i,], 0, gamma.prior.var^0.5, log=TRUE)) 
    mh                  <- exp(mh1 - mh2) 
    if(mh > runif(1)){
      gamma[i+1,]       <- gamma.star; 
      accept[i+1,4]     <- 1
      } else{ 
      gamma[i+1,]       <- gamma[i,]; 
      accept[i+1,4]     <- 0}
    
    if(print.iter == TRUE){if(i%%100 == 0) print(i)} 
  } 
  if(save.u==FALSE){u.save="NA"} 
  list(accept = accept, alpha = alpha, beta = beta, gamma = gamma, phi = phi, theta = theta, u = u.save) 
}

###################################################################################

library(metaBMA)
library(truncnorm)
library(SimDesign)
library(ggplot2)
library(reshape2)

set.seed(4223) 
K                            <- sample(1:dim(df1)[1], replace = FALSE, size = 50000) 
df2                          <- df1[K, ]

spatial.covariates           <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1, 
                                       ymn = 0, ymx = 1, crs = NA) 
spatial.covariates$cell      <- c(1:q) 
spatial.covariates$intercept <- 1 
spatial.covariates$z         <- z 
spatial.covariates$w         <- w 


set.seed(1189) 
ptm <- proc.time() 
chain1 <-fit.model.mcmc_XtypesOfstartingPoints( 
  n.iter             = 5000, ### More iterations
  print.iter         = TRUE, 
  thin.u             = 10, 
  save.u             = TRUE, 
  alpha.start        = alpha, 
  beta.start         = beta, 
  gamma.start        = gamma, 
  theta.start        = theta, ###
  phi.start          = phi,   ###
  alpha.prior.var    = 10^6, 
  beta.prior.var     = 10^6, 
  gamma.prior.var    = 10^6, 
  theta.prior        = 10^6,
  phi.prior          = 10^6, 
  alpha.tune         = c(1/500,1/2500), 
  beta.tune          = 1/50, 
  gamma.tune         = c(1.7/10^5,1/340), 
  theta.tune         = 3,
  phi.tune           = 1/400, 
  y                  = df2$y, 
  X                  = model.matrix(~x-1, data=df2), 
  K                  = K, 
  t.stop             = 20, 
  t.keep             = seq(1, 20, by=1), 
  dt                 = 1, 
  spatial.covariates = spatial.covariates, 
  diffusion.coef     = c(0,1,1,0), 
  growth.coef        = c(0,1,0,1), 
  d                  = d,
  us.fact            = 5) 
proc.time() - ptm

saveRDS(chain1, file = "fittingModel_XtypesOfstartingPoints.rds")

