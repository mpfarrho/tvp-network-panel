library(GIGrvg)
library(MASS)
library(Matrix)
library(bayesm)
library(Rcpp)
library(coda)

rm(list = ls())
source("aux_functions.R")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("threshold_functions.cpp")

# -------------------------------------------------------------------------------------------
# preliminaries
# run <- as.integer(Sys.getenv("SGE_TASK_ID"))
run <- 1
nburn <- 5000
nsave <- 1000
thin <- 10
diag.W <- TRUE # whether diagonal should be included on the W matrices (if FALSE, diag(W)=0)

ntot <- nsave*thin + nburn
save.set <- seq(thin,nsave*thin,thin) + nburn
save.len <- length(save.set)
save.ind <- 0

iplot <- 1000
check.time <- TRUE # indicate progress bar

if(is.null(iplot)){
  toplot <- 0
}else{
  toplot <- seq(1,ntot,by=iplot) # plot online
}

combi <-expand.grid( sample =  "full",
                     pool =    "free", # pool.mu/pool.sig
                     spdep =   "tvp", # spatial dependence
                     tvp =      TRUE, # TVPs
                     dataset = "dTVW",
                     stringsAsFactors = F)

rownames(combi) <- 1:nrow(combi)
sl.run <- combi[run,]

# model specification
sample <- sl.run[["sample"]] # Full versus Subsamples
pool.mu <- sl.run[["pool"]] # c("all","free","shrink")
if(pool.mu == "shrink") pool.sig <- "free" else pool.sig <- sl.run[["pool"]] # c("all","free")
spdep <- sl.run[["spdep"]] # c("none","cons","tvp")
tvp <- sl.run[["tvp"]] # c(FALSE,TRUE)
sl.dataset <- sl.run[["dataset"]]

#alternative model specifications (turned-off): mixture, sv, z-transformation, imputations 
perm.sampler <- TRUE # not permutation, but identification with beta.cons_1
G <- 10 # number of mixture components
common.mean <- FALSE
sv <- FALSE
ztrans <- FALSE
imput <- FALSE

# load specification here
source("00data_script.R")

# Some error messages
if(sv && pool.sig!="free") stop("Cannot pool variances if errors feature stochastic volatility.")
if(tvp && pool.mu=="all") stop("Pooling of coefficients not possible if time-varying parameters are estimated.")

if(sl.dataset == "d1997"){
  sc.fac <- c(100, 0.0005) # prior-scaling for g-prior/upper-bound
}else if(spdep == "cons"){
  sc.fac <- c(100, 0.01) # prior-scaling for g-prior/upper-bound
}else{
  sc.fac <- c(100, 0.005) # prior-scaling for g-prior/upper-bound
}

names(sc.fac) <- c("cons", "tvp")

# ----------------------------
a_sig <- 0.01
b_sig <- 0.01

prior.mat <- rep(0,K)
d0 <- d1 <- 0.6
com.var <- 10

b_draw <- matrix(0,K,N)
bt_draw <- array(0,dim=c(T,K,N))
omega_draw <- matrix(1,K,N)

# prior-scaling (and upper bound for TVP) 
b_ols <- solve(crossprod(Xraw))%*%crossprod(Xraw, rowMeans(y))
XXinv <- diag(solve(crossprod(Xraw)))
sd.shock <- apply(Xraw, 2, sd)[[1]]
XXinv[1] <- XXinv[1]*sd.shock

sc.ub <- XXinv*sc.fac["cons"]
sc.ub.tvp <- c(sc.ub, XXinv*sc.fac["tvp"])
sc.ub.tvp[2*K] <- sc.ub.tvp[2*K]*1e-5
sc.ub.rho <- sc.fac["tvp"]
  
si <- matrix(0,N,G)
if(tvp) mu_draw <- matrix(0, 2*K, G) else mu_draw <- matrix(0, K, G)

if(pool.mu=="shrink"){
  if(tvp){
    b_prior <- cons_draw <- matrix(0,2*K,N)
    V0 <- array(0,dim=c(N,2*K,2*K))
    for(i in 1:N){
      V0[i,,] <- diag(sc.ub.tvp)
    }
    lambda.hs <- nu.hs <- matrix(1,N,2*K)
    tau.hs <- zeta.hs <- rep(1,2*K)
  }else{
    b_prior <- cons_draw <- matrix(0,K,N)
    V0 <- array(0,dim=c(N,K,K))
    for(i in 1:N){
      V0[i,,] <- diag(sc.ub)
    }
    lambda.hs <- nu.hs <- matrix(1,N,K)
    tau.hs <- zeta.hs <- rep(1,K)
  }  
}else if(pool.mu == "mix"){
  # ----------------------------
  # mixture prior
  for(ii in 1:N){
    si[ii,sample(1:G,1)] <- 1
  }
  
  sample.e0 <- TRUE
  e0_j <- e0 <- 0.1
  c_proposal <- 0.5
  accept <- 0
  
  S_b0 <- matrix(1,ncol=G)
  c0 <- 2.5 #+(T-1)/2
  g0 <- 0.5 #+(T-1)/2
  G0 <- 100*g0/c0#/diff(range(y))^2 #100*g0/c0 range(y)^2
  C0 <- 0.01
  
  c_tau <- 0.01
  d_tau <- 0.01
  
  if(tvp){
    b_prior <- cons_draw <- matrix(0,2*K,N)
    V0 <- diag(sc.ub.tvp)
    # common mean
    b0 <- matrix(c(b_ols, 0, 0),2*K)
    B0 <- 2*diag(2*K)
    B0[K+1, K+1] <- sc.ub.tvp[K+1]
    B0[2*K,2*K] <-  sc.ub.tvp[2*K]
  }else{
    b_prior <- cons_draw <- matrix(0,K,N)
    V0 <- diag(sc.ub)
    # common mean
    b0 <- matrix(b_ols,K)
    B0 <- 2*diag(K)
  }
  B0inv <- solve(B0)
  
}else{
  if(tvp){
    b_prior <- cons_draw <- matrix(0,2*K,N)
    V0 <- diag(sc.ub.tvp)
  }else{
    b_prior <- cons_draw <- matrix(0,K,N)
    V0 <- diag(sc.ub)
  }
} 

# error specification
epsilon <- eta <- matrix(0,T,N)
S <- log(matrix(0.1,T,N))

scale_rho <- 0.1
sigma_rho <- 0.001
rhoWt <- list()

tau_draw <- 10
rho_b0 <- 1e-5
rho_b1 <- 0.1

tv.ind <- 1
d0_tau <- 1/2 
s1 <- 0.5

if(spdep%in%c("tvp","cons")){
  rho_draw_c <- 0.2
}else{
  rho_draw_c <- 0
}
rho_draw <- rep(rho_draw_c,T)
rho_pars <- matrix(c(0,0.01))

kappa_draw <- 1
rho_scale <- 0.05; rho_accept <- 0
prev.acc <- FALSE

# ----------------------------
# store objects
omega_store <- b_store <- array(0,c(save.len,K,N))
bt_store <- array(NA,dim=c(save.len,T,K,N))
Sigma_store <- array(NA,dim=c(save.len,T,N))
rho_store <- array(0,dim=c(save.len,T))
sig.rho_store <- matrix(0, save.len, 1)
if(tvp) mu_store <- array(NA,dim=c(save.len,G, 2*K)) else mu_store <- array(NA,dim=c(save.len,G, K)) 
G_store <- matrix(0,save.len, 1)
si_store <- array(NA,c(save.len,N,G))

var.names <- c("mpshock","cons")
ind.names <- colnames(y)

dimnames(omega_store) <- dimnames(b_store) <- list(1:save.len,var.names,ind.names)
dimnames(bt_store) <- list(1:save.len,date.names,var.names,ind.names)
dimnames(Sigma_store) <- list(1:save.len,date.names,ind.names)
dimnames(rho_store) <- list(1:save.len,date.names)
dimnames(si_store) <- list(1:save.len,ind.names,1:G)

# -------------------------------------------------------------------------------------------
# sampling
if(check.time){
  pb <- txtProgressBar(min = 1, max = ntot, style = 3)
}

t.start <- Sys.time()
for(irep in 1:ntot){
  Sy <- matrix(0,T,N)
  if(spdep%in%c("cons","tvp")){
    for(tt in 1:T){
      if(sl.dataset=="dTVW") W <- Wt[sl.W[tt],,]
      Sy[tt,] <- (diag(N)-rho_draw[tt]*W)%*%as.numeric(y[tt,])
    }
    y_ <- Sy
  }else{
    y_ <- Sy <- y
  }
  
  # --------------------------
  # sample mean coefficients
  if(pool.mu=="all"){
    y__ <- as.vector(y_)*as.vector(exp(-S/2))
    X__ <- Xst*as.vector(exp(-S/2))
    
    V_post <- try(solve(crossprod(X__)+diag(1/diag(V0))),silent=TRUE)
    if (is(V_post,"try-error")) V_post <- ginv(crossprod(X__)+diag(1/diag(V0)))
    b_post <- V_post%*%(crossprod(X__,y__))
    
    bi_draw <- try(b_post + t(chol(V_post))%*%rnorm(K),silent=TRUE)
    if (is(bi_draw,"try-error")){bi_draw <- mvrnorm(1,b_post,V_post)}
    b_draw <- matrix(rep(bi_draw,N),ncol=N)
    
    for(i in 1:N) bt_draw[,,i] <- t(matrix(bi_draw, K, T))
    eta <- matrix(as.vector(y_) - Xst %*% bi_draw,T,N)
  }else if(pool.mu=="free"){
    if(tvp){
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        YY <- (y_[,i]-X[,i,]%*%b_draw[,i]) * normalizer
        XX <- X[,i,]%*%diag(omega_draw[,i]) * normalizer
        b_it <- KF_fast(t(as.matrix(YY)),XX,as.matrix(exp(S[,i])),t(matrix(1,K,T)),K,1,T,matrix(0,K,1),diag(K)/10^15)
        
        Xtilde <- X[,i,]*t(b_it)
        XX <- cbind(X[,i,],Xtilde)*normalizer
        YY <- y_[,i] * normalizer
        
        V0inv <- diag(1/diag(V0))
        V_post <- try(solve((crossprod(XX)) + V0inv),silent=T)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(XX) + V0inv))
        b_post <- V_post %*% (crossprod(XX, YY)+ V0inv%*%b_prior[,i])
        b_i <- try(b_post+t(chol(V_post))%*%rnorm(2*K,0,1),silent=T)
        if (is(b_i,"try-error")) b_i <- mvrnorm(1,b_post,V_post)
        
        b_draw[,i] <- b_i[1:K]
        sqrtomega <- b_i[(K+1):(2*K)]
        
        # non-identified element
        uu <- runif(1,0,1)
        if (uu>0.5){
          sqrtomega <- -sqrtomega
          b_it <- -b_it
        }
        omega_draw[,i] <- sqrtomega
        for(tt in 1:T){
          eta[tt,i] <- y_[tt,i]-X[tt,i,]%*%(as.numeric(b_draw[,i])+(sqrtomega*b_it))[,tt]
        }
        bt_draw[,,i] <- t((as.numeric(b_draw[,i])+sqrtomega*b_it))
      }
    }else{
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        XX <- X[,i,] * normalizer
        YY <- y_[,i] * normalizer
        
        V_post <- try(solve(crossprod(XX)+diag(1/diag(V0))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(XX)+diag(1/diag(V0)))
        b_post <- V_post%*%(crossprod(XX,YY)+diag(1/diag(V0)) %*% b_prior[,i])
        
        bi_draw <- try(b_post + t(chol(V_post))%*%rnorm(K),silent=TRUE)
        if (is(bi_draw,"try-error")){bi_draw <- mvrnorm(1,b_post,V_post)}
        
        b_draw[,i] <- bi_draw
        bt_draw[,,i] <- t(matrix(bi_draw, K, T))
        
        eta[,i] <- y_[,i] - (X[,i,] %*% bi_draw)
      }
    }
  }else if(pool.mu=="shrink"){
    # sampling of coefficients
    if(tvp){
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        YY <- (y_[,i]-X[,i,]%*%b_draw[,i]) * normalizer
        XX <- X[,i,]%*%diag(omega_draw[,i]) * normalizer
        b_it <- KF_fast(t(as.matrix(YY)),XX,as.matrix(exp(S[,i])),t(matrix(1,K,T)),K,1,T,matrix(0,K,1),diag(K)*1/10^15)
        
        Xtilde <- X[,i,]*t(b_it)
        XX <- cbind(X[,i,],Xtilde)*normalizer
        YY <- y_[,i] * normalizer
        
        V_post <- try(solve((crossprod(XX)) + diag(1/diag(V0[i,,]))),silent=T)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(XX) + diag(1/diag(V0[i,,]))))
        b_post <- V_post %*% (crossprod(XX, YY)+diag(1/diag(V0[i,,]))%*%b_prior[,i])
        b_i <- try(b_post+t(chol(V_post))%*%rnorm(2*K,0,1),silent=T)
        if (is(b_i,"try-error")) b_i <- mvrnorm(1,b_post,V_post)
        
        b_draw[,i] <- b_i[1:K]
        sqrtomega <- b_i[(K+1):(2*K)]
        
        # non-identified element
        uu <- runif(1,0,1)
        if (uu>0.5){
          sqrtomega <- -sqrtomega
          b_it <- -b_it
        }
        omega_draw[,i] <- sqrtomega
        for(tt in 1:T){
          eta[tt,i] <- y_[tt,i]-X[tt,i,]%*%(as.numeric(b_draw[,i])+(sqrtomega*b_it))[,tt]
        }
        bt_draw[,,i] <- t((as.numeric(b_draw[,i])+sqrtomega*b_it))
      }
    }else{
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        XX <- X[,i,] * normalizer
        YY <- y_[,i] * normalizer
        
        V_post <- try(solve(crossprod(XX)+diag(1/diag(V0[i,,]))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(XX)+diag(1/diag(V0[i,,])))
        b_post <- V_post%*%(crossprod(XX,YY)+diag(1/diag(V0[i,,])) %*% b_prior[,i])
        
        bi_draw <- try(b_post + t(chol(V_post))%*%rnorm(K),silent=TRUE)
        if (is(bi_draw,"try-error")){bi_draw <- mvrnorm(1,b_post,V_post)}
        
        b_draw[,i] <- bi_draw
        bt_draw[,,i] <- t(matrix(bi_draw, K, T))
        eta[,i] <- y_[,i] - (X[,i,] %*% bi_draw)
      }
    }
    
    # shrinkage by variable type
    if(tvp){
      Ktilde <- 2*K
      b_temp <- rbind(b_draw,omega_draw)
    }else{
      Ktilde <- K
      b_temp <- b_draw
    }
    
    for(jj in 1:Ktilde){
      #if(jj!=Ktilde){
        hs_draw <- get.hs(bdraw=b_temp[jj,],lambda.hs=lambda.hs[,jj],nu.hs=nu.hs[,jj],tau.hs=tau.hs[jj],zeta.hs=zeta.hs[jj])
        psi_draw <- hs_draw$psi
        psi_draw[psi_draw > sc.ub.tvp[jj]] <- sc.ub.tvp[[jj]] 
        V0[,jj,jj] <- psi_draw
        lambda.hs[,jj] <- hs_draw$lambda
        nu.hs[,jj] <- hs_draw$nu
        tau.hs[jj] <- hs_draw$tau
        zeta.hs[jj] <- hs_draw$zeta  
      #}else{
      #  V0[,jj,jj] <- 1e-8
      #}
    }
  }else if(pool.mu=="mix"){
    if(tvp){
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        YY <- (y_[,i]-X[,i,]%*%b_draw[,i]) * normalizer
        XX <- X[,i,]%*%diag(omega_draw[,i]) * normalizer
        b_it <- KF_fast(t(as.matrix(YY)),XX,as.matrix(exp(S[,i])),t(matrix(1,K,T)),K,1,T,matrix(0,K,1),diag(K)*1/10^15)
        
        Xtilde <- X[,i,]*t(b_it)
        XX <- cbind(X[,i,],Xtilde)*normalizer
        YY <- y_[,i] * normalizer
        
        V_post <- try(solve((crossprod(XX)) + diag(1/diag(V0))),silent=T)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(XX) + diag(1/diag(V0))))
        b_post <- V_post %*% (crossprod(XX, YY)+diag(1/diag(V0))%*%b_prior[,i])
        
        b_i <- try(b_post+t(chol(V_post))%*%rnorm(2*K,0,1),silent=T)
        if (is(b_i,"try-error")) b_i <- mvrnorm(1,b_post,V_post)
        
        b_draw[,i] <- b_i[1:K]
        sqrtomega <- b_i[(K+1):(2*K)]
        
        # non-identified element
        sign.omega <- sign(sqrtomega)
        sqrtomega <- sqrtomega*sign.omega
        b_it <- b_it*sign.omega
        
        omega_draw[,i] <- sqrtomega
        for(tt in 1:T){
          eta[tt,i] <- y_[tt,i]-X[tt,i,]%*%(as.numeric(b_draw[,i])+(sqrtomega*b_it))[,tt]
        }
        bt_draw[,,i] <- t((as.numeric(b_draw[,i])+sqrtomega*b_it))  
      }
      
      Xhat <- si
      ysi <-  rbind(b_draw,omega_draw)/diag(V0)
      XsiXsi <- crossprod(Xhat)
      Xsiysi <- crossprod(Xhat,t(ysi))
      diagXsi <- matrix(diag(XsiXsi),G,2*K) %*% diag(1/diag(V0))
      
      V_mix <- 1/(matrix(diagXsi,2*K*G,1)+matrix(rep(diag(B0inv), each = G), 2*K*G, 1))
      A_mix <- (as.vector(Xsiysi)+rep(B0inv%*%b0, each = G))*(V_mix)
      A_draw <- A_mix+rnorm(G*2*K,0,sqrt(V_mix))
      mu_draw <- t(matrix(A_draw,G,2*K))
    }else{
      for(i in 1:N){
        normalizer <- exp(-S[,i]/2)
        XX <- X[,i,] * normalizer
        YY <- y_[,i] * normalizer
        
        V_post <- try(solve(crossprod(XX)+diag(1/diag(V0))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(XX)+diag(1/diag(V0)))
        b_post <- V_post%*%(crossprod(XX,YY)+diag(1/diag(V0)) %*% b_prior[,i])
        
        bi_draw <- try(b_post + t(chol(V_post))%*%rnorm(K),silent=TRUE)
        if (is(bi_draw,"try-error")){bi_draw <- mvrnorm(1,b_post,V_post)}
        
        b_draw[,i] <- bi_draw
        bt_draw[,,i] <- t(matrix(bi_draw, K, T))
        eta[,i] <- y_[,i] - (X[,i,] %*% bi_draw)
        Xb[,i] <- X[,i,] %*% bi_draw
      }
      
      Xhat <- si
      ysi <-  b_draw/diag(sqrt(V0))
      XsiXsi <- crossprod(Xhat)
      Xsiysi <- crossprod(Xhat,t(ysi))
      diagXsi <- matrix(diag(XsiXsi),G,K) %*% diag(1/diag(V0))
      
      V_mix <- 1/(matrix(diagXsi,K*G,1)+matrix(rep(diag(B0inv), each = G), K*G, 1))
      A_mix <- (as.vector(Xsiysi)+rep(B0inv%*%b0, each = G))*(V_mix)
      A_draw <- A_mix+rnorm(G*K,0,diag(V_mix))
      mu_draw <- t(matrix(A_draw,G,K))
      if(tvp) mu_draw[(K+1):(2*K),] <- apply(mu_draw[(K+1):(2*K),], 2, abs)
    }
    
    # group intensity parameter
    ek <- apply(si,2,sum) + e0
    zeta <- bayesm::rdirichlet(ek)
    kappa <- zeta/sum(zeta)
    
    # shrinkage parameters
    if(sample.e0){
      a_gam <- 10
      b_gam <- a_gam  # e0~G(a_gam,b_gam*Kmax)
      Gmax <- G
      const <- Gmax * b_gam
      
      e0_j <- e0 # current value
      le0_p <- log(e0_j) + rnorm(1, 0, c_proposal)
      e0_p <- exp(le0_p) # proposal value
      
      kappa[kappa == 0] <- 1e-10
      lalpha1 <- (e0_p - e0_j) * sum(log(kappa)) + lgamma(G * e0_p) - lgamma(G * e0_j) - G *
        (lgamma(e0_p) - lgamma(e0_j)) + (a_gam - 1) * (log(e0_p) - log(e0_j)) - (e0_p - e0_j) * const+ log(e0_p) - log(e0_j)
      alpha1 <- min(exp(lalpha1), 1)
      alpha2 <- runif(1)
      
      if (alpha2 <= alpha1) {
        e0_j <- e0_p
        e0 <- e0_j
        accept <- accept+1
      }
      if (irep<(0.5*nburn)){
        if (accept/irep <0.2) c_proposal <- 0.99*c_proposal
        if (accept/irep >0.4) c_proposal <- 1.01*c_proposal
      }
    }
    
    # sample group indicators 
    probs <- matrix(0,N,G)
    si <- matrix(0,N,G)
    for (jj in 1:N){
      get_probs <- matrix(0,G,1)
      if(tvp){
        b_vec <- rbind(b_draw,omega_draw)[,jj]
      }else{
        b_vec <- b_draw[,jj]
      }
      
      var.vec <- diag(V0)
      var.vec <- kronecker(sqrt(exp(S_b0)), var.vec)
      
      get_probs <- apply(dnorm(b_vec,mu_draw,var.vec,log=TRUE),2,sum)+log(kappa)
      pjj <- (exp(get_probs-max(get_probs)))/sum(exp(get_probs-max(get_probs)))
      probs[jj,] <- pjj
      sample.j <- sample(seq(1,G),1,prob=pjj)
      si[jj,sample.j] <- 1
    }
    
    if(tvp){
      Ktilde <- 2*K
    }else{
      Ktilde <- K
    }
    if (common.mean){
      if(tvp){
        R <- apply(t(rbind(b_draw,omega_draw)), 2, function(x) diff(range(x)))
      }else{
        R <- apply(t(b_draw), 2, function(x) diff(range(x)))
      }
      R[R > 1] <- 1
      
      scale0 <- rep(1,Ktilde)
      for (jj in 1:Ktilde) scale0[[jj]] <- sum((mu_draw[jj,]-b0[[jj]])^2)/R[[jj]]^2
      scale0[scale0 > 10] <- 10
      v0 <- v1 <- 0.1
      aj <- 2*v1
      pk <- v0-G/2
      
      Lambda <- matrix(0,Ktilde,1)
      for (ss in 1:Ktilde) Lambda[ss,1] <- GIGrvg::rgig(1,pk,scale0[[ss]]+1e-25,aj)
      
      if (K>1){
        B0 <- diag((R^2) * as.numeric(Lambda) + 1/10)
        B0inv <- diag(1/diag(B0)) 
      }else{
        B0 <- R^2*Lambda + 1/10
        B0inv <- 1/B0
      }
      if(tvp){
        if(diag(B0)[K+1] > sc.ub.tvp[K+1]) B0[K+1,K+1] <- sc.ub.tvp[K+1]
        B0[Ktilde,Ktilde] <- sc.ub.tvp[Ktilde]
        B0inv <- diag(1/diag(B0)) 
      }
      
      #b0 <- MASS::mvrnorm(1, matrix(0,Ktilde,1), 1/G*B0)
      b0 <- MASS::mvrnorm(1, matrix(rowMeans(mu_draw), Ktilde,1), 1/G*B0)
      if(tvp) b0[(K+1):(2*K)] <- 0
      
    }
    
    if (perm.sampler){
      perm <- sort(mu_draw[1,],index.return=TRUE)$ix#sample(G)
      mu_draw <- mu_draw[,perm,drop=F]
      kappa <- kappa[perm,1,drop = F]
      si <- si[,perm]
      S_b0 <- S_b0[,perm,drop = F]
    }
    
    Sig_common <- matrix(0.1,Ktilde,1)
    for(ss in 1:Ktilde){
      scale1 <- 0
      if(tvp){
        b_temp <- rbind(b_draw,omega_draw)
      }else{
        b_temp <- b_draw        
      }
      for(kk in 1:N){
        sl.mixture <- which(si[kk,]==1)
        scale1 <- scale1+(b_temp[ss,kk]-mu_draw[ss,sl.mixture])^2
      }
      
      sig_q <- 1/rgamma(1,(N+c_tau)/2, (d_tau+scale1)/2)
      
      if(tvp){
        sig_q[sig_q > sc.ub.tvp[ss]] <- sc.ub.tvp[ss] 
        sig_q[sig_q < 0.001*sc.ub.tvp[ss]] <- 0.001*sc.ub.tvp[ss]
      }else{
        sig_q[sig_q > sc.ub[ss]] <- sc.ub[ss]
        sig_q[sig_q < 0.001*sc.ub[ss]] <- 0.001*sc.ub[[ss]]
      }
      Sig_common[ss] <- sig_q
      
    }
    if(tvp) Sig_common[Ktilde] <- 1e-8
    V0 <- diag(as.numeric(Sig_common),Ktilde)
  }
  
  
  fit <- matrix(0,T,N)
  if(tvp){
    for(i in 1:N){
      for(tt in 1:T){
        fit[tt,i] <- X[tt,i,] %*% bt_draw[tt,,i]
      }
    }
  }else{
    for(i in 1:N) fit[,i] <- X[,i,] %*% b_draw[,i]
  }

  # --------------------------
  # sample variances
  if(pool.sig=="free"){
    for(i in 1:N){
      eps <- eta[,i]
      S[,i] <- log(1/rgamma(1,a_sig + T/2,b_sig + crossprod(eps)/2))
    }
  }else if(pool.sig=="all"){
    eps <- as.vector(eta)
    c0p <- 0.01 + (T*N)/2
    C0p <- 0.01 + crossprod(eps)/2
    S[,] <- log(1/rgamma(1,c0p,C0p))
  }else if(pool.sig=="mix"){
    nk <- apply(si,2,sum)
    for(g in 1:G){
      id.sl <- si[,g]==1 #already permuted
      if(nk[g]==0){
        S[,id.sl] <- log(1/rgamma(1,c0,C0))
      }else{
        eps <- as.numeric(eta[,id.sl])
        c0p <- c0 + T*nk[g]/2
        C0p <- C0 + crossprod(eps)/2
        S[,id.sl] <- S_b0[g] <- log(1/rgamma(1,c0p,C0p))
        b_prior[,id.sl] <- mu_draw[,g]
      }
    }
    C0 <- 1/rgamma(1,g0+G*c0,G0+sum(1/exp(S_b0)))
  }
  
  if(pool.mu == "mix" && irep %in% toplot){
     print(paste0("Groups: ", paste0(nk, collapse = ", ")))
     print(paste0("Group mean - coeff 1: ", paste0(round(mu_draw[1,],3), collapse = ", ")))
     print(paste0("V0 - coeffs: ", paste0(round(Sig_common,3), collapse = ", ")))
     print(paste0("B0 - coeffs: ", paste0(round(diag(B0),3), collapse = ", ")))
     print(paste0("Group sig: ", paste0(round(exp(S_b0),3), collapse = ", ")))
  }

  # --------------------------
  # Sample the spatial lag
  if(spdep=="cons"){
    rho_prop <- rnorm(1,rho_draw_c,rho_scale)
    
    if(sl.dataset=="dTVW"){
      logprop <- get.spatial3(y,t=T,Xb=fit,Wt=Wt,sl.W=sl.W,rho=rho_prop,par="rho",S=S) # change here dependent on pooling
      logdraw <- get.spatial3(y,t=T,Xb=fit,Wt=Wt,sl.W=sl.W,rho=rho_draw_c,par="rho",S=S)
    }else{
      logprop <- get.spatial2(y,t=T,Xb=fit,W=W,rho=rho_prop,par="rho",S=S) # change here dependent on pooling
      logdraw <- get.spatial2(y,t=T,Xb=fit,W=W,rho=rho_draw_c,par="rho",S=S)
    }
    accprop <- logprop-logdraw
    
    if(is.nan(accprop)){accprop <- -Inf}
    if(accprop > log(runif(1,0,1))){
      rho_draw_c <- rho_prop
      rho_accept <- rho_accept + 1
    }
    
    if(irep < (nburn/2)){
      if((rho_accept/irep) > 0.3){
        rho_scale <- 1.01*rho_scale
      }
      if((rho_accept/irep) < 0.2){
        rho_scale <- 0.99*rho_scale
      }
    }
    rho_draw[] <- rho_draw_c
    tv.ind <- 0
  }else if(spdep=="tvp"){
    rho0 <- 0
    rho0V <- 0.02
    
    for(tt in 0:T){
      if(tt==0){
        r_sig <- (rho0V*sigma_rho)/(rho0V+sigma_rho)
        r_mu <- r_sig*(rho0/rho0V + rho_draw[1]/sigma_rho)
        r0 <- r_mu + sqrt(r_sig)*rnorm(1)
      }else{
        if(sl.dataset=="dTVW") W <- Wt[sl.W[tt],,]
        normalizer <- exp(-S[tt,]/2)
        if(tt==1){
          r_mu <- (r0 + rho_draw[1])/2
          r_sig <- sigma_rho/2
        }else if(tt==T){
          r_mu <- rho_draw[T-1]
          r_sig <- sigma_rho
        }else{
          r_mu <- (rho_draw[tt-1] + rho_draw[tt+1])/2
          r_sig <- sigma_rho/2
        }
        rho_prop <- rho_draw[tt] + sqrt(scale_rho)*rnorm(1)
        xtt <- (W %*% y[tt,])*normalizer
        ytt <- (y[tt,] - fit[tt,])*normalizer
        
        det_old <- determinant(diag(N)-rho_draw[tt]*W,logarithm = T)$modulus[[1]]
        det_new <- determinant(diag(N)-rho_prop*W,logarithm = T)$modulus[[1]]
        
        liki_old <- det_old - (0.5 * crossprod(ytt-rho_draw[tt]*xtt)) + dnorm(rho_draw[tt],mean=r_mu,sd=sqrt(r_sig),log=T)
        liki_prop <- det_new - (0.5 * crossprod(ytt-rho_prop*xtt)) + dnorm(rho_prop,mean=r_mu,sd=sqrt(r_sig),log=T)
        
        alpha <- liki_prop-liki_old
        if(is.nan(alpha)) alpha <- -Inf
        if(log(runif(1))<alpha) rho_draw[tt] <- rho_prop
      }
    }
    
    # draw variance of state equation for rho
    rho_t <- rho_draw[2:T]
    rho_t1 <- rho_draw[1:(T-1)]
    eps_rho <- (rho_t-rho_t1)^2
    
    # with shrinkage prior
    p0 <- dnorm(sqrt(sigma_rho),0,rho_b0,log=F)+1e-50
    p1 <- dnorm(sqrt(sigma_rho),0,rho_b1,log=F)+1e-50

    alpha <- p1/(p0+p1)
    if(alpha>runif(1)) tv.ind <- 1 else tv.ind <- 0
    scl <- sum(eps_rho)

    # prior on state innovation variances for rho
    rho_pars[2] <- sigma_rho <- 1/rgamma(1, 3 + (T-1)/2, 0.03 + sum(eps_rho)/2)
    if(sigma_rho > sc.ub.rho) rho_pars[2] <- sigma_rho <- sc.ub.rho
  }
  
  # --------------------------
  # store objects
  if(irep %in% save.set){
    save.ind <- save.ind + 1
    omega_store[save.ind,,] <- omega_draw
    b_store[save.ind,,] <- b_draw
    bt_store[save.ind,,,] <- bt_draw
    Sigma_store[save.ind,,] <- exp(S)
    mu_store[save.ind,,] <- mu_draw
    si_store[save.ind,,] <- si
    G_store[save.ind,] <- as.numeric(sum(apply(si,2,mean)>0))
    sig.rho_store[save.ind,] <- sigma_rho
    rho_store[save.ind,] <- rho_draw
  }
  
  if(irep %in% toplot){
    par(mfrow = c(1,1))
    ts.plot(rho_draw,main= paste0("Draw #",irep, ": rho"),ylim=c(0,1))
    # ts.plot(bt_draw[,2,],main= paste0("Draw #",irep, ": beta_0 (interc.)"))
    # ts.plot(bt_draw[,1,],main= paste0("Draw #",irep, ": beta_1"))
  }
  
  if(check.time == TRUE){
    setTxtProgressBar(pb, irep)
  }
}

t.end <- Sys.time()
time.min <- (ts(t.end)-ts(t.start))/60 # Time in minutes

est_obj <- list("spec"=sl.run,"b"=b_store,"bt"=bt_store,"sig"=Sigma_store,"mu"=mu_store,"si"=si_store,"G"=G_store,"rho"=rho_store, "rho_sig" = sig.rho_store, "time"=time.min)
if(spdep!="none"){
  if(sl.dataset=="dTVW"){
    imp_obj <- get.impacts2(b=bt_store,rho=rho_store,Wt=Wt,sl.W=sl.W) 
  }else{
    imp_obj <- get.impacts(b=bt_store,rho=rho_store,W=W) 
  }
}else{
  imp_obj <- list()
}

# Inefficiency factors
rho_ineff <- matrix(NA,T,1)
rho_sig_ineff <- nsave/effectiveSize(sig.rho_store)
b_ineff <-  omega_ineff <- matrix(NA,K,N)
sigma_ineff <- matrix(NA, N,1)

for (jj in 1:T){
  rho_ineff[jj] <- nsave/effectiveSize(rho_store[,jj])
}

for(ii in 1:N){
  sigma_ineff[ii] <- nsave/effectiveSize(Sigma_store[,1,ii])
  for(kk in 1:K){
    b_ineff[kk,ii] <- nsave/effectiveSize(b_store[,kk,ii])
    omega_ineff[kk,ii] <- nsave/effectiveSize(omega_store[,kk,ii])
  }
}
  
ineff_obj <- list("b" = b_ineff, "omega" = omega_ineff, "sig" = sigma_ineff, "rho" = rho_ineff, "rho_sig" = rho_sig_ineff)
save(est_obj,imp_obj, ineff_obj,file=paste0(paste0(sl.run,collapse="_"),".rda"))

# ------------------------------------------------------------------------------------------------------------------
get.post <- function(store,dims){
   mean <- apply(store,2:dims,median)
   lo <- apply(store,2:dims,quantile,0.16)
   hi <- apply(store,2:dims,quantile,0.84)
   lo1 <- apply(store,2:dims,quantile,0.05)
   hi1 <- apply(store,2:dims,quantile,0.95)
   return(list(mean,lo,hi,lo1,hi1))
 }

ht <- get.post(rho_store,dims=2)

pdf(paste0(paste0(sl.run,collapse="_"), "_rho.pdf"), width = 12, height =  8)
ts.plot(cbind(0,ht[[1]],ht[[2]],ht[[3]],ht[[4]],ht[[5]]),
        col=c("red",rep("black",5)),
        lwd=c(2,2,1,1,0.5,0.5),ylab="rho")
dev.off()

bt <- get.post(bt_store,dims=4)
pdf(paste0(paste0(sl.run,collapse="_"), "_bt.pdf"), width = 30, height =  20)
par(mfrow = c(5, 7))
for(ii in 1:dim(bt_store)[[4]]){
  ts.plot(cbind(0,bt[[1]][,,ii],bt[[2]][,,ii],bt[[3]][,,ii]),
          col=c("red",rep("black",6)),
          lwd=c(2,2,2,1,1,1,1), lty = c(rep(1,3), rep(2,4)), ylab= paste0("bt"), main = colnames(y)[[ii]])
}
dev.off()
