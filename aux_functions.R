get.impacts <- function(b,rho,W,quiet=FALSE){
  nsave <- dim(b)[1]
  K <- dim(b)[3]
  N <- dim(b)[4]
  T <- dim(rho)[2]  
  ind.names <- dimnames(b)[[4]]
  date.names <- dimnames(b)[[2]]
  
  imp_mat <- array(NA,dim=c(nsave,T,N,N))
  impdir <- matrix(NA,nsave,T)
  imptot <- matrix(NA,nsave,T)
  impind <- matrix(NA,nsave,T)
  imptotind <- impdirind <- impindind <- array(NA,dim=c(nsave,T,N))
  
  if(!quiet){
    pb <- txtProgressBar(min = 1, max = nsave, style = 3)
  }
  for(irep in 1:nsave){
    rho_d <- rho[irep,]
    b_d <- b[irep,,,]
    for(tt in 1:T){
      imp_mat[irep,tt,,] <- solve(diag(N)-rho_d[tt]*W) %*% diag(b_d[tt,1,])
      impdir[irep,tt] <- mean(diag(imp_mat[irep,tt,,]))
      imptot[irep,tt] <- sum(imp_mat[irep,tt,,])/N
      impind[irep,tt] <- imptot[irep,tt]-impdir[irep,tt]
      
      imptotind[irep,tt,] <- apply(imp_mat[irep,tt,,],c(1),sum)
      impdirind[irep,tt,] <- diag(imp_mat[irep,tt,,])
      impindind[irep,tt,] <- imptotind[irep,tt,]-impdirind[irep,tt,]
    }
    if(!quiet){
      setTxtProgressBar(pb, irep)
    }
  }
  
  imp_post <- apply(imp_mat,c(2,3,4),median)
  dimnames(imp_post) <- list(date.names,ind.names,ind.names)
  return(list("mat"=imp_post,"tot"=imptot,"dir"=impdir,"ind"=impind,
              "imptotind"=imptotind,"impdirind"=impdirind,"impindind"=impindind))
}

get.impacts2 <- function(b,rho,Wt,sl.W,quiet=FALSE){
  nsave <- dim(b)[1]
  K <- dim(b)[3]
  N <- dim(b)[4]
  T <- dim(rho)[2]  
  ind.names <- dimnames(b)[[4]]
  date.names <- dimnames(b)[[2]]
  
  imp_mat <- array(NA,dim=c(nsave,T,N,N))
  impdir <- matrix(NA,nsave,T)
  imptot <- matrix(NA,nsave,T)
  impind <- matrix(NA,nsave,T)
  imptotind <- impdirind <- impindind <- array(NA,dim=c(nsave,T,N))
  
  if(!quiet){
    pb <- txtProgressBar(min = 1, max = nsave, style = 3)
  }
  for(irep in 1:nsave){
    rho_d <- rho[irep,]
    b_d <- b[irep,,,]
    for(tt in 1:T){
      W <- Wt[sl.W[tt],,]
      imp_mat[irep,tt,,] <- solve(diag(N)-rho_d[tt]*W) %*% diag(b_d[tt,1,])
      impdir[irep,tt] <- mean(diag(imp_mat[irep,tt,,]))
      imptot[irep,tt] <- sum(imp_mat[irep,tt,,])/N
      impind[irep,tt] <- imptot[irep,tt]-impdir[irep,tt]
      
      imptotind[irep,tt,] <- apply(imp_mat[irep,tt,,],c(1),sum)
      impdirind[irep,tt,] <- diag(imp_mat[irep,tt,,])
      impindind[irep,tt,] <- imptotind[irep,tt,]-impdirind[irep,tt,]
    }
    if(!quiet){
      setTxtProgressBar(pb, irep)
    }
  }
  
  imp_post <- apply(imp_mat,c(2,3,4),median)
  dimnames(imp_post) <- list(date.names,ind.names,ind.names)
  return(list("mat"=imp_post,"tot"=imptot,"dir"=impdir,"ind"=impind,
              "imptotind"=imptotind,"impdirind"=impdirind,"impindind"=impindind))
}

# --------------------------------------------------------------------------------
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

get.facload <- function(yy,xx,l_sd){
  V_lambda <- try(solve(crossprod(xx)+diag(ncol(xx))*l_sd),silent = T)
  if(is(V_lambda,"try-error")){
    V_lambda <- ginv(crossprod(xx)+diag(ncol(xx))*l_sd)
  }
  
  lambda_mean <- V_lambda%*%(crossprod(xx,yy))
  lambda_draw <- try(lambda_mean+t(chol(V_lambda))%*%rnorm(ncol(xx)))
  if(is(lambda_draw,"try-error")){
    lambda_draw <- mvtnorm::rmvnorm(1,lambda_mean,V_lambda)
  }
  
  return(lambda_draw)
}

# -------------------------------------------------------
getmess.spatial <- function(y,X,W,beta,rho,Lft,Sig,pr){
  if(rho_tvp==TRUE){
    rhoWt <- list()
    for(tt in 1:T){
      rhoWt[[tt]] <- expm(rho_draw[tt]*W)
    }
    IoW <- as.matrix(bdiag(rhoWt))
  }else{
    IoW <- kronecker(diag(T),expm(rho*W))
  }
  Sy <- t(matrix(IoW %*% matrix(t(y),ncol=1),N,T))
  yhat <- matrix(0,T,N)
  for(i in 1:N){yhat[,i] <- cbind(ylag[,-i] %*% W[i,-i], ylag[,i], X[,i,]) %*% b_draw[,i]}
  
  yy <- (Sy - yhat - Lft) * exp(-0.5*Sig)
  logpost <- sum(dnorm(yy,mean=0,sd=1,log=T)) + dnorm(rho,mean=0,sd=pr,log=T)
  return(logpost)
}

get.spatial <- function(y,t,Xb,W,rho,S,par){
  ys <- matrix(t(y),ncol=1)
  IoW <- kronecker(diag(t),rho*W)
  
  Sy <- t(matrix(ys - IoW %*% ys,N,t))
  yhat <- matrix(0,t,N)
  yy <- (Sy - Xb)/sqrt(S)
  
  if(par=="rho"){
    logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      dunif(rho,min=0,max=1,log=T)
  }
  if(par=="kappa"){
    logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      log(1/w)
  }
  return(logpost)
}

get.spatial2 <- function(y,t,Xb,W,rho,S,par){
  # t=T;W=W;rho=rho_prop;par="rho"
  Sy <- y*0
  for(tt in 1:t){
    Sy[tt,] <- (diag(N)-rho*W)%*%as.numeric(y[tt,])
  }
  yy <- (Sy - Xb) * exp(-S/2)
  
  if(par=="rho"){
    # logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
    logpost <- t*log(det(diag(N)-rho*W)) + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      dunif(rho,min=0,max=1,log=T)
  }
  if(par=="kappa"){
    logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      log(1/w)
  }
  return(logpost)
}

get.spatial3 <- function(y,t,Xb,Wt,sl.W,rho,S,par){
  Sy <- y*0
  for(tt in 1:t){
    W <- Wt[sl.W[tt],,]
    Sy[tt,] <- (diag(N)-rho*W)%*%as.numeric(y[tt,])
  }
  yy <- (Sy - Xb) * exp(-S/2)
  
  if(par=="rho"){
    # logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
    logpost <- sum(sl.W==1)*log(det(diag(N)-rho*Wt[1,,])) + 
      sum(sl.W==2)*log(det(diag(N)-rho*Wt[2,,])) + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      dunif(rho,min=0,max=1,log=T)
  }
  if(par=="kappa"){
    logpost <- t*determinant(diag(N)-rho*W,logarithm = T)$modulus[[1]] + 
      sum(dnorm(yy,mean=0,sd=1,log=T)) + 
      log(1/w)
  }
  return(logpost)
}

# get normal gamma shrinkage prior
get.ng <- function(bdraw,psi,pr=0,d0=0.01,d1=0.01,th=0.1){
  k <- NROW(bdraw)
  
  # global shrinkage parameter
  d0_po <- d0 + th*k
  d1_po <- d1 + (th * sum(psi))/2
  lambda2 <- rgamma(1,d0_po,d1_po)
  
  # local shrinkage scalings
  psi <- matrix(0,k,1)
  for(kk in 1:k){
    scale <- ifelse(bdraw[kk]^2<=1e-8, 1e-8, bdraw[kk]^2) # for stability
    psi[kk] <- GIGrvg::rgig(n=1, lambda=th-0.5, chi=scale, psi=lambda2*th)
  }
  return(psi)
}
