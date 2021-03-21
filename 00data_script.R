# --------------------------------------------------------------------------------------------
# load new dataset
source("aux_MPdata.R")
if(sl.dataset=="d1992"){
  load("data/data1992.rda")
}else if(sl.dataset=="d1997"){
  load("data/data1997.rda")
}else if(sl.dataset=="d2002"){
  load("data/data2002.rda")
}else if(sl.dataset=="dTVW"){
  load("data/dataTVW.rda")
}

y <- mkret[date.names,]
sl.ind <- apply(is.na(y),2,sum)<=0
y <- y[,sl.ind]
y[is.na(y)] <- 0

T <- nrow(y)
N <- ncol(y)

ind.names <- colnames(y)
if(sl.dataset=="dTVW"){
  W1 <- t(Wmatr[1,ind.names,ind.names])
  W2 <- t(Wmatr[2,ind.names,ind.names])
  if(!diag.W){
    diag(W1) <- 0
    diag(W2) <- 0
  }
  
  # customer/supplier matrix
  W1 <- W1/apply(W1,1,sum)
  W2 <- W2/apply(W2,1,sum)
  
  # Normalize the spatial weighting matrix (supplier times customer)
  W1 <- t(W1)
  W2 <- t(W2)
  W1 <- W1/apply(W1,1,sum)
  W2 <- W2/apply(W2,1,sum)
  
  Wt <- array(NA,dim=c(2,length(ind.names),length(ind.names)))
  dimnames(Wt) <- list(NULL,ind.names,ind.names)
  Wt[1,,] <- W1
  Wt[2,,] <- W2
  
  sl.W <- (as.numeric(substr(date.names,1,4))>=2002)+1
}else{
  W <- t(Wmatr[ind.names,ind.names])
  if(!diag.W){
    diag(W) <- 0
  }
  # customer/supplier matrix
  W <- W/apply(W,1,sum)
  
  # Normalize the spatial weighting matrix (supplier times customer)
  W <- t(W)
  W <- W/apply(W,1,sum)
}

Xraw <- cbind(apply(shock[date.names,],1,mean,na.rm=T),1)
K <- ncol(Xraw)
if(ztrans){
  y <- apply(y,2,function(x){(x-mean(x))/sd(x)})
  Xraw[,2] <- (Xraw[,2] - mean(Xraw[,2]))/sd(Xraw[,2])
}

Xst <- as.matrix(kronecker(matrix(1,ncol=1,nrow=N),Xraw))
Xb <- y*0
sd.shock <- sd(Xraw[,2])

X <- array(NA,dim=c(T,N,K))
for(i in 1:N){
  X[,i,] <- Xraw
}
dimnames(X) <- list(date.names,ind.names,c("mpshock","cons"))
