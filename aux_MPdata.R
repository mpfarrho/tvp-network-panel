# check stuff from the original dataset
load("data/MPData.rda")

# construct the y's and add indicator for missings
sl.dates <- !Xraw$Date=="2001-01-03" # exclude emergency meeting
if(sample == "subT61_120"){
  sl.dates[1:61] <- FALSE   
}else if(sample == "subT1_60"){
  sl.dates[62:121] <- FALSE 
}

date.names <- Xraw$Date[sl.dates]
y <- as.matrix(Yraw.list$MktRet)
y <- y[sl.dates,]
T <- nrow(y)

sl.nona <- as.numeric(which(apply(is.na(y),2,sum)==0))
N.nona <- length(sl.nona)
y <- y[,sl.nona]
N <- ncol(y)
id.na <- is.na(y)

sd.shock <- sd(Xraw[,2])

# Remove dates and one weirdo observation
Xraw <- as.matrix(Xraw[sl.dates,-1])

# blow up X for all observations
X <- array(NA,dim=c(T,N.nona,2))
for(ii in 1:N.nona)  X[,ii,] <- Xraw

Xst <- as.matrix(kronecker(matrix(1,ncol=1,nrow=N),Xraw))
Xb <- y*0

# get weights matrix
W <- W1992[sl.nona,sl.nona]
diag(W) <- 0 # CHECK HERE FOR DIAGONAL OF W!
W <- W/apply(W,1,sum)
K <- dim(X)[3]

rm(W1992,Xraw,Yraw.list,sl.dates,y,T,sl.nona,N.nona,id.na,sd.shock,X,Xst,Xb,W,K)
