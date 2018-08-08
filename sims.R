## Simulations relating to 
## "A comparison of truncated and time-weighted Plackett-Luce models for probabilistic forecasting of Formula One results"
## 
## daniel.henderson@ncl.ac.uk

source("functions.R")
library(coda)

## load data
load(file="F1.dat")

## 2010-2013
these.years <- seq(2010,2013,1)
Drivers <- NULL
driver.id <- NULL
teams <- NULL
n <- 0
n.i.max <- 0
for(i in 1:length(these.years)){ 
  n.races <- length(F1.res[[these.years[i]-1949]])
  ##n <- length(F1.res[[i]])
  for(j in 1:n.races){
    n.i.max <- max(n.i.max,dim(F1.res[[these.years[i]-1949]][[j]]$res)[1])
    Drivers <- unique(c(Drivers,as.character(F1.res[[these.years[i]-1949]][[j]]$res[,3])))
    driver.id <- unique(c(driver.id,F1.res[[these.years[i]-1949]][[j]]$res[,4]))
    teams <- unique(c(teams,as.character(F1.res[[these.years[i]-1949]][[j]]$res[,5])))
  }
  n <- n+n.races
}

## put data in format for R/rjags
K <- length(Drivers)
x <- matrix(NA,nrow=n,ncol=n.i.max)
dates <- NULL
ind.i <- 0
for(i in 1:length(these.years)){ 
  n.races <- length(F1.res[[these.years[i]-1949]])
  for(j in 1:n.races){
    ind.i <- ind.i+1
##    x[ind.i,1:dim(F1.res[[these.years[i]-1949]][[j]]$res)[1]] <- F1.res[[these.years[i]-1949]][[j]]$res[,4]
    x[ind.i,1:dim(F1.res[[these.years[i]-1949]][[j]]$res)[1]] <- match(F1.res[[these.years[i]-1949]][[j]]$res[,4],driver.id)
    dates[ind.i] <- F1.res[[these.years[i]-1949]][[j]]$date
  }
}
dates <- as.Date(dates,origin="1970-1-1") ## fix the dates

present.date <- as.Date("2013-11-25")

points.actual <- matrix(0,nrow=n,ncol=K)
for(tt in 1:n){
  these.drivers <- x[tt,is.na(x[tt,])==FALSE]
  points.actual[tt,x[tt,is.na(x[tt,])==FALSE]] <- points.calc.2013(1:length(these.drivers))
}


####################################################################
## Sequential analysis - started from draw from prior
####################################################################

## Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.PL <- PL.sequential(x=x,K=K,tau.vec=dates,a=a,b=b,its=its,burn=bur
n,thin=thin))
save(res.PL,file="res.PL") 

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9924
unix.time(res.PL.tw0.9924 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9924,file="res.PL.tw0.9924")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9990
unix.time(res.PL.tw0.9990 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9990,file="res.PL.tw0.9990")

## Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a,file="res.a") 

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9924
unix.time(res.a.tw0.9924 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9924,file="res.a.tw0.9924")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9981
unix.time(res.a.tw0.9981 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9981,file="res.a.tw0.9981")

## truncated Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(10,n)
unix.time(res.PL.t10 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=1,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.t10,file="res.PL.t10") 

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9981
unix.time(res.PL.tw0.9981 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9981,file="res.PL.tw0.9981")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9990
unix.time(res.a.tw0.9990 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9990,file="res.a.tw0.9990")

## truncated Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(14,n)
unix.time(res.PL.t14 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=1,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.t14,file="res.PL.t14") 

## truncated Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(6,n)
unix.time(res.PL.t6 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=1,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.t6,file="res.PL.t6") 

## Plackett-Luce model
a <- rep(10,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.PL.a10 <- PL.sequential(x=x,K=K,tau.vec=dates,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.a10,file="res.PL.a10") 

## Plackett-Luce model
a <- rep(0.1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.PL.a0.1 <- PL.sequential(x=x,K=K,tau.vec=dates,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.a0.1,file="res.PL.a0.1") 

## Plackett-Luce model
## prior hyperparameters from 2009 points (very crude!)
a <- c(26,22,49,84,34.5,1,95,69.5,1,77,17,5,1,1,22,6,32.5,1,1,24,1,3,1,1,1,19,1,1,1,1,1,1,1,48,1,1,1,1,1,1,1,1)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.PL.av <- PL.sequential(x=x,K=K,tau.vec=dates,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.av,file="res.PL.av") 


## Attrition model 
cc <- rep(10,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c10 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c10,file="res.a.c10") 

## Attrition model 
cc <- rep(0.1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c0.1 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c0.1,file="res.a.c0.1") 

## Attrition model 
## prior hyperparameters from 2009 points (very crude!)
cc <- 1/c(26,22,49,84,34.5,1,95,69.5,1,77,17,5,1,1,22,6,32.5,1,1,24,1,3,1,1,1,19,1,1,1,1,1,1,1,48,1,1,1,1,1,1,1,1)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.cv <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.cv,file="res.a.cv") 

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9970
unix.time(res.a.tw0.9970 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9970,file="res.a.tw0.9970")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9960
unix.time(res.a.tw0.9960 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9960,file="res.a.tw0.9960")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9950
unix.time(res.a.tw0.9950 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9950,file="res.a.tw0.9950")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9940
unix.time(res.a.tw0.9940 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9940,file="res.a.tw0.9940")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9772
unix.time(res.a.tw0.9772 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9772,file="res.a.tw0.9772")

## time-weighted Attrition model 
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9500
unix.time(res.a.tw0.9500 <- twattrition.sequential(x=x,K=K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.tw0.9500,file="res.a.tw0.9500")

## Attrition model 
cc <- rep(2,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c2 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c2,file="res.a.c2") 

## Attrition model 
cc <- rep(3,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c3 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c3,file="res.a.c3") 

## Attrition model 
cc <- rep(4,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c4 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c4,file="res.a.c4") 

## Attrition model 
cc <- rep(0.75,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c0.75 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c0.75,file="res.a.c0.75") 

## Attrition model 
cc <- rep(0.5,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c0.5 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c0.5,file="res.a.c0.5") 

## Attrition model 
cc <- rep(0.25,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
unix.time(res.a.c0.25 <- attrition.sequential(x=x,K=K,tau.vec=dates,cc=cc,d=d,its=its,burn=burn,thin=thin))
save(res.a.c0.25,file="res.a.c0.25") 

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9970
unix.time(res.PL.tw0.9970 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9970,file="res.PL.tw0.9970")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9960
unix.time(res.PL.tw0.9960 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9960,file="res.PL.tw0.9960")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9950
unix.time(res.PL.tw0.9950 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9950,file="res.PL.tw0.9950")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9940
unix.time(res.PL.tw0.9940 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9940,file="res.PL.tw0.9940")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9772
unix.time(res.PL.tw0.9772 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9772,file="res.PL.tw0.9772")

## time-weighted Plackett-Luce model
a <- rep(1,K)
b <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
r.vec <- rep(30,n)
xi <- 0.9500
unix.time(res.PL.tw0.9500 <- twtPL.sequential(x=x,K=K,r.vec=r.vec,tau.vec=dates,tau=present.date,xi=xi,a=a,b=b,its=its,burn=burn,thin=thin))
save(res.PL.tw0.9500,file="res.PL.tw0.9500")

## time-weighted attrition model
cc <- rep(1,K)
d <- 1
its <- 10000 
burn <- 100
thin <- 1   
set.seed(41)
xi <- 0.9970
res.a.tw.last <- twattrition.gibbs(x,K,tau.vec=dates,tau=present.date,xi=xi,cc=cc,d=d,gamma.init=rgamma(K,cc,d),its=its,burn=burn,thin=thin)
save(res.a.tw.last,file="res.a.tw.last")


q(save="no")
