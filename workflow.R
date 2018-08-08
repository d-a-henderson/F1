## R workflow for analysing simulations and producing plots/summaries for
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
wins.actual <- matrix(0,nrow=n,ncol=K)
podiums.actual <- matrix(0,nrow=n,ncol=K)
top10s.actual <- matrix(0,nrow=n,ncol=K)
for(tt in 1:n){
  these.drivers <- x[tt,is.na(x[tt,])==FALSE]
  points.actual[tt,x[tt,is.na(x[tt,])==FALSE]] <- points.calc.2013(1:length(these.drivers))
  wins.actual[tt,x[tt,1]] <- 1
  podiums.actual[tt,x[tt,1:3]] <- 1
  top10s.actual[tt,x[tt,1:10]] <- 1
}

##names(which(sapply(iconvlist(),iconv, x=Drivers[34])=="Kimi Räikkönen"))
##iconv(Drivers[34],to="UTF8")
kimi <- Drivers[34]
##Drivers[34] <- "Kimi Räikkönen"
Drivers[34] <- "Kimi Raikkonen"


####################################################################
##read in simulation results 
####################################################################

load(file="res.PL") 
load(file="res.PL.t10") 
load(file="res.PL.t14") 
load(file="res.PL.t6") 
load(file="res.PL.tw0.9924")
load(file="res.PL.tw0.9981")
load(file="res.PL.tw0.9990")
load(file="res.a") 
load(file="res.a.tw0.9924")
load(file="res.a.tw0.9981")
load(file="res.a.tw0.9990")
load(file="res.PL.a10") 
load(file="res.PL.a0.1") 
load(file="res.PL.av") 
load(file="res.a.c0.1") 
load(file="res.a.c0.25") 
load(file="res.a.c0.5") 
load(file="res.a.c0.75") 
load(file="res.a.c2") 
load(file="res.a.c3") 
load(file="res.a.c4") 
load(file="res.a.c10") 
load(file="res.a.cv") 
load(file="res.a.tw0.9500")
load(file="res.a.tw0.9772")
load(file="res.a.tw0.9940")
load(file="res.a.tw0.9950")
load(file="res.a.tw0.9960")
load(file="res.a.tw0.9970")
load(file="res.PL.tw0.9500")
load(file="res.PL.tw0.9772")
load(file="res.PL.tw0.9940")
load(file="res.PL.tw0.9950")
load(file="res.PL.tw0.9960")
load(file="res.PL.tw0.9970")

##############################################################
## prior sensitivity
models <- c("a","a.c2","a.c3","a.c4","a.c0.75")
n.models <- length(models)
model.names <- c("a","a.c2","a.c3","a.c4","a.c0.75")

##############################################################
## Winner, podium and points - log scoring rule

lsws <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsws[,i] <- eval(parse(text=paste("res.",models[i],"$lsws",sep="")))
}

lsps <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsps[,i] <- eval(parse(text=paste("res.",models[i],"$lsps",sep="")))
}

lst10s <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lst10s[,i] <- eval(parse(text=paste("res.",models[i],"$lst10s",sep="")))
}

##############################################################
## sensitivity to prior
postscript("prob-vettel-button-sens.eps",pointsize=18)
par(mfrow=c(2,3))
plot(dates,res.a$pp[,4,1],col=1,type="l",ylim=c(0,0.6),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.a.c2$pp[,4,1],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,4,1],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,4,1],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,4,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),0.6,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.a.c2$pp[,4,2],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,4,2],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,4,2],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,4,2],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.a.c2$pp[,4,3],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,4,3],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,4,3],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,4,3],col=5,lty=5,lwd=2)
legend(as.Date("2012-10-1"),0.4,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,7,1],col=1,type="l",ylim=c(0,0.6),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.a.c2$pp[,7,1],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,7,1],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,7,1],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,7,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),0.6,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.a.c2$pp[,7,2],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,7,2],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,7,2],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,7,2],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.a.c2$pp[,7,3],col=2,lty=2,lwd=2) ## not much difference
lines(dates,res.a.c3$pp[,7,3],col=3,lty=3,lwd=2) ## not much difference
lines(dates,res.a.c4$pp[,7,3],col=4,lty=4,lwd=2) ## not much difference
lines(dates,res.a.c0.75$pp[,7,3],col=5,lty=5,lwd=2)
legend(as.Date("2012-10-1"),0.4,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
par(mfrow=c(1,1))
dev.off()


## head-to-head
postscript("prior-sens.eps",pointsize=18)
par(mfrow=c(1,2))
plot(seq(0,10,0.01),dgamma(seq(0,10,0.01),1,1),type="l",lwd=2,xlab="lambda",ylab="Density",ylim=c(0,2))
lines(seq(0,10,0.01),dgamma(seq(0,10,0.01),2,1),lwd=2,col=2,lty=2)
lines(seq(0,10,0.01),dgamma(seq(0,10,0.01),3,1),lwd=2,col=3,lty=3)
lines(seq(0,10,0.01),dgamma(seq(0,10,0.01),4,1),lwd=2,col=4,lty=4)
lines(seq(0,10,0.01),dgamma(seq(0,10,0.01),0.75,1),lwd=2,col=5,lty=5)
legend(6,2,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
##dev.off()
##postscript("head-to-head-sens.eps",pointsize=18)
##par(mfrow=c(1,1))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),1,1),type="l",lwd=2,xlab="Prior probability of head-to-head win",ylab="Density",ylim=c(0,2.5))
lines(seq(0,1,0.01),dbeta(seq(0,1,0.01),2,2),lwd=2,col=2,lty=2)
lines(seq(0,1,0.01),dbeta(seq(0,1,0.01),3,3),lwd=2,col=3,lty=3)
lines(seq(0,1,0.01),dbeta(seq(0,1,0.01),4,4),lwd=2,col=4,lty=4)
lines(seq(0,1,0.01),dbeta(seq(0,1,0.01),0.75,0.75),lwd=2,col=5,lty=5)
legend(0.65,2.5,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
dev.off()


postscript("racewin-sens.eps",pointsize=18)
par(mfrow=c(2,2))
plot(res.a$pp[58,,1],pch=19,xlab="Driver ID",ylab="Probability of winning",main="2012 Brazilian Grand Prix (25/11/12)",ylim=c(0,0.4))
points(res.a.c0.75$pp[58,,1])
##points(res.a.c2$pp[59,,1],pch=2)
plot(res.a$pp[59,,1],pch=19,xlab="Driver ID",ylab="Probability of winning",main="2013 Australian Grand Prix (17/3/13)",ylim=c(0,0.4))
points(res.a.c0.75$pp[59,,1])
##points(res.a.c2$pp[59,,1],pch=2)
plot(res.a$pp[60,,1],pch=19,xlab="Driver ID",ylab="Probability of winning",main="2013 Malaysian Grand Prix (24/3/13)",ylim=c(0,0.4))
points(res.a.c0.75$pp[60,,1])
par(mfrow=c(1,1))
dev.off()

####################################################################
## log marginal likelihoods

lml.a <- sum(res.a$lpp)
lml.a.c10 <- sum(res.a.c10$lpp)
lml.a.c0.1 <- sum(res.a.c0.1$lpp)
lml.a.cv <- sum(res.a.cv$lpp)
lml.a.c2 <- sum(res.a.c2$lpp)
lml.a.c3 <- sum(res.a.c3$lpp)
lml.a.c4 <- sum(res.a.c4$lpp)
lml.a.c0.25 <- sum(res.a.c0.25$lpp)
lml.a.c0.5 <- sum(res.a.c0.5$lpp)
lml.a.c0.75 <- sum(res.a.c0.75$lpp)

lml <- c(lml.a,lml.a.c2,lml.a.c3,lml.a.c4,lml.a.c10,lml.a.c0.1,lml.a.c0.25,lml.a.c0.5,lml.a.c0.75)
a.seq <- c(1,2,3,4,10,0.1,0.25,0.5,0.75)
models <-  c("a","a.c2","a.c3","a.c4","a.c10","a.c0.1","a.c0.25","a.c0.5","a.c0.75")

n.models <- length(lml)
lpp <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lpp[,i] <- eval(parse(text=paste("res.",models[i],"$lpp",sep="")))
}

## posterior model probability - asssume equal prior
post.mod.prob <- lpp
for(i in 1:n.models){
  sum.ml.rat <- 0
  for(j in 1:n.models){
    sum.ml.rat <- sum.ml.rat+exp(cumsum(lpp[,j])-cumsum(lpp[,i]))
  }
  post.mod.prob[,i] <- 1/sum.ml.rat
}

postscript("lmlprior-sens.eps",pointsize=18)
par(mfrow=c(1,1))
plot(sort(a.seq),lml[order(a.seq)],xlab="a",ylab="Log prior predictive",type="b")
dev.off()

plot(a.seq,lml,xlab="a",ylab="Log prior predictive")
plot(a.seq,exp(lml-max(lml))/sum(exp(lml-max(lml))),xlab="a",ylab="Posterior probability",type="h")

####################################################################
## sensitivity to time-weighting parameter

## focus on these models
models <- c("a","a.tw0.9990","a.tw0.9981", "a.tw0.9970", "a.tw0.9960","a.tw0.9950","a.tw0.9940","a.tw0.9924","a.tw0.9772","a.tw0.9500")
n.models <- length(models)
model.names <- c("a","a.tw9990","a.tw9981", "a.tw9970", "a.tw9960","a.tw9950","a.tw9940","a.tw9924","a.tw9772","a.tw9500")

##############################################################
## Winner, podium and points - log scoring rule

lsws <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsws[,i] <- eval(parse(text=paste("res.",models[i],"$lsws",sep="")))
}

lsps <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsps[,i] <- eval(parse(text=paste("res.",models[i],"$lsps",sep="")))
}

lst10s <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lst10s[,i] <- eval(parse(text=paste("res.",models[i],"$lst10s",sep="")))
}

##############################################################
## sensitivity to time-weighting
postscript("prob-vettel-button-tw.eps",pointsize=18)
par(mfrow=c(2,3))
plot(dates,res.a$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.a.tw0.9981$pp[,4,1],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,4,1],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,4,1],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,4,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.a$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.a.tw0.9981$pp[,4,2],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,4,2],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,4,2],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,4,2],col=5,lty=5,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.a$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.a.tw0.9981$pp[,4,3],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,4,3],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,4,3],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,4,3],col=5,lty=5,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.a$pp[,7,1],col=1,type="l",ylim=c(0,0.6),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.a.tw0.9981$pp[,7,1],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,7,1],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,7,1],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,7,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),0.6,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.a$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.a.tw0.9981$pp[,7,2],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,7,2],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,7,2],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,7,2],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.a$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.a.tw0.9981$pp[,7,3],col=2,lty=2,lwd=2) 
lines(dates,res.a.tw0.9960$pp[,7,3],col=3,lty=3,lwd=2) 
lines(dates,res.a.tw0.9940$pp[,7,3],col=4,lty=4,lwd=2) 
lines(dates,res.a.tw0.9924$pp[,7,3],col=5,lty=5,lwd=2)
legend(as.Date("2012-6-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
par(mfrow=c(1,1))
dev.off()

## collect model probabilities together in one big array
xi.seq <- c(1,0.9990,0.9981,0.9970,0.9960,0.9950,0.9940,0.9924,0.9772,0.9500)
modelprobs <- array(0,c(10,77,42,3))
modelprobs[1,,,] <- res.a$pp
modelprobs[2,,,] <- res.a.tw0.9990$pp
modelprobs[3,,,] <- res.a.tw0.9981$pp
modelprobs[4,,,] <- res.a.tw0.9970$pp
modelprobs[5,,,] <- res.a.tw0.9960$pp
modelprobs[6,,,] <- res.a.tw0.9950$pp
modelprobs[7,,,] <- res.a.tw0.9940$pp
modelprobs[8,,,] <- res.a.tw0.9924$pp
modelprobs[9,,,] <- res.a.tw0.9772$pp
modelprobs[10,,,] <- res.a.tw0.9500$pp
probs.opt <- array(0,dim(res.a$pp))

## winning probs ############################################
## cumulative score
lsws.cum <- apply(lsws,2,cumsum)
## which model is the best so far
this.model <- apply(lsws.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lsws.opt <- lsws[1,this.model[1]]
for(i in 2:length(this.model)){
  lsws.opt <- c(lsws.opt,lsws[i,this.model[i-1]])
}
ts.plot(cumsum(lsws.opt)-lsws.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lsws.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lsws.cum[77,])]

plot(dates,cumsum(lsws.opt)-lsws.cum[,1],type="l") ## comp to non-tw

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,1] <- modelprobs[this.model[i],i,j,1]
  }
}


## podium probs ############################################
## cumulative score
lsps.cum <- apply(lsps,2,cumsum)
## which model is the best so far
this.model <- apply(lsps.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lsps.opt <- lsps[1,this.model[1]]
for(i in 2:length(this.model)){
  lsps.opt <- c(lsps.opt,lsps[i,this.model[i-1]])
}
ts.plot(cumsum(lsps.opt)-lsps.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lsps.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lsps.cum[77,])]

plot(dates,cumsum(lsps.opt)-lsps.cum[,1],type="l") ## comp to non-tw

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,2] <- modelprobs[this.model[i],i,j,2]
  }
}


## points probs ############################################
## cumulative score
lst10s.cum <- apply(lst10s,2,cumsum)
## which model is the best so far
this.model <- apply(lst10s.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lst10s.opt <- lst10s[1,this.model[1]]
for(i in 2:length(this.model)){
  lst10s.opt <- c(lst10s.opt,lst10s[i,this.model[i-1]])
}
ts.plot(cumsum(lst10s.opt)-lst10s.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lst10s.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lst10s.cum[77,])]

plot(dates,cumsum(lst10s.opt)-lst10s.cum[,1],type="l") ## comp to non-tw

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,3] <- modelprobs[this.model[i],i,j,3]
  }
}

par(mfrow=c(2,3))
plot(dates,res.a$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,1],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,2],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,3],col=2,lty=2,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,1],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,1],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,2],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,3],col=2,lty=2,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
par(mfrow=c(1,1))

## log marginal likelihood
lml.a <- sum(res.a$lpp)
lml.a.tw0.9924 <- sum(res.a.tw0.9924$lpp)
lml.a.tw0.9981 <- sum(res.a.tw0.9981$lpp)
lml.a.tw0.9990 <- sum(res.a.tw0.9990$lpp)
lml.a.tw0.9970 <- sum(res.a.tw0.9970$lpp)
lml.a.tw0.9960 <- sum(res.a.tw0.9960$lpp)
lml.a.tw0.9950 <- sum(res.a.tw0.9950$lpp)
lml.a.tw0.9940 <- sum(res.a.tw0.9940$lpp)
lml.a.tw0.9772 <- sum(res.a.tw0.9772$lpp)
lml.a.tw0.9500 <- sum(res.a.tw0.9500$lpp)
lml <- c(lml.a,lml.a.tw0.9990,lml.a.tw0.9981,lml.a.tw0.9970,lml.a.tw0.9960,lml.a.tw0.9950,lml.a.tw0.9940,lml.a.tw0.9924,lml.a.tw0.9772,lml.a.tw0.9500)
plot(sort(xi.seq),lml[order(xi.seq)],type="b",xlab="xi",ylab="Log prior predictive")

## extract individual log prior predictives
n.models <- length(lml)
lpp <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lpp[,i] <- eval(parse(text=paste("res.",models[i],"$lpp",sep="")))
}

## cumulative score
lpp.cum <- apply(lpp,2,cumsum)
## which model is the best so far
this.model <- apply(lpp.cum,1,which.max)
this.model.lpp <- this.model
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lpp.opt <- lpp[1,this.model[1]]
for(i in 2:length(this.model)){
  lpp.opt <- c(lpp.opt,lpp[i,this.model[i-1]])
}
ts.plot(cumsum(lpp.opt)-lpp.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lpp.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lpp.cum[77,])]

plot(dates,cumsum(lpp.opt)-lpp.cum[,1],type="l") ## comp to non-tw

## predictive probabilities based on optimal lml
probs.lml <-  probs.opt
for(i in 2:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.lml[i,j,1] <- modelprobs[this.model[i-1],i,j,1]
    probs.lml[i,j,2] <- modelprobs[this.model[i-1],i,j,2]
    probs.lml[i,j,3] <- modelprobs[this.model[i-1],i,j,3]
  }
}

## log score for optimal lml
lsws.lml <- rep(0,77)
lsps.lml <- rep(0,77)
lst10s.lml <- rep(0,77)
for(i in 1:77){
  lsws.lml[i] <- sum(log(probs.lml[i,x[i,1],1]),log(1-probs.lml[i,x[i,-1][is.na(x[i,-1])==FALSE],1]))
  lsps.lml[i] <- sum(log(probs.lml[i,x[i,1:3],2]),log(1-probs.lml[i,x[i,-(1:3)][is.na(x[i,-(1:3)])==FALSE],2]))
  lst10s.lml[i] <- sum(log(probs.lml[i,x[i,1:10],3]),log(1-probs.lml[i,x[i,-(1:10)][is.na(x[i,-(1:10)])==FALSE],3]))
}

ts.plot(cumsum(lsws.lml)-cumsum(lsws[,1]))
ts.plot(cumsum(lsps.lml)-cumsum(lsps[,1]))
ts.plot(cumsum(lst10s.lml)-cumsum(lst10s[,1]))

##########################################################
## model averaged xi (based on lml)

## posterior model probability - assume equal prior
post.mod.prob <- lpp
for(i in 1:n.models){
  sum.ml.rat <- 0
  for(j in 1:n.models){
    sum.ml.rat <- sum.ml.rat+exp(cumsum(lpp[,j])-cumsum(lpp[,i]))
  }
  post.mod.prob[,i] <- 1/sum.ml.rat
}
ts.plot(post.mod.prob,col=1:10,ylim=c(0,1))


## average xi by posterior probs

plot(dates,post.mod.prob%*%matrix(xi.seq),lwd=2,type="l",ylab="Posterior mean xi")
plot(dates,log(0.5)/log(post.mod.prob%*%matrix(xi.seq)),lwd=2,type="l",ylab="Posterior mean half-life (days)")

## model averaged log prior predictives
lpp.modavg <- lpp[1,1] ## as all the same
for(i in 2:dim(lpp)[1]){  
  lpp.modavg[i] <- log(sum(post.mod.prob[i-1,]*exp(lpp.cum[i,]-max(lpp.cum[i,]))))+max(lpp.cum[i,])
}

## close to optimal if we chose a single value
lpp.modavg[77]-lml

## model averaged predictive probabilities
probs.modavg <-  probs.opt
for(i in 2:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.modavg[i,j,1] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,1])
    probs.modavg[i,j,2] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,2])
    probs.modavg[i,j,3] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,3])
  }
}

par(mfrow=c(2,3))
plot(dates,res.a$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,1],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,1],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,2],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,2],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,3],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,3],col=3,lty=3,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,1],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,1],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,1],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,2],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,2],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.a$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,3],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,3],col=3,lty=3,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw","a.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
par(mfrow=c(1,1))

## log score for model averaging
lsws.modavg <- rep(0,77)
lsps.modavg <- rep(0,77)
lst10s.modavg <- rep(0,77)
for(i in 1:77){
  lsws.modavg[i] <- sum(log(probs.modavg[i,x[i,1],1]),log(1-probs.modavg[i,x[i,-1][is.na(x[i,-1])==FALSE],1]))
  lsps.modavg[i] <- sum(log(probs.modavg[i,x[i,1:3],2]),log(1-probs.modavg[i,x[i,-(1:3)][is.na(x[i,-(1:3)])==FALSE],2]))
  lst10s.modavg[i] <- sum(log(probs.modavg[i,x[i,1:10],3]),log(1-probs.modavg[i,x[i,-(1:10)][is.na(x[i,-(1:10)])==FALSE],3]))
}

## comparison of model averaged and optimal based on lml
postscript("../lscore-comp-xi.eps",pointsize=16)
par(mfrow=c(2,2))
plot(dates,cumsum(lsws.modavg)-cumsum(lsws[,1]),ylab="Log score relative to attrition",ylim=c(-5,25),type="l",lwd=2,main="Winner",xlab="Date")
lines(dates,cumsum(lsws.lml)-cumsum(lsws[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lsws.opt)-cumsum(lsws[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2010-3-1"),25,legend=c("a.tw.ma","a.tw.lpp","a.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,cumsum(lsps.modavg)-cumsum(lsps[,1]),ylab="Log score relative to attrition",ylim=c(-5,25),type="l",lwd=2,main="Top 3",xlab="Date")
lines(dates,cumsum(lsps.lml)-cumsum(lsps[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lsps.opt)-cumsum(lsps[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2010-3-1"),25,legend=c("a.tw.ma","a.tw.lpp","a.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,cumsum(lst10s.modavg)-cumsum(lst10s[,1]),ylab="Log score relative to attrition",ylim=c(-5,25),type="l",lwd=2,main="Top 10",xlab="Date")
lines(dates,cumsum(lst10s.lml)-cumsum(lst10s[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lst10s.opt)-cumsum(lst10s[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2010-3-1"),25,legend=c("a.tw.ma","a.tw.lpp","a.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,lpp.modavg-cumsum(lpp[,1]),ylab="Log score relative to attrition",ylim=c(-5,25),type="l",lwd=2,main="Full finishing order",xlab="Date")
lines(dates,cumsum(lpp.opt)-cumsum(lpp[,1]),col=2,lty=2,lwd=2)
abline(h=0,col=8)
legend(as.Date("2010-3-1"),25,legend=c("a.tw.ma","a.tw.lpp=a.tw.opt"),col=c(1,2),lty=c(1,2),lwd=rep(2,2),cex=0.8)
par(mfrow=c(1,1))
dev.off()

## P(champion) for each driver over time for each model

champprobs <- array(0,c(10,77,42))
for(i in 1:length(Drivers)){
  champprobs[1,,i] <- apply(res.a$cw==i,1,mean)
  champprobs[2,,i] <- apply(res.a.tw0.9990$cw==i,1,mean)
  champprobs[3,,i] <- apply(res.a.tw0.9981$cw==i,1,mean)
  champprobs[4,,i] <- apply(res.a.tw0.9970$cw==i,1,mean)
  champprobs[5,,i] <- apply(res.a.tw0.9960$cw==i,1,mean)
  champprobs[6,,i] <- apply(res.a.tw0.9950$cw==i,1,mean)
  champprobs[7,,i] <- apply(res.a.tw0.9940$cw==i,1,mean)
  champprobs[8,,i] <- apply(res.a.tw0.9924$cw==i,1,mean)
  champprobs[9,,i] <- apply(res.a.tw0.9772$cw==i,1,mean)
  champprobs[10,,i] <- apply(res.a.tw0.9500$cw==i,1,mean)
}

## model averaged and optimal based on lml
champprobs.optlml <- champprobs[1,,]
champprobs.modavg <- champprobs[1,,]
for(i in 2:dim(champprobs.optlml)[1]){
  champprobs.optlml[i,] <- champprobs[this.model.lpp[i-1],i,]
  champprobs.modavg[i,] <- t(post.mod.prob[i-1,])%*%champprobs[,i,]
}

par(mfrow=c(2,1))
par(mar=c(4, 4, 1, 2)+0.1)
plot(1:n,champprobs[1,,4],type="l",xlab="Race",ylab="P(Champion)",ylim=c(0,1.1),main="",lab=c(20,5,7),lwd=2)
lines(1:n,champprobs.optlml[,4],type="l",col=2,lty=2,lwd=2)
lines(1:n,champprobs.modavg[,4],type="l",col=3,lty=3,lwd=2)
lines(1:n,apply(res.PL.t6$cw==4,1,mean),type="l",col=5,lty=5,lwd=2)
text(1:n,rep(0,n),labels=as.character(apply(x==4,1,which.max)),cex=0.4)
abline(v=length(grep("2010",dates))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates)))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates),grep("2012",dates)))+0.5,col="grey",lty=2)
text(9,1.05,"2010")
text(28,1.05,"2011")
text(48,1.05,"2012")
text(69,1.05,"2013")
legend(65,0.55,legend=model.names[c(1,2,3,5)],lty=c(1,2,3,5),col=c(1,2,3,5),lwd=rep(2,3),cex=0.9)
plot(1:n,cumsum(apply(res.a$cw==4,1,mean))-cumsum(apply(res.a$cw==4,1,mean)),type="l",xlab="Race",ylab="Difference from attrition",ylim=c(-4,1.2),main="",lab=c(20,5,7),lwd=2)
lines(1:n,cumsum(champprobs.optlml[,4])-cumsum(apply(res.a$cw==4,1,mean)),type="l",col=2,lty=2,lwd=2)
lines(1:n,cumsum(champprobs.modavg[,4])-cumsum(apply(res.a$cw==4,1,mean)),type="l",col=3,lty=2,lwd=2)
lines(1:n,cumsum(apply(res.PL.t6$cw==4,1,mean))-cumsum(apply(res.a$cw==4,1,mean)),type="l",col=5,lty=5,lwd=2)
abline(v=length(grep("2010",dates))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates)))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates),grep("2012",dates)))+0.5,col="grey",lty=2)
text(9,1.1,"2010")
text(28,1.1,"2011")
text(48,1.1,"2012")
text(69,1.1,"2013")
legend(5,-1.5,legend=model.names[c(1,2,5)],lty=c(1,2,5),col=c(1,2,5),lwd=rep(2,3),cex=0.9)
text(1:n,rep(-4,n),labels=as.character(apply(x==4,1,which.max)),cex=0.4)
par(mar=c(1, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))

## compare to Bookmakers

odds <- read.csv("sky_outright_new.csv",header=FALSE)
## convert to decimal odds
odds[,3] <- odds[,3]+1
Vettel <- odds[odds[,4]=="Sebastian Vettel",]
Alonso <- odds[odds[,4]=="Fernando Alonso",]
Raikkonen <- odds[odds[,4]=="Kimi Raikkonen",]
Hamilton <- odds[odds[,4]=="Kimi Raikkonen",]

par(mfrow=c(2,2))
par(mar=c(4, 4, 1, 2)+0.1)
plot(dates[grep("2013",dates)],apply(res.a$cw==4,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,1),main=Drivers[4],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Vettel[,1]),1/Vettel[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==4,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
lines(dates[grep("2013",dates)],champprobs.optlml[,4][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.25,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(dates[grep("2013",dates)],apply(res.a$cw==1,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,0.6),main=Drivers[1],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Alonso[,1]),1/Alonso[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==1,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
lines(dates[grep("2013",dates)],champprobs.optlml[,1][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.6,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(dates[grep("2013",dates)],apply(res.a$cw==34,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,0.6),main=Drivers[34],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Raikkonen[,1]),1/Raikkonen[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==34,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
lines(dates[grep("2013",dates)],champprobs.optlml[,34][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.6,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),4]),type="s",lwd=2,xlab="",ylab="Drivers' Championship points")
lines(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),1]),col=2,lwd=2,lty=2,type="s")
lines(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),34]),col=4,lwd=2,lty=4,type="s")
legend(as.Date("2013-4-1"),350,legend=Drivers[c(4,1,34)],lty=c(1,2,4),col=c(1,2,4),lwd=c(2,2,2),cex=0.6)
par(mar=c(5, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))

## save important quantities for later
lpp.opt.a <- lpp.opt
lsws.lml.a <- lsws.lml
lsps.lml.a <- lsps.lml
lst10s.lml.a <- lst10s.lml
probs.lml.a <- probs.lml
this.model.lpp.a <- this.model.lpp
champprobs.optlml.a <- champprobs.optlml

###########################################################################
## sensitivity for Plackett-Luce

## focus on these models
models <- c("PL","PL.tw0.9990","PL.tw0.9981", "PL.tw0.9970", "PL.tw0.9960","PL.tw0.9950","PL.tw0.9940","PL.tw0.9924","PL.tw0.9772","PL.tw0.9500")
n.models <- length(models)
model.names <- c("PL","PL.tw9990","PL.tw9981", "PL.tw9970", "PL.tw9960","PL.tw9950","PL.tw9940","PL.tw9924","PL.tw9772","PL.tw9500")

##############################################################
## Winner, podium and points - log scoring rule

lsws <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsws[,i] <- eval(parse(text=paste("res.",models[i],"$lsws",sep="")))
}

lsps <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lsps[,i] <- eval(parse(text=paste("res.",models[i],"$lsps",sep="")))
}

lst10s <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lst10s[,i] <- eval(parse(text=paste("res.",models[i],"$lst10s",sep="")))
}

##############################################################
## sensitivity to time-weighting
par(mfrow=c(2,3))
plot(dates,res.PL$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,4,1],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,4,1],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,4,1],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,4,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.PL$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,4,2],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,4,2],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,4,2],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,4,2],col=5,lty=5,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.PL$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,4,3],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,4,3],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,4,3],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,4,3],col=5,lty=5,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.PL$pp[,7,1],col=1,type="l",ylim=c(0,0.6),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,7,1],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,7,1],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,7,1],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,7,1],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),0.6,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.PL$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,7,2],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,7,2],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,7,2],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,7,2],col=5,lty=5,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
plot(dates,res.PL$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,res.PL.tw0.9981$pp[,7,3],col=2,lty=2,lwd=2) 
lines(dates,res.PL.tw0.9960$pp[,7,3],col=3,lty=3,lwd=2) 
lines(dates,res.PL.tw0.9940$pp[,7,3],col=4,lty=4,lwd=2) 
lines(dates,res.PL.tw0.9924$pp[,7,3],col=5,lty=5,lwd=2)
legend(as.Date("2012-6-1"),0.4,legend=model.names[c(1,3,5,7,8)],col=1:5,lty=1:5,lwd=rep(2,5),cex=0.6)
par(mfrow=c(1,1))

## collect model probabilities together in one big array
xi.seq <- c(1,0.9990,0.9981,0.9970,0.9960,0.9950,0.9940,0.9924,0.9772,0.9500)
modelprobs <- array(0,c(10,77,42,3))
modelprobs[1,,,] <- res.PL$pp
modelprobs[2,,,] <- res.PL.tw0.9990$pp
modelprobs[3,,,] <- res.PL.tw0.9981$pp
modelprobs[4,,,] <- res.PL.tw0.9970$pp
modelprobs[5,,,] <- res.PL.tw0.9960$pp
modelprobs[6,,,] <- res.PL.tw0.9950$pp
modelprobs[7,,,] <- res.PL.tw0.9940$pp
modelprobs[8,,,] <- res.PL.tw0.9924$pp
modelprobs[9,,,] <- res.PL.tw0.9772$pp
modelprobs[10,,,] <- res.PL.tw0.9500$pp
probs.opt <- array(0,dim(res.PL$pp))

## winning probs ############################################
## cumulative score
lsws.cum <- apply(lsws,2,cumsum)
## which model is the best so far
this.model <- apply(lsws.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lsws.opt <- lsws[1,this.model[1]]
for(i in 2:length(this.model)){
  lsws.opt <- c(lsws.opt,lsws[i,this.model[i-1]])
}
ts.plot(cumsum(lsws.opt)-lsws.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lsws.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lsws.cum[77,])]

plot(dates,cumsum(lsws.opt)-lsws.cum[,1],type="l") ## comp to non-tw
##plot(dates,cumsum(lsws.opt)-lsws.cum[,which.max(lsws.cum[77,])],type="l")
##plot(dates,cumsum(lsws.opt)-lsws.cum[,3],type="l")

## adaptive sequential choice is poor compared to a single choice

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,1] <- modelprobs[this.model[i],i,j,1]
  }
}


## podium probs ############################################
## cumulative score
lsps.cum <- apply(lsps,2,cumsum)
## which model is the best so far
this.model <- apply(lsps.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lsps.opt <- lsps[1,this.model[1]]
for(i in 2:length(this.model)){
  lsps.opt <- c(lsps.opt,lsps[i,this.model[i-1]])
}
ts.plot(cumsum(lsps.opt)-lsps.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lsps.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lsps.cum[77,])]

plot(dates,cumsum(lsps.opt)-lsps.cum[,1],type="l") ## comp to non-tw
##plot(dates,cumsum(lsps.opt)-lsps.cum[,which.max(lsps.cum[77,])],type="l")
##plot(dates,cumsum(lsps.opt)-lsps.cum[,3],type="l")

## adaptive sequential choice is poor compared to a single choice

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,2] <- modelprobs[this.model[i],i,j,2]
  }
}


## points probs ############################################
## cumulative score
lst10s.cum <- apply(lst10s,2,cumsum)
## which model is the best so far
this.model <- apply(lst10s.cum,1,which.max)
models[this.model]
plot(dates,xi.seq[this.model],xlab="xi")

## work out the log score for this sequence of models
lst10s.opt <- lst10s[1,this.model[1]]
for(i in 2:length(this.model)){
  lst10s.opt <- c(lst10s.opt,lst10s[i,this.model[i-1]])
}
ts.plot(cumsum(lst10s.opt)-lst10s.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lst10s.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lst10s.cum[77,])]

plot(dates,cumsum(lst10s.opt)-lst10s.cum[,1],type="l") ## comp to non-tw
##plot(dates,cumsum(lst10s.opt)-lst10s.cum[,which.max(lst10s.cum[77,])],type="l")
##plot(dates,cumsum(lst10s.opt)-lst10s.cum[,3],type="l")

## adaptive sequential choice is poor compared to a single choice

for(i in 1:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.opt[i,j,3] <- modelprobs[this.model[i],i,j,3]
  }
}

par(mfrow=c(2,3))
plot(dates,res.PL$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,1],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,2],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,3],col=2,lty=2,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,1],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,1],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,2],col=2,lty=2,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,3],col=2,lty=2,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("a","a.tw"),col=1:2,lty=1:2,lwd=rep(2,2),cex=0.6)
par(mfrow=c(1,1))

## log marginal likelihood
lml.PL <- sum(res.PL$lpp)
lml.PL.tw0.9924 <- sum(res.PL.tw0.9924$lpp)
lml.PL.tw0.9981 <- sum(res.PL.tw0.9981$lpp)
lml.PL.tw0.9990 <- sum(res.PL.tw0.9990$lpp)
lml.PL.tw0.9970 <- sum(res.PL.tw0.9970$lpp)
lml.PL.tw0.9960 <- sum(res.PL.tw0.9960$lpp)
lml.PL.tw0.9950 <- sum(res.PL.tw0.9950$lpp)
lml.PL.tw0.9940 <- sum(res.PL.tw0.9940$lpp)
lml.PL.tw0.9772 <- sum(res.PL.tw0.9772$lpp)
lml.PL.tw0.9500 <- sum(res.PL.tw0.9500$lpp)
lml <- c(lml.PL,lml.PL.tw0.9990,lml.PL.tw0.9981,lml.PL.tw0.9970,lml.PL.tw0.9960,lml.PL.tw0.9950,lml.PL.tw0.9940,lml.PL.tw0.9924,lml.PL.tw0.9772,lml.PL.tw0.9500)
plot(sort(xi.seq),lml[order(xi.seq)],type="b",xlab="xi",ylab="Log prior predictive")

## extract individual log prior predictives
n.models <- length(lml)
lpp <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  lpp[,i] <- eval(parse(text=paste("res.",models[i],"$lpp",sep="")))
}

## cumulative score
lpp.cum <- apply(lpp,2,cumsum)
## which model is the best so far
this.model <- apply(lpp.cum,1,which.max)
this.model.lpp <- this.model
models[this.model]

plot(dates,xi.seq[this.model.lpp],ylab="xi",xlab="Date")

## work out the log score for this sequence of models
lpp.opt <- lpp[1,this.model[1]]
for(i in 2:length(this.model)){
  lpp.opt <- c(lpp.opt,lpp[i,this.model[i-1]])
}
ts.plot(cumsum(lpp.opt)-lpp.cum,col=1:10)

##adaptive choice of $xi$. 
for(i in 1:77){
  plot(xi.seq,lpp.cum[i,])
}

## optimal single choice over all 77 races is 
models[which.max(lpp.cum[77,])]

plot(dates,cumsum(lpp.opt)-lpp.cum[,1],type="l") ## comp to non-tw
##plot(dates,cumsum(lpp.opt)-lpp.cum[,which.max(lpp.cum[77,])],type="l")
##plot(dates,cumsum(lpp.opt)-lpp.cum[,3],type="l")

## the above looks better than our single choice of xi=0.9981

## what are the implications for the other features if we choose the optimal based on log prior predictive at each stage?

## predictive probabilities based on optimal lml

probs.lml <-  probs.opt
for(i in 2:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.lml[i,j,1] <- modelprobs[this.model[i-1],i,j,1]
    probs.lml[i,j,2] <- modelprobs[this.model[i-1],i,j,2]
    probs.lml[i,j,3] <- modelprobs[this.model[i-1],i,j,3]
  }
}

## log score for optimal lml
lsws.lml <- rep(0,77)
lsps.lml <- rep(0,77)
lst10s.lml <- rep(0,77)
for(i in 1:77){
  lsws.lml[i] <- sum(log(probs.lml[i,x[i,1],1]),log(1-probs.lml[i,x[i,-1][is.na(x[i,-1])==FALSE],1]))
  lsps.lml[i] <- sum(log(probs.lml[i,x[i,1:3],2]),log(1-probs.lml[i,x[i,-(1:3)][is.na(x[i,-(1:3)])==FALSE],2]))
  lst10s.lml[i] <- sum(log(probs.lml[i,x[i,1:10],3]),log(1-probs.lml[i,x[i,-(1:10)][is.na(x[i,-(1:10)])==FALSE],3]))
}

ts.plot(cumsum(lsws.lml)-cumsum(lsws[,1]))
ts.plot(cumsum(lsps.lml)-cumsum(lsps[,1]))
ts.plot(cumsum(lst10s.lml)-cumsum(lst10s[,1]))
## on first inspection looks to be not really any better (possibly much worse) 


##########################################################
## model averaged xi (based on lml)

## posterior model probability - assume equal prior
post.mod.prob <- lpp
for(i in 1:n.models){
  sum.ml.rat <- 0
  for(j in 1:n.models){
    sum.ml.rat <- sum.ml.rat+exp(cumsum(lpp[,j])-cumsum(lpp[,i]))
  }
  post.mod.prob[,i] <- 1/sum.ml.rat
}
ts.plot(post.mod.prob,col=1:10,ylim=c(0,1))


## average xi by posterior probs

plot(dates,post.mod.prob%*%matrix(xi.seq),lwd=2,type="l",ylab="xi",xlab="Dates",ylim=c(0.99,1))
plot(dates,log(0.5)/log(post.mod.prob%*%matrix(xi.seq)),lwd=2,type="l",ylab="Posterior mean half-life (days)")


postscript("../postxiopt.eps",pointsize=18)
plot(dates,post.mod.prob%*%matrix(xi.seq),lwd=2,type="l",ylab="xi",xlab="Dates",ylim=c(0.99,1))
points(dates,xi.seq[this.model.lpp])
dev.off()



## model averaged log prior predictives
lpp.modavg <- lpp[1,1] ## as all the same
for(i in 2:dim(lpp)[1]){  
  lpp.modavg[i] <- log(sum(post.mod.prob[i-1,]*exp(lpp.cum[i,]-max(lpp.cum[i,]))))+max(lpp.cum[i,])
}

## close to optimal if we chose a single value
lpp.modavg[77]-lml

## model averaged predictive probabilities

probs.modavg <-  probs.opt
for(i in 2:dim(res.a$pp)[1]){
  for(j in 1:dim(res.a$pp)[2]){
    probs.modavg[i,j,1] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,1])
    probs.modavg[i,j,2] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,2])
    probs.modavg[i,j,3] <- sum(post.mod.prob[i-1,]*modelprobs[,i,j,3])
  }
}

par(mfrow=c(2,3))
plot(dates,res.PL$pp[,4,1],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,1],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,1],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,2],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,2],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,4,3],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,4,3],col=3,lty=3,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,1],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,1],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,1],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,2],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,2],col=3,lty=3,lwd=2) 
legend(as.Date("2010-3-1"),1,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
plot(dates,res.PL$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.opt[,7,3],col=2,lty=2,lwd=2) 
lines(dates,probs.modavg[,7,3],col=3,lty=3,lwd=2) 
legend(as.Date("2012-6-1"),0.4,legend=c("PL","PL.tw","PL.tw.ma"),col=1:3,lty=1:3,lwd=rep(2,2),cex=0.6)
par(mfrow=c(1,1))
## model averaged is close to individual optimal

## log score for model averaging!!
lsws.modavg <- rep(0,77)
lsps.modavg <- rep(0,77)
lst10s.modavg <- rep(0,77)
for(i in 1:77){
  lsws.modavg[i] <- sum(log(probs.modavg[i,x[i,1],1]),log(1-probs.modavg[i,x[i,-1][is.na(x[i,-1])==FALSE],1]))
  lsps.modavg[i] <- sum(log(probs.modavg[i,x[i,1:3],2]),log(1-probs.modavg[i,x[i,-(1:3)][is.na(x[i,-(1:3)])==FALSE],2]))
  lst10s.modavg[i] <- sum(log(probs.modavg[i,x[i,1:10],3]),log(1-probs.modavg[i,x[i,-(1:10)][is.na(x[i,-(1:10)])==FALSE],3]))
}

## comparison of model averaged and optimal based on lml

postscript("../lscore-comp-xi-PL.eps",pointsize=16)
par(mfrow=c(2,2))
plot(dates,cumsum(lsws.modavg)-cumsum(lsws[,1]),ylab="Log score relative to Plackett-Luce",ylim=c(-35,35),type="l",lwd=2,main="Winner",xlab="Date")
lines(dates,cumsum(lsws.lml)-cumsum(lsws[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lsws.opt)-cumsum(lsws[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2012-3-1"),35,legend=c("PL.tw.ma","PL.tw.lpp","PL.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,cumsum(lsps.modavg)-cumsum(lsps[,1]),ylab="Log score relative to Plackett-Luce",ylim=c(-35,35),type="l",lwd=2,main="Top 3",xlab="Date")
lines(dates,cumsum(lsps.lml)-cumsum(lsps[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lsps.opt)-cumsum(lsps[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2012-9-1"),35,legend=c("PL.tw.ma","PL.tw.lpp","PL.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,cumsum(lst10s.modavg)-cumsum(lst10s[,1]),ylab="Log score relative to Plackett-Luce",ylim=c(-35,35),type="l",lwd=2,main="Top 10",xlab="Date")
lines(dates,cumsum(lst10s.lml)-cumsum(lst10s[,1]),col=2,lty=2,lwd=2)
lines(dates,cumsum(lst10s.opt)-cumsum(lst10s[,1]),col=4,lty=4,lwd=2)
abline(h=0,col=8)
legend(as.Date("2012-9-1"),35,legend=c("PL.tw.ma","PL.tw.lpp","PL.tw.opt"),col=c(1,2,4),lty=c(1,2,4),lwd=rep(2,3),cex=0.8)

plot(dates,lpp.modavg-cumsum(lpp[,1]),ylab="Log score relative to Plackett-Luce",ylim=c(-35,35),type="l",lwd=2,main="Full finishing order",xlab="Date")
lines(dates,cumsum(lpp.opt)-cumsum(lpp[,1]),col=2,lty=2,lwd=2)
abline(h=0,col=8)
legend(as.Date("2010-3-1"),-10,legend=c("PL.tw.ma","PL.tw.lpp=PL.tw.opt"),col=c(1,2),lty=c(1,2),lwd=rep(2,2),cex=0.8)
par(mfrow=c(1,1))
dev.off()

lpp.opt.PL <- lpp.opt
lsws.lml.PL <- lsws.lml
lsps.lml.PL <- lsps.lml
lst10s.lml.PL <- lst10s.lml
probs.lml.PL <- probs.lml
this.model.lpp.PL <- this.model.lpp

####################################################################
## Plots for the revised paper

## focus on these models
models <- c("a","a.tw","PL","PL.tw","PL.t6","PL.t10","PL.t14")
n.models <- length(models)
model.names <- c("a","a.tw","PL","PL.tw","PL.t6","PL.t10","PL.t14")

##############################################################
## Winner, podium and points - log scoring rule

lsws <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  if(models[i]=="a.tw"){
     lsws[,i] <- lsws.lml.a 
  }
  else{
    if(models[i]=="PL.tw"){
      lsws[,i] <- lsws.lml.PL
    }
    else{
      lsws[,i] <- eval(parse(text=paste("res.",models[i],"$lsws",sep="")))
    }
  }
}

lsps <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  if(models[i]=="a.tw"){
    lsps[,i] <- lsps.lml.a 
  }
  else{
    if(models[i]=="PL.tw"){
      lsps[,i] <- lsps.lml.PL
    }
    else{
      lsps[,i] <- eval(parse(text=paste("res.",models[i],"$lsps",sep="")))
    }
  }
}

lst10s <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  if(models[i]=="a.tw"){
    lst10s[,i] <- lst10s.lml.a 
  }
  else{
    if(models[i]=="PL.tw"){
      lst10s[,i] <- lst10s.lml.PL
    }
    else{
      lst10s[,i] <- eval(parse(text=paste("res.",models[i],"$lst10s",sep="")))
    }
  }
}

postscript("../prob-vettel-button-noclosedform.eps",pointsize=16)
par(mfrow=c(2,3))
plot(dates,res.a$pp[,4,1],col=1,type="l",ylim=c(0,0.65),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.lml.a[,4,1],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,4,1],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,4,1],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,4,1],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,4,1],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,4,1],col=7,lty=7,lwd=2)
legend(as.Date("2010-3-1"),0.65,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
##text(dates,rnorm(n,0,0.01),labels=as.character(apply(x==4,1,which.max)),cex=0.6)
plot(dates,res.a$pp[,4,2],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,probs.lml.a[,4,2],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,4,2],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,4,2],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,4,2],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,4,2],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,4,2],col=7,lty=7,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,4,3],col=1,type="l",ylim=c(0,1),main=Drivers[4],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,probs.lml.a[,4,3],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,4,3],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,4,3],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,4,3],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,4,3],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,4,3],col=7,lty=7,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)

plot(dates,res.a$pp[,7,1],col=1,type="l",ylim=c(0,0.4),main=Drivers[7],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,probs.lml.a[,7,1],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,7,1],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,7,1],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,7,1],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,7,1],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,7,1],col=7,lty=7,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
##text(dates,rnorm(n,0,0.01),labels=as.character(apply(x==4,1,which.max)),cex=0.6)
plot(dates,res.a$pp[,7,2],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 3 finish",lwd=2)
lines(dates,probs.lml.a[,7,2],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,7,2],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,7,2],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,7,2],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,7,2],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,7,2],col=7,lty=7,lwd=2)
legend(as.Date("2010-3-1"),1,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,res.a$pp[,7,3],col=1,type="l",ylim=c(0,1),main=Drivers[7],xlab="Date",ylab="Probability of a top 10 finish",lwd=2)
lines(dates,probs.lml.a[,7,3],col=2,lty=2,lwd=2)
lines(dates,res.PL$pp[,7,3],type="l",col=3,lty=3,lwd=2)
lines(dates,probs.lml.PL[,7,3],col=4,lty=4,lwd=2)
lines(dates,res.PL.t6$pp[,7,3],col=5,lty=5,lwd=2)
lines(dates,res.PL.t10$pp[,7,3],col=6,lty=6,lwd=2)
lines(dates,res.PL.t14$pp[,7,3],col=7,lty=7,lwd=2)
legend(as.Date("2011-3-1"),0.4,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
par(mfrow=c(1,1))
dev.off()

## comparison of Monte Carlo-based and closed-form probabilities
postscript("../prob-vettel-pl-closedform.eps",pointsize=16)
plot(dates,res.PL$pp[,4,1],col=1,type="l",ylim=c(0,0.2),main=Drivers[4],xlab="Date",ylab="Probability of winning",lwd=2)
lines(dates,res.PL$ppw[,4],col=2,lty=2,lwd=2)
dev.off()


postscript("lscore-comp.eps",pointsize=16)
par(mfrow=c(3,3))
par(mar=c(4, 4, 1, 2)+0.1)
boxplot(lsws,names=model.names,ylab="Log score",main="Winner",las=3)
boxplot(lsps,names=model.names,ylab="Log score",main="Top 3",las=3)
boxplot(lst10s,names=model.names,ylab="Log score",main="Top 10",las=3)
plot(dates,cumsum(lsws[,1])-cumsum(lsws[,3]),type="l",lwd=2,ylim=c(-5,50),xlab="Date",ylab="Log score relative to PL",main="Winner")
for(i in 1:n.models){
  lines(dates,cumsum(lsws[,i])-cumsum(lsws[,3]),type="l",lwd=2,lty=i,col=i)
}
legend(as.Date("2010-3-1"),50,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,cumsum(lsps[,1])-cumsum(lsps[,3]),type="l",lwd=2,ylim=c(-20,120),xlab="Date",ylab="Log score relative to PL",main="Top 3")
for(i in 1:n.models){
  lines(dates,cumsum(lsps[,i])-cumsum(lsps[,3]),type="l",lwd=2,lty=i,col=i)
}
legend(as.Date("2010-3-1"),120,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)
plot(dates,cumsum(lst10s[,1])-cumsum(lst10s[,3]),type="l",lwd=2,ylim=c(-30,120),xlab="Date",ylab="Log score relative to PL",main="Top 10")
for(i in 1:n.models){
  lines(dates,cumsum(lst10s[,i])-cumsum(lst10s[,3]),type="l",lwd=2,lty=i,col=i)
}
legend(as.Date("2010-3-1"),120,legend=model.names,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.6)

plot(dates,cumsum(lsws[,2])-cumsum(lsws[,1]),type="l",lwd=2,ylim=c(-5,15),xlab="Date",ylab="S(a.tw)-S(a)",main="Winner")
abline(0,0,col=8)
plot(dates,cumsum(lsps[,2])-cumsum(lsps[,1]),type="l",lwd=2,ylim=c(-5,15),xlab="Date",ylab="S(a.tw)-S(a)",main="Top 3")
abline(0,0,col=8)
plot(dates,cumsum(lst10s[,2])-cumsum(lst10s[,1]),type="l",lwd=2,ylim=c(-5,15),xlab="Date",ylab="S(a.tw)-S(a)",main="Top 10")
abline(0,0,col=8)
par(mar=c(5, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))
dev.off()

## comparison of observed and expected

observed.stats <- matrix(0,nrow=K,ncol=3)
for(i in 1:K){
  observed.stats[i,1] <- sum(x[,1]==i)
  observed.stats[i,2] <- sum(apply(x[,1:3]==i,1,sum))
  observed.stats[i,3] <- sum(apply(x[,1:10]==i,1,sum))
}

expected.stats <- array(0,c(7,dim(observed.stats)))
expected.stats[1,,] <- apply(res.a$pp,c(2,3),sum)
expected.stats[2,,] <- apply(probs.lml.a,c(2,3),sum)
expected.stats[3,,] <- apply(res.PL$pp,c(2,3),sum)
expected.stats[4,,] <- apply(probs.lml.PL,c(2,3),sum)
expected.stats[5,,] <- apply(res.PL.t6$pp,c(2,3),sum)
expected.stats[6,,] <- apply(res.PL.t10$pp,c(2,3),sum)
expected.stats[7,,] <- apply(res.PL.t14$pp,c(2,3),sum)


this.subset <- c(4,1,3,7,8,5,34,31,2,6)

write(t(t(data.frame(Drivers[this.subset],observed.stats[this.subset,],round(expected.stats[1,this.subset,],2),round(expected.stats[2,this.subset,],2),round(expected.stats[3,this.subset,],2),round(expected.stats[4,this.subset,],2),round(expected.stats[5,this.subset,],2),round(expected.stats[6,this.subset,],2),round(expected.stats[7,this.subset,],2)))),file="obsexp.txt",ncolumns=10)

write(t(t(data.frame(Drivers[this.subset],observed.stats[this.subset,],observed.stats[this.subset,]-round(expected.stats[1,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[2,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[3,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[4,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[5,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[6,this.subset,],1),observed.stats[this.subset,]-round(expected.stats[7,this.subset,],1)))),file="obsminexp.txt",ncolumns=10)


## parameter posteriors
load("res.a.tw.last")
names(res.a.tw.last)

postscript("param-post-tw.eps",pointsize=16)
par(mfrow=c(1,2))
plot(density(res.a.tw.last$gm[,4]),col=1,lwd=2,xlim=c(0,0.4),main="",xlab="lambda")
lines(density(res.a.tw.last$gm[,1]),col=2,lwd=2,lty=2)
lines(density(res.a.tw.last$gm[,34]),col=4,lwd=2,lty=4)
lines(seq(0,0.4,0.001),dgamma(seq(0,0.4,0.001),1,1),col=8)
legend(0.1,40,legend=Drivers[c(4,1,34)],lwd=rep(2,3),lty=c(1,2,4),col=c(1,2,4),cex=0.8)
plot(density(res.a.tw.last$gm[,2]/(res.a.tw.last$gm[,2]+res.a.tw.last$gm[,4])),col=1,lwd=2,xlim=c(0,1),main="",xlab="Probability of head-to-head win")
lines(density(res.a.tw.last$gm[,34]/(res.a.tw.last$gm[,34]+res.a.tw.last$gm[,4])),col=2,lwd=2,lty=2)
lines(density(res.a.tw.last$gm[,34]/(res.a.tw.last$gm[,34]+res.a.tw.last$gm[,1])),col=4,lwd=2,lty=4)
lines(c(0,1),dunif(c(0,1),0,1),col=8)
legend(0,10,legend=c("Vettel vs Alonso","Vettel vs Raikkonen","Alonso vs Raikkonen"),lwd=rep(2,3),lty=c(1,2,4),col=c(1,2,4),cex=0.8)
par(mfrow=c(1,1))
dev.off()

## championship probs
postscript("pchamp-vettel.eps",pointsize=14)
par(mfrow=c(2,1))
par(mar=c(4, 4, 1, 2)+0.1)
plot(1:n,apply(res.a$cw==4,1,mean),type="l",xlab="Race",ylab="P(Champion)",ylim=c(0,1.1),main="",lab=c(20,5,7),lwd=2)
lines(1:n,champprobs.optlml.a[,4],type="l",col=2,lty=2,lwd=2)
lines(1:n,apply(res.PL.t6$cw==4,1,mean),type="l",col=5,lty=5,lwd=2)
text(1:n,rep(0,n),labels=as.character(apply(x==4,1,which.max)),cex=0.4)
abline(v=length(grep("2010",dates))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates)))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates),grep("2012",dates)))+0.5,col="grey",lty=2)
text(9,1.05,"2010")
text(28,1.05,"2011")
text(48,1.05,"2012")
text(69,1.05,"2013")
legend(65,0.55,legend=model.names[c(1,2,5)],lty=c(1,2,5),col=c(1,2,5),lwd=rep(2,3),cex=0.9)
plot(1:n,cumsum(apply(res.a$cw==4,1,mean))-cumsum(apply(res.a$cw==4,1,mean)),type="l",xlab="Race",ylab="Difference from attrition",ylim=c(-4,1.2),main="",lab=c(20,5,7),lwd=2)
lines(1:n,cumsum(champprobs.optlml.a[,4])-cumsum(apply(res.a$cw==4,1,mean)),type="l",col=2,lty=2,lwd=2)
lines(1:n,cumsum(apply(res.PL.t6$cw==4,1,mean))-cumsum(apply(res.a$cw==4,1,mean)),type="l",col=5,lty=5,lwd=2)
abline(v=length(grep("2010",dates))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates)))+0.5,col="grey",lty=2)
abline(v=length(c(grep("2010",dates),grep("2011",dates),grep("2012",dates)))+0.5,col="grey",lty=2)
text(9,1.1,"2010")
text(28,1.1,"2011")
text(48,1.1,"2012")
text(69,1.1,"2013")
legend(5,-1.5,legend=model.names[c(1,2,5)],lty=c(1,2,5),col=c(1,2,5),lwd=rep(2,3),cex=0.9)
text(1:n,rep(-4,n),labels=as.character(apply(x==4,1,which.max)),cex=0.4)
par(mar=c(1, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))
dev.off()

## this has changed to include posteriors
postscript("pchamp-2013.eps",pointsize=18)
par(mfrow=c(2,2))
par(mar=c(4, 4, 1, 2)+0.1)
plot(dates[grep("2013",dates)],apply(res.a$cw==4,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,1),main=Drivers[4],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Vettel[,1]),1/Vettel[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==4,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
##lines(dates[grep("2013",dates)],champprobs.optlml.a[,4][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.25,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(dates[grep("2013",dates)],apply(res.a$cw==1,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,0.6),main=Drivers[1],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Alonso[,1]),1/Alonso[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==1,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
##lines(dates[grep("2013",dates)],champprobs.optlml.a[,1][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.6,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(dates[grep("2013",dates)],apply(res.a$cw==34,1,mean)[grep("2013",dates)],type="l",xlab="",ylab="P(Champion)",ylim=c(0,0.6),main=Drivers[34],lab=c(20,5,7),lwd=2,xlim=as.Date(c("2013-1-1","2013-11-24")),cex.main=0.8)
points(as.Date(Raikkonen[,1]),1/Raikkonen[,3],type="l",col=8)
lines(dates[grep("2013",dates)],apply(res.PL$cw==34,1,mean)[grep("2013",dates)],type="l",col=4,lty=4,lwd=2)
##lines(dates[grep("2013",dates)],champprobs.optlml.a[,34][grep("2013",dates)],type="l",col=3,lty=4,lwd=2)
legend(as.Date("2013-8-1"),0.6,legend=c("Attrition","Plackett-Luce","Sky BET"),lty=c(1,4,1),col=c(1,4,8),lwd=c(2,2,1),cex=0.6)
plot(density(res.a.tw.last$gm[,4]),col=1,lwd=2,xlim=c(0,0.4),main="",xlab="lambda")
lines(density(res.a.tw.last$gm[,1]),col=2,lwd=2,lty=2)
lines(density(res.a.tw.last$gm[,34]),col=4,lwd=2,lty=4)
lines(seq(0,0.4,0.001),dgamma(seq(0,0.4,0.001),1,1),col=8)
legend(0.1,40,legend=Drivers[c(4,1,34)],lwd=rep(2,3),lty=c(1,2,4),col=c(1,2,4),cex=0.8)
##plot(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),4]),type="s",lwd=2,xlab="",ylab="Drivers' Championship points")
##lines(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),1]),col=2,lwd=2,lty=2,type="s")
##lines(dates[grep("2013",dates)],cumsum(points.actual[grep("2013",dates),34]),col=4,lwd=2,lty=4,type="s")
##legend(as.Date("2013-4-1"),350,legend=Drivers[c(4,1,34)],lty=c(1,2,4),col=c(1,2,4),lwd=c(2,2,2),cex=0.6)
par(mar=c(5, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))
dev.off()

## log marginal likelihoods

lml.PL <- sum(res.PL$lpp)
lml.a <- sum(res.a$lpp)
lml.PL.t6 <- sum(res.PL.t6$lpp)
lml.PL.t10 <- sum(res.PL.t10$lpp)
lml.PL.t14 <- sum(res.PL.t14$lpp)
lml <- c(lml.a,sum(lpp.opt.a),lml.PL,sum(lpp.opt.PL),lml.PL.t6,lml.PL.t10,lml.PL.t14)

n.models <- length(lml)
lpp <- matrix(0,nrow=n,ncol=n.models)
for(i in 1:n.models){
  if(models[i]=="a.tw"){
    lpp[,i] <- lpp.opt.a
  }
  else{
    if(models[i]=="PL.tw"){
      lpp[,i] <- lpp.opt.PL
    }
    else{
      lpp[,i] <- eval(parse(text=paste("res.",models[i],"$lpp",sep="")))
    }
  }
}  

## posterior model probability - assume equal prior
post.mod.prob <- lpp
for(i in 1:n.models){
  sum.ml.rat <- 0
  for(j in 1:n.models){
    sum.ml.rat <- sum.ml.rat+exp(cumsum(lpp[,j])-cumsum(lpp[,i]))
  }
  post.mod.prob[,i] <- 1/sum.ml.rat
}


##log prior predictives under various models
postscript("lpp-comp.eps",pointsize=16)
par(mfrow=c(2,2))
par(mar=c(4, 4, 1, 2)+0.1)
boxplot(lpp,names=model.names,ylab="Log prior predictive",main="",las=3)
plot(dates,cumsum(lpp[,1])-cumsum(lpp[,3]),type="l",lwd=2,xlab="Date",ylab="Log Bayes factor relative to PL",ylim=c(-400,400))
for(i in 1:n.models){
  lines(dates,cumsum(lpp[,i])-cumsum(lpp[,3]),type="l",lwd=2,lty=i,col=i)
}
legend(as.Date("2010-5-1"),-100,legend=models,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.5)
plot(dates,cumsum(lpp[,2])-cumsum(lpp[,1]),type="l",lwd=2,xlab="Date",ylab="Log Bayes factor for a.tw vs a",ylim=c(-3,25))
abline(0,0,col=8)
abline(h=5,col=8,lty=4)
abline(h=3,col=8,lty=4)
abline(h=1,col=8,lty=4)
abline(h=-1,col=8,lty=4)
abline(h=-3,col=8,lty=4)
##lines(dates,cumsum(lpp[,3])-cumsum(lpp[,1]),type="l",lwd=2,col=2)
##lines(dates,cumsum(lpp[,4])-cumsum(lpp[,1]),type="l",lwd=2,col=3)
##lines(dates,cumsum(lpp[,5])-cumsum(lpp[,1]),type="l",lwd=2,col=4)
##lines(dates,cumsum(lpp[,6])-cumsum(lpp[,1]),type="l",lwd=2,col=5)
## see table in Kass and Raftery (1995) 2ln(B) > 10 implies very strong evidence!
plot(dates,post.mod.prob[,1],type="l",lwd=2,xlab="Date",ylab="Prior model probabilities",ylim=c(0,1))
for(i in 1:n.models){
  lines(dates,post.mod.prob[,i],type="l",lwd=2,lty=i,col=i)
}
legend(as.Date("2013-1-1"),0.5,legend=models,col=1:n.models,lty=1:n.models,lwd=rep(2,n.models),cex=0.5)
dev.off()		


