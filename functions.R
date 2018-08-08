## Functions relating to 
## "A comparison of truncated and time-weighted Plackett-Luce models for probabilistic forecasting of Formula One results"
## 
## daniel.henderson@ncl.ac.uk


points.calc.2013 <- function(pos)
{
  ## returns points (under 2013 regulations) for a given finishing position
  points <- c(25,18,15,12,10,8,6,4,2,1,0)
  pos <- ifelse(pos>10,11,pos)

  points[pos]
}

rPL <- function(lambda=rep(1,3))
{
  ## generate permuation from Plackett-Luce model for K individuals with
  ## ability parameters lambda

  ## use latent variables to generate the data
  K <- length(lambda)
  Z <- rexp(K,lambda)
  x <- order(Z)
  x
}

likelihood.pl <- function(x,lambda)
{
  ## function for computing Plackett-Luce likelihood
  ## of lambda given single ranking x
  
  ## calculate length of vector x (this is number of drivers in race)
  p <- length(x[is.na(x)==FALSE])
  ## initialise likelihood
  like <- 1
  ## compute likelihood as product of probabilities
  for(j in 1:(p-1)){
    like <- like*lambda[x[j]]/sum(lambda[x[j:p]])
  }
  return(like)
}

loglikelihood.pl <- function(x,lambda)
{
  ## function for computing Plackett-Luce log likelihood
  ## of lambda given single ranking x
  
  ## calculate length of vector x (this is number of drivers in race)
  p <- length(x[is.na(x)==FALSE])
  ## initialise likelihood
  llike <- 0
  ## compute loglikelihood as sum of log probabilities
  for(j in 1:(p-1)){
    llike <- llike + log(lambda[x[j]])-log(sum(lambda[x[j:p]]))
  }
  return(llike)
}


loglike.pl <- function(x,lambda)
{
  ## Plackett-Luce loglikelihood over multiple rankings
  
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  ll.vec <- rep(0,n)
  for(i in 1:n){
    p.i <- length(x[i,is.na(x[i,])==FALSE])
    for(j in 1:(p.i-1)){
      ll.vec[i] <- ll.vec[i]+log(lambda[x[i,j]])-log(sum(lambda[x[i,j:p.i]]))
    }
  }
  return(ll.vec)
}

summary.stats.pl <- function(x,K)
{
  ## compute summary statistics for a given set of rankings data
  ## for the Plackett-Luce model
  
  ## number of multiple comparisons
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- rep(0,n)
  for(i in 1:n){
    p[i] <- length(x[i,is.na(x[i,])==FALSE])
  }
  
  ## identify individual in last place
  ## and count up times for each individual
  last <- rep(0,n)
  participated <- rep(0,K)
  wlast <- rep(0,K)
  for(i in 1:n){
    participated[x[i,1:p[i]]] <- participated[x[i,1:p[i]]]+1
    last[i] <- x[i,p[i]]
    wlast[last[i]] <- wlast[last[i]]+1
  }
  w <- participated-wlast
  
  ## compute delta
  delta <- array(0,c(n,max(p)-1,K))
  for(k in 1:K){
    for(i in 1:n){
      for(j in 1:(p[i]-1)){
        delta[i,j,k] <- sum(k==x[i,j:p[i]])
      }
    }
  }

  return(list(n=n,p=p,w=w,delta=delta))
}

PL.gibbs <- function(x,K,a=rep(1,K),b=1,lambda.init=rep(1,K),its=1000,burn=0,thin=1)
{
  ## Gibbs sampler for Plackett-Luce model
  ## Started: 18/10/13
  ## Updated: 8/5/14
  
  ## preliminary calculations
  ss <- summary.stats.pl(x,K)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- ss$p
  
  ## identify individual in last place
  ## and count up times for each individual
  w <- ss$w
  
  ## compute delta
  delta <- ss$delta

  ## storage for samples
  Y.res <- array(0,c(its,n,max(p)-1))
  lambda.res <- matrix(0,nrow=its,ncol=K)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  lambda.curr <- lambda.init
  Y.curr <- matrix(0,nrow=n,ncol=(max(p)-1))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      ##print(t-burn)
      for(ell in 1:thin){
        ## sample latent variables | everything else
        for(i in 1:n){
          for(j in 1:(p[i]-1)){
            Y.curr[i,j] <- rexp(1,sum(lambda.curr[x[i,j:p[i]]]))
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          lambda.curr[k] <- rgamma(1,a[k]+w[k],b+sum(delta[,,k]*Y.curr))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Lambda.curr <- rgamma(1,sum(a),b)
        lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)
      }
    }
  }
  for(t in 1:its){
    ##print(t)
    for(ell in 1:thin){
      ## sample latent variables | everything else
      for(i in 1:n){
        for(j in 1:(p[i]-1)){
          Y.curr[i,j] <- rexp(1,sum(lambda.curr[x[i,j:p[i]]]))
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        lambda.curr[k] <- rgamma(1,a[k]+w[k],b+sum(delta[,,k]*Y.curr))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Lambda.curr <- rgamma(1,sum(a),b)
      lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)
    }
    
    ## store sampled values
    lambda.res[t,] <- lambda.curr
    Y.res[t,,] <- Y.curr
    
    ## compute observed data loglikelihood
    ll.res[t,] <- loglike.pl(x,lambda.res[t,])
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(lambda.res[t,],a,b,log=TRUE))
    ##print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(Y=Y.res,lambda=lambda.res,ll=ll.res,ljd=ljd.res))
}

psi <- function(t,xi=1)
{
  ## Geometric weighting function

  xi^as.numeric(t)
}

psi.exp <- function(t,xi=0)
{
  ## Exponential weighting function

  exp(-xi * as.numeric(t))
}

summary.stats.twtpl <- function(x,K,r.vec,tau.vec,tau,xi)
{
  ## compute summary statistics for a given set of rankings data
  ## for the time-weighted truncated Plackett-Luce model
  
  ## number of multiple comparisons
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- rep(0,n)
  for(i in 1:n){
    p[i] <- length(x[i,is.na(x[i,])==FALSE])
  }

  ## the rank at which truncation occurs
  r.star <- rep(0,n)
  for(i in 1:n){
    r.star[i] <- min(r.vec[i],p[i]-1)
  }
  
  ## compute time-weighted number of races in which driver k finishes in the top r.star[i]
  w  <- rep(0,K)
  time.weight <- rep(0,n)
  for(k in 1:K){
    for(i in 1:n){
      time.weight[i] <- psi(tau-tau.vec[i],xi)
      for(j in 1:r.star[i]){
        w[k] <- w[k] + ifelse(x[i,j]==k,1,0)*time.weight[i]
      }
    }
  }
  
  ## compute delta
  delta <- array(0,c(n,max(p)-1,K))
  for(k in 1:K){
    for(i in 1:n){
      for(j in 1:(p[i]-1)){
        delta[i,j,k] <- sum(k==x[i,j:p[i]])
      }
    }
  }

  return(list(n=n,p=p,w=w,delta=delta,r.star=r.star,time.weight=time.weight))
}


loglike.twtpl <- function(x,lambda,r.vec,tau.vec,tau,xi)
{
  ## time-weighted truncated Plackett-Luce loglikelihood over multiple rankings
  
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  ll.vec <- rep(0,n)
  time.weight <- rep(0,n)
  for(i in 1:n){
    time.weight[i] <- psi(tau-tau.vec[i],xi)
    p.i <- length(x[i,is.na(x[i,])==FALSE])
    r.i.star <- min(r.vec[i],p.i-1)
    for(j in 1:r.i.star){
      ll.vec[i] <- ll.vec[i]+time.weight[i]*(log(lambda[x[i,j]])-log(sum(lambda[x[i,j:p.i]])))
    }
  }
  return(ll.vec)
}

twtPL.gibbs <- function(x,K,r.vec,tau.vec,tau,xi=1,a=rep(1,K),b=1,lambda.init=rep(1,K),its=1000,burn=0,thin=1)
{
  ## Gibbs sampler for time-weighted truncated Plackett-Luce model
  
  ## preliminary calculations
  ss <- summary.stats.twtpl(x,K,r.vec,tau.vec,tau,xi)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- ss$p
  
  ## the rank at which truncation occurs
  r.star <- ss$r.star

  ## time-weighted number of races in which driver k finishes in top r.star
  w <- ss$w
  
  ## compute delta
  delta <- ss$delta

  ## compute time-weghting
  time.weight <- ss$time.weight

  ## storage for samples
  Y.res <- array(0,c(its,n,max(p)-1))
  lambda.res <- matrix(0,nrow=its,ncol=K)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  lambda.curr <- lambda.init
  Y.curr <- matrix(0,nrow=n,ncol=(max(p)-1))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      ##print(t-burn)
      for(ell in 1:thin){
        ## sample latent variables | everything else
        for(i in 1:n){
          ##for(j in 1:(p[i]-1)){
          for(j in 1:r.star[i]){
            Y.curr[i,j] <- rgamma(1,time.weight[i],sum(lambda.curr[x[i,j:p[i]]]))
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          ##lambda.curr[k] <- rgamma(1,a[k]+w[k],b+sum(delta[,1:r.star[i],k]*Y.curr[,1:r.star[i]]))
          lambda.curr[k] <- rgamma(1,a[k]+w[k],b+sum(delta[,,k]*Y.curr[,]))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Lambda.curr <- rgamma(1,sum(a),b)
        lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)
      }
    }
  }
  for(t in 1:its){
    ##print(t)
    for(ell in 1:thin){
      ## sample latent variables | everything else
      for(i in 1:n){
        ##for(j in 1:(p[i]-1)){
        for(j in 1:r.star[i]){
          Y.curr[i,j] <- rgamma(1,time.weight[i],sum(lambda.curr[x[i,j:p[i]]]))
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        lambda.curr[k] <- rgamma(1,a[k]+w[k],b+sum(delta[,1:r.star[i],k]*Y.curr[,1:r.star[i]]))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Lambda.curr <- rgamma(1,sum(a),b)
      lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)
    }
    
    ## store sampled values
    lambda.res[t,] <- lambda.curr
    Y.res[t,,] <- Y.curr
    
    ## compute observed data time-weighted loglikelihood
    ll.res[t,] <- loglike.twtpl(x,lambda.res[t,],r.vec,tau.vec,tau,xi)
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(lambda.res[t,],a,b,log=TRUE))
    ##print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(Y=Y.res,lambda=lambda.res,ll=ll.res,ljd=ljd.res))
}

PL.EM <- function(x,K,a=rep(1,K),b=1,lambda.init=rep(1,K),its=1000)
{
  ## EM algorithm Gibbs for Plackett-Luce model
  
  ## preliminary calculations
  ss <- summary.stats.pl(x,K)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- ss$p
  
  ## identify individual in last place
  ## and count up times for each individual
  w <- ss$w
  
  ## compute delta
  delta <- ss$delta

  ## storage for samples
  lambda.res <- matrix(0,nrow=its,ncol=K)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  lambda.curr <- lambda.init
  lambda.prop <- lambda.curr
  
  ## Markov chain sampling
  for(t in 1:its){
    print(t)

    ## sample parameters | everything else
    for(k in 1:K){
      dsum <- 0
      for(i in 1:n){
        for(j in 1:(p[i]-1)){
          sum.lambda <- sum(lambda.curr[x[i,j:p[i]]])
          dsum <- dsum + (delta[i,j,k]/sum.lambda)
        }
      }
       
      lambda.prop[k] <- (a[k]+w[k]-1)/(b+dsum)       
      ##lambda.curr[k] <- (a[k]+w[k]-1)/(b+dsum)
    }
    lambda.curr <- lambda.prop
    lambda.prop <- lambda.curr
    
    ## store sampled values
    lambda.res[t,] <- lambda.curr
    
    ## compute observed data loglikelihood
    ll.res[t,] <- loglike.pl(x,lambda.res[t,])
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(lambda.res[t,],a,b,log=TRUE))
    print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(lambda=lambda.res,ll=ll.res,ljd=ljd.res))
}


chib.ml.pl <- function(res,a,b,ss)
{
  ## Chib's estimate of the log marginal likelihood
  ## for the Plackett-Luce model

  ## number of iterations
  its <- dim(res$lambda)[1]

  ## choose a lambda with high joint density
  ##this.val <- which.max(res$ljd)
  ##lambda.max <- res$lambda[this.val,]
  lambda.max <- apply(res$lambda,2,mean)
  
  ## compute log of full conditionals for lambda.max (Rao-Blackwellized)
  log.cd <- rep(0,its)
  for(i in 1:its){
    sum.ldg <- 0
    for(k in 1:K){                                    
      sum.ldg <- sum.ldg+dgamma(lambda.max[k],a[k]+ss$w[k],b+sum(ss$delta[,,k]*res$Y[i,,]),log=TRUE)
    }                    
    log.cd[i] <- sum.ldg
  }
  log.post <- log(mean(exp(log.cd-max(log.cd))))+max(log.cd)  
  ##log.ml.chib <- res$ljd[this.val]-log.post
  log.ml.chib <- sum(loglike.pl(x,lambda.max))+sum(dgamma(lambda.max,a,b,log=TRUE))-log.post
  log.ml.chib
}

chib.ml.twtpl <- function(res,a,b,ss,r.vec,tau.vec,tau,xi)
{
  ## Chib's estimate of the log marginal likelihood
  ## for the time-weighted truncated Plackett-Luce model

  ## number of iterations
  its <- dim(res$lambda)[1]

  ## choose a lambda with high joint density
  ##this.val <- which.max(res$ljd)
  ##lambda.max <- res$lambda[this.val,]
  lambda.max <- apply(res$lambda,2,mean)
  
  ## compute log of full conditionals for lambda.max (Rao-Blackwellized)
  log.cd <- rep(0,its)
  for(i in 1:its){
    sum.ldg <- 0
    for(k in 1:K){                                    
      sum.ldg <- sum.ldg+dgamma(lambda.max[k],a[k]+ss$w[k],b+sum(ss$delta[,,k]*res$Y[i,,]),log=TRUE)
    }                    
    log.cd[i] <- sum.ldg
  }
  log.post <- log(mean(exp(log.cd-max(log.cd))))+max(log.cd)  
  ##log.ml.chib <- res$ljd[this.val]-log.post
  log.ml.chib <- sum(loglike.twtpl(x,lambda.max,r.vec,tau.vec,tau,xi))+sum(dgamma(lambda.max,a,b,log=TRUE))-log.post
  log.ml.chib
}

MC.ml.pl <- function(x,K,N,a,b)
{
  ## Monte Carlo estimate of marginal likelihood under the Plackett-Luce model

  log.ml.mc <- rep(0,N)
  for(i in 1:N){
    log.ml.mc[i] <- sum(loglike.pl(x,rgamma(K,a,b)))
  }
  log.ml.MC <- log(mean(exp(log.ml.mc-max(log.ml.mc))))+max(log.ml.mc)
  log.ml.MC
}

PL.sequential <- function(x=x,K=K,tau.vec=tau.vec,a=a,b=b,its=1000,burn=100,thin=1,zero.adjust=TRUE)
{
  ## sequential updating and prediction under the Plackett-Luce model

  n <- dim(x)[1]
  ##K <- max(x[is.na(x)==FALSE])
  
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
  
  ## set up matrices and vectors to store results
  log.prior.pred <- rep(0,n)
  log.score.winner.sim <- rep(0,n)
  log.score.podium.sim <- rep(0,n)
  log.score.top10.sim <- rep(0,n)
  points.mean <- matrix(0,nrow=n,ncol=K)
  lambda.mean <- matrix(0,nrow=n,ncol=K)
  lambda.upper <- matrix(0,nrow=n,ncol=K)
  lambda.lower <- matrix(0,nrow=n,ncol=K)
  ## array for storing probabilities of i beating j
  win.pairs.prob <- array(0,c(n,K,K))
  ## array for storing predictive probabilities for winner, podium, points
  pred.probs <- array(0,c(n,K,3))
  ## store analytical predictive probs of winning (PL)
  pred.prob.winner <- matrix(0,nrow=n,ncol=K)
  ## matrix to store sampled championship winners over time
  champ.winner <- matrix(0,nrow=n,ncol=its)
  ## matrix to store expected number of wins in a season for each driver
  expected.season.wins <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of podiums in a season for each driver
  expected.season.podiums <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of top 10s in a season for each driver
  expected.season.top10s <- matrix(0,nrow=n,ncol=K)
  
  ## initial values 
  
  ## log prior predictive for first race (verified by simulation!)
  log.prior.pred[1] <- log(1/factorial(length(x[1,is.na(x[1,])==FALSE])))
  ##MC.ml.pl(x[1,],K,10000,a,b)
  ## log score for "winner" of first race
  log.score.winner.sim[1] <- log(1/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-1)*log(1-(1/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "podium finish" of first race
  log.score.podium.sim[1] <- 3*log(3/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-3)*log(1-(3/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "points finish (top 10)" of first race
  log.score.top10.sim[1] <- 10*log(10/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-10)*log(1-(10/length(x[1,is.na(x[1,])==FALSE])))
  ## predicted championship winners after 1st race
  champ.winner[1,] <- sample(x[1,is.na(x[1,])==FALSE],its,replace=TRUE)
  ## Expected season wins, podiums, top10s
  season.races <- grep(strsplit(as.character(tau.vec[1]),"-")[[1]][1],tau.vec)
  length.season <- length(season.races)
  expected.season.wins[1,x[1,is.na(x[1,])==FALSE]] <- length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.podiums[1,x[1,is.na(x[1,])==FALSE]] <- 3*length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.top10s[1,x[1,is.na(x[1,])==FALSE]] <- 10*length.season/length(x[1,is.na(x[1,])==FALSE])
  ## predictive probabilities for winner, podium, points for first race
  these.drivers <- x[1,is.na(x[1,])==FALSE]
  for(ell in 1:length(these.drivers)){
    pred.probs[1,these.drivers[ell],1] <- 1/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],2] <- 3/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],3] <- 10/length(x[1,is.na(x[1,])==FALSE])
    points.mean[1,these.drivers[ell]] <- sum(points.calc.2013(1:10))/length(x[1,is.na(x[1,])==FALSE])
    pred.prob.winner[1,these.drivers[ell]] <- 1/length(x[1,is.na(x[1,])==FALSE])
  }
  
  ## analyse one race at a time
  for(tt in 1:(n-1)) {
    print(tt)
    lambda.init <- rgamma(K,a,b)
    res.PL <- PL.gibbs(x[1:tt,],K,a=a,b=b,lambda.init=lambda.init,its=its,burn=burn,thin=thin)
    
    ## predictive simulations for remaining races - assuming drivers
    ## from race tt+1
    these.drivers <- x[tt+1,is.na(x[tt+1,])==FALSE]
    season.races <- grep(strsplit(as.character(tau.vec[tt]),"-")[[1]][1],tau.vec)
    length.season <- length(season.races)
    remaining.races <- season.races[season.races>tt]
    past.races <- season.races[season.races<=tt]
    num.past.races <- length(past.races)
    num.remain.races <- length(remaining.races)
    if(num.remain.races==0){
      remaining.races <- grep(strsplit(as.character(tau.vec[tt+1]),"-")[[1]][1],tau.vec)
      num.remain.races <- length(remaining.races)
    }

    x.pred.sim <- array(0,c(its,num.remain.races,length(these.drivers)))
    points.sim <- array(0,c(its,num.remain.races,K))
    is.winner.sim <- array(0,c(its,num.remain.races,K))
    is.podium.sim <- array(0,c(its,num.remain.races,K))
    is.top10.sim <- array(0,c(its,num.remain.races,K))
    ##pos.sim <- array(0,c(its,n-tt,K))
    pred.final.points <- matrix(0,nrow=its,ncol=K)
    pred.final.wins <- matrix(0,nrow=its,ncol=K)
    pred.final.podiums <- matrix(0,nrow=its,ncol=K)
    pred.final.top10s <- matrix(0,nrow=its,ncol=K)
    prob.winner <- matrix(0,nrow=its,ncol=K)
    
    ## for each sampled set of parameter values - future race results (in season)
    ## convert to points for each driver
    ## convert to finishing position for each driver (NOT USED)
    for(i in 1:its){
      for(k in 1:num.remain.races){
        x.pred.sim[i,k,] <- these.drivers[rPL(res.PL$lambda[i,x[tt+1,is.na(x[tt+1,])==FALSE]])]
        points.sim[i,k,x.pred.sim[i,k,]] <- points.calc.2013(1:length(these.drivers))
        is.winner.sim[i,k,x.pred.sim[i,k,1]] <- 1
        is.podium.sim[i,k,x.pred.sim[i,k,1:3]] <- 1
        is.top10.sim[i,k,x.pred.sim[i,k,1:10]] <- 1
        ##pos.sim[i,k,x.pred.sim[i,k,]] <- 1:length(these.drivers)
      }
      if(num.past.races==1){
        pred.final.points[i,] <- points.actual[past.races,]+apply(points.sim[i,,],2,sum)
        pred.final.wins[i,] <- wins.actual[past.races,]+apply(is.winner.sim[i,,],2,sum)
        pred.final.podiums[i,] <- podiums.actual[past.races,]+apply(is.podium.sim[i,,],2,sum)
        pred.final.top10s[i,] <- top10s.actual[past.races,]+apply(is.top10.sim[i,,],2,sum)
      }
      else{
        if(num.past.races == (length.season-1)){
          pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+points.sim[i,,]
          pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+is.winner.sim[i,,]
          pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+is.podium.sim[i,,]
          pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+is.top10.sim[i,,]
        }
        else{
          if(num.past.races == length.season){
            pred.final.points[i,] <- apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(is.top10.sim[i,,],2,sum)
          }
          else{
            pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+apply(is.top10.sim[i,,],2,sum)
          }
        }
      }
      ## extra just for PL (as in closed form)
      for(j in 1:length(these.drivers)){
        prob.winner[i,these.drivers[j]] <- res.PL$lambda[i,these.drivers[j]]/sum(res.PL$lambda[i,these.drivers])
      }
    }
    ## closed form probability of winning
    pred.prob.winner[tt+1,] <- apply(prob.winner,2,mean)
    ## simulated championship winner
    champ.winner[tt+1,] <-  apply(pred.final.points,1,which.max)
    expected.season.wins[tt+1,] <- apply(pred.final.wins,2,mean)
    expected.season.podiums[tt+1,] <- apply(pred.final.podiums,2,mean)
    expected.season.top10s[tt+1,] <- apply(pred.final.top10s,2,mean)
    
    ## expected number of points in next race
    points.mean[tt+1,] <- apply(points.sim[,1,],2,mean)
    
    ## store parameter summaries
    lambda.mean[tt,] <- apply(res.PL$lambda,2,mean)
    lambda.upper[tt,] <- apply(res.PL$lambda,2,quantile,0.975)
    lambda.lower[tt,] <- apply(res.PL$lambda,2,quantile,0.025)
    ##res.PL.pi <- res.PL$lambda/rowSums(res.PL$lambda)
    for(i in 1:K){
      for(j in 1:K){
        win.pairs.prob[tt,i,j] <- mean(res.PL$lambda[,i]/(res.PL$lambda[,i]+res.PL$lambda[,j]))
      }
    }
    
    ## log scores for winner,podium, top10 etc ##########
    pred.prob.win <- NULL
    pred.prob.podium <- NULL
    pred.prob.top10 <- NULL
    ## deal with zero probs - like Bayes estimate of binomial proportion
    ## with beta prior
    a.zero <- 0
    b.zero <- 0
    if(zero.adjust){
      a.zero <- 1
      b.zero <- length(these.drivers)-1
    }
    new.denom <- a.zero+b.zero+its
    for(j in 1:length(these.drivers)){
      pred.prob.win[j] <- mean(x.pred.sim[,1,1]==these.drivers[j])
      pred.prob.podium[j] <- mean(apply(x.pred.sim[,1,1:3]==these.drivers[j],1,sum))
      pred.prob.top10[j] <- mean(apply(x.pred.sim[,1,1:10]==these.drivers[j],1,sum))
      pred.prob.win[j] <- (a.zero+(pred.prob.win[j]*its))/new.denom
      pred.prob.podium[j] <- (a.zero+(pred.prob.podium[j]*its))/new.denom
      pred.prob.top10[j] <- (a.zero+(pred.prob.top10[j]*its))/new.denom
      pred.probs[tt+1,these.drivers[j],1] <- pred.prob.win[j]
      pred.probs[tt+1,these.drivers[j],2] <- pred.prob.podium[j]
      pred.probs[tt+1,these.drivers[j],3] <- pred.prob.top10[j]    
    }
    outcome.win <- c(1,rep(0,length(these.drivers)-1))
    owlppw <- outcome.win*log(pred.prob.win) ## deal with 0 probs
    omowlomppw <- (1-outcome.win)*log(1-pred.prob.win)
    log.score.winner.sim[tt+1] <- sum(ifelse(owlppw=="NaN",0,owlppw)+ifelse(omowlomppw=="NaN",0,omowlomppw))  
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    oplppp <- outcome.podium*log(pred.prob.podium) ## deal with 0 probs
    omoplomppp <- (1-outcome.podium)*log(1-pred.prob.podium)
    log.score.podium.sim[tt+1] <- sum(ifelse(oplppp=="NaN",0,oplppp)+ifelse(omoplomppp=="NaN",0,omoplomppp))  
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ot10lppt10 <- outcome.top10*log(pred.prob.top10) ## deal with 0 probs
    omot10lomppt10 <- (1-outcome.top10)*log(1-pred.prob.top10)
    log.score.top10.sim[tt+1] <- sum(ifelse(ot10lppt10=="NaN",0,ot10lppt10)+ifelse(omot10lomppt10=="NaN",0,omot10lomppt10))  
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))

    ##outcome.win <- c(1,rep(0,length(these.drivers)-1))
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    ##outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    ##outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))
##################################
    
    ## compute log prior predictive
    log.pred.prob.all <- rep(0,its)
    for(i in 1:its) {
      log.pred.prob.all[i] <- loglikelihood.pl(x[tt+1,],res.PL$lambda[i,])
    }
    max.lppa <- max(log.pred.prob.all)
    log.prior.pred[tt+1] <- log(mean(exp(log.pred.prob.all-max.lppa)))+max.lppa
  }
  ## last race of the season
  tt = n
  print(tt)
  lambda.init <- rgamma(K,a,b)
  res.PL <- PL.gibbs(x[1:tt,],K,a=a,b=b,lambda.init=lambda.init,its=its,burn=burn,thin=thin)
  lambda.mean[tt,] <- apply(res.PL$lambda,2,mean)
  lambda.upper[tt,] <- apply(res.PL$lambda,2,quantile,0.975)
  lambda.lower[tt,] <- apply(res.PL$lambda,2,quantile,0.025)
  ##res.PL.pi <- res.PL$lambda/rowSums(res.PL$lambda)
  for(i in 1:K){
    for(j in 1:K){
      win.pairs.prob[tt,i,j] <- mean(res.PL$lambda[,i]/(res.PL$lambda[,i]+res.PL$lambda[,j]))
    }
  }
  
##  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,ppw=pred.prob.winner))
  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,ppw=pred.prob.winner,lm=lambda.mean,lu=lambda.upper,ll=lambda.lower,wpp=win.pairs.prob,pm=points.mean,cw=champ.winner,esw=expected.season.wins,esp=expected.season.podiums,est10=expected.season.top10s))
}

twtPL.sequential <- function(x=x,K=K,r.vec=r.vec,tau.vec=tau.vec,tau=tau,xi=xi,a=a,b=b,its=1000,burn=100,thin=1,zero.adjust=TRUE)
{
  ## sequential updating and prediction under the time-weighted truncated Plackett-Luce model

  n <- dim(x)[1]
  ##K <- max(x[is.na(x)==FALSE])
  
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
  
  ## set up matrices and vectors to store results
  log.prior.pred <- rep(0,n)
  log.score.winner.sim <- rep(0,n)
  log.score.podium.sim <- rep(0,n)
  log.score.top10.sim <- rep(0,n)
  points.mean <- matrix(0,nrow=n,ncol=K)
  lambda.mean <- matrix(0,nrow=n,ncol=K)
  lambda.upper <- matrix(0,nrow=n,ncol=K)
  lambda.lower <- matrix(0,nrow=n,ncol=K)
  ## array for storing probabilities of i beating j
  win.pairs.prob <- array(0,c(n,K,K))
  ## array for storing predictive probabilities for winner, podium, points
  pred.probs <- array(0,c(n,K,3))
  ## store analytical predictive probs of winning (PL)
  pred.prob.winner <- matrix(0,nrow=n,ncol=K)
  ## matrix to store sampled championship winners over time
  champ.winner <- matrix(0,nrow=n,ncol=its)
  ## matrix to store expected number of wins in a season for each driver
  expected.season.wins <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of podiums in a season for each driver
  expected.season.podiums <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of top 10s in a season for each driver
  expected.season.top10s <- matrix(0,nrow=n,ncol=K)
  
  ## initial values 
  
  ## log prior predictive for first race (verified by simulation!)
  log.prior.pred[1] <- log(1/factorial(length(x[1,is.na(x[1,])==FALSE])))
  ##MC.ml.pl(x[1,],K,10000,a,b)
  ## log score for "winner" of first race
  log.score.winner.sim[1] <- log(1/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-1)*log(1-(1/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "podium finish" of first race
  log.score.podium.sim[1] <- 3*log(3/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-3)*log(1-(3/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "points finish (top 10)" of first race
  log.score.top10.sim[1] <- 10*log(10/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-10)*log(1-(10/length(x[1,is.na(x[1,])==FALSE])))
  ## predicted championship winners after 1st race
  champ.winner[1,] <- sample(x[1,is.na(x[1,])==FALSE],its,replace=TRUE)
  ## Expected season wins, podiums, top10s
  season.races <- grep(strsplit(as.character(tau.vec[1]),"-")[[1]][1],tau.vec)
  length.season <- length(season.races)
  expected.season.wins[1,x[1,is.na(x[1,])==FALSE]] <- length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.podiums[1,x[1,is.na(x[1,])==FALSE]] <- 3*length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.top10s[1,x[1,is.na(x[1,])==FALSE]] <- 10*length.season/length(x[1,is.na(x[1,])==FALSE])
  ## predictive probabilities for winner, podium, points for first race
  these.drivers <- x[1,is.na(x[1,])==FALSE]
  for(ell in 1:length(these.drivers)){
    pred.probs[1,these.drivers[ell],1] <- 1/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],2] <- 3/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],3] <- 10/length(x[1,is.na(x[1,])==FALSE])
    points.mean[1,these.drivers[ell]] <- sum(points.calc.2013(1:10))/length(x[1,is.na(x[1,])==FALSE])
    pred.prob.winner[1,these.drivers[ell]] <- 1/length(x[1,is.na(x[1,])==FALSE])
  }
  
  ## analyse one race at a time
  for(tt in 1:(n-1)) {
    print(tt)
    lambda.init <- rgamma(K,a,b)
    res.PL <- twtPL.gibbs(x[1:tt,],K,r.vec=r.vec[1:tt],tau.vec=tau.vec[1:tt],tau=tau.vec[tt+1],xi=xi,a=a,b=b,lambda.init=lambda.init,its=its,burn=burn,thin=thin)
    
    ## predictive simulations for remaining races - assuming drivers
    ## from race tt+1
    these.drivers <- x[tt+1,is.na(x[tt+1,])==FALSE]
    season.races <- grep(strsplit(as.character(tau.vec[tt]),"-")[[1]][1],tau.vec)
    length.season <- length(season.races)
    remaining.races <- season.races[season.races>tt]
    past.races <- season.races[season.races<=tt]
    num.past.races <- length(past.races)
    num.remain.races <- length(remaining.races)
    if(num.remain.races==0){
      remaining.races <- grep(strsplit(as.character(tau.vec[tt+1]),"-")[[1]][1],tau.vec)
      num.remain.races <- length(remaining.races)
    }

    x.pred.sim <- array(0,c(its,num.remain.races,length(these.drivers)))
    points.sim <- array(0,c(its,num.remain.races,K))
    is.winner.sim <- array(0,c(its,num.remain.races,K))
    is.podium.sim <- array(0,c(its,num.remain.races,K))
    is.top10.sim <- array(0,c(its,num.remain.races,K))
    ##pos.sim <- array(0,c(its,n-tt,K))
    pred.final.points <- matrix(0,nrow=its,ncol=K)
    pred.final.wins <- matrix(0,nrow=its,ncol=K)
    pred.final.podiums <- matrix(0,nrow=its,ncol=K)
    pred.final.top10s <- matrix(0,nrow=its,ncol=K)
    prob.winner <- matrix(0,nrow=its,ncol=K)
    
    ## for each sampled set of parameter values - future race results
    ## convert to points for each driver
    ## convert to finishing position for each driver (NOT USED)
    for(i in 1:its){
      for(k in 1:num.remain.races){
        x.pred.sim[i,k,] <- these.drivers[rPL(res.PL$lambda[i,x[tt+1,is.na(x[tt+1,])==FALSE]])]
        points.sim[i,k,x.pred.sim[i,k,]] <- points.calc.2013(1:length(these.drivers))
        is.winner.sim[i,k,x.pred.sim[i,k,1]] <- 1
        is.podium.sim[i,k,x.pred.sim[i,k,1:3]] <- 1
        is.top10.sim[i,k,x.pred.sim[i,k,1:10]] <- 1
        ##pos.sim[i,k,x.pred.sim[i,k,]] <- 1:length(these.drivers)
      }
      ## predicted end of season points for each driver
      if(num.past.races==1){
        pred.final.points[i,] <- points.actual[past.races,]+apply(points.sim[i,,],2,sum)
        pred.final.wins[i,] <- wins.actual[past.races,]+apply(is.winner.sim[i,,],2,sum)
        pred.final.podiums[i,] <- podiums.actual[past.races,]+apply(is.podium.sim[i,,],2,sum)
        pred.final.top10s[i,] <- top10s.actual[past.races,]+apply(is.top10.sim[i,,],2,sum)
      }
      else{
        if(num.past.races == (length.season-1)){
          pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+points.sim[i,,]
          pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+is.winner.sim[i,,]
          pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+is.podium.sim[i,,]
          pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+is.top10.sim[i,,]
        }
        else{
          if(num.past.races == length.season){
            pred.final.points[i,] <- apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(is.top10.sim[i,,],2,sum)
          }
          else{
            pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+apply(is.top10.sim[i,,],2,sum)
          }
        }
      }
      ## extra just for PL (as in closed form)
      for(j in 1:length(these.drivers)){
        prob.winner[i,these.drivers[j]] <- res.PL$lambda[i,these.drivers[j]]/sum(res.PL$lambda[i,these.drivers])
      }
    }
    ## closed form probability of winning
    pred.prob.winner[tt+1,] <- apply(prob.winner,2,mean)
    ## simulated championship winner
    champ.winner[tt+1,] <-  apply(pred.final.points,1,which.max)
    expected.season.wins[tt+1,] <- apply(pred.final.wins,2,mean)
    expected.season.podiums[tt+1,] <- apply(pred.final.podiums,2,mean)
    expected.season.top10s[tt+1,] <- apply(pred.final.top10s,2,mean)
    
    ## expected number of points in next race
    points.mean[tt+1,] <- apply(points.sim[,1,],2,mean)
    
    ## store parameter summaries
    lambda.mean[tt,] <- apply(res.PL$lambda,2,mean)
    lambda.upper[tt,] <- apply(res.PL$lambda,2,quantile,0.975)
    lambda.lower[tt,] <- apply(res.PL$lambda,2,quantile,0.025)
    ##res.PL.pi <- res.PL$lambda/rowSums(res.PL$lambda)
    for(i in 1:K){
      for(j in 1:K){
        win.pairs.prob[tt,i,j] <- mean(res.PL$lambda[,i]/(res.PL$lambda[,i]+res.PL$lambda[,j]))
      }
    }
    
    ## log scores for winner,podium, top10 etc ##########
    pred.prob.win <- NULL
    pred.prob.podium <- NULL
    pred.prob.top10 <- NULL
    ## deal with zero probs - like Bayes estimate of binomial proportion
    ## with beta prior
    a.zero <- 0
    b.zero <- 0
    if(zero.adjust){
      a.zero <- 1
      b.zero <- length(these.drivers)-1
    }
    new.denom <- a.zero+b.zero+its
    for(j in 1:length(these.drivers)){
      pred.prob.win[j] <- mean(x.pred.sim[,1,1]==these.drivers[j])
      pred.prob.podium[j] <- mean(apply(x.pred.sim[,1,1:3]==these.drivers[j],1,sum))
      pred.prob.top10[j] <- mean(apply(x.pred.sim[,1,1:10]==these.drivers[j],1,sum))
      pred.prob.win[j] <- (a.zero+(pred.prob.win[j]*its))/new.denom
      pred.prob.podium[j] <- (a.zero+(pred.prob.podium[j]*its))/new.denom
      pred.prob.top10[j] <- (a.zero+(pred.prob.top10[j]*its))/new.denom
      pred.probs[tt+1,these.drivers[j],1] <- pred.prob.win[j]
      pred.probs[tt+1,these.drivers[j],2] <- pred.prob.podium[j]
      pred.probs[tt+1,these.drivers[j],3] <- pred.prob.top10[j]    
    }
    outcome.win <- c(1,rep(0,length(these.drivers)-1))
    owlppw <- outcome.win*log(pred.prob.win) ## deal with 0 probs
    omowlomppw <- (1-outcome.win)*log(1-pred.prob.win)
    log.score.winner.sim[tt+1] <- sum(ifelse(owlppw=="NaN",0,owlppw)+ifelse(omowlomppw=="NaN",0,omowlomppw))  
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    oplppp <- outcome.podium*log(pred.prob.podium) ## deal with 0 probs
    omoplomppp <- (1-outcome.podium)*log(1-pred.prob.podium)
    log.score.podium.sim[tt+1] <- sum(ifelse(oplppp=="NaN",0,oplppp)+ifelse(omoplomppp=="NaN",0,omoplomppp))  
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ot10lppt10 <- outcome.top10*log(pred.prob.top10) ## deal with 0 probs
    omot10lomppt10 <- (1-outcome.top10)*log(1-pred.prob.top10)
    log.score.top10.sim[tt+1] <- sum(ifelse(ot10lppt10=="NaN",0,ot10lppt10)+ifelse(omot10lomppt10=="NaN",0,omot10lomppt10))  
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))

    ##outcome.win <- c(1,rep(0,length(these.drivers)-1))
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    ##outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    ##outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))
##################################
    
    ## compute log prior predictive
    log.pred.prob.all <- rep(0,its)
    for(i in 1:its) {
      log.pred.prob.all[i] <- loglikelihood.pl(x[tt+1,],res.PL$lambda[i,])
    }
    max.lppa <- max(log.pred.prob.all)
    log.prior.pred[tt+1] <- log(mean(exp(log.pred.prob.all-max.lppa)))+max.lppa
  }
  ## last race of the season
  tt = n
  print(tt)
  lambda.init <- rgamma(K,a,b)
  res.PL <- twtPL.gibbs(x[1:tt,],K,r.vec=r.vec[1:tt],tau.vec=tau.vec[1:tt],tau=tau,xi=xi,a=a,b=b,lambda.init=lambda.init,its=its,burn=burn,thin=thin)
  lambda.mean[tt,] <- apply(res.PL$lambda,2,mean)
  lambda.upper[tt,] <- apply(res.PL$lambda,2,quantile,0.975)
  lambda.lower[tt,] <- apply(res.PL$lambda,2,quantile,0.025)
  ##res.PL.pi <- res.PL$lambda/rowSums(res.PL$lambda)
  for(i in 1:K){
    for(j in 1:K){
      win.pairs.prob[tt,i,j] <- mean(res.PL$lambda[,i]/(res.PL$lambda[,i]+res.PL$lambda[,j]))
    }
  }
  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,ppw=pred.prob.winner,lm=lambda.mean,lu=lambda.upper,ll=lambda.lower,wpp=win.pairs.prob,pm=points.mean,cw=champ.winner,esw=expected.season.wins,esp=expected.season.podiums,est10=expected.season.top10s))
}


twtPL.varxi.mcmc <- function(x,K,r.vec,tau.vec,tau,a.xi=1,b.xi=1,a=rep(1,K),b=1,lambda.init=rep(1,K),xi.init=0.9999,sigma.xi=0.01,its=1000,burn=0,thin=1)
{
  ## MCMC sampler for time-weighted truncated Plackett-Luce model
  
  ## preliminary calculations
  ss <- summary.stats.twtpl(x,K,r.vec,tau.vec,tau,xi.init)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  p <- ss$p
  
  ## the rank at which truncation occurs
  r.star <- ss$r.star

  ## time-weighted number of races in which driver k finishes in top r.star
  w.init <- ss$w
  
  ## compute delta
  delta <- ss$delta

  ## compute time-weghting
  time.weight.init <- ss$time.weight

  ## storage for samples
  Y.res <- array(0,c(its,n,max(p)-1))
  lambda.res <- matrix(0,nrow=its,ncol=K)
  xi.res <- rep(0,its)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  lambda.curr <- lambda.init
  xi.curr <- xi.init
  w.curr <- w.init
  time.weight.curr <- time.weight.init
  Y.curr <- matrix(0,nrow=n,ncol=(max(p)-1))
  lprior.xi.curr <- dbeta(xi.curr,a.xi,b.xi,log=TRUE)
  llike.xi.curr <- loglike.twtpl(x,lambda.curr,r.vec,tau.vec,tau,xi.curr)
  theta.curr <- log(xi.curr)-log(1-xi.curr)
  ljacobian.xi.curr <- -theta.curr-2*log(1+exp(-theta.curr))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      print(t-burn)
      for(ell in 1:thin){
        ## sample latent variables | everything else
        llatent.xi.curr <- 0
        for(i in 1:n){
          ##for(j in 1:(p[i]-1)){
          for(j in 1:r.star[i]){
            Y.curr[i,j] <- rgamma(1,time.weight.curr[i],sum(lambda.curr[x[i,j:p[i]]]))
            llatent.xi.curr <- llatent.xi.curr+dgamma(Y.curr[i,j],time.weight.curr[i],sum(lambda.curr[x[i,j:p[i]]]),log=TRUE)
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          lambda.curr[k] <- rgamma(1,a[k]+w.curr[k],b+sum(delta[,1:r.star[i],k]*Y.curr[,1:r.star[i]]))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Lambda.curr <- rgamma(1,sum(a),b)
        lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)

        ##sample xi | everything else
        theta.prop <- rnorm(1,log(xi.curr)-log(1-xi.curr),sigma.xi)
        xi.prop <- 1/(1+exp(-theta.prop))
        w.prop  <- rep(0,K)
        time.weight.prop <- rep(0,n)
        for(k in 1:K){
          for(i in 1:n){
            time.weight.prop[i] <- psi(tau-tau.vec[i],xi.prop)
            for(j in 1:r.star[i]){
              w.prop[k] <- w.prop[k] + ifelse(x[i,j]==k,1,0)*time.weight.prop[i]
            }
          }
        }
        lprior.xi.prop <- dbeta(xi.prop,a.xi,b.xi,log=TRUE)
        llike.xi.prop <- loglike.twtpl(x,lambda.curr,r.vec,tau.vec,tau,xi.prop)
        llatent.xi.prop <- 0
        for(i in 1:n){
          for(j in 1:r.star[i]){
            llatent.xi.prop <- llatent.xi.prop+dgamma(Y.curr[i,j],time.weight.prop[i],sum(lambda.curr[x[i,j:p[i]]]),log=TRUE)
          }
        }
        ljacobian.xi.prop <- -theta.prop-2*log(1+exp(-theta.prop))

        l.acc <- lprior.xi.prop+llike.xi.prop+llatent.xi.prop+ljacobian.xi.prop
        l.acc <- l.acc-(lprior.xi.curr+llike.xi.curr+llatent.xi.curr+ljacobian.xi.curr)

        if(log(runif(1,0,1)) < min(0,l.acc)){

          print("accepted xi")
          xi.curr <- xi.prop
          w.curr <- w.prop
          time.weight.curr <- time.weight.prop
          lprior.xi.curr <- lprior.xi.prop
          llike.xi.curr <- llike.xi.prop
          llatent.xi.curr <- llatent.xi.prop
          ljacobian.xi.curr <- ljacobian.xi.prop
        }
      }
    }
  }
  for(t in 1:its){
    print(t)
    for(ell in 1:thin){
      ## sample latent variables | everything else
      llatent.xi.curr <- 0
      for(i in 1:n){
        ##for(j in 1:(p[i]-1)){
        for(j in 1:r.star[i]){
          Y.curr[i,j] <- rgamma(1,time.weight.curr[i],sum(lambda.curr[x[i,j:p[i]]]))
          llatent.xi.curr <- llatent.xi.curr+dgamma(Y.curr[i,j],time.weight.curr[i],sum(lambda.curr[x[i,j:p[i]]]),log=TRUE)
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        lambda.curr[k] <- rgamma(1,a[k]+w.curr[k],b+sum(delta[,1:r.star[i],k]*Y.curr[,1:r.star[i]]))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Lambda.curr <- rgamma(1,sum(a),b)
      lambda.curr <- Lambda.curr*lambda.curr/sum(lambda.curr)

      ##sample xi | everything else
      theta.prop <- rnorm(1,log(xi.curr)-log(1-xi.curr),sigma.xi)
      xi.prop <- 1/(1+exp(-theta.prop))
      w.prop  <- rep(0,K)
      time.weight.prop <- rep(0,n)
      for(k in 1:K){
        for(i in 1:n){
          time.weight.prop[i] <- psi(tau-tau.vec[i],xi.prop)
          for(j in 1:r.star[i]){
            w.prop[k] <- w.prop[k] + ifelse(x[i,j]==k,1,0)*time.weight.prop[i]
          }
        }
      }
      lprior.xi.prop <- dbeta(xi.prop,a.xi,b.xi,log=TRUE)
      llike.xi.prop <- loglike.twtpl(x,lambda.curr,r.vec,tau.vec,tau,xi.prop)
      llatent.xi.prop <- 0
      for(i in 1:n){
        for(j in 1:r.star[i]){
          llatent.xi.prop <- llatent.xi.prop+dgamma(Y.curr[i,j],time.weight.prop[i],sum(lambda.curr[x[i,j:p[i]]]),log=TRUE)
        }
      }
      ljacobian.xi.prop <- -theta.prop-2*log(1+exp(-theta.prop))
      
      l.acc <- lprior.xi.prop+llike.xi.prop+llatent.xi.prop+ljacobian.xi.prop
      l.acc <- l.acc-(lprior.xi.curr+llike.xi.curr+llatent.xi.curr+ljacobian.xi.curr)
      
      if(log(runif(1,0,1)) < min(0,l.acc)){
        
        print("accepted xi")
        xi.curr <- xi.prop
        w.curr <- w.prop
        time.weight.curr <- time.weight.prop
        lprior.xi.curr <- lprior.xi.prop
        llike.xi.curr <- llike.xi.prop
        llatent.xi.curr <- llatent.xi.prop
        ljacobian.xi.curr <- ljacobian.xi.prop
      }
    }
    
    ## store sampled values
    lambda.res[t,] <- lambda.curr
    xi.res[t] <- xi.curr
    Y.res[t,,] <- Y.curr
    
    ## compute observed data time-weighted loglikelihood
    ll.res[t,] <- loglike.twtpl(x,lambda.res[t,],r.vec,tau.vec,tau,xi.res[t])
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(lambda.res[t,],a,b,log=TRUE))+dbeta(xi.res[t],a.xi,b.xi,log=TRUE)
    print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(Y=Y.res,lambda=lambda.res,xi=xi.res,ll=ll.res,ljd=ljd.res))
}

############################################################################

rattrition <- function(gm=rep(1,3))
{
  ## generate permuation from attrition model for K individuals with
  ## ability parameters gm(="gamma")

  ## use latent variables to generate the data
  K <- length(gm)
  V <- rexp(K,gm)
  x <- rev(order(V))
  x
}

likelihood.attrition <- function(x,gm)
{
  ## function for computing attrition model likelihood
  ## of gamma(=gm) given single ranking x
  
  ## calculate length of vector x (this is number of drivers in race)
  n <- length(x[is.na(x)==FALSE])
  ## initialise likelihood
  like <- 1
  ## compute likelihood as product of probabilities
  for(j in 2:n){
    like <- like*gm[x[j]]/sum(gm[x[1:j]])
  }
  return(like)
}

loglikelihood.attrition <- function(x,gm)
{
  ## function for computing attrition model log likelihood
  ## of gamma(=gm) given single ranking x
  
  ## calculate length of vector x (this is number of drivers in race)
  n <- length(x[is.na(x)==FALSE])
  ## initialise likelihood
  llike <- 0
  ## compute loglikelihood as sum of log probabilities
  for(j in 2:n){
    llike <- llike + log(gm[x[j]])-log(sum(gm[x[1:j]]))
  }
  return(llike)
}

loglike.attrition <- function(x,gm)
{
  ## attrition model loglikelihood over multiple rankings
  
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }

  ll.vec <- rep(0,n)
  for(i in 1:n){
    n.i <- length(x[i,is.na(x[i,])==FALSE])
    for(j in 2:n.i){
      ll.vec[i] <- ll.vec[i]+log(gm[x[i,j]])-log(sum(gm[x[i,1:j]]))
    }
  }
  return(ll.vec)
}

summary.stats.attrition <- function(x,K)
{
  ## compute summary statistics for a given set of rankings data
  ## for the attrition model
  
  ## number of multiple comparisons
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  n.i <- rep(0,n)
  for(i in 1:n){
    n.i[i] <- length(x[i,is.na(x[i,])==FALSE])
  }
  
  ## identify individual in first place
  ## and count up times for each individual
  v  <- rep(0,K)
  for(k in 1:K){
    for(i in 1:n){
      for(j in 2:n.i[i]){
        v[k] <- v[k] + ifelse(x[i,j]==k,1,0)
      }
    }
  }
  
  ## compute beta
  beta <- array(0,c(n,max(n.i)-1,K))
  for(k in 1:K){
    for(i in 1:n){
      for(j in 1:(n.i[i]-1)){
        beta[i,j,k] <- sum(k==x[i,1:(j+1)])
      }
    }
  }

  return(list(n=n,n.i=n.i,v=v,beta=beta))
}

attrition.gibbs <- function(x,K,cc=rep(1,K),d=1,gamma.init=rep(1,K),its=1000,burn=0,thin=1)
{
  ## Gibbs sampler for attrition model
  
  ## preliminary calculations
  ss <- summary.stats.attrition(x,K)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  n.i <- ss$n.i
  
  ## identify individual in last 
  ## and count up times for each individual
  v <- ss$v
  
  ## compute beta
  beta <- ss$beta

  ## storage for samples
  U.res <- array(0,c(its,n,max(n.i)-1))
  gamma.res <- matrix(0,nrow=its,ncol=K)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  gamma.curr <- gamma.init
  U.curr <- matrix(0,nrow=n,ncol=(max(n.i)-1))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      ##print(t-burn)
      for(l in 1:thin){
        ## sample latent variables | everything else
        for(i in 1:n){
          for(j in 1:(n.i[i]-1)){
            U.curr[i,j] <- rexp(1,sum(gamma.curr[x[i,1:(j+1)]]))
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          gamma.curr[k] <- rgamma(1,cc[k]+v[k],d+sum(beta[,,k]*U.curr))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Gamma.curr <- rgamma(1,sum(cc),d)
        gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)
      }
    }
  }
  for(t in 1:its){
    ##print(t)
    for(l in 1:thin){
      ## sample latent variables | everything else
      for(i in 1:n){
        for(j in 1:(n.i[i]-1)){
          U.curr[i,j] <- rexp(1,sum(gamma.curr[x[i,1:(j+1)]]))
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        gamma.curr[k] <- rgamma(1,cc[k]+v[k],d+sum(beta[,,k]*U.curr))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Gamma.curr <- rgamma(1,sum(cc),d)
      gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)
    }
    
    ## store sampled values
    gamma.res[t,] <- gamma.curr
    U.res[t,,] <- U.curr
    
    ## compute observed data loglikelihood
    ll.res[t,] <- loglike.attrition(x,gamma.res[t,])
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(gamma.res[t,],cc,d,log=TRUE))
    ##print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(U=U.res,gm=gamma.res,ll=ll.res,ljd=ljd.res))
}

chib.ml.attrition <- function(res,cc,d,ss)
{
  ## Chib's estimate of the log marginal likelihood
  ## for the attrition model

  ## number of iterations
  its <- dim(res$gm)[1]

  ## choose a gamma with high joint density
  ##this.val <- which.max(res$ljd)
  ##gamma.max <- res$gm[this.val,]
  gamma.max <- apply(res$gm,2,mean)

  ## compute log of full conditionals for gamma.max (Rao-Blackwellized)
  log.cd <- rep(0,its)
  for(i in 1:its){
    sum.ldg <- 0
    for(k in 1:K){                                    
      sum.ldg <- sum.ldg+dgamma(gamma.max[k],cc[k]+ss$v[k],d+sum(ss$beta[,,k]*res$U[i,,]),log=TRUE)
    }                    
    log.cd[i] <- sum.ldg
  }
  log.post <- log(mean(exp(log.cd-max(log.cd))))+max(log.cd)  
  ##log.ml.chib <- res$ljd[this.val]-log.post
  log.ml.chib <- sum(loglike.attrition(x,gamma.max))+sum(dgamma(gamma.max,cc,d,log=TRUE))-log.post
  log.ml.chib
}

MC.ml.attrition <- function(x,K,N,cc,d)
{
  ## Monte Carlo estimate of marginal likelihood uder the attrition model

  log.ml.mc <- rep(0,N)
  for(i in 1:N){
    log.ml.mc[i] <- sum(loglike.attrition(x,rgamma(K,cc,d)))
  }
  log.ml.MC <- log(mean(exp(log.ml.mc-max(log.ml.mc))))+max(log.ml.mc)
  log.ml.MC
}

attrition.sequential.old <- function(x=x,K=K,cc=cc,d=d,its=1000,burn=100,thin=1,zero.adjust=TRUE)
{
  ## sequential updating and prediction under the attrition model

  n <- dim(x)[1]
  ##K <- max(x[is.na(x)==FALSE])
  
  points.actual <- matrix(0,nrow=n,ncol=K)
  for(tt in 1:n){
    these.drivers <- x[tt,is.na(x[tt,])==FALSE]
    points.actual[tt,x[tt,is.na(x[tt,])==FALSE]] <- points.calc.2013(1:length(these.drivers))
  }
  
  ## set up matrices and vectors to store results
  log.prior.pred <- rep(0,n)
  log.score.winner.sim <- rep(0,n)
  log.score.podium.sim <- rep(0,n)
  log.score.top10.sim <- rep(0,n)
  points.mean <- matrix(0,nrow=n,ncol=K)
  gamma.mean <- matrix(0,nrow=n,ncol=K)
  gamma.upper <- matrix(0,nrow=n,ncol=K)
  gamma.lower <- matrix(0,nrow=n,ncol=K)
  ## array for storing probabilities of i beating j
  win.pairs.prob <- array(0,c(n,K,K))
  ## array for storing predictive probabilities for winner, podium, points
  pred.probs <- array(0,c(n,K,3))
  ## matrix to store sampled championship winners over time
  champ.winner <- matrix(0,nrow=n,ncol=its)
  
  ## initial values 
  
  ## log prior predictive for first race (verified by simulation!)
  log.prior.pred[1] <- log(1/factorial(length(x[1,is.na(x[1,])==FALSE])))
  ##MC.ml.attrition(x[1,],K,10000,a,b)
  ## log score for "winner" of first race
  log.score.winner.sim[1] <- log(1/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-1)*log(1-(1/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "podium finish" of first race
  log.score.podium.sim[1] <- 3*log(3/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-3)*log(1-(3/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "points finish (top 10)" of first race
  log.score.top10.sim[1] <- 10*log(10/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-10)*log(1-(10/length(x[1,is.na(x[1,])==FALSE])))
  ## predicted championship winners after 1st race
  champ.winner[1,] <- sample(x[1,is.na(x[1,])==FALSE],its,replace=TRUE)
  ## predictive probabilities for winner, podium, points for first race
  these.drivers <- x[1,is.na(x[1,])==FALSE]
  for(ell in 1:length(these.drivers)){
    pred.probs[1,these.drivers[ell],1] <- 1/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],2] <- 3/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],3] <- 10/length(x[1,is.na(x[1,])==FALSE])
    points.mean[1,these.drivers[ell]] <- sum(points.calc.2013(1:10))/length(x[1,is.na(x[1,])==FALSE])
  }
  
  ## analyse one race at a time
  for(tt in 1:(n-1)) {
    print(tt)
    gamma.init <- rgamma(K,cc,d)
    res.attrition <- attrition.gibbs(x[1:tt,],K,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
    
    ## predictive simulations for remaining races - assuming drivers
    ## from race tt+1
    these.drivers <- x[tt+1,is.na(x[tt+1,])==FALSE]
    x.pred.sim <- array(0,c(its,n-tt,length(these.drivers)))
    points.sim <- array(0,c(its,n-tt,K))
    ##pos.sim <- array(0,c(its,n-tt,K))
    pred.final.points <- matrix(0,nrow=its,ncol=K)
    
    ## for each sampled set of parameter values - future race results
    ## convert to points for each driver
    ## convert to finishing position for each driver (NOT USED)
    for(i in 1:its){
      for(k in 1:(n-tt)){
        x.pred.sim[i,k,] <- these.drivers[rattrition(res.attrition$gm[i,x[tt+1,is.na(x[tt+1,])==FALSE]])]
        points.sim[i,k,x.pred.sim[i,k,]] <- points.calc.2013(1:length(these.drivers))
        ##pos.sim[i,k,x.pred.sim[i,k,]] <- 1:length(these.drivers)
      }
      ## predicted end of season points for each driver
      if(tt ==1){
        pred.final.points[i,] <- points.actual[1,]+apply(points.sim[i,,],2,sum)
      }
      else{
        if(tt == (n-1)){
          pred.final.points[i,] <- apply(points.actual[1:tt,],2,sum)+points.sim[i,,]
        }
        else{
          pred.final.points[i,] <- apply(points.actual[1:tt,],2,sum)+apply(points.sim[i,,],2,sum)
        }
      }
    }
    ## simulated championship winner
    champ.winner[tt+1,] <-  apply(pred.final.points,1,which.max)
    
    ## expected number of points in next race
    points.mean[tt+1,] <- apply(points.sim[,1,],2,mean)
    
    ## store parameter summaries
    gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
    gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
    gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
    ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
    for(i in 1:K){
      for(j in 1:K){
        win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
      }
    }
    
    ## log scores for winner,podium, top10 etc ##########
    pred.prob.win <- NULL
    pred.prob.podium <- NULL
    pred.prob.top10 <- NULL
    ## deal with zero probs - like Bayes estimate of binomial proportion
    ## with beta prior
    a.zero <- 0
    b.zero <- 0
    if(zero.adjust){
      a.zero <- 1
      b.zero <- length(these.drivers)-1
    }
    new.denom <- a.zero+b.zero+its
    for(j in 1:length(these.drivers)){
      pred.prob.win[j] <- mean(x.pred.sim[,1,1]==these.drivers[j])
      pred.prob.podium[j] <- mean(apply(x.pred.sim[,1,1:3]==these.drivers[j],1,sum))
      pred.prob.top10[j] <- mean(apply(x.pred.sim[,1,1:10]==these.drivers[j],1,sum))
      pred.prob.win[j] <- (a.zero+(pred.prob.win[j]*its))/new.denom
      pred.prob.podium[j] <- (a.zero+(pred.prob.podium[j]*its))/new.denom
      pred.prob.top10[j] <- (a.zero+(pred.prob.top10[j]*its))/new.denom
      pred.probs[tt+1,these.drivers[j],1] <- pred.prob.win[j]
      pred.probs[tt+1,these.drivers[j],2] <- pred.prob.podium[j]
      pred.probs[tt+1,these.drivers[j],3] <- pred.prob.top10[j]    
    }

    outcome.win <- c(1,rep(0,length(these.drivers)-1))
    owlppw <- outcome.win*log(pred.prob.win) ## deal with 0 probs
    omowlomppw <- (1-outcome.win)*log(1-pred.prob.win)
    log.score.winner.sim[tt+1] <- sum(ifelse(owlppw=="NaN",0,owlppw)+ifelse(omowlomppw=="NaN",0,omowlomppw))  
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    oplppp <- outcome.podium*log(pred.prob.podium) ## deal with 0 probs
    omoplomppp <- (1-outcome.podium)*log(1-pred.prob.podium)
    log.score.podium.sim[tt+1] <- sum(ifelse(oplppp=="NaN",0,oplppp)+ifelse(omoplomppp=="NaN",0,omoplomppp))  
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ot10lppt10 <- outcome.top10*log(pred.prob.top10) ## deal with 0 probs
    omot10lomppt10 <- (1-outcome.top10)*log(1-pred.prob.top10)
    log.score.top10.sim[tt+1] <- sum(ifelse(ot10lppt10=="NaN",0,ot10lppt10)+ifelse(omot10lomppt10=="NaN",0,omot10lomppt10))  
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))
##################################
    
    ## compute log prior predictive
    log.pred.prob.all <- rep(0,its)
    for(i in 1:its) {
      log.pred.prob.all[i] <- loglikelihood.attrition(x[tt+1,],res.attrition$gm[i,])
    }
    max.lppa <- max(log.pred.prob.all)
    log.prior.pred[tt+1] <- log(mean(exp(log.pred.prob.all-max.lppa)))+max.lppa
  }
  ## last race of the season
  tt = n
  print(tt)
  gamma.init <- rgamma(K,cc,d)
  res.attrition <- attrition.gibbs(x[1:tt,],K,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
  gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
  gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
  gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
  ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
  for(i in 1:K){
    for(j in 1:K){
      win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
    }
  }

  ##return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs))
  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,gm=gamma.mean,gu=gamma.upper,gl=gamma.lower,wpp=win.pairs.prob,pm=points.mean,cw=champ.winner))
}

attrition.sequential <- function(x=x,K=K,tau.vec=tau.vec,cc=cc,d=d,its=1000,burn=100,thin=1,zero.adjust=TRUE)
{
  ## sequential updating and prediction under the attrition model

  n <- dim(x)[1]
  ##K <- max(x[is.na(x)==FALSE])
  
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
  
  ## set up matrices and vectors to store results
  log.prior.pred <- rep(0,n)
  log.score.winner.sim <- rep(0,n)
  log.score.podium.sim <- rep(0,n)
  log.score.top10.sim <- rep(0,n)
  points.mean <- matrix(0,nrow=n,ncol=K)
  gamma.mean <- matrix(0,nrow=n,ncol=K)
  gamma.upper <- matrix(0,nrow=n,ncol=K)
  gamma.lower <- matrix(0,nrow=n,ncol=K)
  Gamma.mean <- rep(0,nrow=n,ncol=K)
  Gamma.upper <- matrix(0,nrow=n,ncol=K)
  Gamma.lower <- matrix(0,nrow=n,ncol=K)
  ## array for storing probabilities of i beating j
  win.pairs.prob <- array(0,c(n,K,K))
  ## array for storing predictive probabilities for winner, podium, points
  pred.probs <- array(0,c(n,K,3))
  ## matrix to store sampled championship winners over time
  champ.winner <- matrix(0,nrow=n,ncol=its)
  ## matrix to store expected number of wins in a season for each driver
  expected.season.wins <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of podiums in a season for each driver
  expected.season.podiums <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of top 10s in a season for each driver
  expected.season.top10s <- matrix(0,nrow=n,ncol=K)
  
  ## initial values 
  
  ## log prior predictive for first race (verified by simulation!)
  log.prior.pred[1] <- log(1/factorial(length(x[1,is.na(x[1,])==FALSE])))
  ##MC.ml.attrition(x[1,],K,10000,a,b)
  ## log score for "winner" of first race
  log.score.winner.sim[1] <- log(1/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-1)*log(1-(1/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "podium finish" of first race
  log.score.podium.sim[1] <- 3*log(3/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-3)*log(1-(3/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "points finish (top 10)" of first race
  log.score.top10.sim[1] <- 10*log(10/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-10)*log(1-(10/length(x[1,is.na(x[1,])==FALSE])))
  ## predicted championship winners after 1st race
  champ.winner[1,] <- sample(x[1,is.na(x[1,])==FALSE],its,replace=TRUE)
  ## Expected season wins, podiums, top10s
  season.races <- grep(strsplit(as.character(tau.vec[1]),"-")[[1]][1],tau.vec)
  length.season <- length(season.races)
  expected.season.wins[1,x[1,is.na(x[1,])==FALSE]] <- length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.podiums[1,x[1,is.na(x[1,])==FALSE]] <- 3*length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.top10s[1,x[1,is.na(x[1,])==FALSE]] <- 10*length.season/length(x[1,is.na(x[1,])==FALSE])
  ## predictive probabilities for winner, podium, points for first race
  these.drivers <- x[1,is.na(x[1,])==FALSE]
  for(ell in 1:length(these.drivers)){
    pred.probs[1,these.drivers[ell],1] <- 1/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],2] <- 3/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],3] <- 10/length(x[1,is.na(x[1,])==FALSE])
    points.mean[1,these.drivers[ell]] <- sum(points.calc.2013(1:10))/length(x[1,is.na(x[1,])==FALSE])
  }
  
  ## analyse one race at a time
  for(tt in 1:(n-1)) {
    print(tt)
    gamma.init <- rgamma(K,cc,d)
    res.attrition <- attrition.gibbs(x[1:tt,],K,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
    
    ## predictive simulations for remaining races - assuming drivers
    ## from race tt+1
    these.drivers <- x[tt+1,is.na(x[tt+1,])==FALSE]
    season.races <- grep(strsplit(as.character(tau.vec[tt]),"-")[[1]][1],tau.vec)
    length.season <- length(season.races)
    remaining.races <- season.races[season.races>tt]
    past.races <- season.races[season.races<=tt]
    num.past.races <- length(past.races)
    num.remain.races <- length(remaining.races)
    if(num.remain.races==0){
      remaining.races <- grep(strsplit(as.character(tau.vec[tt+1]),"-")[[1]][1],tau.vec)
      num.remain.races <- length(remaining.races)
    }
    ##print(c(length.season,num.past.races,num.remain.races))
    ##print(season.races)
    ##print(remaining.races)
    ##print(past.races)
    
    x.pred.sim <- array(0,c(its,num.remain.races,length(these.drivers)))
    points.sim <- array(0,c(its,num.remain.races,K))
    is.winner.sim <- array(0,c(its,num.remain.races,K))
    is.podium.sim <- array(0,c(its,num.remain.races,K))
    is.top10.sim <- array(0,c(its,num.remain.races,K))
    ##pos.sim <- array(0,c(its,n-tt,K))
    pred.final.points <- matrix(0,nrow=its,ncol=K)
    pred.final.wins <- matrix(0,nrow=its,ncol=K)
    pred.final.podiums <- matrix(0,nrow=its,ncol=K)
    pred.final.top10s <- matrix(0,nrow=its,ncol=K)
    
    ## for each sampled set of parameter values - future race results
    ## (in season)
    ## convert to points for each driver convert to
    ## finishing position for each driver (NOT USED)
    for(i in 1:its){
      for(k in 1:num.remain.races){
        x.pred.sim[i,k,] <- these.drivers[rattrition(res.attrition$gm[i,x[tt+1,is.na(x[tt+1,])==FALSE]])]
        points.sim[i,k,x.pred.sim[i,k,]] <- points.calc.2013(1:length(these.drivers))
        is.winner.sim[i,k,x.pred.sim[i,k,1]] <- 1
        is.podium.sim[i,k,x.pred.sim[i,k,1:3]] <- 1
        is.top10.sim[i,k,x.pred.sim[i,k,1:10]] <- 1
        ##pos.sim[i,k,x.pred.sim[i,k,]] <- 1:length(these.drivers)
      }
      ## predicted end of season points etc for each driver
      if(num.past.races==1){
        pred.final.points[i,] <- points.actual[past.races,]+apply(points.sim[i,,],2,sum)
        pred.final.wins[i,] <- wins.actual[past.races,]+apply(is.winner.sim[i,,],2,sum)
        pred.final.podiums[i,] <- podiums.actual[past.races,]+apply(is.podium.sim[i,,],2,sum)
        pred.final.top10s[i,] <- top10s.actual[past.races,]+apply(is.top10.sim[i,,],2,sum)
      }
      else{
        if(num.past.races == (length.season-1)){
          pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+points.sim[i,,]
          pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+is.winner.sim[i,,]
          pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+is.podium.sim[i,,]
          pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+is.top10.sim[i,,]
        }
        else{
          if(num.past.races == length.season){
            pred.final.points[i,] <- apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(is.top10.sim[i,,],2,sum)
          }
          else{
            pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+apply(is.top10.sim[i,,],2,sum)
          }
        }
      }
    }
    ## simulated championship winner
    champ.winner[tt+1,] <-  apply(pred.final.points,1,which.max)
    expected.season.wins[tt+1,] <- apply(pred.final.wins,2,mean)
    expected.season.podiums[tt+1,] <- apply(pred.final.podiums,2,mean)
    expected.season.top10s[tt+1,] <- apply(pred.final.top10s,2,mean)
    
    ## expected number of points in next race
    points.mean[tt+1,] <- apply(points.sim[,1,],2,mean)
    
    ## store parameter summaries
    gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
    gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
    gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
    ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
    for(i in 1:K){
      for(j in 1:K){
        win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
      }
    }
    
    ## log scores for winner,podium, top10 etc ##########
    pred.prob.win <- NULL
    pred.prob.podium <- NULL
    pred.prob.top10 <- NULL
    ## deal with zero probs - like Bayes estimate of binomial proportion
    ## with beta prior
    a.zero <- 0
    b.zero <- 0
    if(zero.adjust){
      a.zero <- 1
      b.zero <- length(these.drivers)-1
    }
    new.denom <- a.zero+b.zero+its
    for(j in 1:length(these.drivers)){
      pred.prob.win[j] <- mean(x.pred.sim[,1,1]==these.drivers[j])
      pred.prob.podium[j] <- mean(apply(x.pred.sim[,1,1:3]==these.drivers[j],1,sum))
      pred.prob.top10[j] <- mean(apply(x.pred.sim[,1,1:10]==these.drivers[j],1,sum))
      pred.prob.win[j] <- (a.zero+(pred.prob.win[j]*its))/new.denom
      pred.prob.podium[j] <- (a.zero+(pred.prob.podium[j]*its))/new.denom
      pred.prob.top10[j] <- (a.zero+(pred.prob.top10[j]*its))/new.denom
      pred.probs[tt+1,these.drivers[j],1] <- pred.prob.win[j]
      pred.probs[tt+1,these.drivers[j],2] <- pred.prob.podium[j]
      pred.probs[tt+1,these.drivers[j],3] <- pred.prob.top10[j]    
    }

    outcome.win <- c(1,rep(0,length(these.drivers)-1))
    owlppw <- outcome.win*log(pred.prob.win) ## deal with 0 probs
    omowlomppw <- (1-outcome.win)*log(1-pred.prob.win)
    log.score.winner.sim[tt+1] <- sum(ifelse(owlppw=="NaN",0,owlppw)+ifelse(omowlomppw=="NaN",0,omowlomppw))  
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    oplppp <- outcome.podium*log(pred.prob.podium) ## deal with 0 probs
    omoplomppp <- (1-outcome.podium)*log(1-pred.prob.podium)
    log.score.podium.sim[tt+1] <- sum(ifelse(oplppp=="NaN",0,oplppp)+ifelse(omoplomppp=="NaN",0,omoplomppp))  
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ot10lppt10 <- outcome.top10*log(pred.prob.top10) ## deal with 0 probs
    omot10lomppt10 <- (1-outcome.top10)*log(1-pred.prob.top10)
    log.score.top10.sim[tt+1] <- sum(ifelse(ot10lppt10=="NaN",0,ot10lppt10)+ifelse(omot10lomppt10=="NaN",0,omot10lomppt10))  
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))
##################################
    
    ## compute log prior predictive
    log.pred.prob.all <- rep(0,its)
    for(i in 1:its) {
      log.pred.prob.all[i] <- loglikelihood.attrition(x[tt+1,],res.attrition$gm[i,])
    }
    max.lppa <- max(log.pred.prob.all)
    log.prior.pred[tt+1] <- log(mean(exp(log.pred.prob.all-max.lppa)))+max.lppa
  }
  ## last race of the season
  tt = n
  print(tt)
  gamma.init <- rgamma(K,cc,d)
  res.attrition <- attrition.gibbs(x[1:tt,],K,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
  gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
  gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
  gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
  ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
  for(i in 1:K){
    for(j in 1:K){
      win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
    }
  }

  ##return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs))
  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,gm=gamma.mean,gu=gamma.upper,gl=gamma.lower,wpp=win.pairs.prob,pm=points.mean,cw=champ.winner,esw=expected.season.wins,esp=expected.season.podiums,est10=expected.season.top10s))
}

summary.stats.twattrition <- function(x,K,tau.vec,tau,xi)
{
  ## compute summary statistics for a given set of rankings data
  ## for the time-weighted attrition model
  
  ## number of multiple comparisons
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  n.i <- rep(0,n)
  for(i in 1:n){
    n.i[i] <- length(x[i,is.na(x[i,])==FALSE])
  }
  
  ## compute time-weighted number of races in which driver k doesn't
  ## finish in first place
  v  <- rep(0,K)
  time.weight <- rep(0,n)
  for(k in 1:K){
    for(i in 1:n){
      time.weight[i] <- psi(tau-tau.vec[i],xi)
      for(j in 2:n.i[i]){
        v[k] <- v[k] + ifelse(x[i,j]==k,1,0)*time.weight[i]
      }
    }
  }
  
  ## compute beta
  beta <- array(0,c(n,max(n.i)-1,K))
  for(k in 1:K){
    for(i in 1:n){
      for(j in 1:(n.i[i]-1)){
        beta[i,j,k] <- sum(k==x[i,1:(j+1)])
      }
    }
  }

  return(list(n=n,n.i=n.i,v=v,beta=beta,time.weight=time.weight))
}

loglike.twattrition <- function(x,gm,tau.vec,tau,xi)
{
  ## time-weighted attrition model loglikelihood over multiple rankings
  
  n <- max(dim(x)[1],1)
  if(n == 1){
    x <- t(matrix(x))
  }

  ll.vec <- rep(0,n)
  time.weight <- rep(0,n)
  for(i in 1:n){
    time.weight[i] <- psi(tau-tau.vec[i],xi)
    n.i <- length(x[i,is.na(x[i,])==FALSE])
    for(j in 2:n.i){
      ll.vec[i] <- ll.vec[i]+time.weight[i]*(log(gm[x[i,j]])-log(sum(gm[x[i,1:j]])))
    }
  }
  return(ll.vec)
}

twattrition.gibbs <- function(x,K,tau.vec,tau,xi=1,cc=rep(1,K),d=1,gamma.init=rep(1,K),its=1000,burn=0,thin=1)
{
  ## Gibbs sampler for time-weighted attrition model
  ## Started: 30/6/14
  ## Updated: 30/6/14
  
  ## preliminary calculations
  ss <- summary.stats.twattrition(x,K,tau.vec,tau,xi)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  n.i <- ss$n.i
  
  ## time-weighted number of races in which driver did not win
  v <- ss$v
  
  ## compute beta
  beta <- ss$beta

  ## compute time-weghting
  time.weight <- ss$time.weight

  ## storage for samples
  U.res <- array(0,c(its,n,max(n.i)-1))
  gamma.res <- matrix(0,nrow=its,ncol=K)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  gamma.curr <- gamma.init
  U.curr <- matrix(0,nrow=n,ncol=(max(n.i)-1))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      ##print(t-burn)
      for(l in 1:thin){
        ## sample latent variables | everything else
        for(i in 1:n){
          for(j in 1:(n.i[i]-1)){
            U.curr[i,j] <- rgamma(1,time.weight[i],sum(gamma.curr[x[i,1:(j+1)]]))
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          gamma.curr[k] <- rgamma(1,cc[k]+v[k],d+sum(beta[,,k]*U.curr))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Gamma.curr <- rgamma(1,sum(cc),d)
        gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)
      }
    }
  }
  for(t in 1:its){
    ##print(t)
    for(l in 1:thin){
      ## sample latent variables | everything else
      for(i in 1:n){
        for(j in 1:(n.i[i]-1)){
          U.curr[i,j] <- rgamma(1,time.weight[i],sum(gamma.curr[x[i,1:(j+1)]]))
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        gamma.curr[k] <- rgamma(1,cc[k]+v[k],d+sum(beta[,,k]*U.curr))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Gamma.curr <- rgamma(1,sum(cc),d)
      gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)
    }
    
    ## store sampled values
    gamma.res[t,] <- gamma.curr
    U.res[t,,] <- U.curr
    
    ## compute observed data loglikelihood
    ll.res[t,] <- loglike.twattrition(x,gamma.res[t,],tau.vec,tau,xi)
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(gamma.res[t,],cc,d,log=TRUE))
    ##print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(U=U.res,gm=gamma.res,ll=ll.res,ljd=ljd.res))
}

twattrition.sequential <- function(x=x,K=K,tau.vec=tau.vec,tau=tau,xi=xi,cc=cc,d=d,its=1000,burn=100,thin=1,zero.adjust=TRUE)
{
  ## sequential updating and prediction under the time-weighted attrition model

  n <- dim(x)[1]
  ##K <- max(x[is.na(x)==FALSE])
  
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
  
  ## set up matrices and vectors to store results
  log.prior.pred <- rep(0,n)
  log.score.winner.sim <- rep(0,n)
  log.score.podium.sim <- rep(0,n)
  log.score.top10.sim <- rep(0,n)
  points.mean <- matrix(0,nrow=n,ncol=K)
  gamma.mean <- matrix(0,nrow=n,ncol=K)
  gamma.upper <- matrix(0,nrow=n,ncol=K)
  gamma.lower <- matrix(0,nrow=n,ncol=K)
  ## array for storing probabilities of i beating j
  win.pairs.prob <- array(0,c(n,K,K))
  ## array for storing predictive probabilities for winner, podium, points
  pred.probs <- array(0,c(n,K,3))
  ## matrix to store sampled championship winners over time
  champ.winner <- matrix(0,nrow=n,ncol=its)
  ## matrix to store expected number of wins in a season for each driver
  expected.season.wins <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of podiums in a season for each driver
  expected.season.podiums <- matrix(0,nrow=n,ncol=K)
  ## matrix to store expected number of top 10s in a season for each driver
  expected.season.top10s <- matrix(0,nrow=n,ncol=K)
  
  ## initial values 
  
  ## log prior predictive for first race (verified by simulation!)
  log.prior.pred[1] <- log(1/factorial(length(x[1,is.na(x[1,])==FALSE])))
  ##MC.ml.attrition(x[1,],K,10000,a,b)
  ## log score for "winner" of first race
  log.score.winner.sim[1] <- log(1/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-1)*log(1-(1/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "podium finish" of first race
  log.score.podium.sim[1] <- 3*log(3/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-3)*log(1-(3/length(x[1,is.na(x[1,])==FALSE])))
  ## log score for "points finish (top 10)" of first race
  log.score.top10.sim[1] <- 10*log(10/length(x[1,is.na(x[1,])==FALSE]))+(length(x[1,is.na(x[1,])==FALSE])-10)*log(1-(10/length(x[1,is.na(x[1,])==FALSE])))
  ## predicted championship winners after 1st race
  champ.winner[1,] <- sample(x[1,is.na(x[1,])==FALSE],its,replace=TRUE)
  ## Expected season wins, podiums, top10s
  season.races <- grep(strsplit(as.character(tau.vec[1]),"-")[[1]][1],tau.vec)
  length.season <- length(season.races)
  expected.season.wins[1,x[1,is.na(x[1,])==FALSE]] <- length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.podiums[1,x[1,is.na(x[1,])==FALSE]] <- 3*length.season/length(x[1,is.na(x[1,])==FALSE])
  expected.season.top10s[1,x[1,is.na(x[1,])==FALSE]] <- 10*length.season/length(x[1,is.na(x[1,])==FALSE])
  ## predictive probabilities for winner, podium, points for first race
  these.drivers <- x[1,is.na(x[1,])==FALSE]
  for(ell in 1:length(these.drivers)){
    pred.probs[1,these.drivers[ell],1] <- 1/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],2] <- 3/length(x[1,is.na(x[1,])==FALSE])
    pred.probs[1,these.drivers[ell],3] <- 10/length(x[1,is.na(x[1,])==FALSE])
    points.mean[1,these.drivers[ell]] <- sum(points.calc.2013(1:10))/length(x[1,is.na(x[1,])==FALSE])
  }
  
  ## analyse one race at a time
  for(tt in 1:(n-1)) {
    print(tt)
    gamma.init <- rgamma(K,cc,d)
    res.attrition <- twattrition.gibbs(x[1:tt,],K,tau.vec=tau.vec[1:tt],tau=tau.vec[tt+1],xi=xi,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
    
    ## predictive simulations for remaining races - assuming drivers
    ## from race tt+1
    these.drivers <- x[tt+1,is.na(x[tt+1,])==FALSE]
    season.races <- grep(strsplit(as.character(tau.vec[tt]),"-")[[1]][1],tau.vec)
    length.season <- length(season.races)
    remaining.races <- season.races[season.races>tt]
    past.races <- season.races[season.races<=tt]
    num.past.races <- length(past.races)
    num.remain.races <- length(remaining.races)
    if(num.remain.races==0){
      remaining.races <- grep(strsplit(as.character(tau.vec[tt+1]),"-")[[1]][1],tau.vec)
      num.remain.races <- length(remaining.races)
    }
    x.pred.sim <- array(0,c(its,num.remain.races,length(these.drivers)))
    points.sim <- array(0,c(its,num.remain.races,K))
    is.winner.sim <- array(0,c(its,num.remain.races,K))
    is.podium.sim <- array(0,c(its,num.remain.races,K))
    is.top10.sim <- array(0,c(its,num.remain.races,K))
    ##pos.sim <- array(0,c(its,n-tt,K))
    pred.final.points <- matrix(0,nrow=its,ncol=K)
    pred.final.wins <- matrix(0,nrow=its,ncol=K)
    pred.final.podiums <- matrix(0,nrow=its,ncol=K)
    pred.final.top10s <- matrix(0,nrow=its,ncol=K)
    
    ## for each sampled set of parameter values - future race results (in season)
    ## convert to points for each driver
    ## convert to finishing position for each driver (NOT USED)
    for(i in 1:its){
      for(k in 1:num.remain.races){
        x.pred.sim[i,k,] <- these.drivers[rattrition(res.attrition$gm[i,x[tt+1,is.na(x[tt+1,])==FALSE]])]
        points.sim[i,k,x.pred.sim[i,k,]] <- points.calc.2013(1:length(these.drivers))
        is.winner.sim[i,k,x.pred.sim[i,k,1]] <- 1
        is.podium.sim[i,k,x.pred.sim[i,k,1:3]] <- 1
        is.top10.sim[i,k,x.pred.sim[i,k,1:10]] <- 1
        ##pos.sim[i,k,x.pred.sim[i,k,]] <- 1:length(these.drivers)
      }
      ## predicted end of season points for each driver
      if(num.past.races==1){
        pred.final.points[i,] <- points.actual[past.races,]+apply(points.sim[i,,],2,sum)
        pred.final.wins[i,] <- wins.actual[past.races,]+apply(is.winner.sim[i,,],2,sum)
        pred.final.podiums[i,] <- podiums.actual[past.races,]+apply(is.podium.sim[i,,],2,sum)
        pred.final.top10s[i,] <- top10s.actual[past.races,]+apply(is.top10.sim[i,,],2,sum)
      }
      else{
        if(num.past.races == (length.season-1)){
          pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+points.sim[i,,]
          pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+is.winner.sim[i,,]
          pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+is.podium.sim[i,,]
          pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+is.top10.sim[i,,]
        }
        else{
          if(num.past.races == length.season){
            pred.final.points[i,] <- apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(is.top10.sim[i,,],2,sum)
          }
          else{
            pred.final.points[i,] <- apply(points.actual[past.races,],2,sum)+apply(points.sim[i,,],2,sum)
            pred.final.wins[i,] <- apply(wins.actual[past.races,],2,sum)+apply(is.winner.sim[i,,],2,sum)
            pred.final.podiums[i,] <- apply(podiums.actual[past.races,],2,sum)+apply(is.podium.sim[i,,],2,sum)
            pred.final.top10s[i,] <- apply(top10s.actual[past.races,],2,sum)+apply(is.top10.sim[i,,],2,sum)
          }
        }
      }
    }
    ## simulated championship winner
    champ.winner[tt+1,] <-  apply(pred.final.points,1,which.max)
    expected.season.wins[tt+1,] <- apply(pred.final.wins,2,mean)
    expected.season.podiums[tt+1,] <- apply(pred.final.podiums,2,mean)
    expected.season.top10s[tt+1,] <- apply(pred.final.top10s,2,mean)
    
    ## expected number of points in next race
    points.mean[tt+1,] <- apply(points.sim[,1,],2,mean)
    
    ## store parameter summaries
    gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
    gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
    gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
    ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
    for(i in 1:K){
      for(j in 1:K){
        win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
      }
    }
    
    ## log scores for winner,podium, top10 etc ##########
    pred.prob.win <- NULL
    pred.prob.podium <- NULL
    pred.prob.top10 <- NULL
    ## deal with zero probs - like Bayes estimate of binomial proportion
    ## with beta prior
    a.zero <- 0
    b.zero <- 0
    if(zero.adjust){
      a.zero <- 1
      b.zero <- length(these.drivers)-1
    }
    new.denom <- a.zero+b.zero+its
    for(j in 1:length(these.drivers)){
      pred.prob.win[j] <- mean(x.pred.sim[,1,1]==these.drivers[j])
      pred.prob.podium[j] <- mean(apply(x.pred.sim[,1,1:3]==these.drivers[j],1,sum))
      pred.prob.top10[j] <- mean(apply(x.pred.sim[,1,1:10]==these.drivers[j],1,sum))
      pred.prob.win[j] <- (a.zero+(pred.prob.win[j]*its))/new.denom
      pred.prob.podium[j] <- (a.zero+(pred.prob.podium[j]*its))/new.denom
      pred.prob.top10[j] <- (a.zero+(pred.prob.top10[j]*its))/new.denom
      pred.probs[tt+1,these.drivers[j],1] <- pred.prob.win[j]
      pred.probs[tt+1,these.drivers[j],2] <- pred.prob.podium[j]
      pred.probs[tt+1,these.drivers[j],3] <- pred.prob.top10[j]    
    }
    outcome.win <- c(1,rep(0,length(these.drivers)-1))
    owlppw <- outcome.win*log(pred.prob.win) ## deal with 0 probs
    omowlomppw <- (1-outcome.win)*log(1-pred.prob.win)
    log.score.winner.sim[tt+1] <- sum(ifelse(owlppw=="NaN",0,owlppw)+ifelse(omowlomppw=="NaN",0,omowlomppw))  
    ##log.score.winner.sim[tt+1] <- sum(outcome.win*log(pred.prob.win) + (1-outcome.win)*log(1-pred.prob.win))
    outcome.podium <- c(rep(1,3),rep(0,length(these.drivers)-3))
    oplppp <- outcome.podium*log(pred.prob.podium) ## deal with 0 probs
    omoplomppp <- (1-outcome.podium)*log(1-pred.prob.podium)
    log.score.podium.sim[tt+1] <- sum(ifelse(oplppp=="NaN",0,oplppp)+ifelse(omoplomppp=="NaN",0,omoplomppp))  
    ##log.score.podium.sim[tt+1] <- sum(outcome.podium*log(pred.prob.podium) + (1-outcome.podium)*log(1-pred.prob.podium))
    outcome.top10 <- c(rep(1,10),rep(0,length(these.drivers)-10))
    ot10lppt10 <- outcome.top10*log(pred.prob.top10) ## deal with 0 probs
    omot10lomppt10 <- (1-outcome.top10)*log(1-pred.prob.top10)
    log.score.top10.sim[tt+1] <- sum(ifelse(ot10lppt10=="NaN",0,ot10lppt10)+ifelse(omot10lomppt10=="NaN",0,omot10lomppt10))  
    ##log.score.top10.sim[tt+1] <- sum(outcome.top10*log(pred.prob.top10) + (1-outcome.top10)*log(1-pred.prob.top10))
##################################
    
    ## compute log prior predictive
    log.pred.prob.all <- rep(0,its)
    for(i in 1:its) {
      log.pred.prob.all[i] <- loglikelihood.attrition(x[tt+1,],res.attrition$gm[i,])
    }
    max.lppa <- max(log.pred.prob.all)
    log.prior.pred[tt+1] <- log(mean(exp(log.pred.prob.all-max.lppa)))+max.lppa
  }
  ## last race of the season
  tt = n
  print(tt)
  gamma.init <- rgamma(K,cc,d)
  res.attrition <- twattrition.gibbs(x[1:tt,],K,tau.vec=tau.vec[1:tt],tau=tau,xi=xi,cc=cc,d=d,gamma.init=gamma.init,its=its,burn=burn,thin=thin)
  gamma.mean[tt,] <- apply(res.attrition$gm,2,mean)
  gamma.upper[tt,] <- apply(res.attrition$gm,2,quantile,0.975)
  gamma.lower[tt,] <- apply(res.attrition$gm,2,quantile,0.025)
  ##res.attrition.omega <- res.attrition$gm/rowSums(res.attrition$gm)
  for(i in 1:K){
    for(j in 1:K){
      win.pairs.prob[tt,i,j] <- mean(res.attrition$gm[,j]/(res.attrition$gm[,i]+res.attrition$gm[,j]))
    }
  }

  return(list(lpp=log.prior.pred,lsws=log.score.winner.sim,lsps=log.score.podium.sim,lst10s=log.score.top10.sim,pp=pred.probs,gm=gamma.mean,gu=gamma.upper,gl=gamma.lower,wpp=win.pairs.prob,pm=points.mean,cw=champ.winner,esw=expected.season.wins,esp=expected.season.podiums,est10=expected.season.top10s))
}


twattrition.varxi.mcmc <- function(x,K,tau.vec,tau,a.xi=1,b.xi=1,cc=rep(1,K),d=1,gamma.init=rep(1,K),xi.init=0.9999,sigma.xi=0.01,its=1000,burn=0,thin=1)
{
  ## MCMC sampler for time-weighted attrition model with xi unknown
  
  ## preliminary calculations
  ss <- summary.stats.twattrition(x,K,tau.vec,tau,xi.init)

  ## number of multiple comparisons
  n <- ss$n
  if(n == 1){
    x <- t(matrix(x))
  }
  
  ## number of individuals in each multiple comparison
  n.i <- ss$n.i
  
  ## time-weighted number of races in which driver did not win
  v.init <- ss$v
  
  ## compute beta
  beta <- ss$beta

  ## compute time-weghting
  time.weight.init <- ss$time.weight

  ## storage for samples
  U.res <- array(0,c(its,n,max(n.i)-1))
  gamma.res <- matrix(0,nrow=its,ncol=K)
  xi.res <- rep(0,its)
  ll.res <- matrix(0,nrow=its,ncol=n)
  ljd.res <- rep(0,its)
  
  ## initial values
  gamma.curr <- gamma.init
  xi.curr <- xi.init
  v.curr <- v.init
  time.weight.curr <- time.weight.init
  U.curr <- matrix(0,nrow=n,ncol=(max(n.i)-1))
  lprior.xi.curr <- dbeta(xi.curr,a.xi,b.xi,log=TRUE)
  llike.xi.curr <- loglike.twattrition(x,gamma.curr,tau.vec,tau,xi.curr)
  theta.curr <- log(xi.curr)-log(1-xi.curr)
  ljacobian.xi.curr <- -theta.curr-2*log(1+exp(-theta.curr))
  
  ## Markov chain sampling
  if(burn>0){
    for(t in 1:burn){
      print(t-burn)
      for(l in 1:thin){
        ## sample latent variables | everything else
        llatent.xi.curr <- 0
        for(i in 1:n){
          for(j in 1:(n.i[i]-1)){
            U.curr[i,j] <- rgamma(1,time.weight.curr[i],sum(gamma.curr[x[i,1:(j+1)]]))
            llatent.xi.curr <- llatent.xi.curr+dgamma(U.curr[i,j],time.weight.curr[i],sum(gamma.curr[x[i,1:(j+1)]]),log=TRUE)
          }
        }
        
        ## sample parameters | everything else
        for(k in 1:K){
          gamma.curr[k] <- rgamma(1,cc[k]+v.curr[k],d+sum(beta[,,k]*U.curr))
        }
        
        ## normalize and rescale (Section 5.1 of Caron and Doucet)
        Gamma.curr <- rgamma(1,sum(cc),d)
        gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)

        ##sample xi | everything else
        theta.prop <- rnorm(1,log(xi.curr)-log(1-xi.curr),sigma.xi)
        xi.prop <- 1/(1+exp(-theta.prop))
        v.prop  <- rep(0,K)
        time.weight.prop <- rep(0,n)
        for(k in 1:K){
          for(i in 1:n){
            time.weight.prop[i] <- psi(tau-tau.vec[i],xi.prop)
            for(j in 2:n.i[i]){
              v.prop[k] <- v.prop[k] + ifelse(x[i,j]==k,1,0)*time.weight.prop[i]
            }
          }
        }
        lprior.xi.prop <- dbeta(xi.prop,a.xi,b.xi,log=TRUE)
        llike.xi.prop <- loglike.twattrition(x,gamma.curr,tau.vec,tau,xi.prop)
        llatent.xi.prop <- 0
        for(i in 1:n){
          for(j in 1:(n.i[i]-1)){
            llatent.xi.prop <- llatent.xi.prop+dgamma(U.curr[i,j],time.weight.prop[i],sum(gamma.curr[x[i,1:(j+1)]]),log=TRUE)
          }
        }
        ljacobian.xi.prop <- -theta.prop-2*log(1+exp(-theta.prop))

        l.acc <- lprior.xi.prop+llike.xi.prop+llatent.xi.prop+ljacobian.xi.prop
        l.acc <- l.acc-(lprior.xi.curr+llike.xi.curr+llatent.xi.curr+ljacobian.xi.curr)

        if(log(runif(1,0,1)) < min(0,l.acc)){

          print("accepted xi")
          xi.curr <- xi.prop
          v.curr <- v.prop
          time.weight.curr <- time.weight.prop
          lprior.xi.curr <- lprior.xi.prop
          llike.xi.curr <- llike.xi.prop
          llatent.xi.curr <- llatent.xi.prop
          ljacobian.xi.curr <- ljacobian.xi.prop
        }
        
      }
    }
  }
  for(t in 1:its){
    print(t)
    for(l in 1:thin){
      ## sample latent variables | everything else
      llatent.xi.curr <- 0
      for(i in 1:n){
        for(j in 1:(n.i[i]-1)){
          U.curr[i,j] <- rgamma(1,time.weight.curr[i],sum(gamma.curr[x[i,1:(j+1)]]))
          llatent.xi.curr <- llatent.xi.curr+dgamma(U.curr[i,j],time.weight.curr[i],sum(gamma.curr[x[i,1:(j+1)]]),log=TRUE)
        }
      }
      
      ## sample parameters | everything else
      for(k in 1:K){
        gamma.curr[k] <- rgamma(1,cc[k]+v.curr[k],d+sum(beta[,,k]*U.curr))
      }
      
      ## normalize and rescale (Section 5.1 of Caron and Doucet)
      Gamma.curr <- rgamma(1,sum(cc),d)
      gamma.curr <- Gamma.curr*gamma.curr/sum(gamma.curr)

      ##sample xi | everything else
      theta.prop <- rnorm(1,log(xi.curr)-log(1-xi.curr),sigma.xi)
      xi.prop <- 1/(1+exp(-theta.prop))
      v.prop  <- rep(0,K)
      time.weight.prop <- rep(0,n)
      for(k in 1:K){
        for(i in 1:n){
          time.weight.prop[i] <- psi(tau-tau.vec[i],xi.prop)
          for(j in 2:n.i[i]){
            v.prop[k] <- v.prop[k] + ifelse(x[i,j]==k,1,0)*time.weight.prop[i]
          }
        }
      }
      lprior.xi.prop <- dbeta(xi.prop,a.xi,b.xi,log=TRUE)
      llike.xi.prop <- loglike.twattrition(x,gamma.curr,tau.vec,tau,xi.prop)
      llatent.xi.prop <- 0
      for(i in 1:n){
        for(j in 1:(n.i[i]-1)){
          llatent.xi.prop <- llatent.xi.prop+dgamma(U.curr[i,j],time.weight.prop[i],sum(gamma.curr[x[i,1:(j+1)]]),log=TRUE)
        }
      }
      ljacobian.xi.prop <- -theta.prop-2*log(1+exp(-theta.prop))
      
      l.acc <- lprior.xi.prop+llike.xi.prop+llatent.xi.prop+ljacobian.xi.prop
      l.acc <- l.acc-(lprior.xi.curr+llike.xi.curr+llatent.xi.curr+ljacobian.xi.curr)
      
      if(log(runif(1,0,1)) < min(0,l.acc)){
        
        print("accepted xi")
        xi.curr <- xi.prop
        v.curr <- v.prop
        time.weight.curr <- time.weight.prop
        lprior.xi.curr <- lprior.xi.prop
        llike.xi.curr <- llike.xi.prop
        llatent.xi.curr <- llatent.xi.prop
        ljacobian.xi.curr <- ljacobian.xi.prop
      }
      
    }
    
    ## store sampled values
    gamma.res[t,] <- gamma.curr
    xi.res[t] <- xi.curr
    U.res[t,,] <- U.curr
    
    ## compute observed data loglikelihood
    ll.res[t,] <- loglike.twattrition(x,gamma.res[t,],tau.vec,tau,xi.res[t])
    ## compute log joint density 
    ljd.res[t] <- sum(ll.res[t,])+sum(dgamma(gamma.res[t,],cc,d,log=TRUE))+dbeta(xi.res[t],a.xi,b.xi,log=TRUE)
    print(c(sum(ll.res[t,]),ljd.res[t]))
  }
  return(list(U=U.res,gm=gamma.res,xi=xi.res,ll=ll.res,ljd=ljd.res))
}

