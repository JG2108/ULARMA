#setwd("C:/Users/santa/OneDrive/Área de Trabalho/Mestrado/Dissertação/Scripts")


simu.ularma_v2 <- function(n,phi=NA,theta=NA,alpha=0.0,link="logit",X=NA,beta=NA)
{
  ar<-NA
  ma<-NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  }
  
  link <- structure(list(link = linktemp,
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  if(any(is.na(X) == TRUE)){
    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F))
    {
      #print("ARMA model")
      # seasonal part
      
      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)
      
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }
      
      return(ts(y[(m+1):(n+m)]))
      
    } # ARMA model
    
    
    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T))
    {
      #print("AR model")
      
      p <- max(ar)
      m <- 2*p
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      eta<- y <- rep(NA,n+m)
      
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + (phi%*%ynew[i-ar])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
      }
      
      return( ts(y[(m+1):(n+m)]))
    } # AR model
    
    
    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F))
    {
      #print("MA model")
      
      q <- max(ma)
      m <- 2*q
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      eta <- y <- error <- rep(0,n+m) # E(error)=0
      
      #print(ma)
      
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }
      
      return( ts(y[(m+1):(n+m)]))
    } # fim MA model
    
  } else {
    
    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F))
    {
      #print("ARMA model")
      # seasonal part
      
      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)
      
      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)
      
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + as.numeric(X[i,]%*%as.matrix(beta)) + as.numeric((phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))) + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
       
      }
      
      return(ts(y[(m+1):(n+m)]))
      
    } # ARMA model
    
    
    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T))
    {
      #print("AR model")
      
      p <- max(ar)
      m <- 2*p
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      eta<- y <- rep(NA,n+m)
      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
      }
      
      return( ts(y[(m+1):(n+m)]))
    } # AR model
    
    
    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F))
    {
      #print("MA model")
      
      q <- max(ma)
      m <- 2*q
      
      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)
      
      eta <- y <- error <- rep(0,n+m) # E(error)=0
      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)
      #print(ma)
      
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }
      
      return( ts(y[(m+1):(n+m)]))
    } # fim MA model
  }
  
  
}

