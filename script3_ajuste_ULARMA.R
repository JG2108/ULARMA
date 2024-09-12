#setwd("C:/Users/santa/OneDrive/Área de Trabalho/Mestrado/Dissertação/Scripts")

ULARMA.fit<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid){
  #y: série modelada
  #ar: número correspondente ao termo autorregressivo
  #ma: número corresposdnente ao termo de médias móveis
  #link: função de ligação utilizada
  #names_phi: ordem autorregressiva
  #names_theta: ordem médias móveis
  #names_beta: nome das variáveis exógenas
  #diag: gráficos (0 - sem gráficos, 1 - gráficos de diagnósticos, 2 - gráficos de diagnósticos e pdf)
  #h1: número de previsões realizadas
  #X: matrix de covariáveis exógenas
  #X_hat: matriz de covariáveis exógenas para previsões
  #resid: define o tipo de resíduo utilizado (1 - padronizado, 2 - deviance, 3 - quartílico)
  
  #pacotes necessários
  require("expint")
  source("erro_padrao_heissiana.R" )
  #definindo o número máximo de iterações para atingir a convergência
  maxit1<-100
  
  #saida/resultado da função
  z <- c()
  
  #ajeitando as saídas
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  #aplicando a função de ligação na resposta
  ynew = linkfun(y)
  ystar = log(y/(1-y))#faz pra comparação com codigo acima??
  
  #Identificando os termos temporais
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  #vetor de previsões
  y_prev <- c(rep(NA,(n+h1)))
  
  # Estimando as covariáveis do termo autorregressivo
  if(any(is.na(ar)==F)){
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)#matriz com a depêndencia temporal de ordem p
    
    for(i in 1:(n-m)){
      P[i,] <- ynew[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  } else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  # matriz de desenho do modelo (com ou sem variáveis exógenas)
  if(any(is.na(X)==T)){# sem variávevis exógenas
    x <- as.matrix(Z)
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    ajuste = lm.fit(x, Ynew)#ajuste do modelo linear para estimar a dependência
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    
  }else{
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n1 - k) * (dlink)^2)# pq aparece o termo dlink no denominador? (aula 6 de MLG)
  }
  
  ######### Estrutura do modelo sem a presença de covariáveis exógenas, somente a estrutura temporal (ARMA)
  ##neste caso é considerado os termos autorregressivos e de médas móveis 
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T)){ 
    reg <- c(mqo, rep(0,q1)) # valores iniciais para os parâmetros
    
    #função de log-verossimilhança
    loglik <- function(z){
      #z é o vetor de parâmetros
      #alpha é constante
      #phi é termos autorregressivos
      #theta é medias móveis
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n){
        eta[i] <- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    } 
    
    names_par <- c("alpha",names_phi,names_theta)
    
    opt <- optim(reg, loglik,
                 method = "BFGS", 
                 control = list(fnscale = -1), 
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+1)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }  
  
  ##neste caso é considerado somente o termo autorregressivo
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T)){
    
    q1<-0
    reg <- c(mqo) # valores iniciais dos parâmetros
    
    #função de log-verossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]
      
      eta<-rep(NA,n)
      
      for(i in (m+1):n){
        eta[i]<-alpha + (phi%*%ynew[i-ar])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }
    
    names_par <- c("alpha",names_phi)
    
    opt <- optim(reg, loglik,
                 method = "BFGS", 
                 control = list(fnscale = -1),
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+1)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    
    z$alpha <- alpha
    z$phi <- phi
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + (phi%*%ynew[i-ar])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
    }
    
  }
  
  ##neste caso é considerado somente o termo médias móveis
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T)){
    p1<-0
    reg <- c(mqo,rep(0,q1)) # valores iniciais dos parâmetros
    
    loglik <- function(z){
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      
      eta <- error <- rep(0,n) # E(error)=0
      
      for(i in (m+1):n){
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<-ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_theta)
    
    opt <- optim(reg, loglik,
                 method = "BFGS", 
                 control = list(fnscale = -1),
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+1)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[2:(q1+1)]
    
    z$alpha <- alpha
    z$theta <- theta
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + (theta%*%errorhat[i-ma])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }
  
  ######### Estrutura do modelo com a presença de covariáveis exógenas
  ##neste caso é considerado os termos temporais no modelo (com covariáveis)
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F)){ 
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1),beta1) # valores iniciais dos parâmetros
    
    #função de logverossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      beta <- z[(p1+q1+2):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n){
        
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] # predictor scale
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    } 
    
    names_par <- c("alpha",names_phi,names_theta,names_beta)
    
    opt <- optim(reg, loglik,
                 method = "BFGS", 
                 control = list(fnscale = -1),
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+1+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    beta <- coef[(p1+q1+2):length(coef)]
    
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }  
  
  ##neste caso é considerado somente os termos autorregressivos no modelo (com covariáveis)
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F)){ 
    q1<-0
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], beta1) # valores iniciais dos parâmetros
    
    #função de logverossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]
      beta <- z[(p1+2):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    } 
  
    names_par <- c("alpha",names_phi,names_beta)
    
    opt <- optim(reg, loglik, 
                 method = "BFGS", 
                 control = list(fnscale = -1),
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+1+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    beta <- coef[(p1+2):length(coef)]
    
    z$alpha <- alpha
    z$phi <- phi
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }
  
  ##neste caso é considerado somente os termos de medias moveis no modelo (com covariáveis)
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F)){
    p1<-0
    beta1<- mqo[(2):length(mqo)]
    reg <- c(mqo[1], rep(0,q1), beta1) # valores iniciais dos parâmetros
    
    #função de log-verossimilhança
    loglik <- function(z){
      alpha <- z[1]
      theta = z[(2):(q1+1)]
      beta <- z[(q1+2):length(z)]
      
      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)
      
      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_theta,names_beta)
    
    opt <- optim(reg, loglik, 
                 method = "BFGS", 
                 control = list(fnscale = -1),
                 hessian = TRUE)
    
    if (opt$conv != 0){
      warning("A FUNÇÃO NÃO CONVERGIU!")
    }else{
      warning("A FUNÇÃO CONVERGIU!")
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+1+ncol(X) )]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[(2):(q1+1)]
    beta <- coef[(q1+2):length(coef)]
    
    z$alpha <- alpha
    z$theta <- theta
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    K <- opt$hessian
    z$K <- K
    
    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)#X[i - ar, ]
    
    
    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
    
  } 
  
  z$serie <- y
  z$ularma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  
  ########################################## residuals
  res1 <- (y-z$fitted) #resíduo ordinário
  
  vary <- z$fitted[(m+1):n]*(((1/z$fitted[(m+1):n])-1)^2*exp((1/z$fitted[(m+1):n])-1) *expint_En(x = ((1/z$fitted[(m+1):n])-1), order = 1)
                    -(1/z$fitted[(m+1):n])+2) - z$fitted[(m+1):n]^2
  
  #resíduo pradronizado
  z$resid1 <- (res1[(m+1):n]/sqrt(vary)) 
  
  #cálculo da deviance
  l_tilde <- log(dUL(y,y))
  l_hat <- log(dUL(y,z$fitted))
  
  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)]<-0
  
  z$l_hat <- l_hat
  
  #resíduo da deviance
  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  #resíduos quartilico
  z$resid3 <- as.vector(qnorm(pUL(y[(m+1):n],z$fitted[(m+1):n])))
  
  
  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3
  
  #matriz de informação observada
  #aqui é feita a decomposição de cholesk para verificar se a matriz é positiva definida, ou seja, inversível (https://www.ime.unicamp.br/~marcia/AlgebraLinear/Arquivos%20PDF/demo_cholesky.pdf)
  Kchol<- tryCatch(chol(z$K), error = function(e) return("error"), 
                   warning = function(o) return("error"))#https://r-lang.com/r-trycatch-function/
  
  if(Kchol[1] == "error"){ 
    z$vcov <- try(solve(K))#https://blog.curso-r.com/posts/2017-04-09-try/
    warning("We have problems with information matrix inversion!")
    
  }else{
    vcov <- try(chol2inv(Kchol))
    z$vcov <- vcov
  }
  
  #desvio padrão
  stderror <- sqrt(diag(z$vcov))
  if(any(is.na(stderror))){
    stderror <- tryCatch(SEfromHessian(opt$hessian), 
                         error = function(e){
                           library(numDeriv)
                           
                           rep(1,length(z$coef))
                         } , 
             warning = function(o) return("error"))
    
  } 
  z$stderror <- stderror
  
  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  
  #medidas comparação
  if(any(is.na(X)==F)){
    z$k<- (p1+q1+1+length(beta))
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+1)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  z$model <- model_presentation
  z$link <- link
  
  if(diag>0){
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="Fitted values",ylab="Observed data",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
         xlab="Fitted values",ylab="Residuals")
    
    densidade<-density(residc)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="ACF",xlab="Lag") 
    
    pacf(residc,ylab="PACF",xlab="Lag") 
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantiles",ylab="Empirical quantiles")
    lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="serie",xlab="tempo")
    lines(z$fitted,col="red")
    
    fim<-end(y)[1]+end(y)[2]/12
    
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    par(mfrow=c(1,1))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
    abline(v=fim,lty=2)
    lines(y)
    
    w1<-5
    h1<-4
    
    if(diag>1)
    {
      png(file = "resid_v_ind.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      png(file = "resid_v_fitted.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
             xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      png(file = "obs_v_fit.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
             xlab="Fitted values",ylab="Observed data",
             xlim=c(0.95*min_y,max_y*1.05),
             ylim=c(0.95*min_y,max_y*1.05))
        lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
      }
      dev.off()
      
      png(file = "resid_density.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        
        plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
      
      png(file = "resid_FAC.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(residc,ylab="ACF",xlab="Lag") 
      }
      dev.off()
      
      png(file = "resid_FACP.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        pacf(residc,ylab="PACF",xlab="Lag")
      }
      dev.off()
      
      png(file = "qq_plot.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(residc, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantiles",ylab="Empirical quantiles")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
      
      png(file = "adjusted.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Serie",xlab="Time")
        lines(z$fitted,col="red")
      }
      dev.off()
      
      
      png(file = "forecast.png",horizontal=F,paper="special",width = 6, height = 4.7,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",lty=2,col="red", ylim=c(min(y),max(y)),ylab="RH",xlab="Times")
        abline(v=fim,lty=2)
        lines(y)
        legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"red"))
      }
      dev.off()
      
    }    
  }  
  
  return(z)
  
  
  
}


















