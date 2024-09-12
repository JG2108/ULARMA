# Created by José Guilehrme (santana.guilherme@outlook.com), Maio/2022
#
# Some informations:
# diag = 0 : without graphs (useful for simulations)
# diag = 1 : plot diagnostic graphs
# diga = 2 : make pdf with diagnostic graphs
#
# (a,b) the interval of data
#
# h : number of forecast steps


ULARMA<- function (yab, ar=NA, ma=NA,link="logit",diag=1,h=NA,X=NA,X_hat=NA,resid=3){ 
  #resid: define o tipo de resíduo utilizado (1 - padronizado, 2 - deviance, 3 - quartílico)
  #yab: variável com distribuição UnitLindley
  #ar: termos autorregressivos
  #ma termos de médias móveis
  #link: função de ligação a ser utilizada ("logit", "probit", "cloglog")
  #diag: gráficos (0 - sem gráficos, 1 - gráficos de diagnósticos, 2 - gráficos de diagnósticos e pdf)
  #h: passos a frente de previsões
  #X: matrix de covariáveis exógenas
  #X_hat: matriz de covariáveis exógenas para previsões
  
  
  source("script1_funcoes_ULARMA.R")
  source("script3_ajuste_ULARMA.R")
  
  #verificando se a variável possui valores dentro do intervalo unitário padrão
  if (min(yab) <= 0 || max(yab) >= 1){
    stop("FORA DO INTERVALO (0,1)!")
  } else{
    y<- yab
  }
  
  #verificando se a variável é uma série temporal
  if(is.ts(y)==T){
    freq<-frequency(y)
  }else stop("Os dados podem ser uma série temporal")
  
  #verificando se o vetor com os termos AR/MA estão completos
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))
  
  #verificando a presença de covariaveis
  if(any(is.na(X))==F) names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  
  #pegando as ordens dos termos temporais
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  
  m <- max(p,q,na.rm=T)
  
  p1 <- length(ar)
  q1 <- length(ma)
  
  #ajeitando as saidas da função de ligação
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)#converte expressão em character
    if (linktemp == "link")
      linktemp <- eval(link)#verifica a expressão
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))    stats <- make.link(linktemp)
  else stop(paste(linktemp, "função de ligação não avaliada, usar \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta, 
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t))) # ver esse ponto aqui
  )
  )
  
  
  fit1 <- ULARMA.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat,resid=resid) # model estimation
  
  return(fit1)
}
