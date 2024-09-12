#setwd("C:/Users/santa/OneDrive/Área de Trabalho/Mestrado/Dissertação/Scripts")
require("expint")

##densidade da ULARMA
dUL <- function(y,mu){
  #padronização em termos da média
  d <- ((1-mu)^2/(mu*(1-y)^3))* exp(-((y*(1-mu))/(mu*(1-y))))
  return(d)
}

##distribuição acumulada
pUL <- function(y,mu){
  p <- 1 - ((mu*y-1)/(y-1)) * exp(-((y*(1-mu))/(mu*(1-y))))
  return(p)
}

##função quantil 
qUL <- function(y,mu){
  q <- ((1/mu) + pracma::lambertWn((1/mu)*(y-1)*exp(-(1/mu))))/(1 + pracma::lambertWn((1/mu)*(y-1)*exp(-(1/mu))))
  return(q)
}

##gerando valores da distribuição
rUL=function(n,mi){
  a=1/mi-1
  require(LindleyR)
  x=rlindley(n, a, mixture=TRUE)
  z=x/(1+x)
  return(z)
}

