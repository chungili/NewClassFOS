FOS = function(x, alpha = 0.0027){
  n = length(x)
  Rx <- sort(x)
  Pn_FOS = function(n){
    if (n > 1851) {0.95
    } else if (n > 1481) {0.90
    } else if (n > 1111) {0.80
    } else if (n > 740)  {0.70
    } else if (n > 370)  {0.50
    } else {0.3}
  }
  limit = function(u){
    n1 = 1+n
    if (0<u && u<=1/n1) {
      Rx[1]+(Rx[2]-Rx[1])*log((n1)*u)
    } else if (1/n1<u && u<n/n1) { 
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
    }
  }
  fun1 = function(r){pbeta(alpha, r, n-r+1)-(1+Pn_FOS(n))/2}
  fun2 = function(r){pbeta(1-alpha, r, n-r+1)-(1-Pn_FOS(n))/2}
  kl = uniroot(fun1, c(0, n+1))$root
  ku = uniroot(fun2, c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  LCL.Q = limit(u=ur)
  UCL.Q = limit(u=us)
  limits = cbind(LCL.Q, UCL.Q)
  return(limits)
}


FOS_ad = function(x, alpha = 0.0027, pn = 0.9){
  n = length(x)
  Rx <- sort(x)
  n1 = n + 1
  m = function(pn){
    if (0.5<pn && pn<=0.7){370}
    else if (0.7<=pn && pn<=0.8){740}
    else if (0.8<pn && pn<=0.95){1111}
  }
  sigma = function(pn){
    if (0.5<pn && pn<=0.7){2}
    else if (0.7<pn && pn<=0.8){1}
    else if (0.8<pn && pn<=0.95){0.5}
  }
  pp0 = function(n){
    if(0<pn && pn<=0.5){1}
    else {pnorm(n/m(pn)-1, mean =0, sd =sigma(pn))}
  }
  limit = function(u){
    if (0<u && u<=1/(pp0(n)*n1)) {
      Rx[1]+(Rx[2]-Rx[1])*log(pp0(n)*n1*u)-(Rx[2]-Rx[1])*(log(pp0(n)*n1*u))^2
    } else if (1/(n1*pp0(n))<u && u<1-1/(n1*pp0(n))) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(pp0(n)*n1*(1-u))+(Rx[n]-Rx[n-1])*(log(pp0(n)*n1*(1-u)))^2
    }
  }
  fun1 = function(r){pbeta(alpha, r, n-r+1)-(1+pn)/2}
  fun2 = function(r){pbeta(1-alpha,r, n-r+1)-(1-pn)/2}
  kl <- uniroot(fun1,c(0, n+1))$root
  ku <- uniroot(fun2,c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  LCL.Q = limit(u=ur)
  UCL.Q = limit(u=us)
  limits = cbind(LCL.Q, UCL.Q)
  return(limits)
}


FOS_3terms <- function(x, alpha = 0.0027, pn = 0.9){
  n = length(x)
  Rx = sort(x)
  n1 = n + 1
  limit = function(u){
    if (0<u && u<=1/(n1)) {
      Rx[1]+(Rx[2]-Rx[1])*log(n1*u)-(Rx[2]-Rx[1])*(log(n1*u))^2+(Rx[2]-Rx[1])*(log(n1*u))^3
    } else if (1/n1<u && u<1-1/n1) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))+(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^2-(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^3
    }
  }
  fun1 = function(r){pbeta(alpha, r, n-r+1)-(1+pn)/2}
  fun2 = function(r){pbeta(1-alpha, r, n-r+1)-(1-pn)/2}
  kl = uniroot(fun1,c(0, n+1))$root
  ku = uniroot(fun2,c(0, n+1))$root
  ur = kl/(n+1)
  us = ku/(n+1)
  LCL = limit(u=ur)
  UCL = limit(u=us)
  limits = cbind(LCL, UCL)
  return(limits)
}
