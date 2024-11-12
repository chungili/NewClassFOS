FineTune = function(Sd = c(0.5, 1, 2), m = c(370, 740, 1111), n, pn, 
                    dist = c("normal", "chisq", "t", "beta", "lnrom"), 
                    alpha = 0.0027, B = 10000){
  n.m = length(m)
  n.Sd = length(Sd)
  FOS = function(x, pn, sd, m, alpha = alpha){
    Rx = sort(x)
    n = length(x)
    n1 = n + 1
    p0 = pnorm(n/m-1, mean = 0, sd = sd)
    Xu = function(u){
      logwl = log(p0*n1*u)
      logwu = log(p0*n1*(1-u) )
      if (0<u && u<=1/(p0*n1)) {
        Rx[1]+(Rx[2]-Rx[1])*logwl-(Rx[2]-Rx[1])*(logwl)^2
      } else if (1/(n1*p0)<u && u<1-1/(n1*p0)) {
        (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
      } else {
        Rx[n]-(Rx[n]-Rx[n-1])*logwu+(Rx[n]-Rx[n-1])*(logwu)^2
      }
    }
    findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
    findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
    kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ur = kl/n1
    us = ku/n1
    LCL = Xu(ur)
    UCL = Xu(us)
    return(list(UCL=UCL, LCL=LCL))
  }
  library(doParallel)
  cl = makeCluster(detectCores()-2)
  registerDoParallel(cl)
  result = foreach( m = rep(m, n.Sd), Sd = rep(Sd, each=n.m), 
                     .combine='cbind') %:%
    foreach (1:B, .combine = c) %dopar% {
      x0 = switch(dist,
                  "normal" = rnorm(n, mean = 0, sd = 1),
                  'chisq' = rchisq(n, df = 4),
                  "t" = rchisq(n, df = 4),
                  "beta" = rbeta(n, shape1 = 9, shape2 = 1),
                  "lnorm" = rlnorm(n, meanlog = 0, sdlog = 1))
      x1 = switch(dist,
                  "normal" = rnorm(n, mean = 0, sd = 1),
                  'chisq' = rchisq(n, df = 4),
                  "t" = rchisq(n, df = 4),
                  "beta" = rbeta(n, shape1 = 9, shape2 = 1),
                  "lnorm" = rlnorm(n, meanlog = 0, sdlog = 1))
      FOS_limits = FOS(x0, pn = pn, sd = Sd, m = m, alpha = alpha)
      mean( x1 > FOS_limits$UCL | x1 < FOS_limits$LCL ) < alpha
    }
  stopCluster(cl)
  pa = matrix( colMeans(result), nrow = n.Sd, ncol = n.m, byrow = T)
  loc = which(pa==min(pa), arr.ind = T)
  return(list(Opt.m = m[loc[2]], Opt.sd = Sd[loc[1]]) )
}

FineTuneWithData = function(Sd = c(0.5, 1, 2), m = c(370, 740, 1111), pn, Train, Test,
                    alpha = 0.0027, B = 10000){
  n = dim(Train)[2]
  n.m = length(m)
  n.Sd = length(Sd)
  FOS = function(x, pn, sd, m, alpha = alpha){
    Rx = sort(x)
    n = length(x)
    n1 = n + 1
    p0 = pnorm(n/m-1, mean = 0, sd = sd)
    Xu = function(u){
      logwl = log(p0*n1*u)
      logwu = log(p0*n1*(1-u) )
      if (0<u && u<=1/(p0*n1)) {
        Rx[1]+(Rx[2]-Rx[1])*logwl-(Rx[2]-Rx[1])*(logwl)^2
      } else if (1/(n1*p0)<u && u<1-1/(n1*p0)) {
        (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
      } else {
        Rx[n]-(Rx[n]-Rx[n-1])*logwu+(Rx[n]-Rx[n-1])*(logwu)^2
      }
    }
    findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
    findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
    kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
    ur = kl/n1
    us = ku/n1
    LCL = Xu(ur)
    UCL = Xu(us)
    return(list(UCL=UCL, LCL=LCL))
  }
  library(doParallel)
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  result <- foreach( m = rep(m, n.Sd), Sd = rep(Sd, each=n.m), 
                     .combine='cbind') %:%
    foreach (1:B, .combine = c) %dopar% {
      locTrain = sample(x = 1:B, size = 1, replace = TRUE)
      locTest = sample(x = 1:B, size = 1, replace = TRUE)
      x0 = Train[locTrain, ]
      x1 = Test[locTest, ]
      FOS_limits <- FOS(x0, pn = pn, sd = Sd, m = m, alpha = alpha)
      mean( x1 > FOS_limits$UCL | x1 < FOS_limits$LCL ) < alpha
    }
  stopCluster(cl)
  pa = matrix( colMeans(result), nrow = n.Sd, ncol = n.m, byrow = T)
  pa[pa-pn < 0]=1
  loc = which(pa==min(pa), arr.ind = T)
  return(list(Opt.m = m[loc[2]], Opt.sd = Sd[loc[1]]) )
}
