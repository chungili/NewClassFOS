packageName <- c("dfphase1","DescTools","moments","tolerance")
for(i in 1:length(packageName)) {
  if(!(packageName[i] %in% rownames(installed.packages()))) {
    install.packages(packageName[i])
  }
}
lapply(packageName, require, character.only = TRUE)
library(dfphase1)
library(DescTools)
library(moments)
library(tolerance)

# Read the raw data
dt = read.csv("uci-secom.csv", header = F)

# There are 591 variables in the dataset.
# In the CSV file, the first column represents time stamps, 
# and the subsequent columns represent the quality characteristics.

prepro = function(x){
  #filters out the missing values
  xfull = x[!is.na(x)]
  xtrim=Trim(xfull, trim = 0.01)
  xfinal = xtrim
  return(xfinal)
}

# The 2nd, 25th, 158th and 190th variables are selected and 
# then named as x1, x2, x3, and x4, respectively.

x = x1 = prepro(dt$V2);x=x[1:61] 
# x= x2 = prepro(dt$V25);x=x[1:379] 
# x= x3 = prepro(dt$V158);x=x[1:751]
# x= x4 = prepro(dt$V190);x=x[1:1536]


variable="X1"
Nsim=n=length(x)
mean=round(mean(x),2)
stdev=round(sd(x),2)
skewness=round(skewness(x),2)
kurtosis=round(kurtosis(x),2)
moments <- c(mean,stdev,skewness,kurtosis);moments

Set.Confid=0.80  # 1-p
z=qnorm((1+Set.Confid)/2);z
alpha=(1-pnorm(3));alpha # alpha= upper 3 sigma;
z.alpha=-qnorm(alpha);z.alpha
par(mfrow=c(1,1))


#################################################################################################################
# 1.class_FOS approach
#################################################################################################################
Simu=n=length(x);n
Rx=sort(x)
n1=n+1

Pn_FOS=function(n){
  if (n>1851) {0.95
  } else if (n>1481) {0.90
  } else if (n>1111) {0.80
  } else if (n>740)  {0.70
  } else if (n>370)  {0.50
  } else {0.3}
}

limit=function(u){
  if (0<u && u<=1/n1) {
    Rx[1]+(Rx[2]-Rx[1])*log((n1)*u)
  } else if (1/n1<u && u<n/n1) {
    (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
  } else {
    Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
  }
}


# ----- Two of one-sided ----------------
fun1=function(r){pbeta(alpha,r,Simu-r+1)-(1+Pn_FOS(n))/2}
fun2=function(r){pbeta(1-alpha,r,Simu-r+1)-(1-Pn_FOS(n))/2}
kl <- uniroot(fun1,c(0,Simu+1))$root;kl
ku <- uniroot(fun2,c(0,Simu+1))$root;ku
ur=kl/(Simu+1);ur
us=ku/(Simu+1);us


limit(ur)
limit(us)
LCL.Q=limit(ur);UCL.Q=limit(us)

cbind(LCL.Q,UCL.Q)
ooc=length(x[x<LCL.Q])+length(x[x>UCL.Q])
ooc.FOS=ooc
ooc.FOS.p=round(100*ooc.FOS/Nsim,2)


#################################################################################################################
# 2.adaptive-FOS approach
#################################################################################################################
n=length(x);n
Rx=sort(x)
n1=n+1
pn=Set.Confid

m=function(pn){
  if (0.5<pn && pn<=0.7){370}
  else if (0.7<=pn && pn<=0.8){740}
  else if (0.8<pn && pn<=0.95){1111}
}

sigma=function(pn){
  if (0.5<pn && pn<=0.7){2}
  else if (0.7<pn && pn<=0.8){1}
  else if (0.8<pn && pn<=0.95){0.5}
}

pp0=function(n){
  if(0<pn && pn<=0.5){1}
  else {pnorm(n/m(pn)-1, mean =0, sd =sigma(pn))}
}

limit=function(u){
  if (0<u && u<=1/(pp0(n)*n1)) {
    Rx[1]+(Rx[2]-Rx[1])*log(pp0(n)*n1*u)-(Rx[2]-Rx[1])*(log(pp0(n)*n1*u))^2
  } else if (1/(n1*pp0(n))<u && u<1-1/(n1*pp0(n))) {
    (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
  } else {
    Rx[n]-(Rx[n]-Rx[n-1])*log(pp0(n)*n1*(1-u))+(Rx[n]-Rx[n-1])*(log(pp0(n)*n1*(1-u)))^2
  }
}


# ----- Two of one-sided ----------------
fun1=function(r){pbeta(alpha,r,Simu-r+1)-(1+Set.Confid)/2}
fun2=function(r){pbeta(1-alpha,r,Simu-r+1)-(1-Set.Confid)/2}
kl <- uniroot(fun1,c(0,Simu+1))$root;kl
ku <- uniroot(fun2,c(0,Simu+1))$root;ku
ur=kl/(Simu+1);ur
us=ku/(Simu+1);us


limit(ur)
limit(us)
LCL.ad_FOS=limit(ur);UCL.ad_FOS=limit(us)

cbind(LCL.ad_FOS,UCL.ad_FOS)
ooc=length(x[x<LCL.ad_FOS])+length(x[x>UCL.ad_FOS])
ooc.ad_FOS=ooc
ooc.ad_FOS.p=round(100*ooc.ad_FOS/Nsim,2)


#################################################################################################################
# 3.3-term FOS approach
#################################################################################################################
n=length(x);n
Rx=sort(x)
n1=n+1
pn=Set.Confid

limit=function(u){
  if (0<u && u<=1/(n1)) {
    Rx[1]+(Rx[2]-Rx[1])*log(n1*u)-(Rx[2]-Rx[1])*(log(n1*u))^2+(Rx[2]-Rx[1])*(log(n1*u))^3
  } else if (1/n1<u && u<1-1/n1) {
    (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
  } else {
    Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))+(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^2-(Rx[n]-Rx[n-1])*(log(n1*(1-u)))^3
  }
}


# ----- Two of one-sided ----------------
fun1=function(r){pbeta(alpha,r,Simu-r+1)-(1+Set.Confid)/2}
fun2=function(r){pbeta(1-alpha,r,Simu-r+1)-(1-Set.Confid)/2}
kl <- uniroot(fun1,c(0,Simu+1))$root;kl
ku <- uniroot(fun2,c(0,Simu+1))$root;ku
ur=kl/(Simu+1);ur
us=ku/(Simu+1);us


limit(ur)
limit(us)
LCL.3term_FOS=limit(ur);UCL.3term_FOS=limit(us)

cbind(LCL.3term_FOS,UCL.3term_FOS)
ooc=length(x[x<LCL.3term_FOS])+length(x[x>UCL.3term_FOS])
ooc.3term_FOS=ooc
ooc.3term_FOS.p=round(100*ooc.3term_FOS/Nsim,2)


#=====================================================================================================================================
# ( Nonparametric CL : Two-Sided limits, Does (2020))
# ===================================================================================================================================

out <- nptol.int(x = x, alpha = 1-Set.Confid, P = 1-2*alpha, side = 2,
                 method = "YM", upper = NULL, lower = NULL)
out1=data.frame(out)
LCL=out1$X2.sided.lower
UCL=out1$X2.sided.upper
length3=out1$X2.sided.upper-out1$X2.sided.lower

center=median(x)

UCL3=UCL.Goed=UCL
LCL3=LCL.Goed=LCL


length3=Band.Goed=UCL-LCL

# =========== OOC =====================================
ooc=length(x[x<LCL.Goed])+length(x[x>UCL.Goed])
ooc.Goed=ooc;
ooc.Goed.p=round(100*ooc.Goed/Nsim,2)


# COMBINED SPC CHARTS

LCL.FOS=LCL.Q
UCL.FOS=UCL.Q
Band.Goed=UCL.Goed-LCL.Goed
Band.FOS=UCL.FOS-LCL.FOS

#HJYANGC
# =========== nonPARAMETRIC  =====================================
LD=min(x,LCL.FOS,LCL.ad_FOS,LCL.3term_FOS)
UD=max(x,UCL.FOS,UCL.ad_FOS,UCL.3term_FOS)
lss=abs(UD-LD)/8
uss=abs(UD-LD)/1.5
plot(x,ylim=c(LD-lss,UD+uss),xlab="sample sequence", ylab="values",
     main=paste0(variable,": mean=",moments[1],", stdev=",moments[2],", skew=",moments[3],", kurt=",moments[4], ", n=",Nsim," (Pn=",round(pn,3),")"),cex.main=1.5,cex = 0.8,pch=20,lty=1)
lines(x)


# class-FOS
abline(h=UCL.FOS,col="brown",lty=3,lwd=3) 
abline(h=LCL.FOS,col="brown",lty=3,lwd=3)

# ad-FOS
abline(h=UCL.ad_FOS,col="red",lty=4,lwd=3) 
abline(h=LCL.ad_FOS,col="red",lty=4,lwd=3)

# 3term-FOS
abline(h=UCL.3term_FOS,col="orange",lty=5,lwd=3) 
abline(h=LCL.3term_FOS,col="orange",lty=5,lwd=3)


legend(x = "topright", legend = c(paste0("classic FOS"," ( Pn=",Pn_FOS(n)," )"),"ad-FOS","3-term FOS"),
       col = c("brown","red","orange"),
       lty = c(3,4,5), cex = 1.2, box.lty = 1,lwd=c(3,3,3))


Pn_FOS=Pn_FOS(n)
CL=rbind(LCL.FOS,UCL.FOS,LCL.ad_FOS,UCL.ad_FOS,LCL.3term_FOS,UCL.3term_FOS)
result=rbind(variable,n,CL,Set.Confid,Pn_FOS);result

write.csv(result,"result_READ.csv")

