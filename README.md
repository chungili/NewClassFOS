# Install the packages

```{r}
packageName <- c("dfphase1", "DescTools", "moments", "tolerance")
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
```

# Data Preprocessing

```{r}
dt = read.csv("uci-secom.csv", header = F)
prepro = function(x){
  #filters out the missing values
  xfull = x[!is.na(x)]
  xtrim = DescTools::Trim(xfull, trim = 0.01)
  xfinal = xtrim
  return(xfinal)
}
x1 = prepro(dt$V2)[1:61]
x2 = prepro(dt$V25)[1:379] 
x3 = prepro(dt$V158)[1:751]
x4 = prepro(dt$V190)[1:1536]
```

# Calcuate EPC Control Limits

In cdoe.R, we provide the implements of the three methods, Class FOS, adaptive FOS, and 3-term FOS, for calculating EPC control limits.
For example, we would like to calcuate the EPC control limits for variable `x1`.

```{r}
EPC_FOS(x1 , alpha = 0.0027)
ad_FOS(x1, alpha = 0.0027, pn = 0.9)
Terms3_FOS(x1, alpha = 0.0027, pn = 0.9)
```
