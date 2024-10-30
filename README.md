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

# aaa
