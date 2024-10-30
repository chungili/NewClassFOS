# Install the packages

```{r}
packageName <- c("dfphase1", "DescTools", "moments", "tolerance")
for(i in 1:length(packageName)) {
  if(!(packageName[i] %in% rownames(installed.packages()))) {
    install.packages(packageName[i])
  }
}
```

# aaa
