# Prokounter
The respository contains the source code for the R package, Prokounter, for making differential richness inferences from 16S microbiome data. 

For more details, please refer: `https://www.biorxiv.org/content/10.1101/2021.11.07.467583v1`

```{r}
#installation
library(devtools)
devtools::install_github("mskb01/prokounter", ref="main") 


#usage
library(prokounter)
?prokounter:::getProkounterTrends()
?prokounter:::getDRforSampleGroups()
?prokounter:::getDRforGenera()
?prokounter:::getDRforTaxaCollection()
```
