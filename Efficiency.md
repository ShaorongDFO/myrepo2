qPCR Efficiency Results
================
Shaorong
2025-05-06

## R Markdown

qPCR processing of gill hypoxia plate 5 for stats

``` r
dilut<- read.csv("dilut_plate5.csv",header=T,row.names=1)  #dilutions as rownames, 3 HK at front columns
#
Log10<- log10(c(0.0016,0.008,0.04,0.2,1))

##1. Efficiency

#stack into columns
eff1<- stack(dilut)
eff1$Lconc<- Log10
eff2<- eff1[,c(2,3,1)]; rm(eff1) #reorder columns
colnames(eff2)[c(1,3)]<- c("assay","ct")

#linear model
library(plyr) #required package
models <- dlply(eff2, "assay", function(df) 
  lm(ct ~ Lconc, data = df))

#information from the models
res<- ldply(models, function(x) {
  r.sq <- summary(x)$r.squared
  slope <- summary(x)$coefficients[2]
  data.frame(r.sq, slope) })

#efficiency calculation
res$efficiency<- 10^ -(1/res$slope) - 1
res<- res[,c(1,4,2,3)]  #reorder

#save results

write.csv(res,"eff_hypoxia_plate5.csv",row.names=F)

#Plots to check
res2<- levels(as.factor(eff2$assay))

for (i in 1:length(res2)){
  plot(ct~ Lconc, eff2[eff2$assay==res2[i],], xlim=c(-3.5,0),pch=16, cex=1.2, main=res2[i])
  mod<- lm(ct~ Lconc, eff2[eff2$assay==res2[i],])
  effi<- 10^ -(1/summary(mod)$coefficients[2]) -1
  abline(mod,col=2)
  text(-0.5, max(eff2[eff2$assay==res2[i],]$ct,na.rm=T),paste("R2=",round(100*summary(mod)$r.squared,0)),cex=1.2 )
  text(-3, min(eff2[eff2$assay==res2[i],]$ct,na.rm=T),paste("eff=",round(effi,2)),cex=1.2 )
}
```

![](Efficiency_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-10.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-11.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-12.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-13.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-14.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-15.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-16.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-17.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-18.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-19.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-20.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-21.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-22.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-23.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-24.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-25.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-26.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-27.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-28.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-29.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-30.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-31.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-32.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-33.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-34.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-35.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-36.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-37.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-38.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-39.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-40.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-41.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-42.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-43.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-44.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-45.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-46.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-47.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-48.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-49.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-50.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-51.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-52.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-53.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-54.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-55.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-56.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-57.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-58.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-59.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-60.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-61.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-62.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-63.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-64.png)<!-- -->![](Efficiency_files/figure-gfm/unnamed-chunk-1-65.png)<!-- -->

## R Markdown

This will create a well-formatted report.
