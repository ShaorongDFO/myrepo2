##qPCR processing of gill hypoxia plate 5 for stats

## setwd to source file location

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

#results

write.csv(res,"eff_hypoxia_plate5.csv",row.names=F)

#Plots to check
res2<- levels(as.factor(eff2$assay))
pdf("hypoxia_plate5_eff.pdf")
for (i in 1:length(res2)){
  plot(ct~ Lconc, eff2[eff2$assay==res2[i],], xlim=c(-3.5,0),pch=16, cex=1.2, main=res2[i])
  mod<- lm(ct~ Lconc, eff2[eff2$assay==res2[i],])
  effi<- 10^ -(1/summary(mod)$coefficients[2]) -1
  abline(mod,col=2)
  text(-0.5, max(eff2[eff2$assay==res2[i],]$ct,na.rm=T),paste("R2=",round(100*summary(mod)$r.squared,0)),cex=1.2 )
  text(-3, min(eff2[eff2$assay==res2[i],]$ct,na.rm=T),paste("eff=",round(effi,2)),cex=1.2 )
}
dev.off()

#Check figures and efficiency before proceeding to #2
##2. Normfinder: find best housekeeping

samp1<- read.csv("raw_plate5.csv",header=T,row.names=1)     #samples as rownames, 3 HK at front columns
samp<- samp1[-c(1:5),] #remove dilutions, now 5

rownames(samp)  #check
samp_hk<- samp[,1:3] #three HK
#replicate efficiency
eff_r<- as.data.frame(replicate(nrow(samp_hk),res$efficiency[c(83, 13, 54)])) #samples x HK assays
colnames(eff_r)<- rownames(samp_hk);rownames(eff_r)<- colnames(samp_hk)
eff_t<- t(eff_r)  #transpose

#minimum
min1<- as.data.frame(apply(samp_hk,2,min))  #min of HK assays (column)
min2<- as.data.frame(replicate(nrow(samp_hk),min1)) #replicate min1 for the samples(row)
rownames(min2)<- colnames(samp_hk);colnames(min2)<- rownames(samp_hk)
min_t<- t(min2);rm(min1,min2) #transpose

#linear transformation
dat_line<- eff_t^(min_t-samp_hk)
#housekeeping
hk_line<- dat_line[,c(1:3)]

##load the NormFinder.R file
#function needs samples in the columns, as column names
#function needs assays in the rows: last row is group
hk_lineT1<- as.data.frame(t(hk_line))  #so transpose
#function needs group in last row  (moved some norm and hypoxia together if single)
hk_lineT<- rbind(hk_lineT1,group=c ('SW18HS','SW18H','SW18ND','SW18N','SW18H','SW18ND','SW18N','SW18H','SW18H','SW18N',
                                   'SW18H','SW18H','SW18HD','SW18N','SW18HD','SW18H','SW18HD','SW18N','SW18HS','SW18ND','SW18N','SW18HD','SW18HS','SW18ND',
                                   'SW18N','SW18H','SW18HD','SW18N','SW18HS','SW18ND','SW18N','SW18H','SW18H','SW18N','SW14N','SW18HS','SW18HS','SW18HS',
                                   'SW18HS','SW18HS','SW18HS','SW18HS','SW18H'))

test_line<- Normfinder(hk_lineT, ctVal=F) #linearized CT
#want smallest stability   
test_line  #78d16.1 and COILP84, 2nd 78d16.1, 3rd COILP84 on its own

##3. delta and calculations 

dilut<- samp1[c(1:5),] #from raw dataframe not dilut (efficiency)
#
dilut_hk<- dilut[,1:3] #three housekeeping
samp_hk<- samp[,1:3]

#mean
dilut_m<- data.frame(avg=dilut_hk[,2])  #COILP84, same as challenges
samp_m<- data.frame(avg=samp_hk[,2]) 
#geometric mean (if three genes)
#dilut_m<- data.frame(avg= exp(rowMeans(log(dilut_hk)[,1:3]))) 
#samp_m<- data.frame(avg= exp(rowMeans(log(samp_hk)[,1:3]))) 

##delta CT: target gene - housekeeping gene
samp_dCT<- samp - samp_m$avg

##delta delta CT: target sample - reference sample (cDNA 1.0) 
#cDNA 1.0 delta CT
dilut_dCT1<- dilut[5,] - dilut_m$avg[5]  #check that row 5 is pool concentration 1
dilut_dCT<- as.data.frame(lapply(dilut_dCT1, rep, nrow(samp_dCT))) #expand
samp_ddCT<- samp_dCT- dilut_dCT

#transformation
dat1<- log2(2^(-samp_ddCT))
dat<- dat1[,-c(1:3)] #remove housekeeping

#Export
write.csv(dat,"gill_plate5.csv")