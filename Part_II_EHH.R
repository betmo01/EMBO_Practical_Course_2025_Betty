#Path
input="/Users/patriciopezo/Desktop/EMBO_course_2025/input/"
out="/Users/patriciopezo/Desktop/EMBO_course_2025/data_process/"

#I. Install the rehh R package
install.packages("rehh")

#II. Load rehh R package
library("rehh")


#III. Use the following files
#Chr2_EDAR_LWK_500K.recode.vcf #(African population)
#Chr2_EDAR_CHS_500K.recode.vcf # (East Asian population)


# IV. What is the profile of ancestral and derived haplotypes of the rs3827760 SNP in AFR and EAS?

#1.	Convert the data to haplohh format
data1<-data2haplohh(hap_file = paste0(input,"Chr2_EDAR_LWK_500K.recode.vcf"), polarize_vcf = F, vcf_reader = "data.table") 
## * Reading input file(s) *
## Using package 'data.table' to read vcf.
## Extracting map information.
## Extracting haplotypes.
## Number of individuals which are 
## Haploid Diploid Triploid, ... : 
## 1 2 
## 0 99 
## * Filtering data *
## Discard markers genotyped on less than 100 % of haplotypes.
## No marker discarded.
## Data consists of 198 haplotypes and 29016 markers.
## Number of mono-, bi-, multi-allelic markers:
## 1 2 
## 21289 7727

data2<-data2haplohh(hap_file = paste0(input, "Chr2_EDAR_CHS_500K.recode.vcf"), polarize_vcf = F, vcf_reader = "data.table")
## * Reading input file(s) *
## Using package 'data.table' to read vcf.
## Extracting map information.
## Extracting haplotypes.
## Number of individuals which are 
## Haploid Diploid Triploid, ... : 
## 1 2 
## 0 105 
## * Filtering data *
## Discard markers genotyped on less than 100 % of haplotypes.
## No marker discarded.
## Data consists of 210 haplotypes and 29016 markers.
## Number of mono-, bi-, multi-allelic markers:
## 1 2 
## 24709 4307


#2.	Calculate the EHH for the candidate SNP (rs3827760) in AFR
ehh_calc_AFR<-calc_ehh(data1,mrk = "rs3827760")
ehh_calc_AFR
## An object of class "ehh"
## [[1]]
## [1] "rs3827760"
## 
## [[2]]
## FREQ_A FREQ_D 
##      1      0 
## 
## [[3]]
##              POSITION      EHH_A EHH_D
## rs552689611 109506332 0.05081270     0
## rs188449710 109506337 0.05081270     0
## rs191146014 109506387 0.05081270     0
## rs562458938 109506482 0.05081270     0
## rs182888230 109506509 0.05081270     0
## rs187445332 109506534 0.05081270     0
## rs113027039 109506633 0.05081270     0
## rs569087080 109506722 0.05081270     0
## rs561388021 109506760 0.05081270     0
## rs557899408 109506776 0.05081270     0

## rs535613872 109527004 0.05004358     0
## rs138769166 109527018 0.05004358     0
## rs201588688 109527031 0.05004358     0
## rs142670672 109527045 0.05004358     0
## rs535840595 109527050 0.05004358     0
## 
## [[4]]
##    IHH_A    IHH_D 
## 1767.692    0.000


#3.	Calculate the EHH for the candidate SNP (rs3827760) in EAS
ehh_calc_EAS<-calc_ehh(data2,mrk = "rs3827760")
ehh_calc_EAS
## An object of class "ehh"
## [[1]]
## [1] "rs3827760"
## 
## [[2]]
##    FREQ_A    FREQ_D 
## 0.0952381 0.9047619 
## 
## [[3]]
##              POSITION      EHH_A      EHH_D
## rs144341651 109411690 0.00000000 0.07234754
## rs546981503 109411706 0.00000000 0.07234754
## rs568352591 109411793 0.00000000 0.07234754
## rs529232887 109411800 0.00000000 0.07234754
## rs550870176 109412087 0.00000000 0.07234754
## rs569104695 109412193 0.00000000 0.07234754
## rs113683672 109412248 0.00000000 0.07234754
## rs79168135  109412332 0.00000000 0.07234754
## rs75277911  109412390 0.00000000 0.07234754
## rs147448762 109412402 0.00000000 0.07234754

## rs573491349 109586154 0.00000000 0.05519354
## rs534276564 109586174 0.00000000 0.05329992
## rs554546487 109586260 0.00000000 0.05329992
## rs574390523 109586275 0.00000000 0.05329992
## rs260694    109586313 0.00000000 0.05329992
## rs563073581 109586323 0.00000000 0.05329992
## rs534812588 109586359 0.00000000 0.05329992
## 
## [[4]]
##     IHH_A     IHH_D 
##  9979.001 55231.782

#4.	Plot EHH around “rs3827760” in AFR
plot(ehh_calc_AFR)

#5.	Plot EHH around “rs3827760” in EAS
plot(ehh_calc_EAS)

#6.	Calculate furcation trees around a candidate SNP in AFR
furcation<-calc_furcation(data1, mrk="rs3827760")
plot(furcation)

#7.	Calculate furcation trees around a candidate SNP in EAS
furcation<-calc_furcation(data2, mrk="rs3827760")
plot(furcation)





#iHS
#1.	Calculate the iHS for all SNPs in the file (~5min) for AFR
AFR<-scan_hh(data1)

#2.	Calculate the iHS for all SNPs in the file (~5min) for AFR
EAS<-scan_hh(data2)

#3.	Check eHH statistics for candidate SNP for AFR
AFR[AFR$POSITION==109513601,]
##           CHR  POSITION FREQ_A FREQ_D NHAPLO_A NHAPLO_D    IHH_A IHH_D      IES
## rs3827760   2 109513601      1      0      198        0 1767.692     0 1767.692
##               INES
## rs3827760 1767.692

#4.	Check eHH statistics for candidate SNP for EAS
EAS[EAS$POSITION==109513601,]
##           CHR  POSITION    FREQ_A    FREQ_D NHAPLO_A NHAPLO_D    IHH_A    IHH_D
## rs3827760   2 109513601 0.0952381 0.9047619       20      190 9979.001 55231.78
##               IES    INES
## rs3827760 43767.3 54737.7

#5.	Estimate the iHS in AFR (use min_maf = 0.02, freqbin = 0.01)
iHS.AFR<-ihh2ihs(AFR, min_maf = 0.02, freqbin = 0.01)

#6.	Estimate the iHS in EAS (use min_maf = 0.02, freqbin = 0.01)
iHS.EAS<-ihh2ihs(EAS, min_maf = 0.02, freqbin = 0.01)

#7.	Check the iHS score for the candidate SNP in AFR
iHS.AFR$ihs[iHS.AFR$ihs$POSITION==109513601,]
## [1] CHR       POSITION  IHS       LOGPVALUE
## <0 rows> (or 0-length row.names)

#8.	Check the iHS score for the candidate SNP in EAS
iHS.EAS$ihs[iHS.EAS$ihs$POSITION==109513601,]
##           CHR  POSITION       IHS LOGPVALUE
## rs3827760   2 109513601 -1.588263 0.9499032

#9.	Plot the iHS score in EAS
plot(iHS.EAS$ihs$POSITION, iHS.EAS$ihs$IHS, col=ifelse(iHS.EAS$ihs$POSITION==109513601, "red", "black"), pch=19)
abline(h=c(2,-2), lty=2 )
abline(v=c(109500000,109605000), col=c("red", "red"), lty=c(2,2), lwd=c(1, 1))


#Make a windows-approach analysis
#1.	Create a function to estimate the mean in sliding windows.
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(abs(data[spots[i]:(spots[i] + window - 1)]),na.rm=TRUE)
  }
  return(result)
}

#2.	Estimate the mean over a window of 50 SNPs with steps of 40 SNPs in EAS.
mean_iHS <- slideFunct(iHS.EAS$ihs$IHS, 50,40)

#3.	Identify the starting position of each window
slidePos <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- data[spots[i]]
  }
  return(result)
}
pos_wind_Eas <- slidePos(iHS.EAS$ihs$POSITION, 50,40)

#4.	Put the position information and average iHS in a table
wind_iHS <- as.data.frame(cbind(pos_wind_Eas, mean_iHS), stringsAsFactors=FALSE)

#5.	Identify the window which contains the candidate SNP
Row_WIND_iHS <- wind_iHS[wind_iHS$pos_wind_Eas<=109513601,]
POS_WIND_iHS<-max(wind_iHS[nrow(Row_WIND_iHS),])
wind_iHS[wind_iHS$pos_wind_Eas==POS_WIND_iHS,]
##    pos_wind_Eas mean_iHS
## 21    109505388 1.146079

#6.	Plot the mean iHS per window
plot(ylim=c(0,1.5), x=wind_iHS[,1], y=wind_iHS[,2], xlab='pos', ylab='iHS windows', pch=20, cex=1.5)
points(x=wind_iHS[which(wind_iHS[,1]==POS_WIND_iHS),1],  y=wind_iHS[which(wind_iHS[,1]==POS_WIND_iHS),2], col='red', cex=2)


#6.	Check the distribution of iHS window in quantiles and check if the candidate SNP is an outlier.
windiHS_distrQT <- quantile(wind_iHS$mean_iHS, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99), na.rm=T)
windiHS_distrQT
##        1%        5%       10%       25%       50%       75%       90%       95% 
## 0.4240205 0.4766102 0.5071474 0.5913277 0.7159912 0.9820098 1.1314710 1.1701018 
##       99% 
## 1.3328222

#7.	Add the cut line for the quartile to the graph
plot(ylim=c(0,1.5), x=wind_iHS[,1], y=wind_iHS[,2], xlab='pos', ylab='iHS windows', pch=20, cex=1.5)
points(x=wind_iHS[which(wind_iHS[,1]==POS_WIND_iHS),1],  y=wind_iHS[which(wind_iHS[,1]==POS_WIND_iHS),2], col='red', cex=2)
abline(h=windiHS_distrQT[[7]], lty=2)

#1.	Calculate the xp-EHH between EAS e AFR
xpEHH.EAS.AFR<-ies2xpehh(EAS,AFR)

## Scan of pop1 contains 29016 markers.
## Scan of pop2 contains 29016 markers.
## Merged data contains 29016 markers.

#2.	Calculate the average xp-EHH per 50 SNP window with 40 SNP steps
mean_xpEHH <- slideFunct(xpEHH.EAS.AFR$XPEHH, 50,40)

#3.	Identify the starting position of each window
pos_wind_Eas <- slidePos(xpEHH.EAS.AFR$POSITION, 50,40)

#4.	Put the position information and average xpEHH in a table
wind_xpEHH <- as.data.frame(cbind(pos_wind_Eas, mean_xpEHH), stringsAsFactors=FALSE)

#5.	Identify the window which contains the candidate SNP
Row_WIND_xpEHH <- wind_xpEHH[wind_xpEHH$pos_wind_Eas<=109513601,]
POS_WIND_xpEHH<-max(wind_xpEHH[nrow(Row_WIND_xpEHH),])
wind_xpEHH[wind_xpEHH$pos_wind_Eas==POS_WIND_xpEHH,]
##     pos_wind_Eas mean_xpEHH
## 344    109512468   1.956037

#6.	Plot the mean xpEHH per window
plot(ylim=c(0,2.05), x=wind_xpEHH[,1], y=wind_xpEHH[,2], xlab='pos', ylab='xpEHH windows', pch=20, cex=1.5)
points(x=wind_xpEHH[which(wind_xpEHH[,1]==POS_WIND_xpEHH),1],  y=wind_xpEHH[which(wind_xpEHH[,1]==POS_WIND_xpEHH),2], col='red', cex=2)

#7.	Check the distribution of xpEHH window in quantiles and check if the candidate SNP is an outlier.
windxpEHH_distrQT <- quantile(wind_xpEHH$mean_xpEHH, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99), na.rm=T)
windxpEHH_distrQT
##         1%         5%        10%        25%        50%        75%        90% 
## 0.04627655 0.10333280 0.16789186 0.42277950 0.88493217 1.23910534 1.47871451 
##        95%        99% 
## 1.67418581 1.99188441

#8.	Add the cut line for the quartile to the graph and outline the candidate gene region
plot(ylim=c(0,2.05), x=wind_xpEHH[,1], y=wind_xpEHH[,2], xlab='pos', ylab='xpEHH windows', pch=20, cex=1.5)
points(x=wind_xpEHH[which(wind_xpEHH[,1]==POS_WIND_xpEHH),1],  y=wind_xpEHH[which(wind_xpEHH[,1]==POS_WIND_xpEHH),2], col='red', cex=2)

abline(h= windxpEHH_distrQT [[8]], lty=2)
abline(v=c(109500000,109605000), col="red")
