install.packages("rehh")
library("rehh")
library(ggplot2)

setwd("~/embo_popgen_2025_Betty/EMBO_Practical_Course_2025_Betty")


#####LOAD DATA#####
AFR_chr2 <- data2haplohh("Chr2_EDAR_LWK_500K.recode.vcf", polarize_vcf = FALSE)
EAS_chr2 <- data2haplohh("Chr2_EDAR_CHS_500K.recode.vcf", polarize_vcf = FALSE)


#####CALCULATE EHH#####
AFR_EHH <- calc_ehh(AFR_chr2, mrk = "rs3827760", polarized = F)
plot(AFR_EHH)

EAS_EHH <- calc_ehh(EAS_chr2, mrk = "rs3827760")
plot(EAS_EHH)

#####CALCULATE FURCATION TREE#####

plot(calc_furcation(AFR_chr2, mrk = "rs3827760"))
plot(calc_furcation(EAS_chr2, mrk = "rs3827760"))

#####CALCULATE iHS#####

AFR_iHS <- scan_hh(AFR_chr2, polarized = T)
EAS_iHS <- scan_hh(EAS_chr2, polarized = T)
  #This time, we put true because we already specified F before and we don't want to "add up" falses

AFR_iHS["rs3827760",]

#CHR  POSITION FREQ_A FREQ_D NHAPLO_A NHAPLO_D    IHH_A IHH_D      IES     INES
#rs3827760   2 109513601      1      0      198        0 1767.692     0 1767.692 1767.692

EAS_iHS["rs3827760",]

#CHR  POSITION    FREQ_A    FREQ_D NHAPLO_A NHAPLO_D    IHH_A    IHH_D     IES    INES
#rs3827760   2 109513601 0.0952381 0.9047619       20      190 9979.001 55231.78 43767.3 54737.7


#####iHS#####

iHS_AFR <- ihh2ihs(AFR_iHS,
        min_maf = 0.02,
        freqbin = 0.01)
iHS_AFR[[1]]["rs3827760",]

#CHR POSITION IHS LOGPVALUE
#NA <NA>       NA  NA        NA

iHS_EAS <- ihh2ihs(EAS_iHS,
                   min_maf = 0.02,
                   freqbin = 0.01)
iHS_EAS[[1]]["rs3827760",]

#CHR  POSITION      IHS LOGPVALUE
#rs3827760   2 109513601 2.630072  2.068711

plot(iHS_EAS[[1]]$POSITION, iHS_EAS[[1]]$IHS, 
     col=ifelse(iHS_EAS[[1]]$POSITION=="109513601", 'red', 'darkolivegreen'),
     xlim = c(109400000,109600000))


#####SNP WINDOWS#####


#Make the windows
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(abs(data[spots[i]:(spots[i] + window - 1)]),na.rm=TRUE)
  }
  return(result)
}

mean_iHS <- slideFunct(iHS_EAS$ihs$IHS, 50, 40)


#Identify the starting pos of each window
slidePos <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- data[spots[i]]
  }
  return(result)
}

pos_wind_EAS <- slidePos(iHS_EAS$ihs$POSITION, 50, 40)

#Merge the data

matrix_wind_EAS <- cbind(pos_wind_EAS,mean_iHS)
matrix_wind_EAS <- as.data.frame(matrix_wind_EAS)

plot(matrix_wind_EAS$pos_wind_EAS, matrix_wind_EAS$mean_iHS)      


#####XP-EHH#####


XPEHH <- ies2xpehh(AFR_iHS,
                   EAS_iHS,
                   "AFR",
                   'EAS')
