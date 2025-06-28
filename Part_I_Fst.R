input="/Users/patriciopezo/Desktop/EMBO_course_2025/input/"
out="/Users/patriciopezo/Desktop/EMBO_course_2025/data_process/"

#I. Read the files with the Fst estimates (AFR_EUR.weir.fst, AFR_EAS.weir.fst and EAS_EUR.weir.fst)
names_header <- c("CHROM","POS","WEIR_AND_COCKERHAM_FST","NUM","DEN")

FST_AFR_EAS <- read.table(paste0(input,"AFR_EAS.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,963 SNPs
FST_AFR_EUR <- read.table(paste0(input,"AFR_EUR.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,964 SNPs
FST_EAS_EUR <- read.table(paste0(input,"EAS_EUR.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,964 SNPs


#II. Remove duplicated positions
FST_AFR_EAS_filter <- FST_AFR_EAS[!duplicated(FST_AFR_EAS$POS),] #582,963 SNPs
FST_AFR_EUR_filter <- FST_AFR_EUR[!duplicated(FST_AFR_EUR$POS),] #582,964 SNPs
FST_EAS_EUR_filter <- FST_EAS_EUR[!duplicated(FST_EAS_EUR$POS),] #582,964 SNPs


#III. Take a look at the weir.fst file

head(FST_AFR_EAS_filter)
head(FST_AFR_EUR_filter)
head(FST_EAS_EUR_filter)


#IV. Exclude NAs position in Fst estimations
FST_AfrEas_data <- FST_AFR_EAS_filter[-which(is.na(FST_AFR_EAS_filter[,3])),] #582,694 SNPs
FST_AfrEur_data <- FST_AFR_EUR_filter[-which(is.na(FST_AFR_EUR_filter[,3])),] #582,363 SNPs
FST_EasEur_data <- FST_EAS_EUR_filter[-which(is.na(FST_EAS_EUR_filter[,3])),] #566,650 SNPs


#V. Overlaping SNPs
#565,779 SNPs

overlap_AfrEas_AfrEur <- FST_AfrEas_data[FST_AfrEas_data$POS %in% FST_AfrEur_data$POS,] 
overlap_AfrEasEur_EasEur <- overlap_AfrEas_AfrEur[overlap_AfrEas_AfrEur$POS %in% FST_EasEur_data$POS,]

FST_AfrEas_data_clean <- FST_AfrEas_data[FST_AfrEas_data$POS %in% overlap_AfrEasEur_EasEur$POS,]
FST_AfrEur_data_clean <- FST_AfrEur_data[FST_AfrEur_data$POS %in% overlap_AfrEasEur_EasEur$POS,]
FST_EasEur_data_clean <- FST_EasEur_data[FST_EasEur_data$POS %in% overlap_AfrEasEur_EasEur$POS,]

#VI. Convert negative values to zero
FST_AfrEas_data_clean[which(FST_AfrEas_data_clean[,3]<0),3] <- 0
FST_AfrEur_data_clean[which(FST_AfrEur_data_clean[,3]<0),3] <- 0
FST_EasEur_data_clean[which(FST_EasEur_data_clean[,3]<0),3] <- 0

#VII. Check if the SNP rs3827760 (pos 109513601) is a candidate for natural selection
#1. Check if the SNP rs3827760, located at position 109513601, is an outlier in the FST distribution for any of the population pairs

POS <- 109513601

FST_AfrEas_data_clean[FST_AfrEas_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST      NUM      DEN
## 266093     2 109513601               0.872881 0.762039 0.873016

FST_AfrEur_data_clean[FST_AfrEur_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST         NUM       DEN
## 266093     2 109513601             0.00997189 0.000108929 0.0109237

FST_EasEur_data_clean[FST_EasEur_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST      NUM      DEN
## 266093     2 109513601               0.859066 0.743056 0.864958


#2. Check which quartile percentile the rs3827760 distribution fall in each analyzed population pair?
FST_AfrEas_distr <- sort(FST_AfrEas_data_clean[,3])
FST_AfrEur_distr <- sort(FST_AfrEur_data_clean[,3])
FST_EasEur_distr <- sort(FST_EasEur_data_clean[,3])

FST_AfrEas_distrQT <- quantile(FST_AfrEas_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))

FST_AfrEas_distrQT
##           1%           5%          10%          25%          50%          75% 
## 0.0000000000 0.0000079206 0.0031241100 0.0262112250 0.1056270000 0.2222135000 
##          90%          95%          99% 
## 0.3675160000 0.4685402500 0.6521344200



FST_AfrEur_distrQT <- quantile(FST_AfrEur_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
FST_AfrEur_distrQT
##          1%          5%         10%         25%         50%         75% 
## 0.000000000 0.000000000 0.002228373 0.018979525 0.081086300 0.186617500 
##         90%         95%         99% 
## 0.308151200 0.395242000 0.551998730



FST_EasEur_distrQT <- quantile(FST_EasEur_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
FST_EasEur_distrQT
##           1%           5%          10%          25%          50%          75% 
## 0.0000000000 0.0000000000 0.0004382527 0.0092101950 0.0476573500 0.1294852500 
##          90%          95%          99% 
## 0.2363264000 0.3175715000 0.4790830000




#3. Ploting FST values in a 10,000 base pair region adjacent to the SNP at position 109513601. Highlight the SNPs that are outliers in the 95th percentile in each population pair.


#Delimit the region of interest to 10000bp adjacent
SNPfrom_BP <- POS - 10000
SNPto_BP <- POS + 10000


#How many SNPs are in this bounded region?

SNPfrom_id_AfrEas <- max(which(FST_AfrEas_data_clean[,2]<=SNPfrom_BP))
SNPto_id_Afr_Eas <- min(which(FST_AfrEas_data_clean[,2]>=SNPto_BP))
length(FST_AfrEas_data_clean[SNPfrom_BP:SNPto_BP, 2])
## [1] 20001


SNPfrom_id_AfrEur <- max(which(FST_AfrEur_data_clean[,2]<=SNPfrom_BP))
SNPto_id_Afr_Eur <- min(which(FST_AfrEur_data_clean[,2]>=SNPto_BP))
length(FST_AfrEur_data_clean[SNPfrom_BP:SNPto_BP, 2])
## [1] 20001


SNPfrom_id_EasEur <- max(which(FST_EasEur_data_clean[,2]<=SNPfrom_BP))
SNPto_id_EasEur <- min(which(FST_EasEur_data_clean[,2]>=SNPto_BP))
length(FST_EasEur_data_clean[SNPfrom_BP:SNPto_BP, 2])
## [1] 20001


#Select from the FST_AfrEur_data the region of interest
FSTdata_SNP_AfrEas <- FST_AfrEas_data_clean[SNPfrom_id_AfrEas:SNPto_id_Afr_Eas, ]
head(FSTdata_SNP_AfrEas)
##        CHROM       POS WEIR_AND_COCKERHAM_FST        NUM       DEN
## 266063     2 109503245              0.0826870 0.00859760 0.1039780
## 266064     2 109503631              0.0747407 0.00718525 0.0961357
## 266066     2 109503778              0.0668349 0.00590110 0.0882937
## 266067     2 109503862              0.1239490 0.02692580 0.2172320
## 266068     2 109504022              0.0448158 0.00616223 0.1375010
## 266069     2 109504287              0.2989630 0.12097500 0.4046480



FSTdata_SNP_AfrEur <- FST_AfrEur_data_clean[SNPfrom_id_AfrEur:SNPto_id_Afr_Eur, ]
FSTdata_SNP_EasEur <- FST_EasEur_data_clean[SNPfrom_id_EasEur:SNPto_id_EasEur, ]



#4. PLOT
#a. AfrEas
plot(ylim=c(0,1), x=FSTdata_SNP_AfrEas[,2], y=FSTdata_SNP_AfrEas[,3], xlab='pos', ylab='FST AFR EAS', pch=20, cex=1.5)

points(x=FSTdata_SNP_AfrEas[which(FSTdata_SNP_AfrEas[,2]==POS),2],  y=FSTdata_SNP_AfrEas[which(FSTdata_SNP_AfrEas[,2]==POS),3], col='blue', cex=2)

abline(h=FST_AfrEas_distrQT[[8]], lty=2)

#b.	AfrEur
plot(ylim=c(0,1), x=FSTdata_SNP_AfrEur[,2], y=FSTdata_SNP_AfrEur[,3], xlab='pos', ylab='FST AFR EUR', pch=20, cex=1.5)

points(x=FSTdata_SNP_AfrEur[which(FSTdata_SNP_AfrEur[,2]==POS),2],  y=FSTdata_SNP_AfrEur[which(FSTdata_SNP_AfrEur[,2]==POS),3], col='blue', cex=2)

abline(h=FST_AfrEur_distrQT[[8]], lty=2)

#c. EasEur
plot(ylim=c(0,1), x=FSTdata_SNP_EasEur[,2], y=FSTdata_SNP_EasEur[,3], xlab='pos', ylab='FST EAS EUR', pch=20, cex=1.5)

points(x=FSTdata_SNP_EasEur[which(FSTdata_SNP_EasEur[,2]==POS),2],  y=FSTdata_SNP_EasEur[which(FSTdata_SNP_EasEur[,2]==POS),3], col='blue', cex=2)

abline(h=FST_EasEur_distrQT[[8]], lty=2)

#VIII. Can the candidate SNP be considered an outlier in all populations? What is the interpretation of this result?
#1. Estimate the p-value for the candidate SNP from the distribution of FST values for each population pair
#	AfrEas
p_value_out_FST_AfrEas <- sum(FST_AfrEas_data_clean$WEIR_AND_COCKERHAM_FST>=FST_AfrEas_data_clean[FST_AfrEas_data_clean$POS==109513601,3])/nrow(FST_AfrEas_data_clean)
p_value_out_FST_AfrEas #0.0004772181


#AfrEur
p_value_out_FST_AfrEur <- sum(FST_AfrEur_data_clean$WEIR_AND_COCKERHAM_FST>=FST_AfrEur_data_clean[FST_AfrEur_data_clean$POS==109513601,3])/nrow(FST_AfrEur_data_clean)
p_value_out_FST_AfrEur #0.8143781


#EurEas
p_value_out_FST_EasEur <- sum(FST_EasEur_data_clean$WEIR_AND_COCKERHAM_FST>=FST_EasEur_data_clean[FST_EasEur_data_clean$POS==109513601,3])/nrow(FST_EasEur_data_clean)
p_value_out_FST_EasEur #5.302424e-06






#Path
input="/Users/patriciopezo/Desktop/EMBO_course_2025/input/"
out="/Users/patriciopezo/Desktop/EMBO_course_2025/data_process/"


#1. Perform PBS test, using EAS as candidate population for selection
#Build the PBS Topology and why is it important to measure the distance
PBS_EAS <- ((-log(1-FST_AfrEas_data_clean$WEIR_AND_COCKERHAM_FST))+(-log(1-FST_EasEur_data_clean$WEIR_AND_COCKERHAM_FST))-(-log(1-FST_AfrEur_data_clean$WEIR_AND_COCKERHAM_FST)))/2


#2. Convert negative PBS values to O
PBS_EAS[which(PBS_EAS<0)] <- 0



#3. Add to the data.table with FST values, a new column with PBS values
fst_pbs<-as.data.frame(cbind(FST_EasEur_data_clean$POS, FST_AfrEas_data_clean$WEIR_AND_COCKERHAM_FST, FST_AfrEur_data_clean$WEIR_AND_COCKERHAM_FST, FST_EasEur_data_clean$WEIR_AND_COCKERHAM_FST, PBS_EAS), stringsAsFactors=FALSE)

head(fst_pbs)
##      V1       V2       V3          V4      PBS_EAS
## 1 10554 0.318827 0.334988 0.00000e+00 0.0000000000
## 2 10560 0.317773 0.333938 0.00000e+00 0.0000000000
## 3 10566 0.315969 0.333938 2.24832e-06 0.0000000000
## 4 10574 0.116785 0.132890 4.42702e-04 0.0000000000
## 5 10587 0.368198 0.355929 0.00000e+00 0.0096164573
## 6 10595 0.205058 0.204892 0.00000e+00 0.0001043992



#4. Check the PBS value for the candidate SNP.
pbs_EDAR<-fst_pbs[fst_pbs$V1==POS,]

pbs_EDAR
##               V1       V2         V3       V4  PBS_EAS
## 259224 109513601 0.872881 0.00997189 0.859066 2.006037




#5.	In which quartile of the distribution does the PBS value for the SNP rs3827760 fall?
PBS_distrQT <- quantile(PBS_EAS, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
PBS_distrQT # estimar os quantils
##         1%         5%        10%        25%        50%        75%        90% 
## 0.00000000 0.00000000 0.00000000 0.00000000 0.02506831 0.10666149 0.23044290 
##        95%        99% 
## 0.33611959 0.60783061



#6.	Plot the PBS values in a 10,000 base pair region adjacent to the SNP at position 109513601. 
#Highlight the SNPs that are outliers in the 95th percentile.
#Select 10000bp adjacent to candidate SNP
SNP_FROM <- POS - 10000
SNP_TO <- POS + 10000

SNPfrom_PBS <- max(which(fst_pbs[,1]<=SNP_FROM))
SNPto_PBS <- min(which(fst_pbs[,1]>=SNP_TO))
length(fst_pbs[SNPfrom_PBS:SNPto_PBS, 1])
## [1] 55




#7. Subset the candidate SNP region
subset_fst_PBS <- fst_pbs[SNPfrom_PBS:SNPto_PBS, ]
head(subset_fst_PBS)
##               V1        V2        V3         V4    PBS_EAS
## 259198 109503245 0.0826870 0.0000000 0.07902310 0.08431343
## 259199 109503631 0.0747407 0.0000000 0.07902310 0.08000079
## 259200 109503778 0.0668349 0.0000000 0.07902310 0.07574673
## 259201 109503862 0.1239490 0.0803297 0.00715918 0.02788793
## 259202 109504022 0.0448158 0.0171244 0.00664388 0.01762220
## 259203 109504287 0.2989630 0.2609700 0.00305409 0.02791831




#8. subset_fst_PBS[subset_fst_PBS$V1==109513601,]
##               V1       V2         V3       V4  PBS_EAS
## 259224 109513601 0.872881 0.00997189 0.859066 2.006037



#9. Plot PBS values
plot(ylim=c(0,2.5), x=subset_fst_PBS[,1], y=subset_fst_PBS[,5], xlab='pos', ylab='PBS', pch=20, cex=1.5)

points(x=subset_fst_PBS[which(subset_fst_PBS$V1==POS),1],  y=subset_fst_PBS[which(subset_fst_PBS$V1==POS),5], col='blue', cex=2)

abline(h=PBS_distrQT[[9]], lty=2)
