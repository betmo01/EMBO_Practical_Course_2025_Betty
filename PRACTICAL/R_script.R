setwd("~/embo_popgen_2025_Betty/EMBO_Practical_Course_2025_Betty")

AFR_EAS_fst <- read.delim("~/embo_popgen_2025_Betty/EMBO_Practical_Course_2025_Betty/AFR_EAS.weir.fst")
AFR_EUR_fst <- read.delim("~/embo_popgen_2025_Betty/EMBO_Practical_Course_2025_Betty/AFR_EUR.weir.fst")
EAS_EUR_fst <- read.delim("~/embo_popgen_2025_Betty/EMBO_Practical_Course_2025_Betty/EAS_EUR.weir.fst")



#####CHECK DUPLICATES#####

for (i in c(1:(length(EAS_EUR_fst$POS)-1))) {
  if (unique(EAS_EUR_fst$POS[[i]])) {
    break
  }
  else  {
    print(i)
  }
}

    #No duplicates

#EAS_EUR_fst_filter <- EAS_EUR_fst[!duplicates(EAS_EUR_fst$POP,)]

#####CHECK DUPLICATES#####

for (i in c(1:(length(AFR_EAS_fst$POS)-1))) {
  if (unique(EAS_EUR_fst$POS[[i]])) {
    break
  }
  else  {
    print(i)
  }
}

#No duplicates

#####CHECK DUPLICATES#####

for (i in c(1:(length(AFR_EUR_fst$POS)-1))) {
  if (unique(EAS_EUR_fst$POS[[i]])) {
    break
  }
  else  {
    print(i)
  }
}

#No duplicates

#####EXCLUDE NAs#####

AFR_EUR_fst_na <- na.omit(AFR_EUR_fst)

AFR_EAS_fst_na <- na.omit(AFR_EAS_fst)

EAS_EUR_fst_na <- na.omit(EAS_EUR_fst)

#####KEEP ONLY OVERLAPING SNP #####

overlap_AFREAS_AFREUR <- AFR_EAS_fst_na[AFR_EAS_fst_na$POS %in% AFR_EUR_fst_na$POS,]
overlap_AFR_EUR_EAS <- overlap_AFREAS_AFREUR[overlap_AFREAS_AFREUR$POS %in% EAS_EUR_fst_na$POS,]

AFR_EAS_clean <- AFR_EAS_fst_na[AFR_EAS_fst_na$POS %in% overlap_AFR_EUR_EAS$POS,]
AFR_EUR_clean <- AFR_EUR_fst_na[AFR_EUR_fst_na$POS %in% overlap_AFR_EUR_EAS$POS,]
EAS_EUR_clean <- EAS_EUR_fst_na[EAS_EUR_fst_na$POS %in% overlap_AFR_EUR_EAS$POS,]


#####NEGATIVE =0 #####


AFR_EAS_clean[AFR_EAS_clean$WEIR_AND_COCKERHAM_FST<0, 'WEIR_AND_COCKERHAM_FST'] <- 0
AFR_EUR_clean[AFR_EUR_clean$WEIR_AND_COCKERHAM_FST<0, 'WEIR_AND_COCKERHAM_FST'] <- 0
EAS_EUR_clean[EAS_EUR_clean$WEIR_AND_COCKERHAM_FST<0, 'WEIR_AND_COCKERHAM_FST'] <- 0


PBS = (((-log(1 - EAS_EUR_clean$WEIR_AND_COCKERHAM_FST)) + (-log(1 - AFR_EAS_clean$WEIR_AND_COCKERHAM_FST)) - (-log(1 - AFR_EUR_clean$WEIR_AND_COCKERHAM_FST)))/2)
PBS
PBS[PBS<0] <- 0

df <- cbind(EAS_EUR_clean$POS, EAS_EUR_clean$WEIR_AND_COCKERHAM_FST,AFR_EUR_clean$WEIR_AND_COCKERHAM_FST,AFR_EAS_clean$WEIR_AND_COCKERHAM_FST,PBS)
colnames(df) <- c("POS","FST_AB", "FST_BC", "FST_AC", "PBS")
df$POS=="3827760"


