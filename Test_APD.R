# FILENAME: Analysis_APD.R

#Last updated 16 May 2024

#Usage: sources other files: 
#RScript Test_APD.R /Users/nikhilbhagwat/Desktop/7_PUBLICATIONS_ONGOING/7.1_APD_MS_2024-01/7.1.0_Data/APD_Neurology_MS_data.csv



cat("\n\n\n\n###############################################################################################\n",
            "SOURCE PACKAGES AND FUNCTIONS\n",
            "###############################################################################################\n\n\n")

arg.vec <- commandArgs(trailingOnly = T) #Defines a vector of arguments. In this case, it is only one argument.  

source('Functions.R') #Source the functions from file in same repository (for now)
	#tidyverse() loaded in Functions.R

source('DFCreationFunction_APD.R') #Source the function taht creates the new dataframe on which analysis is performed
	#lubridate() loaded in DFCreationFunction_APD.R


##LOAD LIBRARIES
library(car) #very useful for Anova(), Levene, etc. Please note that base R anova() and car:Anova() differ for:
	##contrast option already set as c("contr.sum", "contr.poly") for car::Anova(). Also, car::Anova() does type II SS
	##method as default, whereas base R anova() does type I SS method (same results as summary(aov))
library(emmeans) #estimated marginal means of a model
library(ggfortify) #plots for QC model
library(caret)
library(ggplot2) #plot figures with stats
library(ggpubr) #ggscatter
library(ggsurvfit) #survival analysis
library(survival) #survival analysis: survdiff() 
library(performance) #check_normality


# Consider removing if not used: 
# library(multcomp) #multiple comparisons
# library(lmtest) #Breusch-Pagan
# library(ggpmisc) #annotate figures
# library(modelbased) #extract some parameters and estimates from a model
# library(psych)
# library(rcompanion) #Cramer's V
# library(ggstatsplot) #plot figures
# library(fmsb) #radar plot



cat("\n\n\n\n###############################################################################################\n",
            "DATAFRAME LOADING\n",
            "###############################################################################################\n\n\n")

#Loading the QCed dataframe created from the raw data by the function var.func.1

df <- read.func(arg.vec[1], 'df') #arg.vec[1] is the csv that I entered data into. 
df <- var.func.1(df) #the df for analysis is created

write.csv(df, "dataframe.csv") #a copy of object created at the time of submission was kept for reference. Can be shared upon request. 



cat("\n\n\n\n###############################################################################################\n",
            "SUBSET DF FOR ASSUMPTION TESTING\n",
            "###############################################################################################\n\n\n")


##CREATE SOME OF THE SUBSETS USED LATER (MOSTLY FOR OUTLIER IDENTIFICATION)
#Not very efficient - update later if occasion to

CBSdf <- subset(df, DX_APD=="CBS") #for the description of AD+/- cohort
PSPdf <- subset(df, DX_APD=="PSP") #for the description of AD+/- cohort
RTposdf <- subset(df, RTQUIC=="aSyn-SAA positive") #for the main results on RT+/-
RTnegdf <- subset(df, RTQUIC=="aSyn-SAA negative") #for the main results on RT+/-
ADposdf <- subset(df, AD=="AD Positive") #for the description of AD+/- cohort
ADnegdf <- subset(df, AD=="AD Negative") #for the description of AD+/- cohort
APOEposdf <- subset(df, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf <- subset(df, APOEe4=="Negative") #for the description of APOE+/- cohort
YODdf <- subset(df, Early_onset=="Early-onset") #for the description of APOE+/- cohort
LODdf <- subset(df, Early_onset=="Late-onset") #for the description of APOE+/- cohort


# DEFENSE
if (sum(c(nrow(CBSdf), nrow(PSPdf), nrow(RTposdf), nrow(RTnegdf), nrow(ADposdf), nrow(ADnegdf), nrow(ADnegdf), nrow(APOEposdf), nrow(APOEnegdf), nrow(YODdf), nrow(LODdf)) == c(39, 28, 22, 45, 15, 52, 52, 14, 47, 35, 32))!= 11) {
	cat("Warning: potential error when subsetting based on diagnosis, AD, RTQUIC, APOE, or onset.")
}


cat("\n\n\n\n###########################################################################################################\n",
		   "COHORT CHARACTERISTICS\n",
	 	   "##########################################################################################################\n\n")


cat("COMPARISONS OF DEMOGRAPHICS FOR DX: \n")
df %>% count(DX_APD)


cat("\n\n#######################################################################################################\n",
	"                                    COHORT CHARACTERISTICS: NUMERICAL VARIABLES\n",
	    "#######################################################################################################\n")

cat("\n\n#######################################################################################################\n",
	    "########################################              AGE                     #########################\n",
	    "#######################################################################################################\n")


# AGE STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Age) #normal
shapiro.test(PSPdf$Age) #normal
var.test(Age ~ DX_APD, data = df) #homoscedasticity

# AGE STATISTICS: SUMMARY
cat("MEAN AGE AT LP (FOR TABLE): \n")
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))


# AGE STATISTICS: TTEST
t.test <-t.test(df$Age ~ df$DX_APD, var.equal=TRUE)
	if (t.test[3] <= 0.05) {
		cat("There is a significant difference in age at LP between CBS and PSP. p-value:", t.test[3][[1]], "\n")
		cat(t.test[3][[1]])
	} else cat("There is no significant difference in age at LP between CBS and PSP. \n")



cat("\n\n#######################################################################################################\n",
	   "#####################################              EDUCATION                  #########################\n",
	   "#######################################################################################################\n")


# EDUCATION STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Education) #not normal
shapiro.test(PSPdf$Education) #not normal
leveneTest(Education ~ DX_APD, data = df) #homoscedasticity

# EDUCATION STATISTICS: SUMMARY
cat("MEDIAN EDUCATION (FOR TABLE): \n")
df %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))

cat("MEAN EDUCATION (FOR REF ONLY): \n")
df %>% summarize(count=n(), format(round(mean(Education, na.rm=T),2),2), sd=sd(Education, na.rm=T))
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Education, na.rm=T),2),2), sd=sd(Education, na.rm=T))


# EDUCATION STATISTICS: WILCOX
wilcox.test(df$Education ~ df$DX_APD, paired=F)


cat("\n\n#######################################################################################################\n",
	   "################################              ONSET & DURATION                  #######################\n",
	   "#######################################################################################################\n")

# ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Onset_age) #normal
shapiro.test(PSPdf$Onset_age) #normal
leveneTest(Onset_age ~ DX_APD, data = df) #homoscedasticity

# ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))

# ONSET STATISTICS: TTEST
t.test(df$Onset_age ~ df$DX_APD, var.equal=TRUE) 
		

# PARK ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Park_onset) #borderline
shapiro.test(PSPdf$Park_onset) #borderline
leveneTest(Park_onset ~ DX_APD, data = df) #homoscedasticity

# PARK ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# ONSET STATISTICS: TTEST
t.test(df$Park_onset ~ df$DX_APD, var.equal=TRUE) 


# DURATION STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$LP2_Disease_Duration) #not normal
shapiro.test(PSPdf$LP2_Disease_Duration) #not normal
leveneTest(LP2_Disease_Duration ~ DX_APD, data = df) #homoscedasticity but borderline


# DURATION STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))
df %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))


# DURATION STATISTICS: WILCOX
wilcox.test(df$LP2_Disease_Duration ~ df$DX_APD, paired=F)
		


cat("\n\n######################################################################################################\n",
	   "################################              COGNITIVE Z-SCORES                 ######################\n",
	   "#######################################################################################################\n")


# MoCA Z-SCORE STATISTICS: DISTRIBUTION
boxplot(LP2_MOCA_Z.score ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(LP2_MOCA_Z.score ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(LP2_MOCA_Z.score ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(LP2_MOCA_Z.score ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(CBSdf$LP2_MOCA_Z.score) #normal
shapiro.test(PSPdf$LP2_MOCA_Z.score) #not normal
hist(df$LP2_MOCA_Z.score)
hist(CBSdf$LP2_MOCA_Z.score)
hist(PSPdf$LP2_MOCA_Z.score)
leveneTest(LP2_MOCA_Z.score ~ DX_APD, data = df) #heterodasticity


# MoCA Z-SCORE STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_MOCA_Z.score, na.rm=T),2),2), IQR=IQR(LP2_MOCA_Z.score, na.rm=T), min=min(LP2_MOCA_Z.score, na.rm=T), max=max(LP2_MOCA_Z.score, na.rm=T))
df %>% summarize(count=n(), format(round(median(LP2_MOCA_Z.score, na.rm=T),2),2), IQR=IQR(LP2_MOCA_Z.score, na.rm=T), min=min(LP2_MOCA_Z.score, na.rm=T), max=max(LP2_MOCA_Z.score, na.rm=T))

# MoCA Z-SCORE STATISTICS: WILCOX
wilcox.test(df$LP2_MOCA_Z.score ~ df$DX_APD, paired=F) #Since looking at z-score, no need to correct by age for this comparison

	# ALL COGNITIVE SCORES DISTRIBUTION/SUMMARY/STATISTICS: NO NEED IN TABLE 1 (REDUNDANT)
	shapiro.test(CBSdf$LP2_Cognitive_Z.score) #not al
	shapiro.test(PSPdf$LP2_Cognitive_Z.score) #not al
	leveneTest(LP2_Cognitive_Z.score ~ DX_APD, data = df) #heterodasticity
	wilcox.test(df$LP2_Cognitive_Z.score ~ df$DX_APD, paired=F)
	df %>% summarize(count=n(), format(round(median(LP2_Cognitive_Z.score, na.rm=T),2),2), IQR=IQR(LP2_Cognitive_Z.score, na.rm=T), min=min(LP2_Cognitive_Z.score, na.rm=T), max=max(LP2_Cognitive_Z.score, na.rm=T))
	df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_Cognitive_Z.score, na.rm=T),2),2), IQR=IQR(LP2_Cognitive_Z.score, na.rm=T), min=min(LP2_Cognitive_Z.score, na.rm=T), max=max(LP2_Cognitive_Z.score, na.rm=T))


cat("\n\n######################################################################################################\n",
	   "################################              BIOMARKERS: ABETA42                 #####################\n",
	   "#######################################################################################################\n")


# Abeta42 STATISTICS: DISTRIBUTION
# Be aware of the outliers but do not remove from descriptive Table 1 unless very problematic
# If removed, give value in the table
boxplot(abeta_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(abeta_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(abeta_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(abeta_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(CBSdf$logabeta) #normal
shapiro.test(PSPdf$logabeta) #normal
leveneTest(logabeta ~ DX_APD, data = df) #homoscedasticity

# Abeta42 STATISTICS: SUMMARY
df %>% summarize(count=n(), format(round(mean(abeta_2, na.rm=T),2),2), sd=sd(abeta_2, na.rm=T)) #Rounds up the sd for some reason
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(abeta_2, na.rm=T),2),2), sd=sd(abeta_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdf$abeta_2, na.rm=T)
sd(PSPdf$abeta_2, na.rm=T)

# Abeta42 STATISTICS: ANCOVA
t.test(df$logabeta ~ df$DX_APD, var.equal=TRUE) 
aov <- aov(logabeta ~ Age + DX_APD, df) 
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


cat("\n\n######################################################################################################\n",
	   "################################              BIOMARKERS: PTAU181                 #####################\n",
	   "#######################################################################################################\n")

# PTAU181 STATISTICS: DISTRIBUTION
# Be aware of the outliers but do not remove from descriptive Table 1 unless very problematic
# If removed, give value in the table
boxplot(ptau_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ptau_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(ptau_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ptau_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(CBSdf$logptau) #normal
shapiro.test(PSPdf$logptau) #normal
leveneTest(logptau ~ DX_APD, data = df) #homoscedasticity

# PTAU181 STATISTICS: SUMMARY
df %>% summarize(count=n(), format(round(mean(ptau_2, na.rm=T),2),2), sd=sd(ptau_2, na.rm=T)) #Rounds up the sd for some reason
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ptau_2, na.rm=T),2),2), sd=sd(ptau_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdf$ptau_2, na.rm=T)
sd(PSPdf$ptau_2, na.rm=T)

# PTAU181 STATISTICS: ANCOVA
t.test(df$logptau ~ df$DX_APD, var.equal=TRUE) 
aov <- aov(logptau ~ Age + DX_APD, df) 
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


cat("\n\n######################################################################################################\n",
	   "################################              BIOMARKERS: TTAU                 ########################\n",
	   "#######################################################################################################\n")


# T-TAU STATISTICS: DISTRIBUTION
# Be aware of the outliers but do not remove from descriptive Table 1 unless very problematic
# If removed, give value in the table
boxplot(ttau_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ttau_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(ttau_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ttau_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(CBSdf$logttau) #normal
shapiro.test(PSPdf$logttau) #normal
leveneTest(logttau ~ DX_APD, data = df) #homoscedasticity

# TAU STATISTICS: SUMMARY
# Chose mean because ran ANCOVA. 
df %>% summarize(count=n(), format(round(mean(ttau_2, na.rm=T),2),2), sd=sd(ttau_2, na.rm=T)) #Rounds up the sd for some reason
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ttau_2, na.rm=T),2),2), sd=sd(ttau_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdf$ttau_2, na.rm=T)
sd(PSPdf$ttau_2, na.rm=T)

# TAU STATISTICS: ANCOVA
t.test(df$logttau ~ df$DX_APD, var.equal=TRUE) 
aov <- aov(logttau ~ Age + DX_APD, df) 
Anova(aov, type="II") #Compare with type III
	check_normality(aov)

cat("\n\n######################################################################################################\n",
	   "################################              BIOMARKERS: ATI                 #########################\n",
	   "#######################################################################################################\n")


# ATI STATISTICS: DISTRIBUTION
# Be aware of the outliers but do not remove from descriptive Table 1 unless very problematic
# If removed, give value in the table
boxplot(ATI_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ATI_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(ATI_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ATI_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(CBSdf$ATI_2) #non-normal
shapiro.test(PSPdf$ATI_2) #non-normal
leveneTest(ATI_2 ~ DX_APD, data = df) #homoscedasticity
hist(df$ATI_2)

# ATI STATISTICS: SUMMARY
# Chose mean because ran ANCOVA. 
df %>% summarize(count=n(), format(round(mean(ATI_2, na.rm=T),2),2), sd=sd(ATI_2, na.rm=T)) #Rounds up the sd for some reason
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ATI_2, na.rm=T),2),2), sd=sd(ATI_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdf$ATI_2, na.rm=T)
sd(PSPdf$ATI_2, na.rm=T)

# ATI STATISTICS: ANCOVA
wilcox.test(df$ATI_2 ~ df$DX_APD, paired=F)
aov <- aov(ATI_2 ~ Age + DX_APD, df) 
Anova(aov, type="II") #Compare with type III
	check_normality(aov)



cat("\n\n######################################################################################################\n",
	   "################################              BIOMARKERS: NFL                 #########################\n",
	   "#######################################################################################################\n")

# NFL STATISTICS: DISTRIBUTION
# For NFL the outliers are really impactful. Is best to remove them and report their value. However, will only remove outliers of whole dataset and not per diagnosis. 
# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# This is because in FTLD population, the data can be very right-skewed. 
threshold <- min(max(df$logNFL,na.rm=T), as.numeric(quantile(df$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (df$logNFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
threshold
sort(df$logNFL)
subset(df, !(logNFL<threshold))$ID
subset(df, !(logNFL<threshold))$NFL
dfnfl <- subset(df, (logNFL<threshold))


CBSdfnfl <- subset(dfnfl, DX_APD=="CBS")
PSPdfnfl <- subset(dfnfl, DX_APD=="PSP")
shapiro.test(CBSdfnfl$logNFL) #normal
shapiro.test(PSPdfnfl$logNFL) #normal
shapiro.test(dfnfl$logNFL) #normal
leveneTest(logNFL ~ DX_APD, data = dfnfl) #homoscedasticity

# NFL STATISTICS: SUMMARY
dfnfl %>% summarize(count=n(), format(round(mean(NFL_2, na.rm=T),2),2), sd=sd(NFL_2, na.rm=T)) #Rounds up the sd for some reason
dfnfl %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(NFL_2, na.rm=T),2),2), sd=sd(NFL_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdfnfl$NFL_2)
sd(PSPdfnfl$NFL_2)


# NFL STATISTICS: ANCOVA
t.test(dfnfl$logNFL ~ dfnfl$DX_APD, var.equal=TRUE) 
aov <- aov(logNFL ~ Age + DX_APD, dfnfl) 
Anova(aov, type="II")
check_normality(aov)



cat("\n\n#######################################################################################################\n",
	"                                    COHORT CHARACTERISTICS: CATEGORICAL VARIABLES\n",
	    "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "################################              SEX                         #############################\n",
	   "#######################################################################################################\n")

# SEX STATISTICS: SUMMARY
cat("Sex distribution in the dataset is: \n")
df %>% group_by(DX_APD) %>% count(Sex) 


# SEX STATISTICS: PERCENTAGES
totalmatrix <- df %>% count(DX_APD)
sexmatrix <- df %>% group_by(DX_APD) %>% count(Sex) 

cat("Proportion of females in CBS is: \n")
(as.numeric(sexmatrix$n[1]) + as.numeric(sexmatrix$n[3]))/(as.numeric(totalmatrix$n[1]) + as.numeric(totalmatrix$n[2]))
cat("Proportion of females in CBS is:")
as.numeric(sexmatrix$n[1])/as.numeric(totalmatrix$n[1])
cat("Proportion of females in PSP is:")
as.numeric(sexmatrix$n[3])/as.numeric(totalmatrix$n[2])

# SEX STATISTICS: CHI-SQUARE
table(df$Sex, df$DX_APD)
chisq.test(table(df$Sex, df$DX_APD), correct=F)


cat("\n\n######################################################################################################\n",
	   "################################              APOE                         #############################\n",
	   "#######################################################################################################\n")

# APOE STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% count(APOEe4)
table(df$APOEe4, df$DX_APD)


# APOE STATISTICS: CHI-SQUARE
chisq.test(table(df$APOEe4, df$DX_APD), correct=F)


cat("\n\n######################################################################################################\n",
	   "################################              AD                         ##############################\n",
	   "#######################################################################################################\n")

# AD STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% count(AD)
table(df$DX_APD, df$AD)

# AD STATISTICS: CHI-SQUARE
chisq.test(table(df$AD, df$DX_APD), correct=F)



cat("\n\n\n\n###########################################################################################################\n",
		   "ASYN-SAA+ & DEMOGRAPHICS\n",
	 	   "##########################################################################################################\n\n")

cat("\n\n#######################################################################################################\n",
	"                                    COHORT CHARACTERISTICS: ASYN-SAA+ TOTAL NUMBER & DX\n",
	   "#######################################################################################################\n")

df %>% group_by(RTQUIC) %>% count(DX_APD)
sum(df$DX_APD=="CBS" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")
sum(df$DX_APD=="PSP" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")

cat("\n\n#######################################################################################################\n",
	"                                    COHORT CHARACTERISTICS: ASYN-SAA+ TOTAL NUMBER & AGE & ONSET \n",
	   "#######################################################################################################\n")


#BONFERRONI CALCULATION FOR AGE/ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO RTQUIC
#Looking at three similar tests: Difference in mean age at LP, onset, and Parkinsonism onset, in RT-QUIC+ vs RT-QUIC-
0.05/3 #0.01666667


cat("\n\n######################################################################################################\n",
	   "################################              RTQUIC*AGE                  #############################\n",
	   "#######################################################################################################\n")

# RTQUIC*AGE STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Age) #normal
shapiro.test(RTnegdf$Age) #normal
var.test(Age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*AGE STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))

# RTQUIC*AGE STATISTICS: ANCOVA
t.test(df$Age ~ df$RTQUIC, var.equal=TRUE)
aov <- aov(Age ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")



cat("\n\n######################################################################################################\n",
	   "################################              RTQUIC*ONSET                #############################\n",
	   "#######################################################################################################\n")

# RTQUIC*ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Onset_age) #normal
shapiro.test(RTnegdf$Onset_age) #normal
var.test(Onset_age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Onset_age, na.rm=T), sd=sd(Onset_age, na.rm=T))

# RTQUIC*ONSET STATISTICS: ANCOVA
t.test(df$Onset_age ~ df$RTQUIC, var.equal=TRUE)
aov <- aov(Onset_age ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")



cat("Below are the values for the creation of Fig 1.A: \n")
cat("Lower left quadrant: \n")
n1 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and young-onset and AD+: ", n1, "\n")
n2 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD+: ", n2, "\n")
n3 <-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and young-onset and AD+: ", n3, "\n")
n4 <-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD+: ", n4, "\n\n")


cat("Upper left quadrant: \n")
n5 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Early-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and young-onset and AD-: ", n5, "\n")
n6 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Early-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD-: ", n6, "\n")

n7<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and young-onset and AD-: ", n7, "\n")
n8<-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD-: ", n8, "\n\n")

cat("Lower right quadrant: \n")
n9<-nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and late-onset and AD+: ", n9, "\n")
n10<- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and late-onset and AD+: ", n10, "\n")

n11<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and late-onset and AD+: ", n11, "\n")
n12<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Early-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC- and late-onset and AD+: ", n12, "\n\n")

cat("Upper right quadrant: \n")
n13<- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and late-onset and AD-: ", n13, "\n")
n14<-nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and late-onset and AD-: ", n14, "\n")

n15<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Late-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and late-onset and AD-: ", n15, "\n")
n16<-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Late-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC- and late-onset and AD-: ", n16, "\n\n")

if (n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n12+n13+n14+n15+n16 != 67) {
	cat("Issue when calculating the values for Figure 1A. Check that there was no error in condition setup. \n")
	cat(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n12+n13+n14+n15+n16)
}


cat("\n\n######################################################################################################\n",
	   "################################              RTQUIC*PARK_ONSET           #############################\n",
	   "#######################################################################################################\n")

# RTQUIC*PARK_ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Park_onset) #normal
shapiro.test(RTnegdf$Park_onset) #normal
var.test(Park_onset ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*PARK_ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Park_onset, na.rm=T))

# RTQUIC*PARK_ONSET STATISTICS: ANCOVA
t.test(df$Park_onset ~ df$RTQUIC, var.equal=TRUE)
aov <- aov(Park_onset ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")

# cat("\n\n#######################################################################################################\n",
# 	"                                    COHORT CHARACTERISTICS: ASYN-SAA+ TOTAL NUMBER & AGE & ONSET \n",
# 	   "#######################################################################################################\n")
