# FILENAME: Analysis_APD.R

#Last updated May 2024

#Usage: sources other files: 
#RScript Analysis_APD.R /Users/nikhilbhagwat/Desktop/7_PUBLICATIONS_ONGOING/7.1_APD_MS_2024-01/7.1.0_Data/APD_Neurology_MS_data.csv



cat("\n\n\n\n###############################################################################################\n",
            "1. SOURCE PACKAGES AND FUNCTIONS\n",
            "###############################################################################################\n\n\n")

arg.vec <- commandArgs(trailingOnly = T) #Defines a vector of arguments. In this case, it is only one argument.  

source('Functions.R') #Source the functions from file in same repository (for now)
	#tidyverse() loaded in Functions.R

source('DataQC_PROCESSING.R') #Source the function taht creates the new dataframe on which analysis is performed
	#lubridate() loaded in DataQC_PROCESSING.R


## LOAD LIBRARIES
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
library(survminer) #ggsurvplot
library(performance) #check_normality
library(rcompanion) #Cramer's V
library(janitor) #remove empty rows
library(nlme) #weighted least square regression
library(lmtest) #bptest()
library(fmsb) #create_beautiful_radarchart()
library(pscl) #McFadden index
library(lsr) #eta-squared aov()
library(ggsignif) #annotations figures with stats

# Consider removing if not used: 
# library(multcomp) #multiple comparisons
# library(ggpmisc) #annotate figures
# library(modelbased) #extract some parameters and estimates from a model
# library(psych)
# library(ggstatsplot) #plot figures



cat("\n\n\n\n###############################################################################################\n",
            "2. DATAFRAME LOADING\n",
            "###############################################################################################\n\n\n")

#Loading the QCed dataframe created from the raw data by the function var.func.1

df <- read.func(arg.vec[1], 'df') #arg.vec[1] is the csv that I entered data into. 
df <- var.func.1(df) #the df for analysis is created

write.csv(df, "dataframe.csv") #a copy of object created at the time of submission was kept for reference. Can be shared upon request. 



cat("\n\n\n\n###############################################################################################\n",
            "3. SUBSET DF FOR ASSUMPTION TESTING\n",
            "###############################################################################################\n\n\n")


## CREATE SOME OF THE SUBSETS USED LATER (MOSTLY FOR OUTLIER IDENTIFICATION)
CBSdf <- subset(df, DX_APD=="CBS") #for the description of AD+/- cohort
PSPdf <- subset(df, DX_APD=="PSP") #for the description of AD+/- cohort
RTposdf <- subset(df, RTQUIC=="aSyn-SAA positive") #for the main results on RT+/-
RTnegdf <- subset(df, RTQUIC=="aSyn-SAA negative") #for the main results on RT+/-
ADposdf <- subset(df, AD=="AD Positive") #for the description of AD+/- cohort
ADnegdf <- subset(df, AD=="AD Negative") #for the description of AD+/- cohort
APOEposdf <- subset(df, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf <- subset(df, APOEe4=="Negative") #for the description of APOE+/- cohort
YODdf <- subset(df, Early_onset=="Young-onset") #for the description of APOE+/- cohort
LODdf <- subset(df, Early_onset=="Late-onset") #for the description of APOE+/- cohort


# DEFENSE
if (sum(c(nrow(CBSdf), nrow(PSPdf), nrow(RTposdf), nrow(RTnegdf), nrow(ADposdf), nrow(ADnegdf), nrow(ADnegdf), nrow(APOEposdf), nrow(APOEnegdf), nrow(YODdf), nrow(LODdf)) == c(39, 28, 22, 45, 15, 52, 52, 14, 47, 35, 32))!= 11) {
	cat("Warning: potential error when subsetting based on diagnosis, AD, RTQUIC, APOE, or onset.")
}


cat("\n\n\n\n####################################################################################################\n",
		   "4. COHORT CHARACTERISTICS\n",
	 	   "#####################################################################################################\n\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

cat("COMPARISONS OF DEMOGRAPHICS FOR DX: \n")
df %>% count(DX_APD)


cat("\n\n#######################################################################################################\n",
	"                                    4.1. COHORT CHARACTERISTICS: NUMERICAL VARIABLES\n",
	    "#######################################################################################################\n")

cat("\n\n#######################################################################################################\n",
	    "####################################             4.1.1. AGE                   #########################\n",
	    "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

# AGE STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Age) #normal
shapiro.test(PSPdf$Age) #normal
hist(CBSdf$Age)
hist(PSPdf$Age)
var.test(Age ~ DX_APD, data = df) #homoscedasticity

# AGE STATISTICS: SUMMARY
cat("MEAN AGE AT LP (FOR TABLE): \n")
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))


# AGE STATISTICS: STUDENTS TTEST
t.test <-t.test(df$Age ~ df$DX_APD, var.equal=TRUE) #var.equal = Student's t-test is performed 
	if (t.test[3] <= 0.05) {
		cat("There is a significant difference in age at LP between CBS and PSP. p-value:", t.test[3][[1]], "\n")
		cat(t.test[3][[1]])
	} else cat("There is no significant difference in age at LP between CBS and PSP. \n")



cat("\n\n#######################################################################################################\n",
	   "##################################              4.1.2. EDUCATION              #########################\n",
	   "#######################################################################################################\n\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# EDUCATION STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Education) #not normal
shapiro.test(PSPdf$Education) #not normal
leveneTest(Education ~ DX_APD, data = df) #homoscedasticity

# EDUCATION STATISTICS: SUMMARY
cat("MEDIAN EDUCATION (FOR TABLE): \n")
df %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))

# EDUCATION STATISTICS: WILCOXON
wilcox.test(df$Education ~ df$DX_APD, paired=F)


cat("\n\n#######################################################################################################\n",
	   "###########################            4.1.3. ONSET & DURATION                  #######################\n",
	   "#######################################################################################################\n\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

# ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Onset_age) #normal
shapiro.test(PSPdf$Onset_age) #normal
hist(CBSdf$Onset_age)
hist(PSPdf$Onset_age)
leveneTest(Onset_age ~ DX_APD, data = df) #homoscedasticity

# ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))

# ONSET STATISTICS: STUDENTS TTEST
t.test(df$Onset_age ~ df$DX_APD, var.equal=TRUE) 
		

# PARK ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$Park_onset) #borderline
shapiro.test(PSPdf$Park_onset) #borderline
hist(CBSdf$Park_onset)
hist(PSPdf$Park_onset)
leveneTest(Park_onset ~ DX_APD, data = df) #homoscedasticity

# PARK ONSET STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# PARK ONSET STATISTICS: STUDENTS TTEST
t.test(df$Park_onset ~ df$DX_APD, var.equal=TRUE) 


# DURATION STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$LP2_Disease_Duration) #not normal
shapiro.test(PSPdf$LP2_Disease_Duration) #not normal
leveneTest(LP2_Disease_Duration ~ DX_APD, data = df) #homoscedasticity but borderline

# DURATION STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))
df %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))


# DURATION STATISTICS: WILCOXON
wilcox.test(df$LP2_Disease_Duration ~ df$DX_APD, paired=F)
		

cat("\n\n######################################################################################################\n",
	   "#########################              4.1.4. COGNITIVE Z-SCORES                 ######################\n",
	   "#######################################################################################################\n\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

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

# MoCA Z-SCORE STATISTICS: WILCOXON
wilcox.test(df$LP2_MOCA_Z.score ~ df$DX_APD, paired=F) #Since looking at z-score, no need to correct by age for this comparison

	# ALL COGNITIVE SCORES DISTRIBUTION/SUMMARY/STATISTICS: NO NEED IN TABLE 1 (REDUNDANT)
	shapiro.test(CBSdf$LP2_Cognitive_Z.score) #not normal
	shapiro.test(PSPdf$LP2_Cognitive_Z.score) #not normal
	leveneTest(LP2_Cognitive_Z.score ~ DX_APD, data = df) #heterodasticity
	wilcox.test(df$LP2_Cognitive_Z.score ~ df$DX_APD, paired=F)
	df %>% summarize(count=n(), format(round(median(LP2_Cognitive_Z.score, na.rm=T),2),2), IQR=IQR(LP2_Cognitive_Z.score, na.rm=T), min=min(LP2_Cognitive_Z.score, na.rm=T), max=max(LP2_Cognitive_Z.score, na.rm=T))
	df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_Cognitive_Z.score, na.rm=T),2),2), IQR=IQR(LP2_Cognitive_Z.score, na.rm=T), min=min(LP2_Cognitive_Z.score, na.rm=T), max=max(LP2_Cognitive_Z.score, na.rm=T))


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.5. BIOMARKERS: ABETA42                 #####################\n",
	   "#######################################################################################################\n")

## Biomarkers: 
# The goal here is to share data that can be compared to other datasets as reported values for CBS and FTLD can have a wide range. 
# Therefeore, mean of raw values are presented here, even if analyses are on logged values
# Mean and sd are calculated after removing outliers, calculated based on Tukey method (Q1 - 1.5IQR or Q3 + 1.5IQR)
# Exception for NfL which has a strong positive skew so instead the threshold is calculated by 3*IQR (Delaby et al 2020 on similar cohort)
# Outliers are shown in notes of Tables every time. 
# Analyses were done without outlier exclusion as they incorporated multiple variables and were complex models so given the heterogeneity of the dataset, 
## it would lead to groups with small sample sizes being possibly under-represented. Instead, models were run on dataset minus the outlier, 1 outlier at a time. 
# This allowed us to see how different outliers could affect the results of the main model selected. 


cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# Abeta42 STATISTICS: OUTLIERS
# First, display the data points of one group. boxplot() will identify outliers in given df but stripchart function can be funky hence convoluted script. 
# Code below plots (boxplot() + assigns outliers to diagnosis-specific vector)
# Outliers are removed before calculating the values that are presented in table. 
boxplot(abeta_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(abeta_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	CBSvec.outliers <- boxplot(abeta_2 ~ DX_APD, data= CBSdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
boxplot(abeta_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(abeta_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	PSPvec.outliers <- boxplot(abeta_2 ~ DX_APD, data= PSPdf, col = "white")$out 

dfabeta<- df[!(df$DX_APD=="CBS" & (df$abeta_2 %in% CBSvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfabeta<- dfabeta[!(dfabeta$DX_APD=="PSP" & (dfabeta$abeta_2 %in% PSPvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfabeta[, c("DX_APD", "abeta_2")]

if (nrow(dfabeta) == nrow(df)) {cat("No outliers were removed for Abeta42 comparisons \n")} else {cat("The following outliers were removed from the CBS:", CBSvec.outliers," from the PSP:", PSPvec.outliers, " .\n")}

# Abeta42 STATISTICS: DISTRIBUTION
shapiro.test(dfabeta[(dfabeta$DX_APD == "CBS"), ]$logabeta) #normal
shapiro.test(dfabeta[(dfabeta$DX_APD == "PSP"), ]$logabeta) #normal
hist(df$logabeta)
hist(CBSdf$logabeta)
hist(PSPdf$logabeta)
leveneTest(logabeta ~ DX_APD, data = dfabeta) #homoscedasticity

# Abeta42 STATISTICS: SUMMARY
dfabeta %>% summarize(count=n(), format(round(mean(abeta_2, na.rm=T),2),2), sd=sd(abeta_2, na.rm=T)) #Rounds up the sd for some reason
dfabeta %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(abeta_2, na.rm=T),2),2), sd=sd(abeta_2, na.rm=T)) #Rounds up the sd for some reason
sd(dfabeta[(dfabeta$DX_APD == "CBS"), ]$abeta_2, na.rm=T)
sd(PSPdf$abeta_2, na.rm=T)

# Abeta42 STATISTICS: ANCOVA
t.test(df$logabeta ~ df$DX_APD, var.equal=TRUE) #before moving to more complex model, run simple t-test
aov <- aov(logabeta ~ Age + DX_APD, df) 
Anova(aov, type="II") #no interaction in the model so choosing type II
	check_normality(aov)
etaSquared(aov)


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.6. BIOMARKERS: PTAU181                 #####################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# PTAU181 STATISTICS: OUTLIERS
# First, display the data points of one group. boxplot() will identify outliers in given df but stripchart function can be funky hence convoluted script. 
# Code below plots (boxplot() + assigns outliers to diagnosis-specific vector)
# Outliers are removed before calculating the values that are presented in table. 
boxplot(ptau_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ptau_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	CBSvec.outliers <- boxplot(ptau_2 ~ DX_APD, data= CBSdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
boxplot(ptau_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ptau_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	PSPvec.outliers <- boxplot(ptau_2 ~ DX_APD, data= PSPdf, col = "white")$out 

dfptau<- df[!(df$DX_APD=="CBS" & (df$ptau_2 %in% CBSvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfptau<- dfptau[!(dfptau$DX_APD=="PSP" & (dfptau$ptau_2 %in% PSPvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold

if (nrow(dfptau) == nrow(df)) {cat("No outliers were removed for ptau181 comparisons \n")} else {cat("The following outliers were removed from the CBS: ", CBSvec.outliers," from the PSP: ", PSPvec.outliers, " .\n")}


# PTAU181 STATISTICS: DISTRIBUTION
shapiro.test(dfptau[(dfptau$DX_APD == "CBS"), ]$logptau) #normal
shapiro.test(PSPdf$logptau) #normal
hist(df$logptau)
hist(CBSdf$logptau)
hist(PSPdf$logptau)
leveneTest(logptau ~ DX_APD, data = df) #homoscedasticity

# PTAU181 STATISTICS: SUMMARY
dfptau %>% summarize(count=n(), format(round(mean(ptau_2, na.rm=T),2),2), sd=sd(ptau_2, na.rm=T)) #Rounds up the sd for some reason
dfptau %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ptau_2, na.rm=T),2),2), sd=sd(ptau_2, na.rm=T)) #Rounds up the sd for some reason
sd(dfptau[(dfptau$DX_APD == "CBS"), ]$ptau_2, na.rm=T)
sd(PSPdf$ptau_2, na.rm=T)

# PTAU181 STATISTICS: ANCOVA
t.test(dfptau$logptau ~ dfptau$DX_APD, var.equal=TRUE) 
aov <- aov(logptau ~ Age + DX_APD, dfptau) 
Anova(aov, type="II") #no interaction in the model
	check_normality(aov)
etaSquared(aov)


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.7. BIOMARKERS: TTAU                 ########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# T-TAU STATISTICS: OUTLIERS
# First, display the data points of one group. boxplot() will identify outliers in given df but stripchart function can be funky hence convoluted script. 
# Code below plots (boxplot() + assigns outliers to diagnosis-specific vector)
# Outliers are removed before calculating the values that are presented in table. 
boxplot(ttau_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ttau_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	CBSvec.outliers <- boxplot(ttau_2 ~ DX_APD, data= CBSdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
boxplot(ttau_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ttau_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	PSPvec.outliers <- boxplot(ttau_2 ~ DX_APD, data= PSPdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)

dfttau<- df[!(df$DX_APD=="CBS" & (df$ttau_2 %in% CBSvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfttau<- dfttau[!(dfttau$DX_APD=="PSP" & (dfttau$ttau_2 %in% PSPvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold

if (nrow(dfttau) == nrow(df)) {cat("No outliers were removed for t-tau comparisons \n")} else {cat("The following outliers were removed from the CBS: ", CBSvec.outliers," from the PSP: ", PSPvec.outliers, " .\n")}

# T-TAU STATISTICS: DISTRIBUTION
shapiro.test(dfttau[(dfttau$DX_APD == "CBS"), ]$logttau) #normal
shapiro.test(dfttau[(dfttau$DX_APD == "PSP"), ]$logttau) #borderline
hist(df$logttau)
hist(CBSdf$logttau)
hist(PSPdf$logttau)
leveneTest(logttau ~ DX_APD, data = df) #homoscedasticity

# TAU STATISTICS: SUMMARY
# Chose mean because ran ANCOVA. 
dfttau %>% summarize(count=n(), format(round(mean(ttau_2, na.rm=T),2),2), sd=sd(ttau_2, na.rm=T)) #Rounds up the sd for some reason
dfttau %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ttau_2, na.rm=T),2),2), sd=sd(ttau_2, na.rm=T)) #Rounds up the sd for some reason
sd(dfttau[(dfttau$DX_APD == "CBS"), ]$ttau_2, na.rm=T)
sd(dfttau[(dfttau$DX_APD == "PSP"), ]$ttau_2, na.rm=T)

# TAU STATISTICS: ANCOVA
t.test(df$logttau ~ df$DX_APD, var.equal=TRUE) 
aov <- aov(logttau ~ Age + DX_APD, df) 
Anova(aov, type="II") 
	check_normality(aov)
etaSquared(aov)


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.8. BIOMARKERS: ATI                 #########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")


# ATI STATISTICS: OUTLIERS
# First, display the data points of one group. boxplot() will identify outliers in given df but stripchart function can be funky hence convoluted script. 
# Code below plots (boxplot() + assigns outliers to diagnosis-specific vector)
# Outliers are removed before calculating the values that are presented in table. 
boxplot(ATI_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(ATI_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	CBSvec.outliers <- boxplot(ATI_2 ~ DX_APD, data= CBSdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
boxplot(ATI_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(ATI_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	PSPvec.outliers <- boxplot(ATI_2 ~ DX_APD, data= PSPdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)

dfati<- df[!(df$DX_APD=="CBS" & (df$ATI_2 %in% CBSvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfati<- dfati[!(dfati$DX_APD=="PSP" & (dfati$ATI_2 %in% PSPvec.outliers)), ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold

if (nrow(dfati) == nrow(df)) {cat("No outliers were removed for ATI comparisons \n")} else {cat("The following outliers were removed from the CBS: ", CBSvec.outliers," from the PSP: ", PSPvec.outliers, " .\n")}

# ATI STATISTICS: DISTRIBUTION
shapiro.test(CBSdf$ATI_2) #non-normal
shapiro.test(PSPdf$ATI_2) #non-normal
leveneTest(ATI_2 ~ DX_APD, data = df) #homoscedasticity
hist(df$ATI_2)
hist(CBSdf$ATI_2)
hist(PSPdf$ATI_2)


# ATI STATISTICS: SUMMARY
# Chose mean because ran ANCOVA. 
df %>% summarize(count=n(), format(round(mean(ATI_2, na.rm=T),2),2), sd=sd(ATI_2, na.rm=T)) #Rounds up the sd for some reason
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ATI_2, na.rm=T),2),2), sd=sd(ATI_2, na.rm=T)) #Rounds up the sd for some reason
sd(CBSdf$ATI_2, na.rm=T)
sd(PSPdf$ATI_2, na.rm=T)

# ATI STATISTICS: ANCOVA
wilcox.test(df$ATI_2 ~ df$DX_APD, paired=F)
aov <- aov(ATI_2 ~ Age + DX_APD, df) 
Anova(aov, type="II") 
	check_normality(aov)
etaSquared(aov)


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.9. BIOMARKERS: NFL                 #########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# NFL STATISTICS: OUTLIERS
# NfL values are typically right skewed especially in FTLD-related diagnoses. Therefore, different approach for outlier identification was chosen.  
# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# For reference, outliers are added to the notes of Tables. 
boxplot(NFL_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(NFL_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(NFL_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(NFL_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

thresholdCBS <- min(max(CBSdf$NFL_2,na.rm=T), as.numeric(quantile(CBSdf$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdf$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
thresholdPSP <- min(max(PSPdf$NFL_2,na.rm=T), as.numeric(quantile(PSPdf$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (PSPdf$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
cat("Outliers are values above ", thresholdCBS, " in CBS subset. \n")
cat("Outliers are values above ", thresholdPSP, " in PSP subset. \n")

dfnfl<- df[df$DX_APD=="PSP" | df$NFL_2 <= thresholdCBS, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfnfl<- dfnfl[dfnfl$DX_APD=="CBS" | dfnfl$NFL_2 <= thresholdPSP, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are PSP and either have no NFL value or one over the threshold

removed <- setdiff(df, dfnfl) 
cat("Following values were removed for the descriptive stats on NfL: ", removed$NFL_2, "\n")

# NFL STATISTICS: DISTRIBUTION
shapiro.test(dfnfl[dfnfl$DX_APD=="CBS",  ]$logNFL) #normal
shapiro.test(dfnfl[dfnfl$DX_APD=="PSP",  ]$logNFL) #normal
hist(df$logNFL)
hist(CBSdf$logNFL)
hist(PSPdf$logNFL)
leveneTest(logNFL ~ DX_APD, data = dfnfl) #homoscedasticity

# NFL STATISTICS: SUMMARY
dfnfl %>% summarize(count=n(), format(round(mean(NFL_2, na.rm=T),2),2), sd=sd(NFL_2, na.rm=T)) #Rounds up the sd for some reason
dfnfl %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(NFL_2, na.rm=T),2),2), sd=sd(NFL_2, na.rm=T)) #Rounds up the sd for some reason
sd(dfnfl[dfnfl$DX_APD=="CBS",  ]$NFL_2)
sd(dfnfl[dfnfl$DX_APD=="PSP",  ]$NFL_2)

# NFL STATISTICS: ANCOVA
t.test(dfnfl$logNFL ~ dfnfl$DX_APD, var.equal=TRUE) 
aov <- aov(logNFL ~ Age + DX_APD, dfnfl) 
Anova(aov, type="II")
check_normality(aov)
etaSquared(aov)


cat("\n\n#######################################################################################################\n",
	"                                    4.2. COHORT CHARACTERISTICS: CATEGORICAL VARIABLES\n",
	    "#######################################################################################################\n")

# For categorical variables, when an association between 2 variables was tested, the table of contingency was computed,
## and if the expected count was under 5 for any of the 4 cells, the Fisher test was chosen over the Chi-square. 

cat("\n\n######################################################################################################\n",
	   "################################            4.2.1. SEX                         ########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

# SEX STATISTICS: SUMMARY
cat("Sex distribution in the dataset is: \n")
df %>% group_by(DX_APD) %>% count(Sex) 


# SEX STATISTICS: PERCENTAGES
totalmatrix <- df %>% count(DX_APD)
sexmatrix <- df %>% group_by(DX_APD) %>% count(Sex) 

cat("Proportion of females in total is: \n")
(as.numeric(sexmatrix$n[1]) + as.numeric(sexmatrix$n[3]))/(as.numeric(totalmatrix$n[1]) + as.numeric(totalmatrix$n[2]))
cat("Proportion of females in CBS is:")
as.numeric(sexmatrix$n[1])/as.numeric(totalmatrix$n[1])
cat("Proportion of females in PSP is:")
as.numeric(sexmatrix$n[3])/as.numeric(totalmatrix$n[2])

# SEX STATISTICS: CHI-SQUARE
table(df$Sex, df$DX_APD)
chisq.test(table(df$Sex, df$DX_APD), correct=F)


cat("\n\n######################################################################################################\n",
	   "########################              4.2.2. APOE                         #############################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

# APOE STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% count(APOEe4)
table(df$APOEe4, df$DX_APD)


# APOE STATISTICS: CHI-SQUARE
chisq.test(table(df$APOEe4, df$DX_APD), correct=F)


cat("\n\n######################################################################################################\n",
	   "#########################              4.2.3. AD                         ##############################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")
cat("---------------------------   GOES IN TEXT - RESULTS - COHORT CHARACTERISTICS  ----------------------------\n")

# AD STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% count(AD)
table(df$DX_APD, df$AD)

# AD STATISTICS: CHI-SQUARE
chisq.test(table(df$AD, df$DX_APD), correct=F)



cat("\n\n\n\n###################################################################################################\n",
		   "5. ASYN-SAA+ & COHORT CHARACTERISTICS\n",
	 	   "####################################################################################################\n\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")

df %>% group_by(RTQUIC) %>% count(DX_APD)
sum(df$DX_APD=="CBS" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")
sum(df$DX_APD=="PSP" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")


cat("\n\n#######################################################################################################\n",
	"                         5.1. ASYN-SAA+ COHORT CHARACTERISTICS: CATEGORICAL VARIABLES  \n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")


# FOR COMPARISONS WITHIN DIAGNOSIS, GO LOWER to 8.
# SEX STATISTICS: DISTRIBUTION
table(df$Sex, df$RTQUIC)

# SEX STATISTICS: CHI-SQUARE AND FISHER'S TEST
chisq.test(table(df$Sex, df$RTQUIC), correct=F) 

# APOE STATISTICS: DISTRIBUTION
table(df$APOEe4, df$RTQUIC)

# APOE STATISTICS: FISHER'S TEST
fisher.test(table(df$APOEe4, df$RTQUIC)) # Expected count is <5 for one cell

# PARK ONSET STATISTICS: DISTRIBUTION
table(df$Parkinsonian_onset, df$RTQUIC)

# PARK ONSET STATISTICS: CHI-SQUARE AND FISHER'S TEST
chisq.test(table(df$Parkinsonian_onset, df$RTQUIC), correct=F) 


cat("\n\n#######################################################################################################\n",
	"                         5.2. ASYN-SAA+ COHORT CHARACTERISTICS: NUMERICAL VARIABLES  \n",
	   "#######################################################################################################\n")

#BONFERRONI CALCULATION FOR AGE/ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO RTQUIC
#Looking at three similar tests: Difference in mean age at LP, onset, and Parkinsonism onset, in RT-QUIC+ vs RT-QUIC-
0.05/3 #0.01666667


cat("\n\n######################################################################################################\n",
	   "#########################              5.2.1. ASYN-SAA*AGE                #############################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")

# ASYN-SAA*AGE STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Age) #normal
shapiro.test(RTnegdf$Age) #normal
var.test(Age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*AGE STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))

# RTQUIC*AGE STATISTICS: ANCOVA
t.test(df$Age ~ df$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Age ~ DX_APD*RTQUIC, df)
Anova(aov, type="III")
etaSquared(aov)

# p=.007 < .017 so reported bonf-adjusted p-value: .007x3=.02

cat("\n\n######################################################################################################\n",
	   "######################               5.2.2. ASYN-SAA*ONSET                #############################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")

# ASYN-SAA*ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Onset_age) #normal
shapiro.test(RTnegdf$Onset_age) #normal
var.test(Onset_age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*ONSET STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset_age, na.rm=T), sd=sd(Onset_age, na.rm=T))

# RTQUIC*ONSET STATISTICS: ANCOVA
t.test(df$Onset_age ~ df$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Onset_age ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")
etaSquared(aov)

# p >.017 so no need to bonferroni correct as it is ns


cat("\n\n######################################################################################################\n",
	   "################################          5.2.3. ASYN-SAA*PARK_ONSET           ########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")

# RTQUIC*PARK_ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Park_onset) #normal
shapiro.test(RTnegdf$Park_onset) #normal
var.test(Park_onset ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*PARK_ONSET STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Park_onset, na.rm=T))

# RTQUIC*PARK_ONSET STATISTICS: ANCOVA
t.test(df$Park_onset ~ df$RTQUIC, var.equal=TRUE)
aov <- aov(Park_onset ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")
etaSquared(aov)

# p >.017 so no need to bonferroni correct as it is ns

cat("\n\n\n\n###################################################################################################\n",
		   "6. AD+ & COHORT CHARACTERISTICS\n",
	 	   "####################################################################################################\n\n")

cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")
df %>% group_by(DX_APD) %>% count(AD)

cat("\n\n#######################################################################################################\n",
	"                          6.1. AD+ WHOLE COHORT CHARACTERISTICS: CATEGORICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("--------------------------------------   FOR REFERENCE ONLY  --------------------------------------\n")

# SEX STATISTICS: DISTRIBUTION
table(df$AD, df$Sex) 

# SEX STATISTICS: CHI SQUARE
chisq.test(table(df$AD, df$Sex), correct=F)

# PPA STATISTICS: DISTRIBUTION
table(df$AD, df$anyPPA) 

# PPA STATISTICS: CHI SQUARE
fisher.test(table(df$AD, df$anyPPA)) #expected outcome < 5. Significant

# PARK ONSET STATISTICS: DISTRIBUTION
table(df$AD, df$Parkinsonian_onset) 

# PARK ONSET STATISTICS: CHI SQUARE
chisq.test(table(df$AD, df$Parkinsonian_onset), correct=F)


# PARK ONSET STATISTICS: DISTRIBUTION
table(df$AD, df$APOEe4) 

# PARK ONSET STATISTICS: CHI SQUARE
fisher.test(table(df$AD, df$APOEe4)) #expected outcome < 5


cat("\n\n#######################################################################################################\n",
	"                          6.2. AD+ WHOLE COHORT CHARACTERISTICS: NUMERICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("--------------------------------------   FOR REFERENCE ONLY  --------------------------------------\n")

# AGE STATISTICS: DISTRIBUTION
shapiro.test(ADposdf$Age) #normal
shapiro.test(ADnegdf$Age) #normal
leveneTest(Age ~ AD, data = df) #homoscedasticity

# AGE STATISTICS: SUMMARY
df %>% group_by(AD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))

# AGE STATISTICS: TTEST
t.test(df$Onset_age ~ df$AD, var.equal=TRUE) 

# ONSET STATISTICS: DISTRIBUTION
shapiro.test(ADposdf$Onset_age) #normal
shapiro.test(ADnegdf$Onset_age) #normal
leveneTest(Onset_age ~ AD, data = df) #homoscedasticity

# ONSET STATISTICS: SUMMARY
df %>% group_by(AD) %>% summarize(count=n(), format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))

# ONSET STATISTICS: TTEST
t.test(df$Onset_age ~ df$AD, var.equal=TRUE) 
		

# PARK ONSET STATISTICS: DISTRIBUTION
shapiro.test(ADposdf$Park_onset) #borderline
shapiro.test(ADnegdf$Park_onset) #borderline
leveneTest(Park_onset ~ AD, data = df) #homoscedasticity

# PARK ONSET STATISTICS: SUMMARY
df %>% group_by(AD) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# PARK ONSET STATISTICS: WILCOX
wilcox.test(df$Park_onset ~ df$AD, paired=F)


# DURATION STATISTICS: DISTRIBUTION
shapiro.test(ADposdf$LP2_Disease_Duration) #not normal
shapiro.test(ADnegdf$LP2_Disease_Duration) #not normal
leveneTest(LP2_Disease_Duration ~ AD, data = df) #homoscedasticity but borderline


# DURATION STATISTICS: SUMMARY
df %>% group_by(AD) %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))


# DURATION STATISTICS: WILCOX
wilcox.test(df$LP2_Disease_Duration ~ df$AD, paired=F)
		

cat("\n\n#######################################################################################################\n",
	"                          6.3. AD+ CBS COHORT CHARACTERISTICS: CATEGORICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-     ----------------------------------\n")

# SEX STATISTICS: DISTRIBUTION
table(CBSdf$AD, CBSdf$Sex) 

# SEX STATISTICS: CHI SQUARE
chisq.test(table(CBSdf$AD, CBSdf$Sex), correct=F)

# APOE STATISTICS: DISTRIBUTION
table(CBSdf$APOEe4, CBSdf$AD)

# APOE STATISTICS: FISHER'S TEST
fisher.test(table(CBSdf$APOEe4, CBSdf$AD)) # Expected count is <5 for one cell

# LANG ONSET STATISTICS: DISTRIBUTION
table(CBSdf$Language_onset, CBSdf$AD)

# LANG ONSET STATISTICS: FISHER'S TEST
fisher.test(table(CBSdf$Language_onset, CBSdf$AD)) # Expected count is <5 for one cell

# PPA STATISTICS: DISTRIBUTION
table(CBSdf$anyPPA, CBSdf$AD)

# PPA STATISTICS: FISHER'S TEST
fisher.test(table(CBSdf$anyPPA, CBSdf$AD)) # Expected count is <5 for one cell



cat("\n\n#######################################################################################################\n",
	"                          6.4. AD+ CBS COHORT CHARACTERISTICS: NUMERICAL VARIABLES \n",
	"###########################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "#########################        6.4.1. AGE, ONSET, PARK ONSET                #########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-     ----------------------------------\n")

# BONFERRONI CALCULATION FOR AGE/ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO AD 
0.05/3 #0.017

# AGE STATISTICS: DISTRIBUTION
shapiro.test(CBSdf[CBSdf$AD=="AD Positive", ]$Age) #normal
shapiro.test(CBSdf[CBSdf$AD=="AD Negative", ]$Age) #normal
leveneTest(Age ~ AD, CBSdf) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# AGE STATISTICS: SUMMARY
CBSdf %>% summarize(count=n(), mean=format(round(mean(Age, na.rm=T),3),3), sd=format(round(sd(Age, na.rm=T),3),3))
CBSdf %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Age, na.rm=T),3),3), sd=format(round(sd(Age, na.rm=T),3),3))

# AGE STATISTICS: TTEST
t.test(CBSdf$Age ~ CBSdf$AD, var.equal=TRUE)


# ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf[CBSdf$AD=="AD Positive", ]$Onset_age) #normal
shapiro.test(CBSdf[CBSdf$AD=="AD Negative", ]$Onset_age) #normal
leveneTest(Onset_age ~ AD, CBSdf) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# ONSET STATISTICS: SUMMARY
CBSdf %>% summarize(count=n(), mean=format(round(mean(Onset_age, na.rm=T),3),3), sd=format(round(sd(Onset_age, na.rm=T),3),3))
CBSdf %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Onset_age, na.rm=T),3),3), sd=format(round(sd(Onset_age, na.rm=T),3),3))

# ONSET STATISTICS: TTEST
t.test(CBSdf$Onset_age ~ CBSdf$AD, var.equal=TRUE)


# PARK ONSET STATISTICS: DISTRIBUTION
shapiro.test(CBSdf[CBSdf$AD=="AD Positive", ]$Park_onset) #normal
shapiro.test(CBSdf[CBSdf$AD=="AD Negative", ]$Park_onset) #normal
leveneTest(Park_onset ~ AD, CBSdf) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# PARK ONSET STATISTICS: SUMMARY
CBSdf %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))
CBSdf %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))

# PARK ONSET STATISTICS: TTEST
t.test(CBSdf$Park_onset ~ CBSdf$AD, var.equal=TRUE)


cat("\n\n######################################################################################################\n",
	   "#####################################               6.4.2. NFL                #########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-     ----------------------------------\n")

# NFL STATISTICS: DISTRIBUTION
boxplot(NFL_2 ~ AD, data= CBSdf[CBSdf$AD=="AD Positive", ], col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(NFL_2 ~ AD, data = CBSdf[CBSdf$AD=="AD Positive", ], method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(NFL_2 ~ AD, data= CBSdf[CBSdf$AD=="AD Negative", ], col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(NFL_2 ~ AD, data = CBSdf[CBSdf$AD=="AD Negative", ], method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# This is because in FTLD population, the data can be very right-skewed. 
# Although technically median would still be ok to show in the Table even with outliers, but because it may be of intrest to share these extreme values, I prefer having thenm in the Notes.  
thresholdAD <- min(max(CBSdf[CBSdf$AD=="AD Positive", ]$NFL_2,na.rm=T), as.numeric(quantile(CBSdf[CBSdf$AD=="AD Positive", ]$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdf[CBSdf$AD=="AD Positive", ]$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
thresholdnonAD <- min(max(CBSdf[CBSdf$AD=="AD Negative", ]$NFL_2,na.rm=T), as.numeric(quantile(CBSdf[CBSdf$AD=="AD Negative", ]$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdf[CBSdf$AD=="AD Negative", ]$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
cat("Outliers are values above ", thresholdAD, " in CBS-AD+ subset. \n")
cat("Outliers are values above ", thresholdnonAD, " in CBS-AD- subset. \n")

CBSdfnfl<- CBSdf[CBSdf$AD=="AD Negative" | CBSdf$NFL_2 <= thresholdAD, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
CBSdfnfl<- CBSdfnfl[CBSdfnfl$AD=="AD Positive" | CBSdfnfl$NFL_2 <= thresholdnonAD, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are PSP and either have no NFL value or one over the threshold

removed <- setdiff(CBSdf, CBSdfnfl) 
cat("Following values were removed for the descriptive stats on NfL: ", removed$NFL_2, "\n")

shapiro.test(CBSdfnfl[CBSdfnfl$AD=="AD Positive",  ]$logNFL) #non-normal
shapiro.test(CBSdfnfl[CBSdfnfl$AD=="AD Negative",  ]$logNFL) #non-normal
leveneTest(logNFL ~ AD, data = CBSdfnfl) #heterodasticity

# NFL STATISTICS: SUMMARY
CBSdfnfl %>% summarize(count=n(), format(round(median(NFL_2, na.rm=T),2),2), IQR=IQR(NFL_2, na.rm=T), min=min(NFL_2, na.rm=T), max=max(NFL_2, na.rm=T))
CBSdfnfl %>% group_by(AD) %>% summarize(count=n(), format(round(median(NFL_2, na.rm=T),2),2), IQR=IQR(NFL_2, na.rm=T), min=min(NFL_2, na.rm=T), max=max(NFL_2, na.rm=T))
IQR(CBSdfnfl[CBSdfnfl$AD=="AD Positive",  ]$NFL_2,  na.rm=T)
IQR(CBSdfnfl[CBSdfnfl$AD=="AD Negative",  ]$NFL_2,  na.rm=T)


# NFL STATISTICS: WEIGHTED LEAST-SQUARE REGRESSION
t.test(CBSdfnfl$logNFL ~ CBSdfnfl$AD, var.equal=FALSE) #Welch t-test, to deal with the heterodasticity
gls1 <- gls(logNFL ~ Age*AD, CBSdfnfl, weights=varPower())  #to deal with the heterodasticity
gls2 <- gls(logNFL ~ Age + AD, CBSdfnfl, weights=varPower()) 
AIC(gls1, gls2) #gls2 better than gls1 
summary(gls2)


cat("\n\n######################################################################################################\n",
	   "###############################           6.4.3. MOCA Z-SCORES                #########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-     ----------------------------------\n")


# Of note, while the MoCA scores are z-scored to control for the effect of age, that does not mean there is no effect of age when analyzing them. 
# MoCA z-scores are worse in younger subjects, since the latter are compared to a population with no cognitive deficit. So even small changes in MoCA score
## in young subjects trnaslates to very  low z-scores.  
# Here it doesnt matter as we are just presenting the values in the table to summarize the dataset.

# MOCA Z-SCORE STATISTICS: DISTRIBUTION
shapiro.test(CBSdf[CBSdf$AD=="AD Positive", ]$LP2_MOCA_Z.score) #normal
shapiro.test(CBSdf[CBSdf$AD=="AD Negative", ]$LP2_MOCA_Z.score) #borderline
hist(CBSdf[CBSdf$AD=="AD Negative", ]$LP2_MOCA_Z.score) #not normal
leveneTest(LP2_MOCA_Z.score ~ AD, CBSdf) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# MOCA Z-SCORE STATISTICS: SUMMARY
CBSdf %>% summarize(count=n(), format(round(median(LP2_MOCA_Z.score, na.rm=T),2),2), IQR=IQR(LP2_MOCA_Z.score, na.rm=T), min=min(LP2_MOCA_Z.score, na.rm=T), max=max(LP2_MOCA_Z.score, na.rm=T))
CBSdf %>% group_by(AD) %>% summarize(count=n(), format(round(median(LP2_MOCA_Z.score, na.rm=T),2),2), IQR=IQR(LP2_MOCA_Z.score, na.rm=T), min=min(LP2_MOCA_Z.score, na.rm=T), max=max(LP2_MOCA_Z.score, na.rm=T))

# MOCA Z-SCORE STATISTICS: WILCOX
wilcox.test(CBSdf$LP2_MOCA_Z.score ~ CBSdf$AD)


cat("\n\n\n\n###################################################################################################\n",
		   "7. ASYN-SAA+ & AD\n",
	 	   "####################################################################################################\n\n")

cat("\n\n#######################################################################################################\n",
	"                                    7.1. ASYN-SAA+ & AD: CATEGORICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & DEMOGRAPHICS  ----------------------------\n")

sum(df$AD=="AD Positive" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")
sum(df$AD=="AD Positive" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")
sum(df$AD=="AD Negative" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")
sum(df$AD=="AD Negative" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")


cat("\n\n######################################################################################################\n",
	   "############################              7.1.1.FIGURE 1A           		###########################\n",
	   "#######################################################################################################\n")

cat("-----------------------   GOES IN FIGURE 1A: ONSET * AD QUADRANTS   --------------------------------------\n")

cat("Below are the values for the creation of Fig 1.A: \n")
cat("Lower left quadrant: AD+ and young-onset \n")
n1 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Young-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and young-onset and AD+: ", n1, "\n")
n2 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Young-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD+: ", n2, "\n")
n3 <-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Young-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and young-onset and AD+: ", n3, "\n")
n4 <-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Young-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD+: ", n4, "\n\n")


cat("Upper left quadrant: AD- and young-onset\n")
n5 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Young-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and young-onset and AD-: ", n5, "\n")
n6 <- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Young-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD-: ", n6, "\n")
n7<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Young-onset") & (df$AD=="AD Negative") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and young-onset and AD-: ", n7, "\n")
n8<-nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Young-onset") & (df$AD=="AD Negative") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and young-onset and AD-: ", n8, "\n\n")

cat("Lower right quadrant: \n")
n9<-nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC+ and late-onset and AD+: ", n9, "\n")
n10<- nrow(df[(df$RTQUIC=="aSyn-SAA positive") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
cat("PSP subjects who are RTQUIC+ and late-onset and AD+: ", n10, "\n")
n11<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="CBS"),])
cat("CBS subjects who are RTQUIC- and late-onset and AD+: ", n11, "\n")
n12<- nrow(df[(df$RTQUIC=="aSyn-SAA negative") & (df$Early_onset=="Late-onset") & (df$AD=="AD Positive") & (df$DX_APD=="PSP"),])
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


# DEFENSE
if (n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n12+n13+n14+n15+n16 != 67) {
	cat("Issue when calculating the values for Figure 1A. Check that there was no error in condition setup. \n")
	cat(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n12+n13+n14+n15+n16)
}

if (n1+n2+n3+n4 != 9) {cat("Issue in lower left quadrant: AD+ and young-onset. \n", n1, n2, n3, n4, "ie ", n1+n2+n3+n4, "instead of 9\n")}
if (n5+n6+n7+n8 != 26) {cat("Issue in upper left quadrant: AD- and young-onset. \n", n5, n6, n7, n8, "ie ", n5+n6+n7+n8, "instead of 26\n")}
if (n9+n10+n11+n12 != 6) {cat("Issue in lower right quadrant: AD+ and late-onset. \n", n9, n10, n11, n12,"ie ", n9+n10+n11+n12,  "instead of 6\n")}
if (n13+n14+n15+n16 != 26) {cat("Issue in upper right quadrant: AD- and late-onset. \n", n13, n14, n15, n16, "ie ", n13+n14+n15+n16, "instead of 26\n")}


	# If there is an issue, check below the dfs: 
	# df[(df$AD=="AD Positive") & (df$Early_onset=="Young-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")]
	# nrow(df[(df$AD=="AD Positive") & (df$Early_onset=="Young-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])
	# nrow(df[(df$AD=="AD Positive") & (df$Early_onset=="Young-onset") & (df$RTQUIC=="aSyn-SAA positive"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])

	# df[(df$AD=="AD Negative") & (df$Early_onset=="Young-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")]
	# nrow(df[(df$AD=="AD Negative") & (df$Early_onset=="Young-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])
	# nrow(df[(df$AD=="AD Negative") & (df$Early_onset=="Young-onset") & (df$RTQUIC=="aSyn-SAA positive"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])

	# df[(df$AD=="AD Positive") & (df$Early_onset=="Late-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")]
	# nrow(df[(df$AD=="AD Positive") & (df$Early_onset=="Late-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])
	# nrow(df[(df$AD=="AD Positive") & (df$Early_onset=="Late-onset") & (df$RTQUIC=="aSyn-SAA positive"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])

	# df[(df$AD=="AD Negative") & (df$Early_onset=="Late-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")]
	# nrow(df[(df$AD=="AD Negative") & (df$Early_onset=="Late-onset"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])
	# nrow(df[(df$AD=="AD Negative") & (df$Early_onset=="Late-onset") & (df$RTQUIC=="aSyn-SAA positive"), c("ID", "RTQUIC", "DX_APD", "AD", "AD_lifetime_ATHENA", "Onset_age", "Early_onset")])



cat("\n\n######################################################################################################\n",
	   "########################              7.1.2. FREQUENCY DATA           		###########################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & AD+  ------------------------------------\n")
cat("-----------------------   GOES IN FIGURE 1A: ONSET * AD QUADRANTS   --------------------------------------\n")

#BONFERRONI CALCULATION FOR FISHER TESTS USED FOR ASYN/AD RELATIONSHIP IN OVERALL DATASET AND YO DATASET AND LO DATASET
0.05/3 #0.017

# WHOLE COHORT STATISTICS: DISTRIBUTION
table(df$RTQUIC, df$AD)
sum(YODdf$AD=="AD Positive" & YODdf$RTQUIC=="aSyn-SAA positive")
sum(YODdf$RTQUIC=="aSyn-SAA positive")
sum(LODdf$AD=="AD Positive" & LODdf$RTQUIC=="aSyn-SAA positive")
sum(LODdf$RTQUIC=="aSyn-SAA positive")


# WHOLE COHORT STATISTICS: FISHER
fisher.test(table(df$AD, df$RTQUIC)) # Expected count is <5 for one cell. Gives OR
cramerV(table(df$AD, df$RTQUIC)) 
phi(table(df$AD, df$RTQUIC), digits=6)
	
	# NOTES ON CRAMER AND PHI:
	# df =(r - 1) * (c - 1) = 1*1 = 1. At df=1, medium Cramer's V=.3
	# CRAMER: COHEN 1988 REFERS TO AS CRAMER'S PHI. See p. 223 https://www.utstat.toronto.edu/~brunner/oldclass/378f16/readings/CohenPower.pdf
	# PHI: IN 2X2 CASE, SAME AS W IN COHEN 1988 7.2. See p.223 https://www.utstat.toronto.edu/~brunner/oldclass/378f16/readings/CohenPower.pdf

#YOUNG-ONSET STATISTICS+: DISTRIBUTION
table(YODdf$AD, YODdf$RTQUIC)

#YOUNG-ONSET STATISTICS+: FISHER
fisher.test(table(YODdf$AD, YODdf$RTQUIC)) # Expected count is <5 for one cell
cramerV(table(YODdf$AD, YODdf$RTQUIC)) 
phi(table(YODdf$AD, YODdf$RTQUIC), digits=6) 

# p=.015 < .017 so reported bonf-adjusted p-value: .015x3=.045

# LATE-ONSET STATISTICS: DISTRIBUTION
table(LODdf$AD, LODdf$RTQUIC)

# LATE-ONSET STATISTICS: FISHER
fisher.test(table(LODdf$AD, LODdf$RTQUIC)) # Expected count is <5 for one cell
cramerV(table(LODdf$AD, LODdf$RTQUIC)) 
phi(table(LODdf$AD, LODdf$RTQUIC), digits=6) 


cat("\n\n######################################################################################################\n",
	   "########################         7.1.3. FREQUENCY DATA WITHIN DX        	###########################\n",
	   "#######################################################################################################\n")

cat("-----------------------   GOES IN SUPP TEXT - ASYN-SAA+ & AD+ WITHIN DX ----------------------------------\n")

# 1. TRY TO SEE IF IN CBS, AD AND RTQUIC ARE ASSOCIATED. THEN, SEE IF EFFECT OF AGE. 

# CBS ASSOCIATION OF AD+/RTQUIC+ OVERALL
table(CBSdf$RTQUIC, CBSdf$AD)
fisher.test(table(CBSdf$AD, CBSdf$RTQUIC)) # Expected count is <5 for one cell. Gives OR
cramerV(table(CBSdf$AD, CBSdf$RTQUIC))
phi(table(CBSdf$AD, CBSdf$RTQUIC), digits=6)
# There is almost a significant relationship between AD+ and RTQUIC+ in CBS (p<.10)

# CBS ASSOCIATION OF AD+/RTQUIC+ IN YOUNG-ONSET ONLY (25/39 subjects)
nrow(CBSdf[CBSdf$"Early_onset"=="Young-onset", ])
table(CBSdf[CBSdf$"Early_onset"=="Young-onset", ]$RTQUIC, CBSdf[CBSdf$"Early_onset"=="Young-onset", ]$AD)
fisher.test(table(CBSdf[CBSdf$"Early_onset"=="Young-onset", ]$RTQUIC, CBSdf[CBSdf$"Early_onset"=="Young-onset", ]$AD)) # Expected count is <5 for one cell. Gives OR
# There is almost a significant relationship between AD+ and RTQUIC+ in young-onset CBS (p<.10)

nrow(CBSdf[CBSdf$"Early_onset"=="Late-onset", ])
table(CBSdf[CBSdf$"Early_onset"=="Late-onset", ]$RTQUIC, CBSdf[CBSdf$"Early_onset"=="Late-onset", ]$AD)
fisher.test(table(CBSdf[CBSdf$"Early_onset"=="Late-onset", ]$RTQUIC, CBSdf[CBSdf$"Early_onset"=="Late-onset", ]$AD)) # Expected count is <5 for one cell. Gives OR

# From above we can conclude that there is evidence of a potential relationship between AD+ and RTQUIC+ in CBS, but not necessarily affected by age at onset. 
# Hard to make conclusions regarding late-onset CBS as the sample size is low (14, including 5 who are RTQUIC+) but there is no evidence supporting that AD+ and RTQUIC+
## are associated in late-onset CBS.

# CBS ASSOCIATION OF RTQUIC+ WITH AGE
shapiro.test(CBSdf[CBSdf$"RTQUIC"=="aSyn-SAA positive", ]$Onset_age) #normal
shapiro.test(CBSdf[CBSdf$"RTQUIC"=="aSyn-SAA negative", ]$Onset_age)#normal
var.test(Onset_age ~ RTQUIC, data = CBSdf) #homoscedasticity

# CBS RTQUIC*ONSET STATISTICS: ANCOVA
t.test(CBSdf$Onset_age ~ CBSdf$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Onset_age ~ RTQUIC + AD, data=CBSdf)
summary(aov)
etaSquared(aov)

# CBS RTQUIC*AGE STATISTICS: EXTRA
t.test(CBSdf$Age ~ CBSdf$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Age ~ AD + RTQUIC, CBSdf) 
summary(aov)
etaSquared(aov)

# Conclusion for CBS: There is no evidence of relationship between age and RTQUIC status within the CBS diagnosis.
## There IS evidence of a relationship between AD+ and RTQUIC+. 
## These results could be attributable to the low sample size (especially fewer subjects who have late-onset disease. )
## The same is seen when looking directly at Age instead of Onset age. 

# 2. TRY TO SEE IF IN PSP, THERE IS AN EFFECT OF AGE. 

# PSP ASSOCIATION OF RTQUIC+ WITH AGE: 
table(PSPdf$RTQUIC, PSPdf$Early_onset)
fisher.test(table(PSPdf$RTQUIC, PSPdf$Early_onset)) # Expected count is <5 for one cell. Gives OR

shapiro.test(PSPdf[PSPdf$"RTQUIC"=="aSyn-SAA positive", ]$Onset_age) #normal
shapiro.test(PSPdf[PSPdf$"RTQUIC"=="aSyn-SAA negative", ]$Onset_age)#normal
var.test(Onset_age ~ RTQUIC, data = PSPdf) #homoscedasticity

# PSP RTQUIC*ONSET STATISTICS: ANCOVA
t.test(PSPdf$Onset_age ~ PSPdf$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Onset_age ~ RTQUIC, data=PSPdf)
summary(aov)
etaSquared(aov)
t.test(PSPdf$Age ~ PSPdf$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Age ~ AD + RTQUIC, CBSdf) 
etaSquared(aov)

# Conclusion for PSP: There is no evidence of relationship between age and RTQUIC status within the PSP diagnosis.
## However, due to the much smaller size of this group, it is hard to tell whether overall age-related/pathology-specific changes in PSP 
## are resposnible for alpha-syn misfolding (as PSP subjects overall are older)
## Of note, among the 8 RTQUCI+ PSP subjects, only 1 has an onset age <65. Given the aSyn-SAA- is balanced between young-onset and late-osnet, 
## this would suggest that we are simply underpowered. 


cat("\n\n#######################################################################################################\n",
	"                                    7.2. ASYN-SAA+ & AD: NUMERICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "########################              7.2.1. ASYN-SAA+ & ABETA42          	###########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN TABLE 2: MLR MODEL OUTPUT   --------------------------------------\n")
cat("-------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & AD+  ------------------------------------\n")

# Due to relative low sample size and heterogeneity of the cohort (especially for RTQUIC+ or young-onset AD+ subjects) for model complexity,
# preferred approach is to run model on all points then run diagnostic analyses and present different model runs (leaving one outlier out)
boxplot(abeta_2 ~ RTQUIC, data= RTposdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(abeta_2 ~ RTQUIC, data = RTposdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(abeta_2 ~ RTQUIC, data= RTnegdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(abeta_2 ~ RTQUIC, data = RTnegdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

# ABETA STATS: DISTRIBUTION
shapiro.test(RTposdf$logabeta) #normal
shapiro.test(RTnegdf$logabeta)  #normal
shapiro.test(CBSdf$logabeta)  #normal
shapiro.test(PSPdf$logabeta)  #normal
leveneTest(logabeta ~ DX_APD*RTQUIC, df) #homoscedasticity

# ABETA STATS: LINEAR REGRESSION MODEL SELECTION
test1 <- lm(logabeta ~ RTQUIC, df)
test2 <- lm(logabeta ~ RTQUIC  + DX_APD, df)
test3 <- lm(logabeta ~ RTQUIC + Sex, df) #not adding any value to the model
anova(test1, test2) 
anova(test1, test3) 

# Inclusion of Onset as a covariate. Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
ggscatter(df, x = "Onset_age", y = "logabeta", add = "reg.line")+ 
	stat_regline_equation(aes()) 
ggscatter(df, x = "Onset_age", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
stat_regline_equation(aes(color = DX_APD)) 
ggscatter(df, x = "Onset_age", y = "logabeta", color = "RTQUIC", add = "reg.line")+
	stat_regline_equation(aes(color = RTQUIC))
cor.test(df$logabeta, df$Onset_age) #corr
cor.test(df$logabeta, df$Onset_age, method="spearman") #corr
summary(lm(df$logabeta ~ df$Onset_age)) #linear relationship 

## Inclusion of NFL as a covariate. Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
ggscatter(df, x = "NFL_2", y = "logabeta", add = "reg.line")+ 
	stat_regline_equation(aes()) 
ggscatter(df, x = "NFL_2", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
	stat_regline_equation(aes(color = DX_APD)) 
ggscatter(df, x = "NFL_2", y = "logabeta", color = "RTQUIC", add = "reg.line")+
	stat_regline_equation(aes(color = RTQUIC))
cor.test(df$logabeta, df$NFL_2) #corr
cor.test(df$logabeta, df$NFL_2, method="spearman") #corr
summary(lm(df$logabeta ~ df$NFL_2)) #linear relationship 

## Inclusion of other variables as covariate:
summary(lm(df$logabeta ~ df$LP2_Disease_Duration)) #no linear relationship
summary(lm(df$logabeta ~ df$ttau_2)) #no linear relationship
summary(lm(df$logabeta ~ df$ptau_2)) #no linear relationship

## Comparison of complex models
# For biomarkers, only kept logged value on one side of the model equation
stdmlr <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD + scale(NFL_2), df) 
summary(stdmlr)
AIC(stdmlr) 

stdmlr2 <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD, df) 
summary(stdmlr2)
AIC(stdmlr2)

stdmlr3 <- lm(logabeta ~ scale(Onset_age)*RTQUIC, df) 
summary(stdmlr3)
AIC(stdmlr3)

stdmlr4 <- lm(logabeta ~ scale(Onset_age)+ RTQUIC, df) 
summary(stdmlr4)
AIC(stdmlr4)

stdmlr5 <- lm(logabeta ~ scale(Onset_age)*RTQUIC + scale(NFL_2), df) 
summary(stdmlr5)
AIC(stdmlr5)

cat("The model with the lowest AIC is the one including an interaction term of Onset age by aSyn-SAA status as NfL levels. Its AIC is almost identical
	to the one that incorporates Diagnosis. Diagnosis is an extremely important covariate that is relevant clinically. Therefore, we prefer to report
	the model with teh lowest AIC taht still includes diagnosis, model stdmlr. Of note, all models have similar AIC and the interaction term  was
	significant in every model that it was included in. Removing NfL did not affect the results, which is important since NfL outliers were not removed
	prior to running the analysis. \n")

# ABETA STATS: MULTIVARIABLE LINEAR REGRESSION MODEL
stdmlr <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD + scale(NFL_2), df) 
summary(stdmlr)

# ABETA STATS: MULTIVARIABLE LINEAR REGRESSION MODEL DIAGNOSTICS
check_normality(stdmlr) #Ok normality of residuals
autoplot(stdmlr, which = 1:6) #All plots including scale location plot

# Visualize the values for each of the IDs that are indexed on above plots
df[5, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[21, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[27, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[37, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[45, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[65, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
		
# for-loop to test the model without each of these values (ie once without outlier 1, then outlier 2, etc)
vecIDs <- df[c(5,21,27,37,45,65), "ID"] #Create vector of each ID
for (i in vecIDs) {
	test <- subset(df, ID!=i) #for some analysis, need to exclude the potential false negative too
	teststdmlr <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD + scale(NFL_2), test) 
	print(summary(teststdmlr))
	}

plot(stdmlr, which = 3) # 3 = Scale-Location plot. Variance of residuals
bptest(stdmlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
durbinWatsonTest(stdmlr) #Ok autocorrelation of residuals	

# Multicollinearity checks
car::vif(stdmlr) #no multicollinearity at all

# Simple slopes for onset: 
emtrends(stdmlr,pairwise ~  RTQUIC, var="Onset_age") # for intereaction of Onset age by RTQUIC

# Main effects: 
emmeans(stdmlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
emmeans(stdmlr, ~ DX_APD) #adjusted means. Ok because no interaction. 


cat("\n\n######################################################################################################\n",
	   "###########################                 7.2.2. FIGURE 1.B.         	##############################\n",
	   "#######################################################################################################\n")

cat("-----------------------   GOES IN FIGURE 1B: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# FIG1B: ABETA42 over time from the real datapoints

##Create df for the figure based on the actual model
mylist <- list(Onset_age=seq(35,85,by=5), RTQUIC=c("aSyn-SAA positive", "aSyn-SAA negative")) 
emmip(stdmlr, RTQUIC ~ Onset_age, at=mylist, CIs=TRUE)

#Dataset: data are from the model (predicted values)
figdf <- emmip(stdmlr, RTQUIC ~ Onset_age, at=mylist, CIs=TRUE, plotit=FALSE) #
    
label1 <- "paste(F*'(5, 59)'==4.29*', ' ~~ italic(p), \"< .01\"*', ' ~~ italic(R)^2==20.47*'%')" #First annotation of the plot: Model diagnostics. Comma has to be entered as text not in mathematical notation. 
label2 <- "paste(''*italic(p), \" < .05\")" #Second annotation is p-value for the interaction

 # Change name of variable RTQUIC
fig1b <- ggplot(data=figdf, aes(x=Onset_age,y=yvar, color=RTQUIC)) + #yvar is Abeta42 logged 	#Ggplot figure basic layout

	# Add the actual trend/slopes + CI around it
	geom_line() + #DO NOT use show.legend as it will create a duplicate legend if combined with scale_color changes
	geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=RTQUIC), alpha=0.2, show.legend=FALSE) + #show.legend=FALSE to avoid having both ribbon and line in legend

	# Fix legends
	scale_fill_manual(values=cbPalette_RTQUIC) +
	scale_color_manual(values=cbPalette_RTQUIC, name= "ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
	scale_shape_manual(values=c(16,17), name="Diagnosis", breaks=c("CBS", "PSP"), labels=c("CBS", "PSP")) + #Color for the datapoints

	# Add actual data
	geom_point(data=df, aes(x=Onset_age, y=logabeta, shape=DX_APD), size=4) + #Actual datapoints 

	# Annotate:
	stat_regline_equation(label.x=37, label.y=c(7.4, 7.2), size=6, show.legend=FALSE) + #shows the equation for each geom_line. Cannot use other options since it would not be based on the model
			
	annotate("text", size=5, x=70, y=4, label=label1, parse=TRUE)+ #label as defined above. Parse allows for use of mathematical notation
	annotate("text", size=6, x=40, y=7, label=label2, parse=TRUE)+
				

	# Add the labels
	labs(title=expression(bold("Visual representation of the interaction of age at onset by "*alpha*"Syn-SAA")),
	subtitle=expression("From the model: log(A"*beta*"42) ~ Age at onset * "*alpha*"Syn-SAA + Diagnosis + NfL"),
	x="Age at onset (years)",
	y=expression(bold("CSF A"*beta*"42 levels (pg/mL) (log)")))+ #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

	# Aesthetic only
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(plot.subtitle = element_text(size=14, hjust=0.5)) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))

fig1b

ggsave(fig1b, filename = "Fig1b.png", bg= "transparent", width=9, height=10)



cat("\n\n\n\n###################################################################################################\n",
		   "8. ASYN-SAA+ & ALL VARIABLES\n",
	 	   "####################################################################################################\n\n")

cat("\n\n#######################################################################################################\n",
	"                                    8.1. ASYN-SAA+ & NFL \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "###########################              8.1.1. NFL MLR        	#######################################\n",
	   "#######################################################################################################\n")

cat("---------------------------   GOES IN TEXT - RESULTS - ASYN-SAA+ & NFL  ---------------------------------\n")

# Here the variable of interest is NfL, which is why removing the outliers is importabnt. 

# NFL STATISTICS: OUTLIERS
# NfL values are typically right skewed especially in FTLD-related diagnoses. Therefore, different approach for outlier identification was chosen.  
# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# For reference, outliers are added to the notes of Tables. 
boxplot(NFL_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(NFL_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(NFL_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(NFL_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

thresholdCBS <- min(max(CBSdf$NFL_2,na.rm=T), as.numeric(quantile(CBSdf$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdf$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
thresholdPSP <- min(max(PSPdf$NFL_2,na.rm=T), as.numeric(quantile(PSPdf$NFL_2, 0.75, na.rm=T)) + (IQR(na.rm=T, (PSPdf$NFL_2)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
cat("Outliers are values above ", thresholdCBS, " in CBS subset. \n")
cat("Outliers are values above ", thresholdPSP, " in PSP subset. \n")

dfnfl<- df[df$DX_APD=="PSP" | df$NFL_2 <= thresholdCBS, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are CBS and either have no NFL value or one over the threshold
dfnfl<- dfnfl[dfnfl$DX_APD=="CBS" | dfnfl$NFL_2 <= thresholdPSP, ] %>% remove_empty("rows") %>% data.frame() #Removes all subjects who are PSP and either have no NFL value or one over the threshold

removed <- setdiff(df, dfnfl) 
cat("Following values were removed for the descriptive stats on NfL: ", removed$NFL_2, "\n")

# NFL STATISTICS: DISTRIBUTION
shapiro.test(dfnfl[dfnfl$DX_APD=="CBS",  ]$logNFL) #normal
shapiro.test(dfnfl[dfnfl$DX_APD=="PSP",  ]$logNFL) #normal
shapiro.test(dfnfl[dfnfl$RTQUIC=="aSyn-SAA positive",  ]$logNFL) #borderline
shapiro.test(dfnfl[dfnfl$RTQUIC=="aSyn-SAA negative",  ]$logNFL) #normal
shapiro.test(dfnfl[dfnfl$AD=="AD Positive",  ]$logNFL) #normal
shapiro.test(dfnfl[dfnfl$AD=="AD Negative",  ]$logNFL) #borderline
hist(dfnfl[dfnfl$RTQUIC=="aSyn-SAA positive",  ]$logNFL)
hist(dfnfl[dfnfl$AD=="AD Negative",  ]$logNFL)
leveneTest(logNFL ~ DX_APD*AD*RTQUIC, data = dfnfl) #homoscedasticity

# NFL STATS: LINEAR REGRESSION MODEL SELECTION
# Compare models with Ftest
test1 <- lm(logNFL ~ RTQUIC, dfnfl)
test2 <- lm(logNFL ~ RTQUIC + DX_APD, dfnfl) #adding value to the model
test3 <- lm(logNFL ~ RTQUIC + AD, dfnfl) #not adding any value to the model
anova(test1, test2) #models that are nested in each other
anova(test1, test3) 

# Inclusion of other variables as covariate:
summary(lm(dfnfl$logNFL ~ dfnfl$Age)) #no linear relationship
summary(lm(dfnfl$logNFL ~ dfnfl$LP2_Disease_Duration))#no linear relationship 
summary(lm(dfnfl$logNFL ~ dfnfl$abeta_2)) #linear relationship
summary(lm(dfnfl$logNFL ~ dfnfl$ptau_2)) #no linear relationship
summary(lm(dfnfl$logNFL ~ dfnfl$ttau_2)) #no linear relationship

# NFL STATISTICS: LINEAR REGRESSION MODEL
# Based on clear linear relationship between abeta and NfL as seen above, it doesnt make sense to exclude Abeta from the model. 
# Meanwhile, dx is an important covariate,a nd RTQUIC is the IV of interest. So just comparing two models: with or without interaction. 
stdmlr <- lm(logNFL ~ RTQUIC*DX_APD + scale(abeta_2), dfnfl)
summary(stdmlr) 
AIC(stdmlr) #Slightly better AIC

stdmlr2 <- lm(logNFL ~ RTQUIC + DX_APD + scale(abeta_2), dfnfl) 
summary(stdmlr2) 
AIC(stdmlr2)

# AIC is slightly better when inclduing the interaction term between DX and RTQUIC, so will keep it as more informative. 
# Moreover, it seems in model without interaction, the main effect of diagnosis si not signifciant. Suggests that not including the interaction term
## woudl be misleading. 

# NFL STATISTICS: LINEAR REGRESSION MODEL DIAGNOSTICS
check_normality(stdmlr) #Ok normality of residuals
autoplot(stdmlr, which = 1:6) #All plots including scale location plot

# Visualize the values for each of the IDs that are indexed on above plots
dfnfl[8, c("abeta_2", "logabeta", "RTQUIC", "DX_APD", "NFL_2")] 
dfnfl[37, c("abeta_2", "logabeta", "RTQUIC", "DX_APD","NFL_2")] 
dfnfl[55, c("abeta_2", "logabeta", "RTQUIC", "DX_APD","NFL_2")] 
dfnfl[57, c("abeta_2", "logabeta", "RTQUIC", "DX_APD","NFL_2")] 

# for-loop to test the model without each of these values (ie once without outlier 1, then outlier 2, etc)
vecIDs <- dfnfl[c(8, 37, 55, 57), "ID"] #Create vector of each ID. For 1.5*IQR outlier threshold
for (i in vecIDs) {
	test <- subset(dfnfl, ID!=i) #for some analysis, need to exclude the potential false negative too
	teststdmlr <- lm(logNFL ~ RTQUIC*DX_APD + scale(abeta_2), test) 
	print(summary(teststdmlr))
	}

bptest(stdmlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
durbinWatsonTest(stdmlr) #Ok autocorrelation of residuals	

# Multicollinearity checks
car::vif(stdmlr) #no multicollinearity at all

# Main effects: 
emmeans(stdmlr, ~ RTQUIC:DX_APD) #adjusted means.  
pairs(emmeans(stdmlr, ~ RTQUIC:DX_APD)) #adjusted means.  


cat("\n\n######################################################################################################\n",
	   "###########################              8.1.2. FIG 1.D.       	#######################################\n",
	   "#######################################################################################################\n")

# FIG1D: NFL: Different between DX, but not RTQUIC, and linearly related to Abeta42.
# Option 1: Plot Abeta42 by NFL relationship and add colors for DX and shapes for RTQUIC. This is done here but not in paper. 
#Option 2: Plot boxplot between RTQUIC status (see below this figure)

cat("-----------------------  FOR REFERENCE ONLY   ------------------------------------\n")

# Create dataframe usable for plotting (ie kick out the scale())
plotmlr <- lm(logNFL ~ RTQUIC*DX_APD + abeta_2, dfnfl) 
summary(plotmlr) 

label <- "paste(''*italic(p), \" < .05\")" #Second annotation is p-value for the interaction

# General layout of the plot: x and y + stats_smooth for linear relationship
fig1d_ver1 <- ggplot(dfnfl, aes(x=abeta_2, y=logNFL)) + #No need for color or fill 

	# Add actual datapoints
	geom_point(aes(color=RTQUIC, shape=DX_APD), size=4) + 
	stat_smooth(aes(), method="lm", color="red", fill="red", linewidth=0.5, alpha=0.2, level = 0.95) + #No need to do separate lines for RTQUIC diagnosis. Lines are very similar and CIs encompass each other. Tried with geom_smooth too. 

	# Fix legends and set up the appearance of the points
	scale_color_manual(values=cbPalette_RTQUIC, name="ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
	scale_shape_manual(values=c(16,17), name="Diagnosis", breaks=c("CBS", "PSP"), labels=c("CBS", "PSP")) + #Color for the datapoints

	# Fix labs
	labs(title=expression(bold("Linear relationship between "*Alpha*beta*"42 and NfL")),
	subtitle=(""),
	x=expression(bold("CSF A"*beta*"42 levels (pg/mL)")),
	y=expression(bold("CSF NfL levels (pg/mL) (log)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

	# Annotate:
	stat_regline_equation(label.x=120, label.y=8.7, color="red", size=6, show.legend=FALSE) + #shows the equation for each geom_line. Cannot use other options since it would not be based on the model
					
	annotate("text", size=6, x=175, y=8.5, color="red", label=label, parse=TRUE)+ #label as defined above. Parse allows for use of mathematical notation

	# Aesthetic only
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))
		
	fig1d_ver1

ggsave(fig1d_ver1, filename = "Fig1d_ver1.png", bg= "transparent", width=9, height=10)



cat("-----------------------   GOES IN FIGURE 1D: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

#Option 2: Plot boxplot of NfL between RTQUIC/DX. 
# General layout of the plot: boxplot by diagnosis where diagnosis is sig, within diagnosis rtquic groups are not sig

fig1d_ver2 <- ggplot(dfnfl, aes(x=DX_APD, y=logNFL, color=RTQUIC))+ #No need for color or fill 

	#Add actual datapoints
	geom_boxplot() +
	geom_jitter(aes(color=RTQUIC), position=position_jitterdodge()) +
	 
	# Fix legends and set up the appearance of the points
	scale_color_manual(values=cbPalette_RTQUIC, name="ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
	scale_fill_manual(values=cbPalette_RTQUIC, name="ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
	scale_shape_manual(values=c(16,17), name="Diagnosis", breaks=c("CBS", "PSP"), labels=c("CBS", "PSP")) + #Color for the datapoints

	# Fix labs
	labs(title=expression("NfL levels by diagnosis and "*alpha*"Syn-SAA status"),
		subtitle="",
		y=expression(bold("CSF NfL levels (pg/mL) (log)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below


	# Annotate:
	geom_signif(y_position = c(9), xmin = c(1), xmax = c(2), annotation = c("p < .05"), fontface = "italic",
				col="black", size=1.1, textsize=6) +
	
	#Aesthetic only #Needs to be after Annotate() otherwise messes with the size of the asterisk
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14)) 

	fig1d_ver2

ggsave(fig1d_ver2, filename = "Fig1d_ver2.png", bg= "transparent", width=9, height=10)
 


cat("\n\n#######################################################################################################\n",
	"                                    8.2. ASYN-SAA+ & SYMPTOMS \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "##########################         8.2.1. CBS-ONLY NUMERICAL VARIABLES      ###########################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# CBS ONLY STATISTICS: DISTRIBUTION
shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA positive",  ]$Age) #normal
shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA negative",  ]$Age) #normal
var.test(Age ~ RTQUIC, data = CBSdf) #homoscedasticity

shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA positive",  ]$Onset_age) #normal
shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA negative",  ]$Onset_age) #normal
var.test(Onset_age ~ RTQUIC, data = CBSdf) #homoscedasticity

shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA positive",  ]$Park_onset) #normal
shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA negative",  ]$Park_onset) #normal
var.test(Park_onset ~ RTQUIC, data = CBSdf) #homoscedasticity

shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA positive",  ]$Age) #normal
shapiro.test(CBSdf[CBSdf$RTQUIC=="aSyn-SAA negative",  ]$Age) #normal
var.test(Age ~ RTQUIC, data = CBSdf) #homoscedasticity

# CBS ONLY STATISTICS: SUMMARY
CBSdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
CBSdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))
CBSdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# CBS ONLY STATISTICS: T-TEST
t.test(CBSdf$Age ~ CBSdf$RTQUIC, var.equal=TRUE) 
t.test(CBSdf$Onset_age ~ CBSdf$RTQUIC, var.equal=TRUE) 
t.test(CBSdf$Park_onset ~ CBSdf$RTQUIC, var.equal=TRUE) 


cat("\n\n######################################################################################################\n",
	   "##########################         8.2.2. CBS-ONLY CATEGORICAL VARIABLES      #########################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# CBS ONLY STATISTICS: SUMMARY
CBSdf %>% group_by(RTQUIC) %>% count(Sex)
CBSdf %>% group_by(RTQUIC) %>% count(APOEe4)
CBSdf %>% group_by(RTQUIC) %>% count(Parkinsonian_onset)
CBSdf %>% group_by(RTQUIC) %>% count(Tremor_binary)
CBSdf %>% group_by(RTQUIC) %>% count(RestTremor) #Very rare symptom (<3 in both groups) so no stats comparison within CBS (ie CBS-aSyn+ vs CBS-aSyn-). 
CBSdf %>% group_by(RTQUIC) %>% count(LimbRigidity)
CBSdf %>% group_by(RTQUIC) %>% count(Slowness_binary)
CBSdf %>% group_by(RTQUIC) %>% count(LP2_falls_PI)
CBSdf %>% group_by(RTQUIC) %>% count(LP2_gait)
CBSdf %>% group_by(RTQUIC) %>% count(RBD_binary) #Very rare symptom (<3 in both groups) so no stats comparison within CBS (ie CBS-aSyn+ vs CBS-aSyn-). 
CBSdf %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
CBSdf %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
CBSdf %>% group_by(RTQUIC) %>% count(Constipation_binary)
CBSdf %>% group_by(RTQUIC) %>% count(Sexual_binary) #Very rare symptom (<3 in both groups) so no stats comparison within CBS (ie CBS-aSyn+ vs CBS-aSyn-). 
CBSdf %>% group_by(RTQUIC) %>% count(Urinary_binary)
CBSdf %>% group_by(RTQUIC) %>% count(Bowel_binary)
CBSdf %>% group_by(RTQUIC) %>% count(anyPPA)

# CBS ONLY STATISTICS: CHISQUARE OR FISHER TEST
if (chooseX2.func(table(CBSdf$Sex, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Sex, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Sex, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$APOEe4, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$APOEe4, CBSdf$RTQUIC))
} else fisher.test(CBSdf$APOEe4, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Parkinsonian_onset, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Parkinsonian_onset, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Parkinsonian_onset, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Tremor_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Tremor_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Tremor_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$LimbRigidity, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$LimbRigidity, CBSdf$RTQUIC))
} else fisher.test(CBSdf$LimbRigidity, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Slowness_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Slowness_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Slowness_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$LP2_falls_PI, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$LP2_falls_PI, CBSdf$RTQUIC))
} else fisher.test(CBSdf$LP2_falls_PI, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$LP2_gait, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$LP2_gait, CBSdf$RTQUIC))
} else fisher.test(CBSdf$LP2_gait, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Lifetime_Dopa_responder_true, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Lifetime_Dopa_responder_true, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Lifetime_Dopa_responder_true, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Lifetime_VisualHallucinations_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Lifetime_VisualHallucinations_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Lifetime_VisualHallucinations_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Constipation_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Constipation_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Constipation_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Urinary_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Urinary_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Urinary_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$Bowel_binary, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$Bowel_binary, CBSdf$RTQUIC))
} else fisher.test(CBSdf$Bowel_binary, CBSdf$RTQUIC)

if (chooseX2.func(table(CBSdf$anyPPA, CBSdf$RTQUIC)) == "chisquare") {
	chisq.test(table(CBSdf$anyPPA, CBSdf$RTQUIC))
} else fisher.test(CBSdf$anyPPA, CBSdf$RTQUIC)


cat("\n\n######################################################################################################\n",
	   "##########################         8.2.3. PSP-ONLY NUMERICAL VARIABLES      ###########################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# PSP ONLY STATISTICS: DISTRIBUTION
shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA positive",  ]$Age) #normal
shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA negative",  ]$Age) #normal
var.test(Age ~ RTQUIC, data = PSPdf) #homoscedasticity

shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA positive",  ]$Onset_age) #normal
shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA negative",  ]$Onset_age) #normal
var.test(Onset_age ~ RTQUIC, data = PSPdf) #homoscedasticity

shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA positive",  ]$Park_onset) #normal
shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA negative",  ]$Park_onset) #normal
var.test(Park_onset ~ RTQUIC, data = PSPdf) #homoscedasticity

shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA positive",  ]$Age) #normal
shapiro.test(PSPdf[PSPdf$RTQUIC=="aSyn-SAA negative",  ]$Age) #normal
var.test(Age ~ RTQUIC, data = PSPdf) #homoscedasticity

# PSP ONLY STATISTICS: SUMMARY
PSPdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
PSPdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset_age, na.rm=T),2),2), sd=sd(Onset_age, na.rm=T))
PSPdf %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# PSP ONLY STATISTICS: T-TEST
t.test(PSPdf$Age ~ PSPdf$RTQUIC, var.equal=TRUE) 
t.test(PSPdf$Onset_age ~ PSPdf$RTQUIC, var.equal=TRUE) 
t.test(PSPdf$Park_onset ~ PSPdf$RTQUIC, var.equal=TRUE) 

cat("\n\n######################################################################################################\n",
	   "##########################         8.2.4. PSP-ONLY CATEGORICAL VARIABLES      #########################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# PSP ONLY STATISTICS: SUMMARY
PSPdf %>% group_by(RTQUIC) %>% count(Sex)
PSPdf %>% group_by(RTQUIC) %>% count(APOEe4)
PSPdf %>% group_by(RTQUIC) %>% count(Parkinsonian_onset)
PSPdf %>% group_by(RTQUIC) %>% count(Tremor_binary)
PSPdf %>% group_by(RTQUIC) %>% count(RestTremor) #Very rare symptom (<3 in both groups) so no stats comparison within PSP (ie PSP-aSyn+ vs PSP-aSyn-). 
PSPdf %>% group_by(RTQUIC) %>% count(LimbRigidity)
PSPdf %>% group_by(RTQUIC) %>% count(Slowness_binary)
PSPdf %>% group_by(RTQUIC) %>% count(LP2_falls_PI)
PSPdf %>% group_by(RTQUIC) %>% count(LP2_gait)
PSPdf %>% group_by(RTQUIC) %>% count(RBD_binary)
PSPdf %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
PSPdf %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
PSPdf %>% group_by(RTQUIC) %>% count(Constipation_binary)
PSPdf %>% group_by(RTQUIC) %>% count(Sexual_binary) #Very rare symptom (<3 in both groups) so no stats comparison within PSP (ie PSP-aSyn+ vs PSP-aSyn-). 
PSPdf %>% group_by(RTQUIC) %>% count(Urinary_binary)
PSPdf %>% group_by(RTQUIC) %>% count(Bowel_binary)

# PSP ONLY STATISTICS: CHISQUARE OR FISHER TEST
if (chooseX2.func(table(PSPdf$Sex, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Sex, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Sex, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$APOEe4, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$APOEe4, PSPdf$RTQUIC))
} else fisher.test(PSPdf$APOEe4, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Parkinsonian_onset, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Parkinsonian_onset, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Parkinsonian_onset, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Tremor_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Tremor_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Tremor_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$LimbRigidity, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$LimbRigidity, PSPdf$RTQUIC))
} else fisher.test(PSPdf$LimbRigidity, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Slowness_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Slowness_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Slowness_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$LP2_falls_PI, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$LP2_falls_PI, PSPdf$RTQUIC))
} else fisher.test(PSPdf$LP2_falls_PI, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$LP2_gait, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$LP2_gait, PSPdf$RTQUIC))
} else fisher.test(PSPdf$LP2_gait, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$RBD_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$RBD_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$RBD_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Lifetime_Dopa_responder_true, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Lifetime_Dopa_responder_true, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Lifetime_Dopa_responder_true, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Lifetime_VisualHallucinations_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Lifetime_VisualHallucinations_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Lifetime_VisualHallucinations_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Constipation_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Constipation_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Constipation_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Urinary_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Urinary_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Urinary_binary, PSPdf$RTQUIC)

if (chooseX2.func(table(PSPdf$Bowel_binary, PSPdf$RTQUIC)) == "chisquare") {
	chisq.test(table(PSPdf$Bowel_binary, PSPdf$RTQUIC))
} else fisher.test(PSPdf$Bowel_binary, PSPdf$RTQUIC)


cat("\n\n######################################################################################################\n",
	   "##########################         8.2.5. BOTH DX NUMERICAL VARIABLES      ############################\n",
	   "#######################################################################################################\n")

cat("-------------------------            GET FROM SECTIONS 5.2.   		  ------------------------------------\n")

cat("\n\n######################################################################################################\n",
	   "##########################         8.2.6. BOTH DX CATEGORICAL VARIABLES    ############################\n",
	   "#######################################################################################################\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

# BOTH DX STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% count(Sex)
df %>% group_by(RTQUIC) %>% count(APOEe4)
df %>% group_by(RTQUIC) %>% count(Parkinsonian_onset)
df %>% group_by(RTQUIC) %>% count(Tremor_binary)
df %>% group_by(RTQUIC) %>% count(RestTremor) 
df %>% group_by(RTQUIC) %>% count(LimbRigidity)
df %>% group_by(RTQUIC) %>% count(Slowness_binary)
df %>% group_by(RTQUIC) %>% count(LP2_falls_PI)
df %>% group_by(RTQUIC) %>% count(LP2_gait)
df %>% group_by(RTQUIC) %>% count(RBD_binary)
df %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
df %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
df %>% group_by(RTQUIC) %>% count(Constipation_binary)
df %>% group_by(RTQUIC) %>% count(Sexual_binary) 
df %>% group_by(RTQUIC) %>% count(Urinary_binary)
df %>% group_by(RTQUIC) %>% count(Bowel_binary)


# BOTH DX STATISTICS: CHISQUARE OR FISHER TEST
if (chooseX2.func(table(df$Sex, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Sex, df$RTQUIC))
} else fisher.test(df$Sex, df$RTQUIC)

if (chooseX2.func(table(df$APOEe4, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$APOEe4, df$RTQUIC))
} else fisher.test(df$APOEe4, df$RTQUIC)

if (chooseX2.func(table(df$Parkinsonian_onset, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Parkinsonian_onset, df$RTQUIC))
} else fisher.test(df$Parkinsonian_onset, df$RTQUIC)

if (chooseX2.func(table(df$Tremor_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Tremor_binary, df$RTQUIC))
} else fisher.test(df$Tremor_binary, df$RTQUIC)

if (chooseX2.func(table(df$LimbRigidity, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$LimbRigidity, df$RTQUIC))
} else fisher.test(df$LimbRigidity, df$RTQUIC)

if (chooseX2.func(table(df$Slowness_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Slowness_binary, df$RTQUIC))
} else fisher.test(df$Slowness_binary, df$RTQUIC)

if (chooseX2.func(table(df$LP2_falls_PI, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$LP2_falls_PI, df$RTQUIC))
} else fisher.test(df$LP2_falls_PI, df$RTQUIC)

if (chooseX2.func(table(df$LP2_gait, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$LP2_gait, df$RTQUIC))
} else fisher.test(df$LP2_gait, df$RTQUIC)

if (chooseX2.func(table(df$RBD_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$RBD_binary, df$RTQUIC))
} else fisher.test(df$RBD_binary, df$RTQUIC)

if (chooseX2.func(table(df$Lifetime_Dopa_responder_true, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Lifetime_Dopa_responder_true, df$RTQUIC))
} else fisher.test(df$Lifetime_Dopa_responder_true, df$RTQUIC)

if (chooseX2.func(table(df$Lifetime_VisualHallucinations_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Lifetime_VisualHallucinations_binary, df$RTQUIC))
} else fisher.test(df$Lifetime_VisualHallucinations_binary, df$RTQUIC)

if (chooseX2.func(table(df$Constipation_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Constipation_binary, df$RTQUIC))
} else fisher.test(df$Constipation_binary, df$RTQUIC)

if (chooseX2.func(table(df$Sexual_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Sexual_binary, df$RTQUIC))
} else fisher.test(df$Sexual_binary, df$RTQUIC)

if (chooseX2.func(table(df$Urinary_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Urinary_binary, df$RTQUIC))
} else fisher.test(df$Urinary_binary, df$RTQUIC)

if (chooseX2.func(table(df$Bowel_binary, df$RTQUIC)) == "chisquare") {
	chisq.test(table(df$Bowel_binary, df$RTQUIC))
} else fisher.test(df$Bowel_binary, df$RTQUIC)


cat("\n\n######################################################################################################\n",
	   "###########################              8.2.7. RADAR PLOTS      ######################################\n",
	   "#######################################################################################################\n")

cat("-----------------------   GOES IN FIGURE 1C: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")


# START WITH CBS first because I'd like to overlay two alpha syn on the same diagnosis plot. 
List = list() #Create an empty list to which you will assign the values from the for-loop. By creating this list outside of the main for-loop, you allow for the data to be entered under different entries which means you can call all the values you need. 

rt.value <- c("aSyn-SAA positive", "aSyn-SAA negative")
for (rt in rt.value) { #for-loop that tests each value of the variable

rtCBSdf <- CBSdf[CBSdf$RTQUIC== rt, ] #Within subset of CBS for eg, look for subset of RT+

print(rt)

#SUMMARIZE THE CATEGORICAL VARIABLES INTO ONE SINGLE COUNT OF "YES"
Tremor_perc <- rtCBSdf %>% summarise(Tremor_count= sum(Tremor_binary == "Yes")) %>% mutate(Tremor_count= (as.numeric(Tremor_count)/nrow(rtCBSdf)*100))
RestTremor_perc <- rtCBSdf %>% summarise(RestTremor_count= sum(RestTremor == "Yes")) %>% mutate(RestTremor_count= (as.numeric(RestTremor_count)/nrow(rtCBSdf)*100))
LimbRigidity_perc <- rtCBSdf %>%  summarise(LimbRigidity_count= sum(LimbRigidity == "Yes")) %>% mutate(LimbRigidity_count= (as.numeric(LimbRigidity_count)/nrow(rtCBSdf)*100))
Slowness_perc <- rtCBSdf %>% summarise(Slowness_count= sum(Slowness_binary == "Yes")) %>% mutate(Slowness_count= (as.numeric(Slowness_count)/nrow(rtCBSdf)*100))
Apraxia_perc <- rtCBSdf %>% summarise(Apraxia_count= sum(LP2_apraxia == "Yes")) %>% mutate(Apraxia_count= (as.numeric(Apraxia_count)/nrow(rtCBSdf)*100))
Gait_perc <- rtCBSdf %>%  summarise(Gait_count= sum(LP2_gait == "Yes")) %>% mutate(Gait_count= (as.numeric(Gait_count)/nrow(rtCBSdf)*100))
FallsPI_perc <- rtCBSdf %>%  summarise(FallsPI_count= sum(LP2_falls_PI == "Yes")) %>% mutate(FallsPI_count= (as.numeric(FallsPI_count)/nrow(rtCBSdf)*100))

#Assign these counts to a variable
Tremor_perc <- Tremor_perc[, 1]
RestTremor_perc <- RestTremor_perc[, 1]
LimbRigidity_perc <- LimbRigidity_perc[, 1]
Slowness_perc <- Slowness_perc[, 1]
Apraxia_perc <- Apraxia_perc[, 1]
Gait_perc <- Gait_perc[, 1]
FallsPI_perc <- FallsPI_perc[, 1]

#Save each value in a list which grows with each iteration of the for-loop (just need to be mindful of the values you are entering here)
List[[length(List)+1]] = c(Tremor_perc, RestTremor_perc, LimbRigidity_perc, Slowness_perc, Apraxia_perc, Gait_perc, FallsPI_perc)

} #end of loop

# To understand the structure of List: ##[[1]] is CBS asyn - [[2]] is CBS asyn + [[3]] is PSP asyn -[[4]] is PSP asyn +					
print(List[[1]])
print(List[[2]])

radar.rtCBSdf <- data.frame(row.names = c("aSyn-SAA negative", "aSyn-SAA positive"),
     Tremor = c(List[[1]][[1]], List[[2]][[1]]),
     Rest = c(List[[1]][[2]], List[[2]][[2]]),
     Limb = c(List[[1]][[3]], List[[2]][[3]]),
     Slowness = c(List[[1]][[4]], List[[2]][[4]]),
     Apraxia = c(List[[1]][[5]], List[[2]][[5]]),
     Gait = c(List[[1]][[6]], List[[2]][[6]]),
     Falls = c(List[[1]][[7]], List[[2]][[7]]))

max_mindf <- data.frame(Tremor = c(100, 0), Rest = c(100, 0), Limb = c(100, 0), Slowness = c(100, 0), Apraxia = c(100, 0), Gait = c(100, 0), Falls = c(100, 0))
rownames(max_mindf) <- c("Max", "Min")

# Bind the variable ranges to the data
radar.rtCBSdf <- rbind(max_mindf, radar.rtCBSdf)
radar.rtCBSdf

#Rename some variables for presentation purposes
colnames(radar.rtCBSdf)[which(names(radar.rtCBSdf) == "Rest")] <- "Rest tremor"
colnames(radar.rtCBSdf)[which(names(radar.rtCBSdf) == "Limb")] <- "Limb rigidity"
colnames(radar.rtCBSdf)[which(names(radar.rtCBSdf) == "Falls")] <- "Falls & instability"
colnames(radar.rtCBSdf)[which(names(radar.rtCBSdf) == "Gait")] <- "Gait \n p<0.1"

# Create the radar charts
png(filename ="Fig1c_CBS.png")

create_beautiful_radarchart(data= radar.rtCBSdf, color= cbPalette_RTQUIC, vlcex=1, plty=1, title="Motor symptoms in CBS")
		# Add an horizontal legend
		# legend(-0.65, -1.2, legend=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-")), horiz=TRUE, bty= "o", pch= 15 , col= cbPalette_RTQUIC, text.col= "black", cex= 1, pt.cex= 1.5)


#NOW PSP. 

List = list() #Create an empty list to which you will assign the values from the for-loop. By creating this list outside of the main for-loop, you allow for the data to be entered under different entries which means you can call all the values you need. 

rt.value <- c("aSyn-SAA positive", "aSyn-SAA negative")
for (rt in rt.value) { #for-loop that tests each value of the variable

	rtPSPdf <- PSPdf[PSPdf$RTQUIC== rt, ] #Within subset of CBS for eg, look for subset of RT+

print(rt)

#SUMMARIZE THE CATEGORICAL VARIABLES INTO ONE SINGLE COUNT OF "YES"
Tremor_perc <- rtPSPdf %>% summarise(Tremor_count= sum(Tremor_binary == "Yes")) %>% mutate(Tremor_count= (as.numeric(Tremor_count)/nrow(rtPSPdf)*100))
RestTremor_perc <- rtPSPdf %>% summarise(RestTremor_count= sum(RestTremor == "Yes")) %>% mutate(RestTremor_count= (as.numeric(RestTremor_count)/nrow(rtPSPdf)*100))
LimbRigidity_perc <- rtPSPdf %>%  summarise(LimbRigidity_count= sum(LimbRigidity == "Yes")) %>% mutate(LimbRigidity_count= (as.numeric(LimbRigidity_count)/nrow(rtPSPdf)*100))
AxialRigidity_perc <- rtPSPdf %>%  summarise(AxialRigidity_count= sum(AxialRigidity == "Yes")) %>% mutate(AxialRigidity_count= (as.numeric(AxialRigidity_count)/nrow(rtPSPdf)*100))
Slowness_perc <- rtPSPdf %>% summarise(Slowness_count= sum(Slowness_binary == "Yes")) %>% mutate(Slowness_count= (as.numeric(Slowness_count)/nrow(rtPSPdf)*100))
OM_perc <- rtPSPdf %>% summarise(OM_count= sum(VerticalOM == "Yes")) %>% mutate(OM_count= (as.numeric(OM_count)/nrow(rtPSPdf)*100))
Gait_perc <- rtPSPdf %>%  summarise(Gait_count= sum(LP2_gait == "Yes")) %>% mutate(Gait_count= (as.numeric(Gait_count)/nrow(rtPSPdf)*100))
FallsPI_perc <- rtPSPdf %>%  summarise(FallsPI_count= sum(LP2_falls_PI == "Yes")) %>% mutate(FallsPI_count= (as.numeric(FallsPI_count)/nrow(rtPSPdf)*100))

#Assign these counts to a variable
Tremor_perc <- Tremor_perc[, 1]
RestTremor_perc <- RestTremor_perc[, 1]
LimbRigidity_perc <- LimbRigidity_perc[, 1]
AxialRigidity_perc <- AxialRigidity_perc[, 1]
Slowness_perc <- Slowness_perc[, 1]
OM_perc <- OM_perc[, 1]
Gait_perc <- Gait_perc[, 1]
FallsPI_perc <- FallsPI_perc[, 1]

# #Save each value in a list which grows with each iteration of the for-loop (just need to be mindful of the values you are entering here)
List[[length(List)+1]] = c(Tremor_perc, RestTremor_perc, LimbRigidity_perc, AxialRigidity_perc, Slowness_perc, OM_perc, Gait_perc, FallsPI_perc)


} #end of loop

# To understand the structure of List: ##[[1]] is CBS asyn - [[2]] is CBS asyn + [[3]] is PSP asyn -[[4]] is PSP asyn +					
print(List[[1]])
print(List[[2]])

radar.rtPSPdf <- data.frame(row.names = c("aSyn-SAA negative", "aSyn-SAA positive"),
     Tremor = c(List[[1]][[1]], List[[2]][[1]]),
     Rest = c(List[[1]][[2]], List[[2]][[2]]),
     Limb = c(List[[1]][[3]], List[[2]][[3]]),
     Axial = c(List[[1]][[4]], List[[2]][[4]]),
     Slowness = c(List[[1]][[5]], List[[2]][[5]]),
     OM = c(List[[1]][[6]], List[[2]][[6]]),
     Gait = c(List[[1]][[7]], List[[2]][[7]]),
     Falls = c(List[[1]][[8]], List[[2]][[8]]))

max_mindf <- data.frame(Tremor= c(100, 0), Rest= c(100, 0), Limb= c(100, 0), Axial=c(100, 0), Slowness = c(100, 0), OM= c(100, 0), Gait= c(100, 0), Falls= c(100, 0))
rownames(max_mindf) <- c("Max", "Min")

# Bind the variable ranges to the data
radar.rtPSPdf <- rbind(max_mindf, radar.rtPSPdf)
radar.rtPSPdf

#Rename some variables for presentation purposes
colnames(radar.rtPSPdf)[which(names(radar.rtPSPdf) == "Rest")] <- "Rest tremor"
colnames(radar.rtPSPdf)[which(names(radar.rtPSPdf) == "Limb")] <- "Limb rigidity \n p<0.1"
colnames(radar.rtPSPdf)[which(names(radar.rtPSPdf) == "Axial")] <- "Axial rigidity"
colnames(radar.rtPSPdf)[which(names(radar.rtPSPdf) == "OM")] <- "Oculomotor"
colnames(radar.rtPSPdf)[which(names(radar.rtPSPdf) == "Falls")] <- "Falls & instability"


# Create the radar charts
png(filename ="Fig1c_PSP.png")

create_beautiful_radarchart(data= radar.rtPSPdf, color= cbPalette_RTQUIC, vlcex=1, plty=1, title="Motor symptoms in PSP")
		# Add an horizontal legend
		# legend(-0.65, -1.2, legend=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-")), horiz=TRUE, bty= "o", pch= 15 , col= cbPalette_RTQUIC, text.col= "black", cex= 1, pt.cex= 1.5)

dev.off()#Not needed?


cat("\n\n\n\n###################################################################################################\n",
		   "9. BINARY LOGISTIC REGRESSION\n",
	 	   "####################################################################################################\n\n")

cat("-------------------------   GOES IN TEXT - RESULTS - LOGISTIC REGRESSION  --------------------------------\n")

# DEFENSE
if (sum(df$RTQUIC_BLR ==1) != 22) {
	cat("There is an issue with the binarization of RTQUIC variable for binary logistic regression \n")
}


# BLR STATS: LOGISTIC REGRESSION MODEL SELECTION
# Model is selected based on previous analyses (logabeta*Onset) + the simple comparisons in Supp material (Gait/RBD_binary).
# Model was run with and without AD/DX_APD status and all values remained very similar. 

dfblr<- df[-which(is.na(df$logNFL)), ] #Remove NAs in df for the predictor NfL (anyway would have been automatically excluded)

blr <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial"(logit))
summary(blr)
AIC(blr) #lowest AIC but Dx and NFL are not significant. So should consider excluding them. 

blr2 <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age) + scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial"(logit))
AIC(blr2) #Slightly better to keep interaction term

blr3 <- glm(RTQUIC_BLR ~  scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial"(logit))
AIC(blr3) #Slightly better with DX 

blr4 <- glm(RTQUIC_BLR ~  DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait, data= dfblr, family = "binomial"(logit))
AIC(blr4) #Slightly better with DX 

blr5 <- glm(RTQUIC_BLR ~  DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait, data= dfblr, family = "binomial"(logit))
AIC(blr5) #Slightly better with NfL

blr6 <- glm(RTQUIC_BLR ~  DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + scale(logNFL), data= dfblr, family = "binomial"(logit))
AIC(blr6) #Better with gait

blr7 <- glm(RTQUIC_BLR ~  DX_APD + scale(Onset_age)*scale(logabeta) + LP2_gait + scale(logNFL), data= dfblr, family = "binomial"(logit))
AIC(blr7) #Better with RBD


#IMPORTANT NOTE: As mentioned in Results, we knwo that the generalizability of this model is likely to be limited due to the relatively low
## sample size consideirng hte number of predictors. Therefore, the model is mostly meant to give additional information. For ex, RBD is a rare
## symptom so it would be very unlikely that the OR would remain stable across multiple iterations of the model. Similarly, we preferred to keep
## outliers for the run, but below the model without these outliers is also shown. 

# BLR STATS: LOGISTIC REGRESSION MODEL DIAGNOSTICS
pscl::pR2(blr)["McFadden"]

# COLLINEARITY CHECK WITH VARIABLE INFLATION FACTOR
# First look at Onset_age and logged Abeta42 separately due to interaction
blronset <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial") #Check individual relationship of Onset with RTQUIC
	car::vif(blronset)
blrabeta <- glm(RTQUIC_BLR ~ DX_APD + scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial") #Check individual relationship of Abeta with RTQUIC
	car::vif(blrabeta)
blr_nointerac <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age)+ scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfblr, family = "binomial") #Check individual relationship of Abeta with RTQUIC
	car::vif(blr_nointerac)

# With interaction and scale
car::vif(blr)

# AUTOCORRELATION RESIDUALS
durbinWatsonTest(blr)

# Variable importance
caret::varImp(blr) #Similar, so all variables are of comparable importance in this model. 

# LINEARITY BETWEEN LOGIT AND IVS
# In logistic regression, we assume the relationship is linear on the logit scale. This is assessed with component-plus-residual plots. The “component” is the values of a variable multiplied by its estimated coefficient (meaning that each predictor has its own component vector), and the “residual” is the working residuals, a type of residuals in generalized linear models.
# https://sscc.wisc.edu/sscc/pubs/RegDiag-R/logistic-regression.html#log_lin
crPlots(blr2) #without the interaction

dfblr |> #Good 
  mutate(residuals = coef(blr)[3]*Onset_age + residuals(blr, type="working")) |> 
  ggplot(aes(x = Onset_age, y = residuals)) +
  geom_point() +
  geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
  geom_smooth(se = F)

dfblr |> #Not linear
  mutate(residuals = coef(blr)[4]*logabeta + residuals(blr, type="working")) |> 
  ggplot(aes(x = logabeta, y = residuals)) +
  geom_point() +
  geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
  geom_smooth(se = F)

dfblr |> #Not linear
  mutate(residuals = coef(blr)[7]*logNFL + residuals(blr, type="working")) |> 
  ggplot(aes(x = scale(logNFL), y = residuals)) +
  geom_point() +
  geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
  geom_smooth(se = F)

# These plots indicate that there are important outliers. Since this analysis is just exploratory and underpowered to begin with, the results below are presented for 
## informative purposes to evaluate the robustness of findings.


# TESTING THE MODEL WITHOUT OUTLIERS OR DIFFERENT ITERATIONS - MODEL SELECTION CONTINUED
# As noted in the manuscript, the binary logistic regression modeling is underpowered due to the heterogeneity of hte sample. 

# REMOVING ABETA OUTLIERS (EVEN NOT IDENTIFIED BY TUKEY)
dfabetablr <- dfabeta[!(dfabeta$logabeta <5), ] #Remove the low Abeta outlier
dfabetablr<- dfabetablr[-which(is.na(dfabetablr$logNFL)), ] #Remove NAs in df for the predictor NfL (anyway would have been automatically excluded)

blrtest <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait + scale(logNFL), data= dfabetablr, family = "binomial"(logit))
summary(blrtest)
AIC(blrtest) #lowest AIC but Dx and NFL are not significant. So should consider excluding them. 

dfabetablr |> 
 	mutate(residuals = coef(blrtest)[3]*Onset_age + residuals(blrtest, type="working")) |> 
 	ggplot(aes(x = Onset_age, y = residuals)) +
  	geom_point() +
  	geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
  	geom_smooth(se = F)

dfabetablr |> 
	mutate(residuals = coef(blrtest)[4]*logabeta + residuals(blrtest, type="working")) |> 
	ggplot(aes(x = logabeta, y = residuals)) +
	geom_point() +
	geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
	geom_smooth(se = F)

dfabetablr |> 
	mutate(residuals = coef(blrtest)[7]*logNFL + residuals(blrtest, type="working")) |> 
	ggplot(aes(x = scale(logNFL), y = residuals)) +
	geom_point() +
	geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
	geom_smooth(se = F)

removed <- setdiff(dfblr, dfabetablr) 
removed[, c("RTQUIC", "AD", "Early_onset", "logabeta")]

cat("The following IDs are removed based on Abeta42 values for the corrected logistic regression model with better diagnostics: ", removed$logabeta, "\n")
cat("The issue is, both are aSyn-SAA+ (bringing down total number of aSyn-SAA+ to 20 instead of 22 (aSyn-SAA+ is the outcome measure for this model),
	but 1 especially is one of the 8 aSyn-SAA+/young-onset, 7 aSyn-SAA+/AD+, and 5 aSyn-SAA+/AD+. We already know from Fig 1A that these combinations are rare 
	but they are of interest in our cohort due to 1. hypothesis 2. results of frequency analyses (categorical tests on AD+ vs Young-onset, so removing
	a subject that meets all these criteria would in itself be an issue. \n")
cat("Removing extreme Abeta42 values (even not identified by Tukey method) improves assumption testing output. It removes important information from the model though as it reduces significantly groups that were already underrepresetned but of interest.
In this iteration, Onset age, Gait, RBD, and NfL are significantly predictors of aSyn-SAA+, but not Abeta 42. ")

# REMOVING VARIABLES:
# NFL and DX can be removed as they are not significant in blr model. However: DX is clinically really importnat hence its inclusion. 
# NFL is significant in some iterations of the model, such as the one above - it could be a significant covariate/confounding factor
# that "competes" with Abeta. 
blrtest2 <- glm(RTQUIC_BLR ~ scale(Onset_age)*scale(logabeta) + RBD_binary + LP2_gait, data= dfblr, family = "binomial"(logit))
summary(blrtest2)
AIC(blrtest2) #lowest AIC but Dx and NFL are not significant. So should consider excluding them. 

dfblr |> 
	mutate(residuals = coef(blrtest2)[3]*Onset_age + residuals(blrtest2, type="working")) |> 
	ggplot(aes(x = Onset_age, y = residuals)) +
	geom_point() +
	geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
	geom_smooth(se = F)

dfblr |> 
	mutate(residuals = coef(blrtest2)[4]*logabeta + residuals(blrtest2, type="working")) |> 
	ggplot(aes(x = logabeta, y = residuals)) +
	geom_point() +
	geom_smooth(color = "red", method = "lm", linetype = 2, se = F) +
	geom_smooth(se = F)

cat("Removing variables overall improves assumption testing output, except at very high age at onset(>80). 
In this iteration, Onset age, Gait, and RBD are significantly predictors of aSyn-SAA+. Onset age by Abeta42 is significant at a 90% confidence level.
It is the best compromise between reducing violation of assumptions and representativity of the sample. It also reduces the issue of overfitting.")

# Conclusion: The logit linearity assumption was violated so we examined multiple possibilities. 
# Transforming data was not appropriate due to the complexity of the model.
# Removing outliers was promising based on graph: Removing NfL outliers (not shown) fixed the assumption for NfL but not Abeta42. 
# Removing outliers of Abeta42 based on Tukey method (eg using dfabeta dataset) did not fix the assumption for either variable (not shown).
# Removing the lowest value of Abeta42 (<5 logged Abeta42 value) in addition to the positive outliers of Abeta42 (not of NfL) reduced violation
# of the linearity assumption. However, the low Abeta42 subject is not an outlier based on Tukey method (ie not unreasonably low Abeta42 compared to other AD+/diagnostic group/RT-QUIC group. 
## On the other hand, this low Abeta42 subject is one of few subjects who are RTQUIC+ (rarer outcome, moreover outcome to be predicted), AD+ (<20% of whole dataset), and young-onset RTQUIC+ (8 people total, of which 5 are AD+). 
# The last approach tested was removing variables that were not significant or of interest from the final model. In that model, the whole dataset could be including without violating the linearity assumption (except at very high age at onset, which is easier to interpret). 

cat("\n")

cat("For the reasons above, no further efforts towards prediction are performed. OR will be given for scale. Visualization of the model is kept here for interpreation
	but not included in the manuscript. Influence of outliers should be highlighted, but to remove them entirely would also be inappropriate. \n")
cat("Consistently significant effect of: Age at onset; Gait; and RBD.\n")

# MODEL COEFFICIENTS (FOR SCALE ONLY)
coef(blr)
exp(coef(blr)) 
cbind(coef(blr),odds_ratio=exp(coef(blr)),exp(confint(blr, level=0.95))) #it says 2.5% and 97.5% because these are the two borders to end up with the cnetral 95% of your distribution

coef(blrtest2)
exp(coef(blrtest2)) 
cbind(coef(blrtest2),odds_ratio=exp(coef(blrtest2)),exp(confint(blrtest2, level=0.95))) #it says 2.5% and 97.5% because these are the two borders to end up with the cnetral 95% of your distribution

#Example of interpretation:
# exp(4.88278) #RBD. 81.60169. Odds of someone with RBD being RTquic+ increases by 8000%. IE x80. 
# exp(-2.85663) #Gait. 0.09253855. Odds of someone with gait issues being RTQUIC+ decreases by 10% compared to someone without. 
##### exp(0.15122) #Age: 1.16. #The value indicates that as age increase by one more unit, then the odds of being SAA+ increases by 16%



cat("\n\n\n\n###################################################################################################\n",
		   "10. RTQUIC PARAMETERS SUPP ANALYSES\n",
	 	   "####################################################################################################\n\n")

cat("\n\n#######################################################################################################\n",
	"                                    10.1. RTQUIC PARAMETERS SUPP: LAG \n",
	   "#######################################################################################################\n")

cat("-------------------------------------------------   GOES IN REVIEW  ---------------------------------------\n")


cat("Lag hours is the #hours required to reach max ThT. It makes the most sense to think about it as data suited for survival analysis. \n")
cat("For that purpose, we are censoring the subjects who never reached positivity (RTQUIC negative) as max ThT is not something that exists in this case 
(they plateau early on in the analysis as their curve never rises - so comparign them to RT-QUIC+ subjects would not make sense. ) \n")

## 1. We want to have a column: Negative vs Positive
## 2. We want to have a column: lag hours max or 40

# LAG STATISTICS: DISTRIBUTION
hist(df$RTQUIC_survival_hours)  
hist(RTposdf$RTQUIC_survival_hours)  
shapiro.test(df[RTposdf$DX_APD =="CBS", ]$RTQUIC_survival_hours) #nonnormal
shapiro.test(df[RTposdf$DX_APD =="PSP", ]$RTQUIC_survival_hours) #nonnormal
shapiro.test(df[RTposdf$Early_onset =="Young-onset", ]$RTQUIC_survival_hours) #nonnormal
shapiro.test(df[RTposdf$Early_onset =="Late-onset", ]$RTQUIC_survival_hours) #nonnormal
shapiro.test(df[RTposdf$AD =="AD Positive", ]$RTQUIC_survival_hours) #nonnormal
shapiro.test(df[RTposdf$AD =="AD Negative", ]$RTQUIC_survival_hours) #nonnormal
leveneTest(RTQUIC_survival_hours ~ DX_APD, data = RTposdf) #homoscedasticity  
leveneTest(RTQUIC_survival_hours ~ Early_onset, data = RTposdf) #homoscedasticity  
leveneTest(RTQUIC_survival_hours ~ AD, data = RTposdf) #homoscedasticity  


# LAG STATISTICS: SUMMARY
RTposdf %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))
RTposdf %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))

RTposdf %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))
RTposdf %>% group_by(AD) %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))

RTposdf %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))
RTposdf %>% group_by(Early_onset) %>% summarize(count=n(), format(round(median(RTQUIC_survival_hours, na.rm=T),2),2), IQR=IQR(RTQUIC_survival_hours, na.rm=T), min=min(RTQUIC_survival_hours, na.rm=T), max=max(RTQUIC_survival_hours, na.rm=T))


# LAG STATISTICS: CORRELATIONS
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$Onset_age, method="spearman") #correlated
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$Age, method="spearman") 
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logabeta, method="spearman") 
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logNFL, method="spearman")
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logptau, method="spearman") #correlated
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logttau, method="spearman")  
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$ATI_2, method="spearman") 


# LAG STATISTICS: CORRELATION PLOTS
# Age at onset
ggscatter(RTposdf, x="Onset_age", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with age at onset (Spearman)",
		xlab="Age at onset (years)", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="Onset_age", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="Age at onset (years)", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

# Abeta42
ggscatter(RTposdf, x="logabeta", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with Abeta42 levels (Spearman)", 
		xlab="Logged Abeta42 (pg/mL)", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="logabeta", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="Logged Abeta42 (pg/mL)", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

# Ptau181
ggscatter(RTposdf, x="logptau", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with ptau181 levels (Spearman)", 
		xlab="Logged ptau181 (pg/mL)", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="logptau", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="Logged ptau181 (pg/mL)", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

# T-tau
ggscatter(RTposdf, x="logttau", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with t-tau levels (Spearman)", 
		xlab="Logged t-tau (pg/mL)", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="logttau", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="Logged t-tau (pg/mL)", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

# NfL
ggscatter(RTposdf, x="logNFL", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with t-tau levels (Spearman)", 
		xlab="Logged NfL (pg/mL)", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="logNFL", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="Logged NfL (pg/mL)", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

# ATI
ggscatter(RTposdf, x="ATI_2", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21,
		palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with t-tau levels (Spearman)", 
		xlab="ATI", ylab = "Lag (hours)") +
		scale_fill_discrete(name = 'Diagnosis',labels = c("CBS", "PSP"))

ggscatter(RTposdf, x="ATI_2", y= "RTQUIC_survival_hours", fill="DX_APD",
		xlab="ATI", ylab = "Lag (hours)",
		add = "reg.line", conf.int = TRUE) + 
		stat_cor(aes(color = DX_APD, label = paste("R^2~=,~")), show.legend = FALSE) +
		scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))


# LAG STATISTICS: LOG-RANK TEST
summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df), times = 48)
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df) #p=.05

summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df), times = 48)
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df) 

summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df), times = 48)
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df) 


# LAG STATISTICS: KAPLAN-MEIER CURVE FOR WHOLE DATASET
# Uses survfit2() in ggsurvfit() package:
survfit2(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df)  %>% 
    ggsurvfit() + add_confidence_interval() + labs(x="Lag (hours)", "Survival probability") #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

survfit2(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df)  %>% 
    ggsurvfit() + add_confidence_interval() + labs(x="Lag (hours)", "Survival probability") #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below


survfit2(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df)  %>% 
    ggsurvfit() + add_confidence_interval() + labs(x="Lag (hours)", "Survival probability") #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

# LAG STATISTICS: ANALYSES ON ONLY THE RT-QUIC SUBJECTS
# Fit a Cox proportional hazards model
# Redundant with the fisher analyses already shown above
fit.coxph <- coxph(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset + AD, data = df)
ggforest(fit.coxph, data = df) # Redundant with previous analyses using chisqaure and fisher approach in the manuscript - but does suggest potential difference between AD+ and AD-.


cat("\n\n#######################################################################################################\n",
	"                                    10.2. RTQUIC PARAMETERS SUPP: THT MAX \n",
	   "#######################################################################################################\n")

cat("THT max is the max fluorescent signal reached after 48 hours of monitoring of the assay. \n")

# THT STATISTICS: DISTRIBUTION
hist(RTposdf$ThTmax)  
shapiro.test(df[RTposdf$DX_APD =="CBS", ]$ThTmax) #normal
shapiro.test(df[RTposdf$DX_APD =="PSP", ]$ThTmax) #normal
shapiro.test(df[RTposdf$Early_onset =="Young-onset", ]$ThTmax) #normal
shapiro.test(df[RTposdf$Early_onset =="Late-onset", ]$ThTmax) #normal
shapiro.test(df[RTposdf$AD =="AD Positive", ]$ThTmax) #normal
shapiro.test(df[RTposdf$AD =="AD Negative", ]$ThTmax) #normal
leveneTest(ThTmax ~ DX_APD, data = RTposdf) #homoscedasticity  
leveneTest(ThTmax ~ Early_onset, data = RTposdf) #homoscedasticity  
leveneTest(ThTmax ~ AD, data = RTposdf) #homoscedasticity  


# THT STATISTICS: SUMMARY
RTposdf %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% group_by(AD) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% group_by(Early_onset) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))

# THT STATISTICS: TTEST
t.test(RTposdf$ThTmax ~ RTposdf$DX_APD, var.equal=TRUE) 
t.test(RTposdf$ThTmax ~ RTposdf$Early_onset, var.equal=TRUE) 
t.test(RTposdf$ThTmax ~ RTposdf$AD, var.equal=TRUE) 

# THT STATISTICS: BARPLOTS
ggplot(RTposdf, aes(x=DX_APD, y=ThTmax, color=DX_APD))+ #No need for color or fill 

	#Add actual datapoints
	geom_boxplot() +
	geom_jitter(aes(color=DX_APD), position=position_jitterdodge()) +
	 
	# Fix legends and set up the appearance of the points
	scale_color_manual(values=cbPalette_DX_APD, name="Diagnosis", breaks=c("CBS", "PSP"), labels=c("CBS", "PSP")) + #expression allows you to add greek letters

	# Fix labs
	labs(title=expression("ThT max by diagnosis"),
		subtitle="",
		x="", y=expression(bold("ThT max (arbitrary fluorescence units)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below
	
	#Aesthetic only #Needs to be after Annotate() otherwise messes with the size of the asterisk
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))  


ggplot(RTposdf, aes(x=AD, y=ThTmax, color=AD))+ #No need for color or fill 

	#Add actual datapoints
	geom_boxplot() +
	geom_jitter(aes(color=AD), position=position_jitterdodge()) +
	 
	# Fix legends and set up the appearance of the points
	scale_color_manual(values=c("red", "green"), name="AD status", breaks=c("AD Positive", "AD Negative"), labels=c("AD Positive", "AD Negative")) + #expression allows you to add greek letters

	# Fix labs
	labs(title=expression("ThT max by AD status"),
		subtitle="",
		x="", y=expression(bold("ThT max (arbitrary fluorescence units)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below
	
	#Aesthetic only #Needs to be after Annotate() otherwise messes with the size of the asterisk
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))  

ggplot(RTposdf, aes(x=Early_onset, y=ThTmax, color=Early_onset))+ #No need for color or fill 

	#Add actual datapoints
	geom_boxplot() +
	geom_jitter(aes(color=Early_onset), position=position_jitterdodge()) +
	 
	# Fix legends and set up the appearance of the points
	scale_color_manual(values=c("hotpink1", "royalblue1"), name="Onset", breaks=c("Young-onset", "Late-onset"), labels=c("Young-onset", "Late-onset")) + #expression allows you to add greek letters

	# Fix labs
	labs(title=expression("ThT max by type of onset"),
		subtitle="",
		x="", y=expression(bold("ThT max (arbitrary fluorescence units)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below
	
	#Aesthetic only #Needs to be after Annotate() otherwise messes with the size of the asterisk
	theme_classic() +
	theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))  


