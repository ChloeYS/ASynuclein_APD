# FILENAME: Analysis_APD.R

#Last updated 16 May 2024

#Usage: sources other files: 
#RScript Test_APD.R /Users/nikhilbhagwat/Desktop/7_PUBLICATIONS_ONGOING/7.1_APD_MS_2024-01/7.1.0_Data/APD_Neurology_MS_data.csv



cat("\n\n\n\n###############################################################################################\n",
            "1. SOURCE PACKAGES AND FUNCTIONS\n",
            "###############################################################################################\n\n\n")

arg.vec <- commandArgs(trailingOnly = T) #Defines a vector of arguments. In this case, it is only one argument.  

source('Functions.R') #Source the functions from file in same repository (for now)
	#tidyverse() loaded in Functions.R

source('DFCreationFunction_APD.R') #Source the function taht creates the new dataframe on which analysis is performed
	#lubridate() loaded in DFCreationFunction_APD.R


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
library(lmtest) #Breusch-Pagan



# Consider removing if not used: 
# library(multcomp) #multiple comparisons
# library(ggpmisc) #annotate figures
# library(modelbased) #extract some parameters and estimates from a model
# library(psych)
# library(ggstatsplot) #plot figures
# library(fmsb) #radar plot



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

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

#Biomarkers: if an outlier is removed, give value in the table

# Abeta42 STATISTICS: OUTLIERS
# First, display the data points of one group. boxplot() will identify outliers in given df but stripchart function can be funky hence convoluted script. 
# Code below plots (boxplot() + assigns outliers to diagnosis-specific vector)
# Outliers are removed before calculating the values that are presented in table. 
boxplot(abeta_2 ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	stripchart(abeta_2 ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	CBSvec.outliers <- boxplot(abeta_2 ~ DX_APD, data= CBSdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
boxplot(abeta_2 ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	stripchart(abeta_2 ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	PSPvec.outliers <- boxplot(abeta_2 ~ DX_APD, data= PSPdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)

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
t.test(df$logabeta ~ df$DX_APD, var.equal=TRUE) 
aov <- aov(logabeta ~ Age + DX_APD, df) 
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


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
	PSPvec.outliers <- boxplot(ptau_2 ~ DX_APD, data= PSPdf, col = "white")$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)

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
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


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
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


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
Anova(aov, type="II") #Compare with type III
	check_normality(aov)


cat("\n\n######################################################################################################\n",
	   "#########################              4.1.9. BIOMARKERS: NFL                 #########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")

# NFL STATISTICS: OUTLIERS
# For NFL the outliers are really impactful. Is best to remove them and report their value. However, will only remove outliers of whole dataset and not per diagnosis. 
# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# This is because in FTLD population, the data can be very right-skewed - so prefer to have representative data/tolerant threshold in the Demographics table. 
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


cat("\n\n#######################################################################################################\n",
	"                                    4.2. COHORT CHARACTERISTICS: CATEGORICAL VARIABLES\n",
	    "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "################################            4.2.1. SEX                         ########################\n",
	   "#######################################################################################################\n")

cat("------------------------------------   GOES IN TABLE 1: DEMOGRAPHICS   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")
cat("-------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-       ------------------------------------\n")

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

# AD STATISTICS: SUMMARY
df %>% group_by(DX_APD) %>% count(AD)
table(df$DX_APD, df$AD)

# AD STATISTICS: CHI-SQUARE
chisq.test(table(df$AD, df$DX_APD), correct=F)




cat("\n\n\n\n###################################################################################################\n",
		   "5. ASYN-SAA+ & COHORT CHARACTERISTICS\n",
	 	   "####################################################################################################\n\n")

cat("-------------------------   GOES IN eTABLE 1: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

df %>% group_by(RTQUIC) %>% count(DX_APD)
sum(df$DX_APD=="CBS" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")
sum(df$DX_APD=="PSP" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")

cat("\n\n#######################################################################################################\n",
	"                         5.1. ASYN-SAA+ COHORT CHARACTERISTICS: CATEGORICAL VARIABLES  \n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")

# SEX STATISTICS: DISTRIBUTION
table(df$Sex, df$RTQUIC)
table(CBSdf$Sex, CBSdf$RTQUIC)
table(PSPdf$Sex, PSPdf$RTQUIC)

# SEX STATISTICS: CHI-SQUARE AND FISHER'S TEST
chisq.test(table(df$Sex, df$RTQUIC), correct=F) 
chisq.test(table(CBSdf$Sex, CBSdf$RTQUIC), correct=F) 
fisher.test(table(PSPdf$APOEe4, PSPdf$RTQUIC)) # Expected count is <5 for one cell

# APOE STATISTICS: DISTRIBUTION
table(df$APOEe4, df$RTQUIC)
table(CBSdf$APOEe4, CBSdf$RTQUIC)
table(PSPdf$APOEe4, PSPdf$RTQUIC)

# APOE STATISTICS: FISHER'S TEST
fisher.test(table(df$APOEe4, df$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(CBSdf$APOEe4, CBSdf$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(PSPdf$APOEe4, PSPdf$RTQUIC)) # Expected count is <5 for one cell

# PARK ONSET STATISTICS: DISTRIBUTION
table(df$Parkinsonian_onset, df$RTQUIC)
table(CBSdf$Parkinsonian_onset, CBSdf$RTQUIC)
table(PSPdf$Parkinsonian_onset, PSPdf$RTQUIC)

# PARK ONSET STATISTICS: CHI-SQUARE AND FISHER'S TEST
chisq.test(table(df$Parkinsonian_onset, df$RTQUIC), correct=F) 
fisher.test(table(CBSdf$Parkinsonian_onset, CBSdf$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(PSPdf$Parkinsonian_onset, PSPdf$RTQUIC)) # Expected count is <5 for one cell


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

# ASYN-SAA*AGE STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Age) #normal
shapiro.test(RTnegdf$Age) #normal
var.test(Age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*AGE STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))
CBSdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))
PSPdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))

# RTQUIC*AGE STATISTICS: ANCOVA
t.test(df$Age ~ df$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Age ~ DX_APD*RTQUIC, df)
Anova(aov, type="III")

# p=.007 < .017 so reported bonf-adjusted p-value: .007x3=.02

cat("\n\n######################################################################################################\n",
	   "######################               5.2.2. ASYN-SAA*ONSET                #############################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")

# ASYN-SAA*ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Onset_age) #normal
shapiro.test(RTnegdf$Onset_age) #normal
var.test(Onset_age ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*ONSET STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset_age, na.rm=T), sd=sd(Onset_age, na.rm=T))
CBSdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset_age, na.rm=T), sd=sd(Age, na.rm=T))
PSPdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset_age, na.rm=T), sd=sd(Age, na.rm=T))

# RTQUIC*ONSET STATISTICS: ANCOVA
t.test(df$Onset_age ~ df$RTQUIC, var.equal=TRUE) #What is reported in eTable 1
aov <- aov(Onset_age ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")

# p >.017 so no need to bonferroni correct as it is ns


cat("\n\n######################################################################################################\n",
	   "################################          5.2.3. ASYN-SAA*PARK_ONSET           ########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 1: ASYN-SAA- vs ASYN-SAA+ ----------------------------------\n")

# RTQUIC*PARK_ONSET STATISTICS: DISTRIBUTION
shapiro.test(RTposdf$Park_onset) #normal
shapiro.test(RTnegdf$Park_onset) #normal
var.test(Park_onset ~ RTQUIC, data = df) #homoscedasticity


# RTQUIC*PARK_ONSET STATISTICS: SUMMARY
df %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Park_onset, na.rm=T))
CBSdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Age, na.rm=T))
PSPdf %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Age, na.rm=T))

# RTQUIC*PARK_ONSET STATISTICS: ANCOVA
t.test(df$Park_onset ~ df$RTQUIC, var.equal=TRUE)
aov <- aov(Park_onset ~ DX_APD*RTQUIC, df)  
Anova(aov, type="III")


# p >.017 so no need to bonferroni correct as it is ns

cat("\n\n\n\n###################################################################################################\n",
		   "6. AD+ & COHORT CHARACTERISTICS\n",
	 	   "####################################################################################################\n\n")

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
	   "#######################################################################################################\n")

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

# For NFL the outliers are really impactful. Is best to remove them and report their value. However, will only remove outliers of whole dataset and not per diagnosis. 
# Remove outliers over full dataset, but at a tolerant threshold (Q3+3*IQR instead of 1.5 IQR. Reference for this is: https://www.nature.com/articles/s41598-020-66090-x 
# This is because in FTLD population, the data can be very right-skewed. 
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
t.test(CBSdfnfl$logNFL ~ CBSdfnfl$AD, var.equal=FALSE) #Welch t-test
gls1 <- gls(logNFL ~ Age*AD, CBSdfnfl, weights=varPower()) 
gls2 <- gls(logNFL ~ Age + AD, CBSdfnfl, weights=varPower()) 
AIC(gls1, gls2) #gls2 better than gls1 
summary(gls2)


cat("\n\n######################################################################################################\n",
	   "###############################           6.4.3. MOCA Z-SCORES                #########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN eTABLE 2: CBS-AD+ vs CBS-AD-     ----------------------------------\n")

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

sum(df$AD=="AD Positive" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")
sum(df$AD=="AD Positive" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")
sum(df$AD=="AD Negative" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="M")
sum(df$AD=="AD Negative" & df$RTQUIC=="aSyn-SAA positive" & df$Sex=="F")


cat("\n\n######################################################################################################\n",
	   "############################              7.1.1.FIGURE 1A           		###########################\n",
	   "#######################################################################################################\n")

cat("-----------------------   GOES IN FIGURE 1A: ASYN-SAA+ vs ASYN-SAA-   ------------------------------------\n")

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
cramerV(table(LODdf$AD, LODdf$RTQUIC)) # Expected count is <5 for one cell
phi(table(LODdf$AD, LODdf$RTQUIC), digits=6) # Expected count is <5 for one cell


cat("\n\n#######################################################################################################\n",
	"                                    7.2. ASYN-SAA+ & AD: NUMERICAL VARIABLES \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "########################              7.2.1. ASYN-SAA+ & ABETA42          	###########################\n",
	   "#######################################################################################################\n")

cat("----------------------------   GOES IN TABLE 2: MLR MODEL OUTPUT   --------------------------------------\n")

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
leveneTest(logabeta ~ DX_APD*RTQUIC, df) #homo

# ABETA STATS: LINEAR REGRESSION MODEL SELECTION
test1 <- lm(logabeta ~ RTQUIC, df)
test2 <- lm(logabeta ~ RTQUIC  + DX_APD, df)
test3 <- lm(logabeta ~ RTQUIC + Sex, df) #not adding any value to the model
anova(test1, test2) 
anova(test1, test3) 

## Inclusion of Onset as a covariate. Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
ggscatter(df, x = "Onset_age", y = "logabeta", add = "reg.line")+ 
	stat_regline_equation(aes()) 
ggscatter(df, x = "Onset_age", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
stat_regline_equation(aes(color = DX_APD)) 
ggscatter(df, x = "Onset_age", y = "logabeta", color = "RTQUIC", add = "reg.line")+
	stat_regline_equation(aes(color = RTQUIC))
cor.test(df$logabeta, df$Onset_age) #corr
summary(lm(df$logabeta ~ df$Onset_age)) #linear relationship 

## Inclusion of NFL as a covariate. Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
ggscatter(df, x = "NFL_2", y = "logabeta", add = "reg.line")+ 
	stat_regline_equation(aes()) 
ggscatter(df, x = "NFL_2", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
	stat_regline_equation(aes(color = DX_APD)) 
ggscatter(df, x = "NFL_2", y = "logabeta", color = "RTQUIC", add = "reg.line")+
	stat_regline_equation(aes(color = RTQUIC))
cor.test(df$logabeta, df$NFL_2) #corr
summary(lm(df$logabeta ~ df$NFL_2)) #linear relationship 

## Inclusion of other variables as covariate:
summary(lm(df$logabeta ~ df$LP2_Disease_Duration)) #no linear relationship
summary(lm(df$logabeta ~ df$ttau_2)) #no linear relationship
summary(lm(df$logabeta ~ df$ptau_2)) #no linear relationship

# ABETA STATS: MULTIVARIABLE LINEAR REGRESSION MODEL
mlr <- lm(logabeta ~ Onset_age*RTQUIC + DX_APD + NFL_2, df) 
summary(mlr)

stdmlr <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD + scale(NFL_2), df) 
summary(stdmlr)

## Diagnostics of the model run
check_normality(stdmlr) #Ok normality of residuals
autoplot(stdmlr, which = 1:6) #All plots including scale location plot

# Visualize the values for each of the IDs that are indexed on above plots
df[5, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[21, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[27, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[37, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[45, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
df[65, c("DX_APD", "abeta_2", "logabeta", "Onset_age", "RTQUIC", "NFL_2")] 
		
## for-loop to test the model without each of these values (ie once without outlier 1, then outlier 2, etc)
vecIDs <- df[c(5,21,27,37,45,65), "ID"] #Create vector of each ID
for (i in vecIDs) {
	test <- subset(df, ID!=i) #for some analysis, need to exclude the potential false negative too
	teststdmlr <- lm(logabeta ~ scale(Onset_age)*RTQUIC + DX_APD + scale(NFL_2), test) 
	print(summary(teststdmlr))
	}

plot(stdmlr, which = 3) # 3 = Scale-Location plot. Variance of residuals
bptest(stdmlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
durbinWatsonTest(stdmlr) #Ok autocorrelation of residuals	

## Multicollinearity checks
car::vif(stdmlr) #no multicollinearity at all

## Simple slopes for onset: 
emtrends(stdmlr,pairwise ~  RTQUIC, var="Onset_age")

## Main effects: 
emmeans(stdmlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
emmeans(stdmlr, ~ DX_APD) #adjusted means. Ok because no interaction. 


cat("\n\n######################################################################################################\n",
	   "###########################                 7.2.2. FIGURE 1.B.         	##############################\n",
	   "#######################################################################################################\n")

# FIG1B: ABETA42 over time from the real datapoints

##Create df for the figure based on the actual model
mylist <- list(Onset_age=seq(35,85,by=5), RTQUIC=c("aSyn-SAA positive", "aSyn-SAA negative")) 
emmip(stdmlr, RTQUIC ~ Onset_age, at=mylist, CIs=TRUE)

#Dataset: data are from the model (predicted values)
figdf <- emmip(stdmlr, RTQUIC ~ Onset_age, at=mylist, CIs=TRUE, plotit=FALSE) #
    
label1 <- "paste(F*'(5, 59)'==4.29*', ' ~~ italic(p), \"< .01\"*', ' ~~ italic(R)^2==20.47*'%')" #First annotation of the plot: Model diagnostics. Comma has to be entered as text not in mathematical notation. 
label2 <- "paste(''*italic(p), \" < .05\")" #Second annotation is p-value for the interaction

 # Change name of variable RTQUIC
fig1b_ver1 <- ggplot(data=figdf, aes(x=Onset_age,y=yvar, color=RTQUIC)) + #yvar is Abeta42 logged 	#Ggplot figure basic layout

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

fig1b_ver1

ggsave(fig1b_ver1, filename = "Fig1b_ver1.png", bg= "transparent", width=9, height=10)


cat("\n\n\n\n###################################################################################################\n",
		   "8. ASYN-SAA+ & ALL VARIABLES\n",
	 	   "####################################################################################################\n\n")

cat("\n\n#######################################################################################################\n",
	"                                    8.1. ASYN-SAA+ & NFL \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "###########################              8.1.1. NFL MLR        	#######################################\n",
	   "#######################################################################################################\n")


# NFL STATS: DISTRIBUTION
boxplot(logNFL ~ DX_APD, data= df)$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
boxplot(logNFL ~ DX_APD, data= df)$out #[1,] lower whisker, [3,] median, [5,] upper whisker
boxplot(logNFL ~ RTQUIC, data= df)$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
boxplot(logNFL ~ DX_APD, data= df)$out #[1,] lower whisker, [3,] median, [5,] upper whisker
vec.outliers <- boxplot(logNFL ~ DXRTQUIC, data= df)$out #list of outlier (Tukey method 1.5 IQR in each RTQUIC_DX combination)
dfnfl<- df[!(df$logNFL %in% vec.outliers), ] #remove any subjects whose logNFL is exactly equal to the value of these outliers
dfnfl <- dfnfl %>% drop_na("logNFL")


IQR1.5value <- as.numeric(quantile(RTposdf2$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTposdf2$logNFL)*1.5)) #reports the value of first point that is above Q3+ IQR*3 (3 is very tolerant threshold)
maxvalue <- max(RTposdf2$logNFL,na.rm=T)
threshold <- min(IQR1.5value, maxvalue)
df2[!((df2$logNFL<threshold & df2$RTQUIC=="aSyn-SAA positive")| (df2$RTQUIC=="aSyn-SAA negative")), c("ID", "NFL", "logNFL")] 
df2nfl <- subset(df2, (logNFL<threshold & RTQUIC=="aSyn-SAA positive") | (RTQUIC=="aSyn-SAA negative"))

boxplot <- boxplot(logNFL ~ DX_APD, data= RTnegdf2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	stripchart(logNFL ~ DX_APD, data = RTnegdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
	threshold <- min(max(RTnegdf2$logNFL,na.rm=T), as.numeric(quantile(RTnegdf2$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTnegdf2$logNFL)*1.5))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
df2nfl[!((df2nfl$logNFL<threshold & df2nfl$RTQUIC=="aSyn-SAA negative")| (df2nfl$RTQUIC=="aSyn-SAA positive")), c("ID", "NFL", "logNFL")] 
	df2nfl <- subset(df2nfl, (logNFL<threshold & RTQUIC=="aSyn-SAA negative") | (RTQUIC=="aSyn-SAA positive"))

#REMOVE OUTLIERS BASED ON DIAGNOSIS
threshold <- min(max(RTposdf2$logNFL,na.rm=T), as.numeric(quantile(RTposdf2$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTposdf2$logNFL)*1.5))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
df2[!((df2$logNFL<threshold & df2$RTQUIC=="aSyn-SAA positive")| (df2$RTQUIC=="aSyn-SAA negative")), c("ID", "NFL", "logNFL")] 
df2nfl <- subset(df2, (logNFL<threshold & RTQUIC=="aSyn-SAA positive") | (RTQUIC=="aSyn-SAA negative"))

nrow(df2nfl)
df2nfl[, c("ID", "NFL", "logNFL")]
sort(df2nfl$NFL)

boxplot(NFL ~ RTQUIC, data= df2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	stripchart(NFL ~ RTQUIC, data = df2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(NFL ~ RTQUIC, data= df2nfl, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	stripchart(NFL ~ RTQUIC, data = df2nfl, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

boxplot(logNFL ~ RTQUIC, data= df2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	stripchart(logNFL ~ RTQUIC, data = df2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
boxplot(logNFL ~ RTQUIC, data= df2nfl, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	stripchart(logNFL ~ RTQUIC, data = df2nfl, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

shapiro.test(RTposdfnfl$logNFL) #al
shapiro.test(RTnegdfnfl$logNFL)  #al
shapiro.test(CBSdfnfl$logNFL)  #al
shapiro.test(PSPdfnfl$logNFL)  #al
shapiro.test(ADposdfnfl$logNFL)  #al
shapiro.test(ADnegdfnfl$logNFL)  #al
leveneTest(logNFL ~ DX_APD*RTQUIC, dfnfl) #homo

# NFL STATS: LINEAR REGRESSION MODEL SELECTION

# Compare models with Ftest
		test1 <- lm(logNFL ~ RTQUIC, df2nfl)
		test2 <- lm(logNFL ~ DX_APD, df2nfl)
		test3 <- lm(logNFL ~ Lifetime_AD_binary, df2nfl)
		test4 <- lm(logNFL ~ RTQUIC + DX_APD, df2nfl) #adding value to the model
		test5 <- lm(logNFL ~ RTQUIC + Lifetime_AD_binary, df2nfl) #not adding any value to the model
		anova(test1, test4) 
		anova(test1, test5) 

		# Inclusion of other variables as covariate:
		summary(lm(df2nfl$logNFL ~ df2nfl$Age)) #no linear relationship
		summary(lm(df2nfl$logNFL ~ df2nfl$LP2_Disease_Duration))#no linear relationship 
		summary(lm(df2nfl$logNFL ~ df2nfl$abeta)) #linear relationship
		summary(lm(df2nfl$logNFL ~ df2nfl$ptau)) #no linear relationship
		summary(lm(df2nfl$logNFL ~ df2nfl$ttau)) #no linear relationship

# NFL STATS: LINEAR REGRESSION MODEL

stdmlr <- lm(logNFL ~ RTQUIC*DX_APD + scale(abeta), df2nfl) 
summary(stdmlr) 

# Diagnostics of the model run
	check_normality(stdmlr) #Ok normality of residuals
	autoplot(stdmlr, which = 1:6) #All plots including scale location plot

		# Visualize the values for each of the IDs that are indexed on above plots
		df2nfl[27, c("ID", "abeta", "logabeta", "RTQUIC", "NFL")] 
		df2nfl[34, c("ID", "abeta", "logabeta", "RTQUIC", "NFL")] 
		df2nfl[57, c("ID", "abeta", "logabeta", "RTQUIC", "NFL")] 
		
		for-loop to test the model without each of these values (ie once without outlier 1, then outlier 2, etc)
			vecIDs <- df2nfl[c(8,34,37,57), "ID"] #Create vector of each ID. For 3*IQR outlier threshold
		vecIDs <- df2nfl[c(27,34,53), "ID"] #Create vector of each ID. For 1.5*IQR outlier threshold
		for (i in vecIDs) {
			test <- subset(df2nfl, ID!=i) #for some analysis, need to exclude the potential false negative too
			teststdmlr <- lm(logNFL ~ RTQUIC*DX_APD + scale(abeta), test) 
			print(summary(teststdmlr))
		}

	bptest(stdmlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
	durbinWatsonTest(stdmlr) #Ok autocorrelation of residuals	

		# Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
		df2nfl[(cooks.distance(stdmlr))>0.25, ]$ID #no subject >0.25
		cooksD <- cooks.distance(stdmlr)
		n <- nrow(df2nfl) #Anothr threshold for Cook's data = 4/n so lower. In this case, too many high Cook. 
		plot(cooksD, main = "Cooks Distance for Influential Obs") #Plot Cook's Distance with a horizontal line at 4/n to see which observations exceed this thresdhold. There is a cluster of high Cook subjects. 
			abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line

	# Multicollinearity checks
	car::vif(stdmlr) #no multicollinearity at all

	# Main effects: 
	emmeans(stdmlr, ~ RTQUIC:DX_APD) #adjusted means.  
    pairs(emmeans(stdmlr, ~ RTQUIC:DX_APD)) #adjusted means.  

cat("\n\n######################################################################################################\n",
	   "###########################              8.1.2. FIG 1.D.       	#######################################\n",
	   "#######################################################################################################\n")

# FIG1D: NFL: not different between DX, RTQUIC, but linearly related to Abeta42. Option 1: Plot Abeta42 by NFL relationship and add colors for DX and shapes for RTQUIC. 
#Option 2: Plot boxplot between RTQUIC status (see below this figure)

	# Create dataframe usable for plotting (ie kick out the scale())
	plotmlr <- lm(logNFL ~ RTQUIC + DX_APD + abeta, df2nfl) 
	summary(plotmlr) 

	label <- "paste(''*italic(p), \" < .05\")" #Second annotation is p-value for the interaction

	# General layout of the plot: x and y + stats_smooth for linear relationship
	fig1d_ver1 <- ggplot(df2nfl, aes(x=abeta, y=logNFL)) + #No need for color or fill 

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


# FIG1Dversion2: NFL: not different between DX, RTQUIC, show boxplots only.
#Option 2: Plot boxplot between RTQUIC status. 

# General layout of the plot: boxplot by diagnosis where diagnosis is sig, within diagnosis rtquic groups are not sig

	fig1d_ver2 <- ggplot(df2, aes(x=DX_APD, y=logNFL, color=RTQUIC))+ #No need for color or fill 

				#Add actual datapoints
				geom_boxplot() +
				geom_jitter(aes(color=RTQUIC), position=position_jitterdodge()) +
 
				# Fix legends and set up the appearance of the points
				scale_color_manual(values=cbPalette_RTQUIC, name="ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
				scale_fill_manual(values=cbPalette_RTQUIC, name="ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters
				scale_shape_manual(values=c(16,17), name="Diagnosis", breaks=c("CBS", "PSP"), labels=c("CBS", "PSP")) + #Color for the datapoints

				# Fix labs
				labs(title="ASyn-SAA does not affect NfL levels",
					subtitle="",
					# x=expression(bold("AD status")),
					y=expression(bold("CSF NfL levels (pg/mL) (log)"))) + #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

				# Annotate:
				stat_regline_equation(label.x=120, label.y=8.7, color="red", size=6, show.legend=FALSE) + #shows the equation for each geom_line. Cannot use other options since it would not be based on the model
				
				# annotate("text", size=6, x=175, y=8.5, color="red", label=label, parse=TRUE)+ #label as defined above. Parse allows for use of mathematical notation

				#Aesthetic only
				theme_classic() +
				theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
				theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
				theme(legend.title = element_text(face="bold", size=16), legend.text= element_text(size=14))
	
	fig1d_ver2

	ggsave(fig1d_ver1, filename = "Fig1d_ver1.png", bg= "transparent", width=9, height=10)

boxplot1<- boxplot(logNFL ~ RTQUIC, data= CBSdf2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
		stripchart(logNFL ~ RTQUIC, data = CBSdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		
boxplot2<-boxplot(logNFL ~ RTQUIC, data= PSPdf2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
		stripchart(logNFL ~ RTQUIC, data = PSPdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

boxplot1$out
boxplot2$out



cat("\n\n#######################################################################################################\n",
	"                                    8.2. ASYN-SAA+ & SYMPTOMS \n",
	   "#######################################################################################################\n")

cat("\n\n######################################################################################################\n",
	   "###########################              8.2.1. CBS-ONLY       #######################################\n",
	   "#######################################################################################################\n")


# DO CBS FIRST: COUNTS + ANALYSES
#################################
RTposCBSdf2 <- subset(CBSdf2, RTQUIC=="aSyn-SAA positive")
RTnegCBSdf2 <- subset(CBSdf2, RTQUIC=="aSyn-SAA negative")

CBSdf2 %>% group_by(RTQUIC) %>% count(Sex)
CBSdf2 %>% group_by(RTQUIC) %>% count(AD_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(APOEe4)

CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# STATISTICAL COMPARISONS FOR AGES AND SEX
table(CBSdf2$Sex, CBSdf2$RTQUIC)
chisq.test(table(CBSdf2$Sex, CBSdf2$RTQUIC), correct=F)

table(CBSdf2$AD_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$AD_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell

table(CBSdf2$APOEe4, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$APOEe4, CBSdf2$RTQUIC)) # Expected count is <5 for one cell


shapiro.test(RTposCBSdf2$Age) #al
shapiro.test(RTnegCBSdf2$Age) #al
var.test(Age ~ RTQUIC, data = CBSdf2) #homo
t.test(CBSdf2$Age ~ CBSdf2$RTQUIC, var.equal=TRUE) 

shapiro.test(RTposCBSdf2$Onset) #al
shapiro.test(RTnegCBSdf2$Onset) #al
var.test(Onset ~ RTQUIC, data = CBSdf2) #homo
t.test(CBSdf2$Onset ~ CBSdf2$RTQUIC, var.equal=TRUE) 

shapiro.test(RTposCBSdf2$Park_onset) #al
shapiro.test(RTnegCBSdf2$Park_onset) #al
var.test(Park_onset ~ RTQUIC, data = CBSdf2) #homo
	# wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
t.test(CBSdf2$Park_onset ~ CBSdf2$RTQUIC, var.equal=TRUE) 

CBSdf2 %>% group_by(RTQUIC) %>% count(Tremor_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(RestTremor)
CBSdf2 %>% group_by(RTQUIC) %>% count(LimbRigidity)
CBSdf2 %>% group_by(RTQUIC) %>% count(Slowness_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Falls_PI)
CBSdf2 %>% group_by(RTQUIC) %>% count(Gait)
CBSdf2 %>% group_by(RTQUIC) %>% count(RBD_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
## CBSdf2 %>% group_by(RTQUIC) %>% count(Anosmia_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Constipation_binary)
## CBSdf2 %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
CBSdf2 %>% group_by(RTQUIC) %>% count(Sexual_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Urinary_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Bowel_binary)
CBSdf2 %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

table(CBSdf2$anyPPA, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$anyPPA, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Tremor_binary, CBSdf2$RTQUIC)
chisq.test(table(CBSdf2$Tremor_binary, CBSdf2$RTQUIC), correct=F)
table(CBSdf2$RestTremor, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$RestTremor, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$LimbRigidity, CBSdf2$RTQUIC)
chisq.test(table(CBSdf2$LimbRigidity, CBSdf2$RTQUIC), correct=F)
table(CBSdf2$Slowness_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Slowness_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Falls_PI, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Falls_PI, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Gait, CBSdf2$RTQUIC)
chisq.test(table(CBSdf2$Gait, CBSdf2$RTQUIC), correct=F)
table(CBSdf2$RBD_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$RBD_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Lifetime_Dopa_responder_true, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Lifetime_Dopa_responder_true, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Lifetime_VisualHallucinations_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Lifetime_VisualHallucinations_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Constipation_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Constipation_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Urinary_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Urinary_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
table(CBSdf2$Bowel_binary, CBSdf2$RTQUIC)
fisher.test(table(CBSdf2$Bowel_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell


cat("\n\n######################################################################################################\n",
	   "###########################              8.2.2. PSP-ONLY       #######################################\n",
	   "#######################################################################################################\n")

# DO PSP SECOND: COUNTS + ANALYSES
################################
RTposPSPdf2 <- subset(PSPdf2, RTQUIC=="aSyn-SAA positive")
RTnegPSPdf2 <- subset(PSPdf2, RTQUIC=="aSyn-SAA negative")

PSPdf2 %>% group_by(RTQUIC) %>% count(Sex)
PSPdf2 %>% group_by(RTQUIC) %>% count(AD_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(APOEe4)
PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))




# STATISTICAL COMPARISONS FOR AGES AND SEX
table(PSPdf2$Sex, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Sex, PSPdf2$RTQUIC)) # Expected count is <5 for one cell

table(PSPdf2$AD_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$AD_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell

table(PSPdf2$APOEe4, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$APOEe4, PSPdf2$RTQUIC)) # Expected count is <5 for one cell

shapiro.test(RTposPSPdf2$Age) #al
shapiro.test(RTnegPSPdf2$Age) #al
var.test(Age ~ RTQUIC, data = PSPdf2) #homo
t.test(PSPdf2$Age ~ PSPdf2$RTQUIC, var.equal=TRUE) 

shapiro.test(RTposPSPdf2$Onset) #al
shapiro.test(RTnegPSPdf2$Onset) #al
var.test(Onset ~ RTQUIC, data = PSPdf2) #homo
t.test(PSPdf2$Onset ~ PSPdf2$RTQUIC, var.equal=TRUE) 

shapiro.test(RTposPSPdf2$Park_onset) #al
shapiro.test(RTnegPSPdf2$Park_onset) #al
var.test(Park_onset ~ RTQUIC, data = PSPdf2) #homo
# # wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
t.test(PSPdf2$Park_onset ~ PSPdf2$RTQUIC, var.equal=TRUE) 


PSPdf2 %>% group_by(RTQUIC) %>% count(anyPPA)
PSPdf2 %>% group_by(RTQUIC) %>% count(Tremor_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(RestTremor)
PSPdf2 %>% group_by(RTQUIC) %>% count(LimbRigidity)
PSPdf2 %>% group_by(RTQUIC) %>% count(Slowness_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Falls_PI)
PSPdf2 %>% group_by(RTQUIC) %>% count(Gait)
PSPdf2 %>% group_by(RTQUIC) %>% count(RBD_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
PSPdf2 %>% group_by(RTQUIC) %>% count(Anosmia_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Constipation_binary)
## PSPdf2 %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
PSPdf2 %>% group_by(RTQUIC) %>% count(Sexual_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Urinary_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Bowel_binary)
PSPdf2 %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

table(PSPdf2$Tremor_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Tremor_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$LimbRigidity, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$LimbRigidity, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Slowness_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Slowness_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Falls_PI, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Falls_PI, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Gait, PSPdf2$RTQUIC)
chisq.test(table(PSPdf2$Gait, PSPdf2$RTQUIC), correct=F)
table(PSPdf2$RBD_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$RBD_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Lifetime_Dopa_responder_true, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Lifetime_Dopa_responder_true, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Lifetime_VisualHallucinations_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Lifetime_VisualHallucinations_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Constipation_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Constipation_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Orthostatism_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Orthostatism_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Urinary_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Urinary_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
table(PSPdf2$Bowel_binary, PSPdf2$RTQUIC)
fisher.test(table(PSPdf2$Bowel_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell


# COMPARISONS IGNORING DX
t.test(df2$Age ~ df2$RTQUIC, var.equal=TRUE) 
t.test(df2$Onset ~ df2$RTQUIC, var.equal=TRUE) 
t.test(df2$Park_onset ~ df2$RTQUIC, var.equal=TRUE) 
table(df2$anyPPA, df2$RTQUIC) #Fisher
table(df2$Lifetime_Dopa_responder_true, df2$RTQUIC) #Fisher
table(df2$AD_binary, df2$RTQUIC) #Fisher
table(df2$APOEe4, df2$RTQUIC) #Fisher
table(df2$Tremor_binary, df2$RTQUIC) #Chisquare
table(df2$RestTremor, df2$RTQUIC) #Fisher
table(df2$LimbRigidity, df2$RTQUIC) #Chisquare
table(df2$Slowness_binary, df2$RTQUIC) #Fisher
table(df2$Falls_PI, df2$RTQUIC) #Fisher
table(df2$Gait, df2$RTQUIC) #Chisquare
table(df2$RBD_binary, df2$RTQUIC) #Fisher
table(df2$Slowness_binary, df2$RTQUIC) #Fisher
table(df2$Lifetime_VisualHallucinations_binary, df2$RTQUIC) #Fisher
table(df2$Constipation, df2$RTQUIC) #Fisher
table(df2$Sexual, df2$RTQUIC) #Fisher
table(df2$Orthostatism_binary, df2$RTQUIC) #Fisher
table(df2$Urinary, df2$RTQUIC) #Chisquare
table(df2$Bowel, df2$RTQUIC) #Fisher
table(df2$Thermoregulatory_binary, df2$RTQUIC) #Fisher

fisher.test(table(df2$anyPPA, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Lifetime_Dopa_responder_true, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$AD_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$APOEe4, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$RestTremor, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Slowness_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$RBD_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Constipation, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Sexual, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Orthostatism_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Bowel, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Thermoregulatory_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Lifetime_VisualHallucinations_binary, df2$RTQUIC)) # Expected count is <5 for one cell
fisher.test(table(df2$Falls_PI, df2$RTQUIC)) # Expected count is <5 for one cell
chisq.test(table(df2$Tremor_binary, df2$RTQUIC), correct=F)
chisq.test(table(df2$LimbRigidity, df2$RTQUIC), correct=F)
chisq.test(table(df2$Gait, df2$RTQUIC), correct=F)
chisq.test(table(df2$Urinary, df2$RTQUIC), correct=F)

cat("\n\n######################################################################################################\n",
	   "###########################              8.2.3. RADAR PLOTS      ######################################\n",
	   "#######################################################################################################\n")

# START WITH DIAGNOSIS. CBS first because I'd like to overlay two alpha syn on the same diagnosis plot. 

List = list() #Create an empty list to which you will assign the values from the for-loop. By creating this list outside of the main for-loop, you allow for the data to be entered under different entries which means you can call all the values you need. 

rt.value <- c("aSyn-SAA positive", "aSyn-SAA negative")
for (rt in rt.value) { #for-loop that tests each value of the variable

	CBSdf2.rt <- CBSdf2[CBSdf2$RTQUIC== rt, ] #Within subset of CBS for eg, look for subset of RT+

print(rt)

#SUMMARIZE THE CATEGORICAL VARIABLES INTO ONE SINGLE COUNT OF "YES"
Tremor_perc <- CBSdf2.rt %>% summarise(Tremor_count= sum(Tremor_binary == "Yes")) %>% mutate(Tremor_count= (as.numeric(Tremor_count)/nrow(CBSdf2.rt)*100))
RestTremor_perc <- CBSdf2.rt %>% summarise(RestTremor_count= sum(RestTremor == "Yes")) %>% mutate(RestTremor_count= (as.numeric(RestTremor_count)/nrow(CBSdf2.rt)*100))
LimbRigidity_perc <- CBSdf2.rt %>%  summarise(LimbRigidity_count= sum(LimbRigidity == "Yes")) %>% mutate(LimbRigidity_count= (as.numeric(LimbRigidity_count)/nrow(CBSdf2.rt)*100))
Slowness_perc <- CBSdf2.rt %>% summarise(Slowness_count= sum(Slowness_binary == "Yes")) %>% mutate(Slowness_count= (as.numeric(Slowness_count)/nrow(CBSdf2.rt)*100))
Apraxia_perc <- CBSdf2.rt %>% summarise(Apraxia_count= sum(Apraxia == "Yes")) %>% mutate(Apraxia_count= (as.numeric(Apraxia_count)/nrow(CBSdf2.rt)*100))
Gait_perc <- CBSdf2.rt %>%  summarise(Gait_count= sum(Gait == "Yes")) %>% mutate(Gait_count= (as.numeric(Gait_count)/nrow(CBSdf2.rt)*100))
FallsPI_perc <- CBSdf2.rt %>%  summarise(FallsPI_count= sum(Falls_PI == "Yes")) %>% mutate(FallsPI_count= (as.numeric(FallsPI_count)/nrow(CBSdf2.rt)*100))

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

radar.CBSdf2 <- data.frame(row.names = c("aSyn-SAA negative", "aSyn-SAA positive"),
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
radar.CBSdf2 <- rbind(max_mindf, radar.CBSdf2)
radar.CBSdf2

#Rename some variables for presentation purposes
colnames(radar.CBSdf2)[which(names(radar.CBSdf2) == "Rest")] <- "Rest tremor"
colnames(radar.CBSdf2)[which(names(radar.CBSdf2) == "Limb")] <- "Limb rigidity"
colnames(radar.CBSdf2)[which(names(radar.CBSdf2) == "Falls")] <- "Falls & instability"
colnames(radar.CBSdf2)[which(names(radar.CBSdf2) == "Gait")] <- "Gait \n p<0.1"


# Create the radar charts
png(filename ="Fig1c_ver1.png")

create_beautiful_radarchart(data= radar.CBSdf2, color= cbPalette_RTQUIC, vlcex=1, plty=1, title="Motor symptoms in CBS")
		# Add an horizontal legend
		# legend(-0.65, -1.2, legend=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-")), horiz=TRUE, bty= "o", pch= 15 , col= cbPalette_RTQUIC, text.col= "black", cex= 1, pt.cex= 1.5)

# dev.off()#Not needed?


#NOW PSP. 

#START WITH DIAGNOSIS. CBS first because I'd like to overlay two alpha syn on the same diagnosis plot. 
List = list() #Create an empty list to which you will assign the values from the for-loop. By creating this list outside of the main for-loop, you allow for the data to be entered under different entries which means you can call all the values you need. 

rt.value <- c("aSyn-SAA positive", "aSyn-SAA negative")
for (rt in rt.value) { #for-loop that tests each value of the variable

	PSPdf2.rt <- PSPdf2[PSPdf2$RTQUIC== rt, ] #Within subset of CBS for eg, look for subset of RT+

print(rt)

#SUMMARIZE THE CATEGORICAL VARIABLES INTO ONE SINGLE COUNT OF "YES"
Tremor_perc <- PSPdf2.rt %>% summarise(Tremor_count= sum(Tremor_binary == "Yes")) %>% mutate(Tremor_count= (as.numeric(Tremor_count)/nrow(PSPdf2.rt)*100))
RestTremor_perc <- PSPdf2.rt %>% summarise(RestTremor_count= sum(RestTremor == "Yes")) %>% mutate(RestTremor_count= (as.numeric(RestTremor_count)/nrow(PSPdf2.rt)*100))
LimbRigidity_perc <- PSPdf2.rt %>%  summarise(LimbRigidity_count= sum(LimbRigidity == "Yes")) %>% mutate(LimbRigidity_count= (as.numeric(LimbRigidity_count)/nrow(PSPdf2.rt)*100))
AxialRigidity_perc <- PSPdf2.rt %>%  summarise(AxialRigidity_count= sum(AxialRigidity == "Yes")) %>% mutate(AxialRigidity_count= (as.numeric(AxialRigidity_count)/nrow(PSPdf2.rt)*100))
Slowness_perc <- PSPdf2.rt %>% summarise(Slowness_count= sum(Slowness_binary == "Yes")) %>% mutate(Slowness_count= (as.numeric(Slowness_count)/nrow(PSPdf2.rt)*100))
OM_perc <- PSPdf2.rt %>% summarise(OM_count= sum(VerticalOM == "Yes")) %>% mutate(OM_count= (as.numeric(OM_count)/nrow(PSPdf2.rt)*100))
Gait_perc <- PSPdf2.rt %>%  summarise(Gait_count= sum(Gait == "Yes")) %>% mutate(Gait_count= (as.numeric(Gait_count)/nrow(PSPdf2.rt)*100))
FallsPI_perc <- PSPdf2.rt %>%  summarise(FallsPI_count= sum(Falls_PI == "Yes")) %>% mutate(FallsPI_count= (as.numeric(FallsPI_count)/nrow(PSPdf2.rt)*100))

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

radar.PSPdf2 <- data.frame(row.names = c("aSyn-SAA negative", "aSyn-SAA positive"),
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
radar.PSPdf2 <- rbind(max_mindf, radar.PSPdf2)
radar.PSPdf2

#Rename some variables for presentation purposes
colnames(radar.PSPdf2)[which(names(radar.PSPdf2) == "Rest")] <- "Rest tremor"
colnames(radar.PSPdf2)[which(names(radar.PSPdf2) == "Limb")] <- "Limb rigidity \n p<0.1"
colnames(radar.PSPdf2)[which(names(radar.PSPdf2) == "Axial")] <- "Axial rigidity"
colnames(radar.PSPdf2)[which(names(radar.PSPdf2) == "OM")] <- "Oculomotor"
colnames(radar.PSPdf2)[which(names(radar.PSPdf2) == "Falls")] <- "Falls & instability"


# Create the radar charts
png(filename ="Fig1d_ver1.png")

create_beautiful_radarchart(data= radar.PSPdf2, color= cbPalette_RTQUIC, vlcex=1, plty=1, title="Motor symptoms in PSP")
		# Add an horizontal legend
		# legend(-0.65, -1.2, legend=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-")), horiz=TRUE, bty= "o", pch= 15 , col= cbPalette_RTQUIC, text.col= "black", cex= 1, pt.cex= 1.5)



cat("\n\n\n\n###################################################################################################\n",
		   "9. BINARY LOGISTIC REGRESSION\n",
	 	   "####################################################################################################\n\n")

# DEFENSE
if (sum(df$RTQUIC_BLR ==1) != 22) {
	cat("There is an issue with the binarization of RTQUIC variable for binary logistic regression \n")
}


# Model is selected based on previous analyses (logabeta*Onset) + the simple comparisons in Supp material (Gait/RBD_binary).
# Model was run with and without AD/DX_APD status and all values remained very similar. 
blr <- glm(RTQUIC_BLR ~ DX_APD + scale(Onset_age)*scale(logabeta) + RBD_binary + Gait_2, data= df2, family = "binomial")
summary(blr)
				
## Model fit
pscl::pR2(blr)["McFadden"]


### Multicollinearity checks
car::vif(blr)
	blronset <- glm(RTQUIC_BLR ~ DX_APD + Onset_age + RBD_binary + Gait_2, data= df, family = "binomial") #Check individual relationship of Onset with RTQUIC
	car::vif(blronset)
blrabeta <- glm(RTQUIC_BLR ~ DX_APD + logabeta + RBD_binary + Gait, data= df, family = "binomial") #Check individual relationship of Abeta with RTQUIC
	car::vif(blrabeta)

df <- df %>% mutate(centredbeta = logabeta- mean(logabeta)) %>%
			mutate(centredonset = Onset_age- mean(Onset_age)) %>% data.frame()

centredblr <- glm(RTQUIC_BLR ~ DX_APD + centredbeta*centredonset + RBD_binary + Gait_2, data= df, family = "binomial")
summary(centredblr)
car::vif(centredblr) #Now all VIF are low, indicating no problem of multicolllinearity. 


df[(cooks.distance(blr))>0.06, ]$ID #identifies subjects with high Cook's distance data. 
test <- subset(df, ID!="T258") 
blrtest <- glm(RTQUIC_BLR ~  DX_APD + Onset_age*logabeta + RBD_binary + Gait, data= test, family = "binomial")
summary(blrtest)

#Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
df[(cooks.distance(blr))>0.25, ]$ID #two subjects >0.25
cooksD <- cooks.distance(blr)
n <- nrow(df) #Another threshold for Cook,s data = 4/n so lower. In this case, too many high Cook. 
plot(cooksD, main = "Cooks Distance for Influential Obs") #Two subjects are clearly different than others. However, 1 is one of few RBD+ subjects & 1 is one of few AD+/RT+: hard to exclude them as analysis pointless then.
	abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line
	df %>% group_by(RTQUIC_BLR) %>% summarize(mean=format(round(mean(Onset_age, na.rm=T),3),3), sd=format(round(sd(Onset_age, na.rm=T),3),3))
	df %>% group_by(RTQUIC_BLR) %>% summarize(mean=format(round(mean(logabeta, na.rm=T),3),3), sd=format(round(sd(logabeta, na.rm=T),3),3))

#https://www.statology.org/dfbetas-in-r/
dfbetas <- as.data.frame(dfbetas(blr))
thresh <- 2/sqrt(nrow(df))
	
par(mfrow=c(2,1))
plot(dfbetas$Gait, type='h')
abline(h = thresh, lty = 2)
abline(h = -thresh, lty = 2)
plot(dfbetas$RBD_binary, type='h')
abline(h = thresh, lty = 2)
abline(h = -thresh, lty = 2)


#Sample size: https://www.statology.org/assumptions-of-logistic-regression/
	
# Check linear relationship between logit and explanatory variables: box-tidwell test
logodds <- blr$linear.predictors #if (any(x1 <= 0)) stop("the variables to be transformed must have only positive values")
boxTidwell(logodds ~ df$DX_APD)
boxTidwell(logodds ~ df$RBD_binary)
boxTidwell(logodds ~ df$Gait_2)

# Coefficients
caret::varImp(blr)	

# Odds ratio
coef(blr)
exp(coef(blr)) 
cbind(coef(blr),odds_ratio=exp(coef(blr)),exp(confint(blr, level=0.95))) #it says 2.5% and 97.5% because these are the two borders to end up with the cnetral 95% of your distribution

#Example of interpretation:
# exp(4.88278) #RBD. 81.60169. Odds of someone with RBD being RTquic+ increases by 8000%. IE x80. 
# exp(-2.85663) #Gait. 0.09253855. Odds of someone with gait issues being RTQUIC+ decreases by 10% compared to someone without. 
##### exp(0.15122) #Age: 1.16. #The value indicates that as age increase by one more unit, then the odds of being SAA+ increases by 16%



##############################################	 LAG HOURS	  		###########################################################
##############################################################################################################################

cat("Lag hours is the #hours required to reach threshold for positivity. It makes the most sense to think about it as data suited for survival analysis. \n")
cat("For that purpose, we are censoring the subjects who never reached positivity (RTQUIC negative) \n")
df %>% count(RTQUIC) 



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
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logptau, method="spearman")
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$logttau, method="spearman")  
cor.test(RTposdf$RTQUIC_survival_hours, RTposdf$ATI_2, method="spearman") 

# LAG STATISTICS: CORRELATIONS
ggscatter(RTposdf, x="Onset_age", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with age at onset",
		xlab="Age at onset (years)", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="Onset_age", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))


ggscatter(RTposdf, x="logabeta", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with logged amyloid beta 42",
		xlab="Age at onset (years)", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="logabeta", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

ggscatter(RTposdf, x="logptau", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with logged ptau181",
		xlab="Age at onset (years)", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="logptau", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))
  

ggscatter(RTposdf, x="logttau", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with logged t-tau",
		xlab="Age at onset (years)", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="logttau", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

ggscatter(RTposdf, x="logNFL", y= "RTQUIC_survival_hours", fill="DX_APD",https://www.facebook.com/photo?fbid=10159531068156876&set=pcb.3597966813754060
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with logged NFL",
		xlab="NFL (pg/mL)", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="logNFL", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))

ggscatter(RTposdf, x="ATI_2", y= "RTQUIC_survival_hours", fill="DX_APD",
		size=5, shape=21, palette = c(CBS= "#56B4E9", PSP = "#CC79A7"),
		add = "reg.line", cor.coef = TRUE, cor.method = "spearman", 
		title="Correlation of lag with ATI",
		xlab="ATI_2", ylab = "Lag (hours)") 
ggscatter(RTposdf, x="ATI_2", y= "RTQUIC_survival_hours", fill="DX_APD", add = "reg.line", conf.int = TRUE) +
	stat_cor(aes(color = DX_APD, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_fill_manual(values=cbPalette_DX_APD, name="Diagnosis", labels=c("CBS", "PSP"))


# # LAG STATISTICS: LOG-RANK TEST
# s1 <- survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df) #uses the survival() package
# s1 <- survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df) #uses the survival() package
# s2 <- survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df) #uses the survival() package
# s3 <- survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df) #uses the survival() package
# ggsurvplot(s1, data = df, pval = TRUE)
# ggsurvplot(s2, data = df, pval = TRUE)
# ggsurvplot(s3, data = df, pval = TRUE)


# survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df) #p=.05
# summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df), times = 24)


# # LAG STATISTICS: KAPLAN-MEIER CURVE FOR WHOLE DATASET
# # Uses survfit2() in ggsurvfit() package:
# survfit2(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df)  %>% 
#     ggsurvfit() + add_confidence_interval() + labs(x="Lag (hours)", "Survival probability") #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

# # LAG STATISTICS: ANALYSES ON ONLY THE RT-QUIC SUBJECTS
# # Fit a Cox proportional hazards model
# # Redundant with the fisher analyses already shown above
# fit.coxph <- coxph(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset + AD, data = df)
# ggforest(fit.coxph, data = df)




############################################		THT FLUO	  		#######################################################
##############################################################################################################################

# cat("THT max is the max fluorescent signal reached after 48 hours of monitoring of the assay. \n")
# STATISTICS: SUMMARY

# THT STATISTICS: DISTRIBUTION
hist(RTposdf$ThTmax)  
shapiro.test(df[RTposdf$DX_APD =="CBS", ]$ThTmax) #nonnormal
shapiro.test(df[RTposdf$DX_APD =="PSP", ]$ThTmax) #nonnormal
shapiro.test(df[RTposdf$Early_onset =="Young-onset", ]$ThTmax) #nonnormal
shapiro.test(df[RTposdf$Early_onset =="Late-onset", ]$ThTmax) #nonnormal
shapiro.test(df[RTposdf$AD =="AD Positive", ]$ThTmax) #nonnormal
shapiro.test(df[RTposdf$AD =="AD Negative", ]$ThTmax) #nonnormal
leveneTest(ThTmax ~ DX_APD, data = RTposdf) #homoscedasticity  
leveneTest(ThTmax ~ Early_onset, data = RTposdf) #homoscedasticity  
leveneTest(ThTmax ~ AD, data = RTposdf) #homoscedasticity  


# THT STATISTICS: SUMMARY
# ONSET STATISTICS: SUMMARY
RTposdf %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% group_by(AD) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% group_by(Early_onset) %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))
RTposdf %>% summarize(count=n(), format(round(mean(ThTmax, na.rm=T),2),2), sd=sd(ThTmax, na.rm=T))

# ONSET STATISTICS: TTEST
t.test(RTposdf$ThTmax ~ RTposdf$DX_APD, var.equal=TRUE) 


# THT STATISTICS: ANCOVA
aov <- aov(ThTmax ~  DX_APD + Onset_age + AD, RTposdf) 
Anova(aov, type="II") #Compare with type III
check_normality(aov)



