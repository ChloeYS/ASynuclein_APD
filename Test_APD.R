# FILENAME: Analysis_APD.R

#Last updated 16 May 2024

#Usage: sources other files: 
#RScript Test_APD.R /Users/nikhilbhagwat/Desktop/7_PUBLICATIONS_ONGOING/7.1_APD_MS_2024-01/7.1.0_Data/APD_Neurology_MS_data.csv



###############################################################################################################################
#SOURCE PACKAGES AND FUNCTIONS
###############################################################################################################################

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


# Consider removing if not used: 
# library(multcomp) #multiple comparisons
# library(lmtest) #Breusch-Pagan
# library(ggpmisc) #annotate figures
# library(modelbased) #extract some parameters and estimates from a model
# library(psych)
# library(performance) #check_normality
# library(rcompanion) #Cramer's V
# library(ggstatsplot) #plot figures
# library(fmsb) #radar plot



###############################################################################################################################
#DATAFRAME LOADING
###############################################################################################################################

#Loading the QCed dataframe created from the raw data by the function var.func.1

df <- read.func(arg.vec[1], 'df') #arg.vec[1] is the csv that I entered data into. 
df <- var.func.1(df) #the df for analysis is created

write.csv(df, "dataframe.csv") #a copy of object created at the time of submission was kept for reference. Can be shared upon request. 



###############################################################################################################################
#CREATE SUBSETS FOR ASSUMPTION TESTING 
###############################################################################################################################

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


###############################################################################################################################
												# COHORT CHARACTERISTICS
###############################################################################################################################

cat("1. COMPARISONS OF DEMOGRAPHICS FOR DX: \n")
df %>% count(DX_APD)

##############################################################################################################################


#############################################			SEX	  		###########################################################
###############################################################################################################################

# STATISTICS: SUMMARY
cat("1.1. TOTAL NUMBER + SEX: \n")
cat("Prior to subject exclusion, the total number of subjects in the dataset is: \n")
df %>% count(DX_APD)
cat("Sex distribution in the dataset is: \n")
df %>% group_by(DX_APD) %>% count(Sex) 


# STATISTICS: PERCENTAGES
totalmatrix <- df %>% count(DX_APD)
sexmatrix <- df %>% group_by(DX_APD) %>% count(Sex) 

cat("Proportion of females in CBS is: \n")
(as.numeric(sexmatrix$n[1]) + as.numeric(sexmatrix$n[3]))/(as.numeric(totalmatrix$n[1]) + as.numeric(totalmatrix$n[2]))
cat("Proportion of females in CBS is:")
as.numeric(sexmatrix$n[1])/as.numeric(totalmatrix$n[1])
cat("Proportion of females in PSP is:")
as.numeric(sexmatrix$n[3])/as.numeric(totalmatrix$n[2])

# STATISTICS: CHI-SQUARE
table(df$Sex, df$DX_APD)
chisq.test(table(df$Sex, df$DX_APD), correct=F)



#############################################			AGE	  		###########################################################
###############################################################################################################################


# STATISTICS: SUMMARY
cat("1.2. AGE AT LP: \n")
df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
df %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))


# STATISTICS: TTEST
shapiro.test(CBSdf$Age) #normal
shapiro.test(PSPdf$Age) #normal
var.test(Age ~ DX_APD, data = df) #homoscedasticity
t.test <-t.test(df$Age ~ df$DX_APD, var.equal=TRUE)
	if (t.test[3] <= 0.05) {
		cat("There is a significant difference in age at LP between CBS and PSP. p-value:", t.test[3][[1]], "\n")
		cat(t.test[3][[1]])
	} else cat("There is no significant difference in age at LP between CBS and PSP. \n")


#############################################		EDUCATION	  		#######################################################
###############################################################################################################################

## ADD LATER ALL THE EDUCATION, BIOMARKERS, ETC. 




###############################################	 LAG HOURS	  		###########################################################
###############################################################################################################################

cat("Lag hours is the #hours required to reach threshold for positivity. It makes the most sense to think about it as data suited for survival analysis. \n")
cat("For that purpose, we are censoring the subjects who never reached positivity (RTQUIC negative) \n")
df %>% count(RTQUIC) 


# STATISTICS: SUMMARY

## 1. We want to have a column: Negative vs Positive
## 2. We want to have a column: lag hours max or 40

Surv(df$RTQUIC_survival_hours, df$RTQUIC_survived)
s1 <- survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df) #uses the survival() package



# STATISTICS: KAPLAN-MEIER CURVE FOR WHOLE DATASET
# Uses survfit2() in ggsurvfit() package:
survfit2(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df)  %>% 
    ggsurvfit() + add_confidence_interval() +
  add_risktable() #Summary of events vs nonevents

summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df), times = 24)
summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df), times = 13)
summary(survfit(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ 1, data = df), times = 48)



# STATISTICS: LOG-RANK TESTS FOR GROUP COMPARISONS
# Uses survival() package
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df)
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df)
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df) #p=.05
cat("Differences in overall survival (ie never crossing threshold) between the subjects who have an early vs late-onset are seen. \n")

# SANITY CHECK:
survdiff(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ RTQUIC, data = df)


# # STATISTICS: COX REGRESSION
# ## CHECK ASSUMPTIONS. NOT SURE IF APPROPRIATE. 
# # Uses survival() package
# coxph(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ DX_APD, data = df)
# coxph(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ AD, data = df)
# coxph(Surv(RTQUIC_survival_hours, RTQUIC_survived) ~ Early_onset, data = df)


#############################################		THT FLUO	  		#######################################################
###############################################################################################################################

cat("THT max is the max fluorescent signal reached after 48 hours of monitoring of the assay. \n")
# STATISTICS: SUMMARY




