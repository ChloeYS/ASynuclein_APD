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
ADposdf <- subset(df, AD=="AD Positive") #for the description of AD+/- cohort
ADnegdf <- subset(df, AD=="AD Negative") #for the description of AD+/- cohort
APOEposdf <- subset(df, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf <- subset(df, APOEe4=="Negative") #for the description of APOE+/- cohort
YODdf <- subset(df, Early_onset=="Early-onset") #for the description of APOE+/- cohort
LODdf <- subset(df, Early_onset=="Late-onset") #for the description of APOE+/- cohort


# DEFENSE
c(nrow(CBSdf), nrow(PSPdf), nrow(ADposdf), nrow(ADnegdf), nrow(ADnegdf), nrow(APOEposdf), nrow(APOEnegdf), nrow(YODdf), nrow(LODdf)) == c(39, 28, 15, 52, 52, 14, 47, 35, 32)

###############################################################################################################################
#CREATE SUBSETS FOR ASSUMPTION TESTING
###############################################################################################################################

##CREATE SOME OF THE SUBSETS USED LATER (MOSTLY FOR OUTLIER IDENTIFICATION)
#Not very efficient - update later if occasion to
