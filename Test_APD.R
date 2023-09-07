#Test_APDs.R
#In the order of the manuscript (approx) to be able to find stats more easily as we read it. 

##LOAD THE INPUT DATA + RSCRIPTS
arg.vec <- commandArgs(trailingOnly = T) #Defines a vector of arguments. In this case, it is only one argument.  

source('Functions.R') #Source the functions. 
source('DFCreationFunction1.R') #Source the function taht creates the new dataframe on which analysis is performed. 

##CREATE DATAFRAME
df <- read.func(arg.vec[1], 'df') #arg.vec[1] is the csv that I entered data into. 
df <- var.func.1(df) #the df for analysis is created

write.csv(df, "dataframe.csv") #a copy is kept for reference. 

##LOAD LIBRARIES
library(car) #very useful for Anova(), Levene, etc. Please note that base R anova() and car:Anova() differ for:
	##contrast option already set as c("contr.sum", "contr.poly") for car::Anova(). Also, car::Anova() does type II SS
	##method as default, whereas base R anova() does type I SS method (same results as summary(aov))
library(emmeans) #estimated marginal means of a model
library(modelbased) #extract some parameters and estimates from a model
library(ggfortify) #plots for QC model
library(performance) #check_normality
library(lmtest) #Breusch-Pagan
library(rcompanion) #Cramer's V
library(caret)
library(psych)
library(ggplot2) #plot figures with stats
library(ggpubr) #ggscatter
library(ggstatsplot) #plot figures
library(ggpmisc) #annotate figures
# library(multcomp) #multiple comparisons



###############################################################################################################################
													#SETUP OF THE DATA
###############################################################################################################################


#EXCLUDE SUBJECTS
df1 <- subset(df, ID!="T102, T192") #for some analyses, I want to have all the CBS-AD including 2308, but not the subjects with less certain DX (T192)
df2 <- subset(df, ID!="2308_v1" & ID!="T102, T192") #for some analysis, need to exclude the potential false negative too


##CREATE SOME OF THE SUBSETS USED LATER (MOSTLY FOR OUTLIER IDENTIFICATION)
##CREATE A FUNCTION FOR THIS LATER
#df
CBSdf <- subset(df, DX_APD=="CBS") #for the description of AD+/- cohort
PSPdf <- subset(df, DX_APD=="PSP") #for the description of AD+/- cohort
ADposdf <- subset(df, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdf <- subset(df, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdf <- subset(df, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf <- subset(df, APOEe4=="Negative") #for the description of APOE+/- cohort

#df1
RTposdf1 <- subset(df1, RTQUIC=="aSyn-SAA positive") #for the main results on RT+/-
RTnegdf1 <- subset(df1, RTQUIC=="aSyn-SAA negative") #for the main results on RT+/-
CBSdf1 <- subset(df1, DX_APD=="CBS") #for the main results on RT+/-
PSPdf1 <- subset(df1, DX_APD=="PSP") #for the main results on RT+/-
ADposdf1 <- subset(df1, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdf1 <- subset(df1, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdf1 <- subset(df1, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf1 <- subset(df1, APOEe4=="Negative") #for the description of APOE+/- cohort

#df2
RTposdf2 <- subset(df2, RTQUIC=="aSyn-SAA positive") #for the main results on RT+/-
RTnegdf2 <- subset(df2, RTQUIC=="aSyn-SAA negative") #for the main results on RT+/-
CBSdf2 <- subset(df2, DX_APD=="CBS") #for the main results on RT+/-
PSPdf2 <- subset(df2, DX_APD=="PSP") #for the main results on RT+/-
ADposdf2 <- subset(df2, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdf2 <- subset(df2, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdf2 <- subset(df2, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf2 <- subset(df2, APOEe4=="Negative") #for the description of APOE+/- cohort

	##SCALE THE DATA FOR MLR/BLR
	# procdf1 <- preProcess(as.data.frame(df1), method=c("range")) #All numerical data are put in a range of 0-1
	# normdf1 <- predict(procdf1, as.data.frame(df1))

	# procdf2 <- preProcess(as.data.frame(df2), method=c("range")) #All numerical data are put in a range of 0-1
	# normdf2 <- predict(procdf2, as.data.frame(df2))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #colorblind-friendly palette. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette_RTQUIC <- c("#E69F00", "#999999") #colorblind-friendly palette. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette


###############################################################################################################################
														#METHODS
###############################################################################################################################

#######################			1. DEMOGRAPHICS ON INITIAL DATASET (IE: NO EXCLUSIONS)		  #################################
###############################################################################################################################



# cat("1. COMPARISONS OF DEMOGRAPHICS FOR DX PRIOR TO SUBJECT EXCLUSION: \n")
#############################################################################

# # 1.1. TXT: TOTAL NUMBER + SEX%
# cat("1.1. TOTAL NUMBER + SEX: \n")
# cat("Prior to subject exclusion, the total number of subjects in the dataset is: \n")
# df %>% count(DX_APD)

# cat("Prior to subject exclusion, sex distribution in the dataset is: \n")
# df %>% group_by(DX_APD) %>% count(Sex) 

# # Calculation of %
# totalmatrix <- df %>% count(DX_APD)
# sexmatrix <- df %>% group_by(DX_APD) %>% count(Sex) 

# cat("Prior to subject exclusion, proportion of females in CBS is: \n")
# (as.numeric(sexmatrix$n[1]) + as.numeric(sexmatrix$n[3]))/(as.numeric(totalmatrix$n[1]) + as.numeric(totalmatrix$n[2]))
# cat("Prior to subject exclusion, proportion of females in CBS is:")
# as.numeric(sexmatrix$n[1])/as.numeric(totalmatrix$n[1])
# cat("Prior to subject exclusion, proportion of females in PSP is:")
# as.numeric(sexmatrix$n[3])/as.numeric(totalmatrix$n[2])


# # 1.2. TXT: AGE AT LP
# cat("1.2. AGE AT LP: \n")
	# shapiro.test(CBSdf$Age) #normal
	# shapiro.test(PSPdf$Age) #normal
	# var.test(Age ~ DX_APD, data = df) #homoscedasticity
	# t.test <-t.test(df$Age ~ df$DX_APD, var.equal=TRUE)
	# if (t.test[3] <= 0.05) {
		# cat("Prior to subject exclusion,there is a significant difference in age at LP between CBS and PSP. p-value:", t.test[3][[1]], "\n")
		# cat(t.test[3][[1]])
	# } else cat("Prior to subject exclusion,there is no significant difference in age at LP between CBS and PSP. \n")
# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
# df %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))



	##### TABLE1 #####
	# df %>% count(DX_APD)

	#TABLE1: AGE AT LP
	# t.test(df$Age ~ df$DX_APD, var.equal=TRUE) 
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# df %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))


	#TABLE1: SEX
	# df %>% group_by(DX_APD) %>% count(Sex)
	# table(df$Sex, df$DX_APD)
	# chisq.test(table(df$Sex, df$DX_APD), correct=F)


	#TABLE1: EDUCATION
		# shapiro.test(CBSdf$Education) #not normal
		# shapiro.test(PSPdf$Education) #not normal
		# leveneTest(Education ~ DX_APD, data = df) #homoscedasticity
	# wilcox.test(df$Education ~ df$DX_APD, paired=F)
		# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Education, na.rm=T),2),2), sd=sd(Education, na.rm=T))
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))
		# df %>% summarize(count=n(), format(round(mean(Education, na.rm=T),2),2), sd=sd(Education, na.rm=T))
	# df %>% summarize(count=n(), format(round(median(Education, na.rm=T),2),2), IQR=IQR(Education, na.rm=T), min=min(Education, na.rm=T), max=max(Education, na.rm=T))


	#TABLE1: AGE AT ONSET
	# t.test(df$Onset ~ df$DX_APD, var.equal=TRUE) 
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
	# df %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
		

	#TABLE1: AGE AT PARKINSONISM ONSET
	# t.test(df$Park_onset ~ df$DX_APD, var.equal=TRUE) 
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
	# df %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))


	#TABLE1: DURATION OF DISEASE BY LP
		# shapiro.test(CBSdf$LP2_Disease_Duration) #not normal
		# shapiro.test(PSPdf$LP2_Disease_Duration) #not normal
		# leveneTest(LP2_Disease_Duration ~ DX_APD, data = df) #homoscedasticity
	# wilcox.test(df$LP2_Disease_Duration ~ df$DX_APD, paired=F)
		# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(LP2_Disease_Duration, na.rm=T),2),2), sd=sd(LP2_Disease_Duration, na.rm=T))
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))
		# df %>% summarize(count=n(), format(round(mean(LP2_Disease_Duration, na.rm=T),2),2), sd=sd(LP2_Disease_Duration, na.rm=T))
	# df %>% summarize(count=n(), format(round(median(LP2_Disease_Duration, na.rm=T),2),2), IQR=IQR(LP2_Disease_Duration, na.rm=T), min=min(LP2_Disease_Duration, na.rm=T), max=max(LP2_Disease_Duration, na.rm=T))


	#MoCA Z-SCORE
		# boxplot(MOCA_Z ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
		# 	stripchart(ATI ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot(ATI ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
		# 	stripchart(ATI ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

		# shapiro.test(CBSdf$MOCA_Z) #normal
		# shapiro.test(PSPdf$MOCA_Z) #not normal
		# hist(df$MOCA_Z)
		# hist(CBSdf$MOCA_Z)
		# hist(PSPdf$MOCA_Z)
		# leveneTest(MOCA_Z ~ DX_APD, data = df) #heterodasticity
	# wilcox.test(df$MOCA_Z ~ df$DX_APD, paired=F) #Since looking at z-score, no need to correct by age for this comparison
		# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(MOCA_Z, na.rm=T),2),2), IQR=IQR(MOCA_Z, na.rm=T), min=min(MOCA_Z, na.rm=T), max=max(MOCA_Z, na.rm=T))
		# df %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	# df %>% summarize(count=n(), format(round(median(MOCA_Z, na.rm=T),2),2), IQR=IQR(MOCA_Z, na.rm=T), min=min(MOCA_Z, na.rm=T), max=max(MOCA_Z, na.rm=T))


		#Z-SCORE for ALL COGNITIVE SCORES: NO NEED IN TABLE 1 (REDUNDANT)
			# shapiro.test(CBSdf$Cognitive_Z) #not al
			# shapiro.test(PSPdf$Cognitive_Z) #not al
			# leveneTest(Cognitive_Z ~ DX_APD, data = df) #heterodasticity
		# wilcox.test(df$Cognitive_Z ~ df$DX_APD, paired=F)
		# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(Cognitive_Z, na.rm=T),2),2), IQR=IQR(Cognitive_Z, na.rm=T), min=min(Cognitive_Z, na.rm=T), max=max(Cognitive_Z, na.rm=T))
		# df %>% summarize(count=n(), format(round(median(Cognitive_Z, na.rm=T),2),2), IQR=IQR(Cognitive_Z, na.rm=T), min=min(Cognitive_Z, na.rm=T), max=max(Cognitive_Z, na.rm=T))


	#TABLE1: APOE
	# df %>% group_by(DX_APD) %>% count(APOEe4)
	# table(df$APOEe4, df$DX_APD)
	# chisq.test(table(df$APOEe4, df$DX_APD), correct=F)


	#TABLE1: AD
	# df %>% group_by(DX_APD) %>% count(Lifetime_AD_binary)
	# table(df$Lifetime_AD_binary, df$DX_APD)
	# chisq.test(table(df$Lifetime_AD, df$DX_APD), correct=F)


	# TABLE1: ABETA42
		#Be aware of the outliers but do not remove unless very problematic
		# boxplot(abeta ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
			# stripchart(abeta ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot(abeta ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
			# stripchart(abeta ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

		###Then check ality and heterodasticity
		# shapiro.test(CBSdf$logabeta) #normal
		# shapiro.test(PSPdf$logabeta) #normal
		# leveneTest(logabeta ~ DX_APD, data = df) #homoscedasticity

		# t.test(df$logabeta ~ df$DX_APD, var.equal=TRUE) 
	# aov <- aov(logabeta ~ Age + DX_APD, df) 
	# Anova(aov, type="II") #Compare with type III
		# check_normality(aov)
	# df %>% summarize(count=n(), format(round(mean(abeta, na.rm=T),2),2), sd=sd(abeta, na.rm=T)) #Rounds up the sd for some reason
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(abeta, na.rm=T),2),2), sd=sd(abeta, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdf$abeta)
	# sd(PSPdf$abeta)



	#TABLE1: PTAU181
		#Be aware of the outliers but do not remove unless very problematic
		# boxplot(ptau ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
			# stripchart(ptau ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot(ptau ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
			# stripchart(ptau ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

		###Then check ality and heterodasticity
		# shapiro.test(CBSdf$logptau) #normal
		# shapiro.test(PSPdf$logptau) #normal
		# leveneTest(logptau ~ DX_APD, data = df) #homoscedasticity

		# t.test(df$logptau ~ df$DX_APD, var.equal=TRUE) 
	# aov <- aov(logptau ~ Age + DX_APD, df) 
	# Anova(aov, type="II") #Compare with type III
		# check_normality(aov)
	# df %>% summarize(count=n(), format(round(mean(ptau, na.rm=T),2),2), sd=sd(ptau, na.rm=T)) #Rounds up the sd for some reason
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ptau, na.rm=T),2),2), sd=sd(ptau, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdf$ptau)
	# sd(PSPdf$ptau)



	# TABLE1: TOTAL TAU
		#Be aware of the outliers but do not remove unless very problematic
		# boxplot(ttau ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
			# stripchart(ttau ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot(ttau ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
			# stripchart(ttau ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

		# Then check normality and heterodasticity
		# shapiro.test(CBSdf$logttau) #normal
		# shapiro.test(PSPdf$logttau) #normal
		# leveneTest(logttau ~ DX_APD, data = df) #homoscedasticity

		# t.test(df$logttau ~ df$DX_APD, var.equal=TRUE) 
	# aov <- aov(logttau ~ Age + DX_APD, df) 
	# Anova(aov, type="II") #Compare with type III
		# check_normality(aov)
	# df %>% summarize(count=n(), format(round(mean(ttau, na.rm=T),2),2), sd=sd(ttau, na.rm=T)) #Rounds up the sd for some reason
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ttau, na.rm=T),2),2), sd=sd(ttau, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdf$ttau, na.rm=T)
	# sd(PSPdf$ttau)


	# TABLE1: ATI
		#Be aware of the outliers but do not remove unless very problematic
		# boxplot(ATI ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
			# stripchart(ATI ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot(ATI ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
			# stripchart(ATI ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	 	# Then check ality and heterodasticity
		# shapiro.test(CBSdf$ATI) #non normal
		# shapiro.test(PSPdf$ATI) #non normal
		# shapiro.test(df$ATI) #nonnormal
		# hist(df$ATI)
		# leveneTest(ATI ~ DX_APD, data = df) #homoscedasticity

		# wilcox.test(df$ATI ~ df$DX_APD, paired=F)
	# aov <- aov(ATI ~ Age + DX_APD, df) 
	# Anova(aov, type="II") #Compare with type III
		# check_normality(aov) #Check because the ATI was bimodal including within CBS
	# df %>% summarize(count=n(), format(round(mean(ATI, na.rm=T),2),2), sd=sd(ATI, na.rm=T)) #Rounds up the sd for some reason
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ATI, na.rm=T),2),2), sd=sd(ATI, na.rm=T)) #Rounds up the sd for some reason


	# TABLE1: NFL
	#For NFL the outliers are really impactful. Is best to remove them and report their value. However, will only remove outliers of whole dataset and not per diagnosis. 
	#Remove outliers over full dataset
	# threshold <- min(max(df$logNFL,na.rm=T), as.numeric(quantile(df$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (df$logNFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
			# threshold
			# sort(df$logNFL)
			# subset(df, !(logNFL<threshold))$ID
			# subset(df, !(logNFL<threshold))$NFL
	# dfnfl <- subset(df, (logNFL<threshold))

		# Then check ality and heterodasticity
		# CBSdfnfl <- subset(dfnfl, DX_APD=="CBS")
		# PSPdfnfl <- subset(dfnfl, DX_APD=="PSP")
		# shapiro.test(CBSdfnfl$logNFL) #normal
		# shapiro.test(PSPdfnfl$logNFL) #normal
		# shapiro.test(dfnfl$logNFL) #normal
		# leveneTest(logNFL ~ DX_APD, data = dfnfl) #homoscedasticity

		# t.test(dfnfl$logNFL ~ dfnfl$DX_APD, var.equal=TRUE) 
	# aov <- aov(logNFL ~ Age + DX_APD, dfnfl) 
	# Anova(aov, type="II")
		# check_normality(aov)
	# dfnfl %>% summarize(count=n(), format(round(mean(NFL, na.rm=T),2),2), sd=sd(NFL, na.rm=T)) #Rounds up the sd for some reason
	# dfnfl %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(NFL, na.rm=T),2),2), sd=sd(NFL, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdfnfl$NFL)
	# sd(PSPdfnfl$NFL)




###############################################################################################################################
														#RESULTS
##############################################################################################################################


##################################		ASYN SAA POSITIVITY AND DEMOGRAPHICS		 ###########################################
################################################################################################################################


#COMPARISONS OF DEMOGRAPHICS FOR RTQUIC:
######################################## 

#TXT: TOTAL NUMBER RTQUIC + SEX% WITHIN DIAGNOSIS
	# df %>% group_by(RTQUIC) %>% count(DX_APD)
# df2 %>% group_by(RTQUIC) %>% count(DX_APD)
# sum(df2$DX_APD=="CBS" & df2$RTQUIC=="SAA positive" & df2$Sex=="M")
# sum(df2$DX_APD=="PSP" & df2$RTQUIC=="SAA positive" & df2$Sex=="M")

	#BONFERRONI CALCULATION FOR AGE/ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO RTQUIC
	# 0.05/3 #0.01666667

#TXT: AGE
#No fundamental change by adding APOEe4 (which has interaction with diagnosis)
	# shapiro.test(RTposdf2$Age) #al
	# shapiro.test(RTnegdf2$Age) #al
	# var.test(Age ~ RTQUIC, data = df2) #homoscedasticity
	# t.test(df2$Age ~ df2$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Age ~ DX_APD*RTQUIC, df2)  
	# Anova(aov2, type="III")
# aov <- aov(Age ~ DX_APD + RTQUIC, df2)
# Anova(aov, type="II")
# df2 %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))

#TXT: AGE OF ONSET 
	# shapiro.test(RTposdf2$Onset) #al
	# shapiro.test(RTnegdf2$Onset) #al
	# var.test(Onset ~ RTQUIC, data = df2) #homoscedasticity
	# t.test(df2$Onset ~ df2$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Onset ~ DX_APD*RTQUIC, df2)  
	# Anova(aov2, type="III")
# aov <- aov(Onset ~ DX_APD + RTQUIC, df2)
# Anova(aov, type="II")
# df2 %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))
# df2 %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))


#TXT: AGE AT PARKINSONISM ONSET 
	# shapiro.test(RTposdf2$Park_onset) #al
	# shapiro.test(RTnegdf2$Park_onset) #al
	# var.test(Park_onset ~ RTQUIC, data = df2) #homoscedasticity
	# t.test(df2$Park_onset ~ df2$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Park_onset ~ DX_APD*RTQUIC, df2)  
	# Anova(aov2, type="III")
# aov <- aov(Park_onset ~ DX_APD + RTQUIC, df2)
# Anova(aov, type="II")
# df2 %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Park_onset, na.rm=T))
# df2 %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))


#TXT: TYPE OF ONSET
	# fisher.test(table(df2$Parkinsonian_onset, df2$RTQUIC))
# table(df2$Parkinsonian_onset, df2$RTQUIC)
	# fisher.test(table(df2$Cognitive_onset, df2$RTQUIC))
	# fisher.test(table(df2$Language_onset, df2$RTQUIC))

#TXT: SEX
# fisher.test(table(df2$Sex, df2$RTQUIC)) # Expected count is <5 for one cell

# TXT: SEX
# fisher.test(table(df2$APOEe4, df2$RTQUIC)) # Expected count is <5 for one cell




##################################			CHARACTERIZATION OF THE AD+ GROUP		  ###########################################
##################################################################################################################################


#COMPARISONS OF DEMOGRAPHICS FOR AD:
####################################


#TXT: PREVALENCE OF DX BY AD GROUP
# Here remove CBS-HIV as looking at DX #USE LIFETIME_AD
# chisq.test(table(df1$Lifetime_AD_binary, df1$DX_APD), correct=F)
# table(df1$Lifetime_AD_binary, df1$DX_APD) 
# sum(df1$Lifetime_AD_binary=="AD Positive" & df1$DX_APD=="CBS" & df1$Sex=="M")
# sum(df1$Lifetime_AD_binary=="AD Positive" & df1$DX_APD=="PSP" & df1$Sex=="M")


#TXT: PREVALENCE OF APOEe4 BY AD GROUP
#Here interested in AD and vars independent of DX so keeping CBS-HIV
# df1 %>% group_by(Lifetime_AD_binary) %>% count(APOEe4)
	# chisq.test(table(df1$Lifetime_AD_binary, df1$APOEe4), correct=T)

	
	#BONFERRONI CALCULATION FOR ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO AD 
	# 0.05/2 #0.025

#TXT: DIFFERENCES IN AGE OF ONSET BETWEEN AD+ AND AD-
#Here remove CBS-HIV as looking at DX: df1 dataset. 
	# shapiro.test(CBSdf1$Onset) #al
	# shapiro.test(PSPdf1$Onset) #al
	# shapiro.test(ADposdf1$Onset) #al
	# shapiro.test(ADnegdf1$Onset) #al
	# shapiro.test(APOEposdf1$Onset) #al
	# shapiro.test(APOEnegdf1$Onset) #al
	# leveneTest(Onset ~ DX_APD*APOEe4*Lifetime_AD_binary, df1) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

	# Compare models with Ftest. Cannot include APOE as there is missing data. 
	# test1 <- lm(Onset ~ Lifetime_AD_binary, df1)
	# test2 <- lm(Onset ~ DX_APD, df1)
	# test3 <- lm(Onset ~ DX_APD + Lifetime_AD_binary, df1)
	# anova(test1, test3)
	# anova(test2, test3) #dx_APD by itself is better thna with AD
	# df1apoe <- df1[!is.na(df1$APOEe4), ]
	# test1 <- lm(Onset ~ APOEe4, df1apoe)
	# test2 <- lm(Onset ~ DX_APD, df1apoe)
	# test3 <- lm(Onset ~ DX_APD +APOEe4, df1apoe)
	# anova(test1, test3) 
	# anova(test2, test3) #dx_APD by itself is better thna with APOEe4. Also tested with Sex.

	#For the sake of clarity and because the SME analysis does not indicate APOEe4 to be a very important factor, report lm2 instead. 
	# aov <- aov(Onset ~ DX_APD*APOEe4 + Lifetime_AD_binary, df1) #The model reported as we are interested to describe the effect (here non-effect) of AD
	# Anova(aov, type="III") #Since it includes interaction, report III
	# emmeans(aov, ~ DX_APD:APOEe4, opt.digits=T) #opt.digits: if F, base R getOption setting used. However, this may be more precision than justified using SE. 
	# pairs(emmeans(aov,~ DX_APD:APOEe4))

# aov2 <- aov(Onset ~ DX_APD + Lifetime_AD_binary, df1) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(aov2)
# Anova(aov2, type="II") 
	# shapiro.test(residuals(aov2)) #Careful as not al
	# df1 %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
	# df1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))


#TXT: DIFFERENCES IN AGE OF ONSET AT PARKINSONISM BETWEEN AD+ AND AD-
#Here remove CBS-HIV as looking at DX: df1 dataset. 
	# shapiro.test(CBSdf1$Park_onset) #al
	# shapiro.test(PSPdf1$Park_onset) #al
	# shapiro.test(ADposdf1$Park_onset) #al
	# shapiro.test(ADnegdf1$Park_onset) #al
	# shapiro.test(APOEposdf1$Park_onset) #al
	# shapiro.test(APOEnegdf1$Park_onset) #al
	# leveneTest(Park_onset ~ DX_APD*APOEe4*Lifetime_AD_binary, df1) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# aov2 <- aov(Park_onset ~ DX_APD + Lifetime_AD_binary, df1) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(aov2)
# Anova(aov2, type="II") 
	# shapiro.test(residuals(aov2))
	# df1 %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
	# df1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))



#TXT: DIFFERENCES IN COGNITIVE SCORES BETWEEN AD+ and AD- 
# Here remove CBS-HIV as looking at DX
	# boxplot(MOCA_Z ~ DX_APD, data= CBSdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = CBSdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ DX_APD, data= PSPdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = PSPdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADposdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADposdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADnegdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADnegdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEposdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEposdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEnegdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEnegdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	# Then check ality of my group for AD and DX groups, and homoscedasticity
	# shapiro.test(CBSdf1$MOCA_Z) 
	# shapiro.test(PSPdf1$MOCA_Z) #nonal
	# shapiro.test(ADposdf1$MOCA_Z)
	# shapiro.test(ADnegdf1$MOCA_Z) #nonal
	# shapiro.test(APOEposdf1$MOCA_Z)
	# shapiro.test(APOEnegdf1$MOCA_Z) #nonal
	# leveneTest(MOCA_Z ~ DX_APD*APOEe4*Lifetime_AD_binary, df1) #Homoscedasticity. 

		#Observing some bimodality in the z-scores. Explore:
		# hist(df1$MOCA_Z)
		# testdf <- subset(df1, Onset<=65)
		# hist(testdf$MOCA_Z)
		# testdf <- subset(df1, Onset>65)
		# hist(testdf$MOCA_Z) #Cannot tell source of bimodality.
		# Seems it is just some coincidence + artefact of the nature of z-scores (the left skew is al for a demented cohort)
		# As a result, preference for MLR.

	# #Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# cor.test(df1$MOCA_Z, df1$Age) #correlated
	# summary(lm(df1$MOCA_Z ~ df1$Age)) #linear relationship
	# ggscatter(df1, x = "Age", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# 	stat_regline_equation(aes()) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# 	stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
  	# 	stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "Lifetime_AD_binary", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
  	# 	stat_regline_equation(aes(color = Lifetime_AD_binary))
	##NOTE: Interaction of Age with either DX or AD on Cognitive z-score (CBS/AD+ show much more steep differences in cognitive impairment
			#relative to demographics depending on age. IE young AD+ has much worse CI than a young AD-, but an old AD+ and an old AD- are comparable.
			#This points to AD+ subjects being much more affected at young age.  

	#Inclusion of Duration as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# cor.test(df1$MOCA_Z, df1$LP2_Disease_Duration) #no correlation
	# summary(lm(df1$MOCA_Z ~ df1$LP2_Disease_Duration)) #no linear relationship

	##Compare models with Ftest
	# test1 <- lm(MOCA_Z ~ DX_APD, df1)
	# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, df1)
	# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD, df1) #both additions are beneficial so model should have AD and DX
	# anova(test1, test3)
	# anova(test2, test3)

	# #Inclusion of APOE in model:
	# df1apoe <- df1[!is.na(df1$APOEe4), ]
	# test1 <- lm(MOCA_Z ~ APOEe4, df1apoe)
	# test2 <- lm(MOCA_Z ~ DX_APD, df1apoe)
	# test3 <- lm(MOCA_Z ~ DX_APD +APOEe4, df1apoe)
	# anova(test1, test3) 
	# anova(test2, test3) #DX by itself is better than with APOEe4

# Do not include interaction of DX_APD with AD because too few PSP-AD subjects. 
#Include the interaction of age with AD though. 
# mlr <- lm(MOCA_Z ~ DX_APD + Lifetime_AD_binary*Age, df1) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(mlr)
	# shapiro.test(residuals(mlr))
# emmeans(mlr, ~ Lifetime_AD_binary, opt.digits=T) #opt.digits: if F, base R getOption setting used. However, this may be more precision than justified using SE. 
	
	## ALL COGNITIVE SCORES FOR CONFIRMATION
	# mlr <- lm(Cognitive_Z ~ DX_APD + Lifetime_AD_binary*Age, df1)
	# summary(mlr)



###############################			ASYN SAA POSITIVITY AND ASSOCIATION WITH AD  	######################################
##############################################################################################################################


# DESCRIPTION OF AD/RTQUIC RELATIONSHIP WITH FREQUENCY DATA
###########################################################


#TXT: TOTAL NUMBER RTQUIC + SEX% WITHIN AD
	# df %>% count(RTQUIC)
	# df2 %>% count(RTQUIC)
# table(df2$RTQUIC, df2$AD_binary)
# df2 %>% count(AD_binary)
# sum(df2$AD_binary=="AD Positive" & df2$RTQUIC=="SAA positive" & df2$Sex=="M")
# sum(df2$AD_binary=="AD Positive" & df2$RTQUIC=="SAA positive" & df2$Sex=="F")
# sum(df2$AD_binary=="AD Negative" & df2$RTQUIC=="SAA positive" & df2$Sex=="M")
# sum(df2$AD_binary=="AD Negative" & df2$RTQUIC=="SAA positive" & df2$Sex=="F")

	#BONFERRONI CALCULATION FOR FISHER TESTS USED FOR ASYN/AD RELATIONSHIP IN OVERALL DATASET AND YO DATASET 
	# 0.05/2 #0.025

#TXT: RTQUIC CHISQUARE WITH AD
# Compare proportions between AD+ and AD- for RT positivity
	# table(df2$AD_binary, df2$RTQUIC)
	# chisq.test(table(df2$AD_binary, df2$RTQUIC), correct=T) #correct: Yates'continuity correction
# fisher.test(table(df2$AD_binary, df2$RTQUIC)) # Expected count is <5 for one cell. Could use one-way but kind of weak to do so. (alternative="greater")
	# cramerV(table(df2$AD_binary, df2$RTQUIC))
# phi(table(df2$AD_binary, df2$RTQUIC), digits=6)


#TXT: RTQUIC CHISQUARE WITH AD IN YOD ONLY
#Analysis of AD*RTQUIC association in the young cohort
# df2young <- subset(df2, 65>=Onset)
# nrow(df2young)
# fisher.test(table(df2young$AD_binary, df2young$RTQUIC)) # Expected count is <5 for one cell
	# cramerV(table(df2young$AD_binary, df2young$RTQUIC))
	# table(df2young$AD_binary, df2young$RTQUIC)
# phi(table(df2young$AD_binary, df2young$RTQUIC), digits=6)
	# df2young$Onset

# #TXT: RTQUIC CHISQUARE WITH AD IN LOD ONLY
# df2old <- subset(df2, Onset>65)
# nrow(df2old)
# fisher.test(table(df2old$AD_binary, df2old$RTQUIC)) # Expected count is <5 for one cell
	# cramerV(table(df2old$AD_binary, df2old$RTQUIC))
# phi(table(df2old$AD_binary, df2old$RTQUIC), digits=6)
	# df2old$Onset


#TXT: PROPORTION OF AD- RTQUIC+ SUBJECTS in YOD and LOD
# sum(df2$AD_binary=="AD Positive" & df2$RTQUIC=="SAA positive" & 65>=df2$Onset)
# sum(df2$RTQUIC=="SAA positive" & 65>=df2$Onset)
# sum(df2$AD_binary=="AD Positive" & df2$RTQUIC=="SAA positive" & df2$Onset>65)
# sum(df2$RTQUIC=="SAA positive" & df2$Onset>65)



# DESCRIPTION OF AD/RTQUIC RELATIONSHIP WITH ABETA42
#####################################################

#TXT: ABETA COMPARISON BETWEEN ASYN+ VS ASYN- 

#Use preprocessed dataset (alized so all cont variables fall within range)
	# boxplot(logabeta ~ RTQUIC, data= RTposdf2, col = "white")$out #identify outliers: in AD positive df, none. 
		# stripchart(logabeta ~ RTQUIC, data = RTposdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(logabeta ~ RTQUIC, data= RTnegdf2, col = "white")$out #create vector with outlier values (from AD-negative group)
		# stripchart(logabeta ~ RTQUIC, data = RTnegdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# df2abeta <- df2 #No positive outlier
	# df2$abeta
	# nrow(df2abeta)
	# head(df2abeta)

	# RTposdf2abeta<- subset(df2abeta, RTQUIC=="SAA positive")
	# RTnegdf2abeta <- subset(df2abeta, RTQUIC=="SAA negative")
	# CBSdf2abeta <- subset(df2abeta, DX_APD=="CBS")
	# PSPdf2abeta <- subset(df2abeta, DX_APD=="PSP")

	# shapiro.test(RTposdf2abeta$logabeta) #normal
	# shapiro.test(RTnegdf2abeta$logabeta)  #normal
	# shapiro.test(CBSdf2abeta$logabeta)  #normal
	# shapiro.test(PSPdf2abeta$logabeta)  #normal
	# leveneTest(logabeta ~ DX_APD*RTQUIC, df2abeta) #homo

	# Compare models with Ftest
	# test1 <- lm(logabeta ~ RTQUIC, df2abeta)
	# test2 <- lm(logabeta ~ RTQUIC  + DX_APD, df2abeta)
	# test3 <- lm(logabeta ~ RTQUIC + Sex, df2abeta) #not adding any value to the model
	# anova(test1, test2) 
	# anova(test1, test3) 

	# Inclusion of Onset as a covariate:
	# Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(df2abeta, x = "Onset", y = "logabeta", add = "reg.line")+ 
  		# stat_regline_equation(aes()) 
	# ggscatter(df2abeta, x = "Onset", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
  		# stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(df2abeta, x = "Onset", y = "logabeta", color = "RTQUIC", add = "reg.line")+
  		# stat_regline_equation(aes(color = RTQUIC))
	# cor.test(df2abeta$logabeta, df2abeta$Onset) #corr
	# summary(lm(df2abeta$logabeta ~ df2abeta$Onset)) #linear relationship 

	# Inclusion of NFL as a covariate:
	# Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(df2abeta, x = "NFL", y = "logabeta", add = "reg.line")+ 
  		# stat_regline_equation(aes()) 
	# ggscatter(df2abeta, x = "NFL", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
  		# stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(df2abeta, x = "NFL", y = "logabeta", color = "RTQUIC", add = "reg.line")+
  		# stat_regline_equation(aes(color = RTQUIC))
	# cor.test(df2abeta$logabeta, df2abeta$NFL) #corr
	# summary(lm(df2abeta$logabeta ~ df2abeta$NFL)) #linear relationship 

	# Inclusion of other variables as covariate:
	# summary(lm(df2abeta$logabeta ~ df2abeta$LP2_Disease_Duration)) #no linear relationship
	# summary(lm(df2abeta$logabeta ~ df2abeta$ttau)) #no linear relationship
	# summary(lm(df2abeta$logabeta ~ df2abeta$ptau)) #no linear relationship


#Model
# mlr <- lm(logabeta ~ Onset*RTQUIC + DX_APD + NFL, df2abeta) 
# summary(mlr)


	# Diagnostics of the model run
	# check_normality(mlr) #Ok normality of residuals
	# autoplot(mlr, which = 1:6) #All plots including scale location plot
		# df2abeta[5, c(1, 30, 190, 193,195, 232, 301)] #Does not change fundamentally results (p of 0.0516 for interaction term with this addition)
		# df2abeta[21, c(1, 30, 190, 193, 195, 232, 301)] #Does not change fundamentally results
		# df2abeta[27, c(1, 30, 190, 193,195, 232, 301)] #Does not change fundamentally results
		# df2abeta[37, c(1, 30, 190,193, 195, 232, 301)] #Does not change fundamentally results
		# df2abeta[45, c(1, 30, 190, 193,195, 232, 301)] #Does not change fundamentally results
		# df2abeta[65, c(1, 30, 190, 193,195, 232, 301)] #Does not change fundamentally results
		# test <- subset(df2abeta, ID!="T417, 2214") #for some analysis, need to exclude the potential false negative too
		# mlr <- lm(logabeta ~ Onset*RTQUIC + DX_APD + NFL, test) 
		# summary(mlr) 

		# plot(mlr, which = 3) # 3 = Scale-Location plot. Variance of residuals
	# bptest(mlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
	# durbinWatsonTest(mlr) #Ok autocorrelation of residuals	

		#Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
		# df2[(cooks.distance(mlr))>0.25, ]$ID #no subject >0.25
		# cooksD <- cooks.distance(mlr)
		# n <- nrow(df2abeta) #Anothr threshold for Cook's data = 4/n so lower. In this case, too many high Cook. 
		# plot(cooksD, main = "Cooks Distance for Influential Obs") #Plot Cook's Distance with a horizontal line at 4/n to see which observations exceed this thresdhold. There is a cluster of high Cook subjects. 
			# abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line

	#Multicollinearity checks
	# car::vif(mlr)
	# mlrinteracless <- lm(logabeta ~ Onset + RTQUIC + DX_APD + NFL,df2abeta)
		# car::vif(mlrinteracless)

# Simple slopes for onset: 
# emtrends(mlr,pairwise ~  RTQUIC, var="Onset")

# Main effects: 
# emmeans(mlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
# emmeans(mlr, ~ DX_APD) #adjusted means. Ok because no interaction. 


# MODEL SCALING THE CONTINUOUS VARIABLES FOR 
# stdmlr <- lm(logabeta ~ scale(Onset)*RTQUIC + DX_APD + scale(NFL), df2abeta) 
# summary(stdmlr)


	# Diagnostics of the model run
	# check_normality(stdmlr) #Ok normality of residuals
	# autoplot(stdmlr, which = 1:6) #All plots including scale location plot

		# #Visualize the values for each of the IDs that are indexed on above plots
		# df2abeta[5, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		# df2abeta[21, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		# df2abeta[27, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		# df2abeta[37, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		# df2abeta[45, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		# df2abeta[65, c("ID", "abeta", "logabeta", "Onset", "RTQUIC", "NFL")] 
		
		#for-loop to test the model without each of these values (ie once without outlier 1, then outlier 2, etc)
		# vecIDs <- df2abeta[c(5,21,27,37,45,65), "ID"] #Create vector of each ID
		# for (i in vecIDs) {
			# test <- subset(df2abeta, ID!=i) #for some analysis, need to exclude the potential false negative too
			# stdmlr <- lm(logabeta ~ scale(Onset)*RTQUIC + DX_APD + scale(NFL), test) 
			# print(summary(stdmlr))
		# }

		# plot(stdmlr, which = 3) # 3 = Scale-Location plot. Variance of residuals
	# bptest(stdmlr) #Ok variance of residulas. Breusch-Pagan test for heterodasticity.
	# durbinWatsonTest(stdmlr) #Ok autocorrelation of residuals	

		#Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
		# df2abeta[(cooks.distance(stdmlr))>0.25, ]$ID #no subject >0.25
		# cooksD <- cooks.distance(stdmlr)
		# n <- nrow(df2abeta) #Anothr threshold for Cook's data = 4/n so lower. In this case, too many high Cook. 
		# plot(cooksD, main = "Cooks Distance for Influential Obs") #Plot Cook's Distance with a horizontal line at 4/n to see which observations exceed this thresdhold. There is a cluster of high Cook subjects. 
			# abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line

	# Multicollinearity checks
	# car::vif(stdmlr) #no multicollinearity at all

# Simple slopes for onset: 
# emtrends(stdmlr,pairwise ~  RTQUIC, var="Onset")

#Main effects: 
# emmeans(stdmlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
# emmeans(stdmlr, ~ DX_APD) #adjusted means. Ok because no interaction. 

	# FIG1A: ABETA42 over time from the real datapoints
	##Create df for the figure based on the actual model
	# mylist <- list(Onset=seq(35,85,by=5), RTQUIC=c("aSyn-SAA positive", "aSyn-SAA negative")) 
		# emmip(stdmlr, RTQUIC ~ Onset, at=mylist, CIs=TRUE)

	#Dataset: data are from the model (predicted values)
	# figdf <- emmip(stdmlr, RTQUIC ~ Onset, at=mylist, CIs=TRUE, plotit=FALSE) #
    
    # label1 <- "paste(F*'(5, 59)'==4.29*', ' ~~ italic(P)-value, \" = .018\")" #First annotation of the plot: Model diagnostics. Comma has to be entered as text not in mathematical notation. 
    # label2 <- "paste(italic(R)^2==20.47*'%')" #Second annotation is the Rsquare of teh model. To enter as character, use *'character'. There cannot be parentheses inside. 

    #Change name of variable RTQUIC
	# fig1a_ver1 <- ggplot(data=figdf, aes(x=Onset,y=yvar, color=RTQUIC)) + #yvar is Abeta42 logged 	#Ggplot figure basic layout


					# Add the actual trend/slopes + CI around it
					# geom_line() + #DO NOT use show.legend as it will create a duplicate legend if combined with scale_color changes
					# geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=RTQUIC), alpha=0.2, show.legend=FALSE) + #show.legend=FALSE to avoid having both ribbon and line in legend

					#Fix legends
					# scale_fill_manual(values=cbPalette_RTQUIC) +
					# scale_color_manual(values=cbPalette_RTQUIC, name= "ASyn-SAA status", breaks=c("aSyn-SAA positive", "aSyn-SAA negative"), labels=c(expression(alpha*"Syn-SAA+"),expression(alpha*"Syn-SAA-"))) + #expression allows you to add greek letters

					#Add actual data
					# geom_point(data=df2abeta, aes(x=Onset, y=logabeta)) + #Actual datapoints 

					#Annotate:
					# stat_regline_equation(label.x=35, label.y=c(7.3, 7.2), show.legend=FALSE) + #shows the equation for each geom_line. Cannot use other options since it would not be based on the model
			
				  	# annotate("text", size=5, x=48, y=7.7, label=label1, parse=TRUE)+ #label as defined above. Parse allows for use of mathematical notation
					# annotate("text", size=5, x=40, y=7.6, label=label2, parse=TRUE)+
				

					#Add the labels
					# labs(title=expression(bold("Visual representation of the interaction of age at onset by "*alpha*"Syn-SAA")),
						 # subtitle=expression("From the model: log(A"*beta*"42) ~ Age at onset * "*alpha*"Syn-SAA + Diagnosis + NfL"),
						 # x="Age at onset (years)",
						 # y=expression(bold("CSF A"*beta*"42 levels (pg/mL) (log)")))+ #bold() is required otherwise the Y axis will not be bold in spite of element_text specification below

					#Aesthetic only
					# theme_classic() +
					# theme(plot.title = element_text(size=16, hjust=0.5, face="bold")) +
					# theme(plot.subtitle = element_text(size=14, hjust=0.5)) +
			    	# theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
			    	# theme(legend.title = element_text(face="bold", size=14), legend.text= element_text(size=13))
	fig1a_ver1

	# ggsave(fig1a_ver1, filename = "Fig1a_ver1.png", bg= "transparent", width=9, height=10)

	#figure with scatterplot overlaid
	# fig1a_ver2 <- fig1a_ver1 +
	# 				geom_point(data=df2abeta, aes(x=Onset, y=logabeta))

	# fig1a_ver2



###################################	BINARY LOGISTIC REGRESSION INCLUDING CLINICAL ###########################################
##############################################################################################################################

#TXT: BLR
# df2 <- df2 %>% mutate(RTQUIC = case_when(RTQUIC == "SAA positive" ~ 1, RTQUIC == "SAA negative" ~ 0)) %>% data.frame() 
	# sum(df2$RTQUIC ==1)

# library(boot)
# library(perm)

#Model is selected based on previous analyses (logabeta*Onset), simple comparisons in Supp material (Gait/RBD_binary).
#Model was run with and without AD/DX_APD status and all values remained very similar. 
# blr <- glm(RTQUIC ~ DX_APD + scale(Onset)*scale(logabeta) + RBD_binary + Gait, data= df2, family = "binomial")
# summary(blr)
				
#Model fit
# pscl::pR2(blr)["McFadden"]

	#Model diagnostics and assumptions
	#Multicollinearity checks
	# car::vif(blr)
	# blronset <- glm(RTQUIC ~ DX_APD + Onset + RBD_binary + Gait, data= df2, family = "binomial") #Check individual relationship of Onset with RTQUIC
		# car::vif(blronset)
	# blrabeta <- glm(RTQUIC ~ DX_APD + logabeta + RBD_binary + Gait, data= df2, family = "binomial") #Check individual relationship of Abeta with RTQUIC
		# car::vif(blrabeta)
	# df2 <- df2 %>%
        			# mutate(centredbeta = logabeta- mean(logabeta)) %>%
					# mutate(centredonset = Onset- mean(Onset)) %>%
        			# data.frame()
    # centredblr <- glm(RTQUIC ~ DX_APD + centredbeta*centredonset + RBD_binary + Gait, data= df2, family = "binomial")
    # summary(centredblr)
    # car::vif(centredblr) #Now all VIF are low, indicating no problem of multicolllinearity. 

    #Diagnostics
	# df2[(cooks.distance(blr))>0.06, ]$ID #identifies subjects with high Cook's distance data. 
	# test <- subset(df2, ID!="T258") #for some analysis, need to exclude the potential false negative too
	# blr2 <- glm(RTQUIC ~  DX_APD + Onset*logabeta + RBD_binary + Gait, data= test, family = "binomial")
	# summary(blr2)

	#Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
	# df2[(cooks.distance(blr))>0.25, ]$ID #two subjects >0.25
	# cooksD <- cooks.distance(blr)
	# n <- nrow(df2abeta) #Another threshold for Cook,s data = 4/n so lower. In this case, too many high Cook. 
	# plot(cooksD, main = "Cooks Distance for Influential Obs") #Two subjects are clearly different than others. However, 1 is one of few RBD+ subjects & 1 is one of few AD+/RT+: hard to exclude them as analysis pointless then.
	# 	abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line
		# df2 %>% group_by(RTQUIC) %>% summarize(mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
		# df2 %>% group_by(RTQUIC) %>% summarize(mean=format(round(mean(logabeta, na.rm=T),3),3), sd=format(round(sd(logabeta, na.rm=T),3),3))

	#Dfbetas
	#https://www.statology.org/dfbetas-in-r/
	# dfbetas <- as.data.frame(dfbetas(blr))
	# thresh <- 2/sqrt(nrow(df2))
	
	# par(mfrow=c(2,1))
	# plot(dfbetas$Gait, type='h')
	# abline(h = thresh, lty = 2)
	# abline(h = -thresh, lty = 2)
	# plot(dfbetas$RBD_binary, type='h')
	# abline(h = thresh, lty = 2)
	# abline(h = -thresh, lty = 2)


	#Sample size: https://www.statology.org/assumptions-of-logistic-regression/
	
	#Check linear relationship between logit and explanatory variables: box-tidwell test
	# logodds <- blr$linear.predictors #if (any(x1 <= 0)) stop("the variables to be transformed must have only positive values")
	# boxTidwell(logodds ~ df2$DX_APD)
	# boxTidwell(logodds ~ df2$RBD_binary)
	# boxTidwell(logodds ~ df2$Gait)

#Coefficients
# caret::varImp(blr)	

# #Odds ratio
# coef(blr)
# exp(coef(blr)) 
# cbind(coef(blr),odds_ratio=exp(coef(blr)),exp(confint(blr, level=0.95))) #it says 2.5% and 97.5% because these are the two borders to end up with the cnetral 95% of your distribution

	#Example of interpretation:
	# exp(4.88278) #RBD. 81.60169. Odds of someone with RBD being RTquic+ increases by 8000%. IE x80. 
	# exp(-2.85663) #Gait. 0.09253855. Odds of someone with gait issues being RTQUIC+ decreases by 10% compared to someone without. 
	##### exp(0.15122) #Age: 1.16. #The value indicates that as age increase by one more unit, then the odds of being SAA+ increases by 16%




###############################################	ASYN SAA POSITIVITY AND NFL ##################################################
##############################################################################################################################

#TXT: NFL COMPARISON BETWEEN ASYN+ VS ASYN- 

		# Even if it is just a descriptive analysis, remove outlier to better represent the data
		# boxplot <- boxplot(logNFL ~ DX_APD, data= RTposdf2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
		# 	# stripchart(logNFL ~ DX_APD, data = RTposdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
		# threshold <- min(max(RTposdf2$logNFL,na.rm=T), as.numeric(quantile(RTposdf2$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTposdf2$logNFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
		# df2[!((df2$logNFL<threshold & df2$RTQUIC=="SAA positive")| (df2$RTQUIC=="SAA negative")), c("ID", "NFL")] 
		# df2nfl <- subset(df2, (logNFL<threshold & RTQUIC=="SAA positive") | (RTQUIC=="SAA negative"))

		# boxplot <- boxplot(logNFL ~ DX_APD, data= RTnegdf2, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
			# stripchart(logNFL ~ DX_APD, data = RTnegdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
		# boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
		# threshold <- min(max(RTnegdf2$logNFL,na.rm=T), as.numeric(quantile(RTnegdf2$logNFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTnegdf2$logNFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
		# df2nfl[!((df2nfl$logNFL<threshold & df2nfl$RTQUIC=="SAA negative")| (df2nfl$RTQUIC=="SAA positive")), c("ID", "NFL")] 
		# df2nfl <- subset(df2nfl, (logNFL<threshold & RTQUIC=="SAA negative") | (RTQUIC=="SAA positive"))

		# nrow(df2nfl)

		#RESCALE NFL
		# procdf2nfl <- preProcess(as.data.frame(df2nfl), method=c("range"))
		# df2nfl <- predict(procdf2nfl, as.data.frame(df2nfl))

		# RTposdf2nfl <- subset(df2nfl, RTQUIC=="SAA positive") #for the main results on RT+/-
		# RTnegdf2nfl <- subset(df2nfl, RTQUIC=="SAA negative") #for the main results on RT+/-
		# CBSdf2nfl <- subset(df2nfl, DX_APD=="CBS") #for the main results on RT+/-
		# PSPdf2nfl <- subset(df2nfl, DX_APD=="PSP") #for the main results on RT+/-
		# ADposdf2nfl <- subset(df2nfl, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
		# ADnegdf2nfl <- subset(df2nfl, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort

		# shapiro.test(RTposdf2nfl$logNFL) #al
		# shapiro.test(RTnegdf2nfl$logNFL)  #al
		# shapiro.test(CBSdf2nfl$logNFL)  #al
		# shapiro.test(PSPdf2nfl$logNFL)  #al
		# shapiro.test(ADposdf2nfl$logNFL)  #al
		# shapiro.test(ADnegdf2nfl$logNFL)  #al
		# leveneTest(logNFL ~ DX_APD*RTQUIC, df2nfl) #homo

		# Compare models with Ftest
		# test1 <- lm(logNFL ~ RTQUIC, df2nfl)
		# test2 <- lm(logNFL ~ DX_APD, df2nfl)
		# test3 <- lm(logNFL ~ RTQUIC + DX_APD, df2nfl) #not adding any value to the model
		# anova(test1, test3) 
		# anova(test2, test3) 

		# Inclusion of other variables as covariate:
		# summary(lm(df2nfl$logNFL ~ df2nfl$Age)) #no linear relationship
		# summary(lm(df2nfl$logNFL ~ df2nfl$LP2_Disease_Duration))#no linear relationship 
		# summary(lm(df2nfl$logNFL ~ df2nfl$abeta)) #linear relationship
		# summary(lm(df2nfl$logNFL ~ df2nfl$ptau)) #no linear relationship
		# summary(lm(df2nfl$logNFL ~ df2nfl$ttau)) #no linear relationship


#Model
# mlr <- lm(logNFL ~ RTQUIC + DX_APD + abeta, df2nfl) 
# summary(mlr) 


	# Diagnostics of the model run
	# check_ality(mlr) #Ok
	# durbinWatsonTest(mlr) #Ok

		#Examine Cook's distance: not as importnat given dataset characteristics (sample size limited in terms of number of subjects both AD+ and RT+: cannot exclude easily subjects)
		# df2[(cooks.distance(mlr))>0.25, ]$ID #no subject >0.25
		# cooksD <- cooks.distance(mlr)
		# n <- nrow(df2abeta) #Anothr threshold for Cook,s data = 4/n so lower. In this case, too many high Cook. 
		# plot(cooksD, main = "Cooks Distance for Influential Obs") #Plot Cook's Distance with a horizontal line at 4/n to see which observations exceed this thresdhold. There is a cluster of high Cook subjects. 
		# 	abline(h = 4/67, lty = 2, col = "steelblue") # add cutoff line

#Simple slopes for onset: 
# emtrends(mlr,pairwise ~  RTQUIC, var="Onset")


#Main effects: 
# emmeans(mlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
# emmeans(mlr, ~ DX_APD) #adjusted means. Ok because no interaction. 



############################	TABLE 2 ASYN POSITIVITY AND MOTOR/AUTONOMIC/PD SYMPTOMS   ####################################
#############################################################################################################################


	#TABLE S2: CLINICAL COMPARISONS BETWEEN CBS AND PSP ACROSS ASYN SAAs

	# df2 %>% count(Sex)
	# df2 %>% count(AD_binary)
	# df2 %>% count(APOEe4)

	# df2 %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# df2 %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
	# df2 %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

	# DO CBS FIRST: COUNTS + ANALYSES
	#################################
	# RTposCBSdf2 <- subset(CBSdf2, RTQUIC=="SAA positive")
	# RTnegCBSdf2 <- subset(CBSdf2, RTQUIC=="SAA negative")

	# CBSdf2 %>% group_by(RTQUIC) %>% count(Sex)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(AD_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(APOEe4)

	# CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
	# CBSdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

	# STATISTICAL COMPARISONS FOR AGES AND SEX
	# table(CBSdf2$Sex, CBSdf2$RTQUIC)
	# chisq.test(table(CBSdf2$Sex, CBSdf2$RTQUIC), correct=F)

	# table(CBSdf2$AD_binary, CBSdf2$RTQUIC)
	# fisher.test(table(CBSdf2$AD_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell

	# table(CBSdf2$APOEe4, CBSdf2$RTQUIC)
	# fisher.test(table(CBSdf2$APOEe4, CBSdf2$RTQUIC)) # Expected count is <5 for one cell


		# shapiro.test(RTposCBSdf2$Age) #al
		# shapiro.test(RTnegCBSdf2$Age) #al
		# var.test(Age ~ RTQUIC, data = CBSdf2) #homo
	# t.test(CBSdf2$Age ~ CBSdf2$RTQUIC, var.equal=TRUE) 

		# shapiro.test(RTposCBSdf2$Onset) #al
		# shapiro.test(RTnegCBSdf2$Onset) #al
		# var.test(Onset ~ RTQUIC, data = CBSdf2) #homo
	# t.test(CBSdf2$Onset ~ CBSdf2$RTQUIC, var.equal=TRUE) 

		# shapiro.test(RTposCBSdf2$Park_onset) #al
		# shapiro.test(RTnegCBSdf2$Park_onset) #al
		# var.test(Park_onset ~ RTQUIC, data = CBSdf2) #homo
		# # wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
	# t.test(CBSdf2$Park_onset ~ CBSdf2$RTQUIC, var.equal=TRUE) 

	# CBSdf2 %>% group_by(RTQUIC) %>% count(Tremor_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(RestTremor)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(LimbRigidity)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Slowness_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Falls_PI)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Gait)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(RBD_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
	# ## CBSdf2 %>% group_by(RTQUIC) %>% count(Anosmia_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Constipation_binary)
	# ## CBSdf2 %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Sexual_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Urinary_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Bowel_binary)
	# CBSdf2 %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

# table(CBSdf2$anyPPA, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$anyPPA, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Tremor_binary, CBSdf2$RTQUIC)
# chisq.test(table(CBSdf2$Tremor_binary, CBSdf2$RTQUIC), correct=F)
# table(CBSdf2$RestTremor, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$RestTremor, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$LimbRigidity, CBSdf2$RTQUIC)
# chisq.test(table(CBSdf2$LimbRigidity, CBSdf2$RTQUIC), correct=F)
# table(CBSdf2$Slowness_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Slowness_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Falls_PI, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Falls_PI, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Gait, CBSdf2$RTQUIC)
# chisq.test(table(CBSdf2$Gait, CBSdf2$RTQUIC), correct=F)
# table(CBSdf2$RBD_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$RBD_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Lifetime_Dopa_responder_true, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Lifetime_Dopa_responder_true, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Lifetime_VisualHallucinations_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Lifetime_VisualHallucinations_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Constipation_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Constipation_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Urinary_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Urinary_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdf2$Bowel_binary, CBSdf2$RTQUIC)
# fisher.test(table(CBSdf2$Bowel_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell


	# DO PSP SECOND: COUNTS + ANALYSES
	################################
	# RTposPSPdf2 <- subset(PSPdf2, RTQUIC=="SAA positive")
	# RTnegPSPdf2 <- subset(PSPdf2, RTQUIC=="SAA negative")

	# PSPdf2 %>% group_by(RTQUIC) %>% count(Sex)
	# PSPdf2 %>% group_by(RTQUIC) %>% count(AD_binary)
	# PSPdf2 %>% group_by(RTQUIC) %>% count(APOEe4)
	# PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
	# PSPdf2 %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))



# STATISTICAL COMPARISONS FOR AGES AND SEX
# table(PSPdf2$Sex, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Sex, PSPdf2$RTQUIC)) # Expected count is <5 for one cell

# table(PSPdf2$AD_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$AD_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell

# table(PSPdf2$APOEe4, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$APOEe4, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
#
# shapiro.test(RTposPSPdf2$Age) #al
# shapiro.test(RTnegPSPdf2$Age) #al
# var.test(Age ~ RTQUIC, data = PSPdf2) #homo
# t.test(PSPdf2$Age ~ PSPdf2$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposPSPdf2$Onset) #al
# shapiro.test(RTnegPSPdf2$Onset) #al
# var.test(Onset ~ RTQUIC, data = PSPdf2) #homo
# t.test(PSPdf2$Onset ~ PSPdf2$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposPSPdf2$Park_onset) #al
# shapiro.test(RTnegPSPdf2$Park_onset) #al
# var.test(Park_onset ~ RTQUIC, data = PSPdf2) #homo
# # # wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
# t.test(PSPdf2$Park_onset ~ PSPdf2$RTQUIC, var.equal=TRUE) 


# PSPdf2 %>% group_by(RTQUIC) %>% count(anyPPA)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Tremor_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(RestTremor)
# PSPdf2 %>% group_by(RTQUIC) %>% count(LimbRigidity)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Slowness_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Falls_PI)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Gait)
# PSPdf2 %>% group_by(RTQUIC) %>% count(RBD_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Anosmia_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Constipation_binary)
# ## PSPdf2 %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
# PSPdf2 %>% group_by(RTQUIC) %>% count(Sexual_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Urinary_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Bowel_binary)
# PSPdf2 %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

# table(PSPdf2$Tremor_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Tremor_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$LimbRigidity, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$LimbRigidity, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Slowness_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Slowness_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Falls_PI, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Falls_PI, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Gait, PSPdf2$RTQUIC)
# chisq.test(table(PSPdf2$Gait, PSPdf2$RTQUIC), correct=F)
# table(PSPdf2$RBD_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$RBD_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Lifetime_Dopa_responder_true, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Lifetime_Dopa_responder_true, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Lifetime_VisualHallucinations_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Lifetime_VisualHallucinations_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Constipation_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Constipation_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Orthostatism_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Orthostatism_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Urinary_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Urinary_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdf2$Bowel_binary, PSPdf2$RTQUIC)
# fisher.test(table(PSPdf2$Bowel_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell



	# COMPARISONS IGNORING DX
	# t.test(df2$Age ~ df2$RTQUIC, var.equal=TRUE) 
	# t.test(df2$Onset ~ df2$RTQUIC, var.equal=TRUE) 
	# t.test(df2$Park_onset ~ df2$RTQUIC, var.equal=TRUE) 
	# table(df2$anyPPA, df2$RTQUIC) #Fisher
	# table(df2$Lifetime_Dopa_responder_true, df2$RTQUIC) #Fisher
	# table(df2$AD_binary, df2$RTQUIC) #Fisher
	# table(df2$APOEe4, df2$RTQUIC) #Fisher
	# table(df2$Tremor_binary, df2$RTQUIC) #Chisquare
	# table(df2$RestTremor, df2$RTQUIC) #Fisher
	# table(df2$LimbRigidity, df2$RTQUIC) #Chisquare
	# table(df2$Slowness_binary, df2$RTQUIC) #Fisher
	# table(df2$Falls_PI, df2$RTQUIC) #Fisher
	# table(df2$Gait, df2$RTQUIC) #Chisquare
	# table(df2$RBD_binary, df2$RTQUIC) #Fisher
	# table(df2$Slowness_binary, df2$RTQUIC) #Fisher
	# table(df2$Lifetime_VisualHallucinations_binary, df2$RTQUIC) #Fisher
	# table(df2$Constipation, df2$RTQUIC) #Fisher
	# table(df2$Sexual, df2$RTQUIC) #Fisher
	# table(df2$Orthostatism_binary, df2$RTQUIC) #Fisher
	# table(df2$Urinary, df2$RTQUIC) #Chisquare
	# table(df2$Bowel, df2$RTQUIC) #Fisher
	# table(df2$Thermoregulatory_binary, df2$RTQUIC) #Fisher

# fisher.test(table(df2$anyPPA, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Lifetime_Dopa_responder_true, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$AD_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$APOEe4, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$RestTremor, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Slowness_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$RBD_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Constipation, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Sexual, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Orthostatism_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Bowel, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Thermoregulatory_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Lifetime_VisualHallucinations_binary, df2$RTQUIC)) # Expected count is <5 for one cell
# fisher.test(table(df2$Falls_PI, df2$RTQUIC)) # Expected count is <5 for one cell
# chisq.test(table(df2$Tremor_binary, df2$RTQUIC), correct=F)
# chisq.test(table(df2$LimbRigidity, df2$RTQUIC), correct=F)
# chisq.test(table(df2$Gait, df2$RTQUIC), correct=F)
# chisq.test(table(df2$Urinary, df2$RTQUIC), correct=F)


###################################	SURVIVAL ANALYSES ON THE SEVERE MOTOR DISEASE ###########################################
##############################################################################################################################
# library(survival)
# library(lubridate)
# library(ggsurvfit)
# library(gtsummary)
# library(tidycmprsk)
# library(condsurv)

# df2$Severe_Motor_onset_status_12 <- as.numeric(df2$Severe_Motor_onset_status_12)
# CBSdf2$Severe_Motor_onset_status_12 <- as.numeric(CBSdf2$Severe_Motor_onset_status_12)
# Surv(df2$Severe_onset_diff, df2$Severe_Motor_onset_status_12==2)
# surv1 <- survfit(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ 1, data = df2)
# str(surv1)

# summary(survfit(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ 1, data = df2), times = 15) #average surv time for 5 years

# coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ RTQUIC, data = CBSdf2)
# coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ DX_APD, data = df2)

# surv<-coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ RTQUIC, data = df2)
# summary(surv)


####################################################	MISCELLANEOUS ########################################################
##############################################################################################################################

# #Exploratory: remove other influential points
# lev <- hat(model.matrix(mlr1)) #The leverage of an obs. measures its ability to move the regression model, ie the amount by which the predicted value would change if the obs. was shifted by one unit
# plot(lev)
# df2tauratio[lev > mean(lev)*3, c("ID", "NFL")]
# cooks <- cooks.distance(mlr1)
# plot(cooks)
# 	abline(h = 4/nrow(df2tauratio), lty = 2, col = "steelblue") # add cutoff line
# df2tauratio[cooks > (4/nrow(df2tauratio)), c("ID", "tauratio", "ptau", "ttau")]
# # vec <- df2tauratio[cooks > (4/nrow(df1nflDX)), "ID"]
# as.numeric(names(cooks)[(cooks > (4/nrow(df1nfl)))])
# testdf <- subset(df1nflDX, ID!=vec[1] & ID!=vec[2] & ID!=vec[3] & ID!=vec[4] & ID!=vec[5]) 

# CBSdf2ptau <- subset(df2ptau, DX_APD=="CBS")
# shapiro.test(RTposdf2ptau$logptau) #al
# shapiro.test(RTnegdf2ptau$logptau) #al
# # Run aov model
# # https://statisticsbyjim.com/anova/ancova/
# acovmodel <- aov(logptau ~ AD + DX_APD + RTQUIC + abeta, df2ptau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(acovmodel)
# Anova(acovmodel, type="II") #Since it is an ancova, report III
# Anova(acovmodel, type="III")
# AIC(acovmodel)
# 	df2ptau %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	df2ptau %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	#https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
# 	emmeans(acovmodel, ~ RTQUIC: abeta) #adjusted means
# 	emmeans(acovmodel, ~ DX_APD:abeta) #adjusted means
# acovmodel <- aov(logptau ~ AD + RTQUIC + abeta, CBSdf2ptau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(acovmodel)
# Anova(acovmodel, type="II") #Since it is an ancova, report III
# Anova(acovmodel, type="III")
# AIC(acovmodel)
# 	df2ptau %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	emmeans(acovmodel, ~ RTQUIC:abeta) #adjusted means

# #Exploratory: remove other influential points
# lev <- hat(model.matrix(aovmodel)) #The leverage of an obs. measures its ability to move the regression model, ie the amount by which the predicted value would change if the obs. was shifted by one unit
# plot(lev)
# df1nflDX[lev > mean(lev)*3, c("ID", "NFL")]
# cooks <- cooks.distance(aovmodel)
# plot(cooks)
# 	abline(h = 4/nrow(df1nflDX), lty = 2, col = "steelblue") # add cutoff line
# df1nflDX[cooks > (4/nrow(df1nflDX)), c("ID", "NFL")]
# vec <- df1nflDX[cooks > (4/nrow(df1nflDX)), "ID"]
# # as.numeric(names(cooks)[(cooks > (4/nrow(df1nfl)))])
# testdf <- subset(df1nflDX, ID!=vec[1] & ID!=vec[2] & ID!=vec[3] & ID!=vec[4] & ID!=vec[5]) 
# #Identify and remove outliers in Ttau: here doing it directly on logged data to avoid unnecessary removal of data given it is al enough that log procedure fixes issue. 
# boxplot(RTposdf2$logttau ~ RTposdf2$RTQUIC)$out #no outlier
# 	stripchart(logttau ~ RTQUIC, data = RTposdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# vec <- boxplot(RTnegdf2$logttau ~ RTnegdf2$RTQUIC)$out #2 highest values are outliers
# 	stripchart(logttau ~ RTQUIC, data = RTnegdf2, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# df2ttau <- subset(df2, logttau!=vec[1]) #remove outlier in logged data
# boxplot(df2ptau$logttau ~ df2ttau$RTQUIC)$out
# 	stripchart(logttau ~ RTQUIC, data = df2ttau, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# RTposdf2ttau <- subset(df2ttau, RTQUIC=="SAA positive")
# RTnegdf2ttau <- subset(df2ttau, RTQUIC=="SAA negative")
# shapiro.test(RTposdf2ttau$logttau) #al
# shapiro.test(RTnegdf2ttau$logttau) #al
# #Run aov model
# aovmodel <- aov(logttau ~ AD + DX_APD + RTQUIC + Age + abeta + LP2_Disease_Duration, df2ttau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(aovmodel)
# AIC(aovmodel)
# Anova(aovmodel, type='II') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# df2ttau %>% group_by(RTQUIC) %>% summarize(count=n(), mean= format(round(mean(ttau, na.rm=T),3),3), sd= format(round(sd(ttau, na.rm=T),3),3))
# df2ttau %>% group_by(DX_APD) %>% summarize(count=n(), mean= format(round(mean(ttau, na.rm=T),3),3), sd= format(round(sd(ttau, na.rm=T),3),3))

# vec <- boxplot(df2$abeta ~ df2$RTQUIC)$out
# df2abeta <- subset(df2, abeta!=vec[1] & abeta!=vec[2]) #remove outlier in raw data
# boxplot(df2abeta$logabeta ~ df2abeta$RTQUIC)$out
# 	aovmodel <- aov(logabeta ~ AD + DX_APD + RTQUIC + Age + ptau + ttau, df2abeta) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# 	summary(aovmodel)
# 	Anova(aovmodel, type='II') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# 	df2abeta %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(abeta, na.rm=T), sd=sd(abeta, na.rm=T))
# 	df2abeta %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(abeta, na.rm=T), sd=sd(abeta, na.rm=T))
# # wilcox.test(df2$abeta ~ df2$RTQUIC)
# #Ordinal logistic regression?
# boxplot(abeta ~ RTQUIC, data = df2abeta, col = "white")
# stripchart(abeta ~ RTQUIC,
#            data = df2abeta,
#            method = "jitter",
#            pch = 19,
#            col = 2:4,
#            vertical = TRUE,
#            add = TRUE)

# https://www.scribbr.com/statistics/anova-in-r/






###############################################################################################################################
														#SUPPLEMENTARY MATERIAL
###############################################################################################################################



#####################################			CBS-AD+ vs CBS-AD-	 			###############################################
###############################################################################################################################



	##SUPPLEMENTARY MATERIAL
	# Data are not ally distributed so be careful with residuals
	# Using lm to compare Rsquare to that of test models below (there is no difference with aov)
	# mlr <- lm(Cognitive_Z ~ Lifetime_AD_binary + DX_APD + Age, df1)
	# mlr2 <- lm(Cognitive_Z ~Age*Lifetime_AD_binary + DX_APD, df1) #Technically best model but since interaction remains ns I think for sake of simplicity it is okay to select mlr instead. 
	# mlr3 <- lm(Cognitive_Z ~Age*DX_APD + Lifetime_AD_binary, df1)
	# summary(mlr) 
	# summary(mlr2) 
	# summary(mlr3) 

	# ##Diagnostics of the model run
	# residuals <- residuals(mlr)
	# hist(residuals)
	# qq(residuals)
	# shapiro.test(residuals)
	# check_ality(mlr) #Very non al: model not fittint too well the data
	# durbinWatsonTest(mlr) #Check the residuals are independent (multiple regression assumption)

	##Means to report for model mlr: report EMMs + SE 
	## df1 %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Cognitive_Z, na.rm=T),3),3), sd=format(round(sd(Cognitive_Z, na.rm=T),3),3))
	# emmeans(mlr, ~ AD, opt.digits=T) #adjusted means
	## emmeans(mlr, ~ DX_APD, opt.digits=T) #adjusted means


# MULT LINEAR REGRESSION  TO COMPARE MOCA Z SCORES BETWEEN DX AND AD WITH DURATION AS COVARIATE
# #Here remove CBS-HIV as looking at DX
	# boxplot(MOCA_Z ~ DX_APD, data= CBSdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = CBSdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ DX_APD, data= PSPdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = PSPdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADposdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADposdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADnegdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADnegdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEposdf1, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEposdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEnegdf1, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEnegdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	# #Then check ality of my group for AD and DX groups, and homoscedasticity
	# shapiro.test(CBSdf1$MOCA_Z)
	# shapiro.test(PSPdf1$MOCA_Z) #not al
	# shapiro.test(ADposdf1$MOCA_Z)
	# shapiro.test(ADnegdf1$MOCA_Z) #not al
	# shapiro.test(APOEposdf1$MOCA_Z)
	# shapiro.test(APOEnegdf1$MOCA_Z) #not al
	# leveneTest(MOCA_Z ~ DX_APD*APOEe4*Lifetime_AD_binary, df1) #Homoscedasticity. 

	#Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
	#   stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(df1, x = "Age", y = "MOCA_Z", color = "Lifetime_AD_binary", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
	#   stat_regline_equation(aes(color = Lifetime_AD_binary))
	# cor.test(df1$MOCA_Z, df1$Age)
	# summary(lm(df1$MOCA_Z ~ df1$Age)) #linear relationship
	##NOTE: possible interaction of Age with either DX or AD on Cognitive z-score (CBS/AD+ show much more steep differences in cognitive impairment
			#relative to demographics depending on age. IE young AD+ has much worse CI than a young AD-, but an old AD+ and an old AD- are comparable.
			#This points to AD+ subjects being much more affected at young age. Not sure if makes sense to include in this model). 

	#Inclusion of Duration as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(df1, x = "LP2_Disease_Duration", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(df1, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(df1, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
	#   stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(df1, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "AD", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
	#   stat_regline_equation(aes(color = AD))
	# cor.test(df1$MOCA_Z, df1$LP2_Disease_Duration) #No relationship

	##Compare models with Ftest
	# test1 <- lm(MOCA_Z ~ DX_APD, df1)
	# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, df1)
	# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD, df1) #both additions are beneficial so model should have AD and DX
	# anova(test1, test3)
	# anova(test2, test3)

	# Inclusion of APOE in model:
	# df1apoe <- df1[!is.na(df1$APOEe4), ]
	# test1 <- lm(MOCA_Z ~ APOEe4, df1apoe)
	# test2 <- lm(MOCA_Z ~ DX_APD, df1apoe)
	# test3 <- lm(MOCA_Z ~ DX_APD +APOEe4, df1apoe)
	# anova(test1, test3) 
	# anova(test2, test3) #DX by itself is better than with APOEe4

# Using lm to compare Rsquare to that of test models below (there is no difference with aov)
# mlr <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD + Age, df1)
	# mlr2 <- lm(MOCA_Z ~Age*Lifetime_AD_binary + DX_APD, df1) #Technically best model but since interaction remains ns I think for sake of simplicity it is okay to select mlr instead. 
	# mlr3 <- lm(MOCA_Z ~Age*DX_APD + Lifetime_AD_binary, df1)
# summary(mlr) #type I (no interaction so order does not matter)
	# summary(mlr2) 
	# summary(mlr3) 

# # ##Diagnostics of the model run
# hist(residuals)
# qq(residuals)
# shapiro.test(residuals(mlr))
# check_ality(mlr) #Very non al: model not fittint too well the data
# durbinWatsonTest(mlr) #Check the residuals are independent (multiple regression assumption)

##Means to report for model mlr: report EMMs + SE 
## df1 %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Cognitive_Z, na.rm=T),3),3), sd=format(round(sd(Cognitive_Z, na.rm=T),3),3))
# emmeans(mlr, ~ Lifetime_AD_binary, opt.digits=T) #adjusted means
## emmeans(mlr, ~ DX_APD, opt.digits=T) #adjusted means





	##### TABLES2 #####
	#Remove CBS-HIV patient
	
	# CBSdf1 %>% count(Lifetime_AD_binary)

	# #TABLE1: AGE
	# 	shapiro.test(ADposCBSdf1$Age) #al
	# 	shapiro.test(ADnegCBSdf1$Age) #al
	# 	var.test(Age ~ Lifetime_AD_binary, data = CBSdf1) #homoscedasticity
	# t.test(CBSdf1$Age ~ CBSdf1$Lifetime_AD_binary, var.equal=TRUE)
	# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# CBSdf1 %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))

	# #TABLE1:SEX
	# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% count(Sex)
	# table(CBSdf1$Sex, CBSdf1$Lifetime_AD_binary)
	# chisq.test(table(CBSdf1$Sex, CBSdf1$Lifetime_AD_binary), correct=F)


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			# shapiro.test(ADposCBSdf1$Onset) #al
			# shapiro.test(ADnegCBSdf1$Onset) #al
			# var.test(Onset ~ Lifetime_AD_binary, data = CBSdf1) #homoscedasticity
		# t.test(CBSdf1$Onset ~ CBSdf1$Lifetime_AD_binary, var.equal=TRUE) 
		# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
		# CBSdf1 %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
		
		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		# 	shapiro.test(ADposCBSdf1$Park_onset) #al
		# 	shapiro.test(ADnegCBSdf1$Park_onset) #al
		# 	var.test(Park_onset ~ Lifetime_AD_binary, data = CBSdf1) #homoscedasticity
		# t.test(CBSdf1$Park_onset ~ CBSdf1$Lifetime_AD_binary, var.equal=TRUE) 
		# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
		# CBSdf1 %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			# shapiro.test(ADposCBSdf1$Cognitive_Z) #al
			# shapiro.test(ADnegCBSdf1$Cognitive_Z) #nonal
			# leveneTest(Cognitive_Z ~ Lifetime_AD_binary, data = CBSdf1) #homoscedasticity
		# wilcox.test(CBSdf1$Cognitive_Z ~ CBSdf1$Lifetime_AD_binary, paired=F)
		# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Cognitive_Z, na.rm=T),2),2), sd=sd(Cognitive_Z, na.rm=T))
		# CBSdf1 %>% summarize(count=n(), format(round(mean(Cognitive_Z, na.rm=T),2),2), sd=sd(Cognitive_Z, na.rm=T))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% count(APOEe4)
		# table(CBSdf1$APOEe4, CBSdf1$Lifetime_AD_binary)
		# fisher.test(table(CBSdf1$APOEe4, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell


			###SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
	###COMPARISONS OF BIOMARKERS FOR AD:
	##################################################

			##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		##Look only in CBS: should be no interaction of APOEe4 but possible effect of AD
		# CBSdf1onsetapoe <- CBSdf1onset[!is.na(CBSdf1onset$APOEe4), ]
			# test1 <- lm(Onset ~ Lifetime_AD_binary, CBSdf1onsetapoe)
			# test2 <- lm(Onset ~ APOEe4, CBSdf1onsetapoe)
			# test3 <- lm(Onset ~ APOEe4+ Lifetime_AD_binary, CBSdf1onsetapoe)
			# anova(test1, test3)
			# anova(test2, test3) 
			# ADposCBSdf1onset <- subset(CBSdf1onset, Lifetime_AD_binary=="AD Positive")
			# ADnegCBSdf1onset <- subset(CBSdf1onset, Lifetime_AD_binary=="AD Negative")
			# shapiro.test(ADposCBSdf1onset$Onset) #al
			# shapiro.test(ADnegCBSdf1onset$Onset) #al
			# var.test(Onset ~ Lifetime_AD_binary, data = CBSdf1onset) #homoscedasticity
		# t.test(CBSdf1onset$Onset ~ CBSdf1onset$Lifetime_AD_binary, var.equal=TRUE)
		# CBSdf1onset %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
		# CBSdf1onset %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		##Look only in CBS: should be no interaction of APOEe4 but possible effect of AD
		# CBSdf1onsetapoe <- CBSdf1onset[!is.na(CBSdf1onset$APOEe4), ]
		# test1 <- lm(Park_onset ~ Lifetime_AD_binary, CBSdf1onsetapoe)
		# test2 <- lm(Park_onset ~ APOEe4, CBSdf1onsetapoe)
		# test3 <- lm(Park_onset ~ APOEe4+ Lifetime_AD_binary, CBSdf1onsetapoe)
		# anova(test1, test3)
		# anova(test2, test3) 
		# ADposCBSdf1onset <- subset(CBSdf1onset, Lifetime_AD_binary=="AD Positive")
		# ADnegCBSdf1onset <- subset(CBSdf1onset, Lifetime_AD_binary=="AD Negative")
		# shapiro.test(ADposCBSdf1onset$Park_onset) #al
		# shapiro.test(ADnegCBSdf1onset$Park_onset) #al
		# var.test(Park_onset ~ Lifetime_AD_binary, data = CBSdf1onset) #homoscedasticity
		# t.test(CBSdf1onset$Park_onset ~ CBSdf1onset$Lifetime_AD_binary, var.equal=TRUE)
		# CBSdf1onset %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))
		# CBSdf1onset %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))


	##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
	# table(CBSdf1$Parkinsonian_onset, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Parkinsonian_onset, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$Cognitive_onset, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Cognitive_onset, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$Language_onset, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Language_onset, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$anyPPA, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$anyPPA, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$Tremor_binary, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Tremor_binary, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$Myoclonus, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Myoclonus, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdf1$Language_onset, CBSdf1$Lifetime_AD_binary)
	# fisher.test(table(CBSdf1$Language_onset, CBSdf1$Lifetime_AD_binary)) # Expected count is <5 for one cell



			##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			#Already got rid of outliers earlier. 
			# CBSdf1apoe <- CBSdf1[!is.na(CBSdf1$APOEe4), ]
			# test1 <- lm(MOCA_Z ~ APOEe4, CBSdf1apoe)
			# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, CBSdf1apoe)
			# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary +APOEe4, CBSdf1apoe)
			# anova(test1, test3) 
			# anova(test2, test3) #DX by itself is better than with APOEe
			# shapiro.test(ADposCBSdf1$MOCA_Z) #al
			# shapiro.test(ADnegCBSdf1$MOCA_Z) #al
			# leveneTest(MOCA_Z ~ Lifetime_AD_binary, data = CBSdf1) #homoscedasticity
			# t.test(CBSdf1$MOCA_Z ~ CBSdf1$Lifetime_AD_binary, var.equal=TRUE) 
			# CBSdf1 %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
			# CBSdf1 %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	

	# #MULT LINEAR REGRESSION TO COMPARE NFL LEVELS
	# #For NFL, another approach is chosen, as the data are expected to be very right-skewed. Threshold for outlier is changed.

	# Remove outliers for AD
	#First identify and remove outlier from the AD+ group only. 
	# boxplot <- boxplot(NFL ~ Lifetime_AD_binary, data= ADposCBSdf1, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	# 	stripchart(NFL ~ Lifetime_AD_binary, data = ADposCBSdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# #boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
	# threshold <- min(max(ADposCBSdf1$NFL,na.rm=T), as.numeric(quantile(ADposCBSdf1$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (ADposCBSdf1$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
	# #df1[!((df1$NFL<threshold & df1$Lifetime_AD_binary=="AD Positive")| (df1$Lifetime_AD_binary=="AD Negative")), c("ID", "NFL")] 
	# CBSdf1nfl <- subset(CBSdf1, (NFL<threshold & Lifetime_AD_binary=="AD Positive") | (Lifetime_AD_binary=="AD Negative")) #make sure you don't accidentally erase the AD negative here, or filter their values

	# ##Then identify and remove outlier from the AD- group only
	# boxplot <- boxplot(NFL ~ Lifetime_AD_binary, data= ADnegCBSdf1, col = "white") #Same for AD- group
	# 	stripchart(NFL ~ Lifetime_AD_binary, data = ADnegCBSdf1, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# #boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
	# threshold<- min(max(ADnegCBSdf1$NFL,na.rm=T), as.numeric(quantile(ADnegCBSdf1$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (ADnegCBSdf1$NFL)*3))) 
	# #df1nflAD[!((df1nflAD$NFL<threshold & df1nflAD$Lifetime_AD_binary=="AD Negative")| (df1nflAD$Lifetime_AD_binary=="AD Positive")), c("ID", "NFL")] 
	# CBSdf1nfl <- subset(CBSdf1nfl, (NFL<threshold & Lifetime_AD_binary=="AD Negative") | (Lifetime_AD_binary=="AD Positive")) 


	#Compare models withand without AD after removing the AD outliers with the Ftest
	# test1 <- lm(logNFL ~ Lifetime_AD_binary, CBSdf1nfl)
	# test2 <- lm(logNFL ~ Sex, CBSdf1nfl)
	# test3 <- lm(logNFL ~ Lifetime_AD_binary + Sex, CBSdf1nfl)
	# anova(test1, test3)
	# anova(test2, test3)

	# #Oneway seems better but let's perform F-test to have full picture. 
	# testdf <- CBSdf1nfl %>% drop_na("APOEe4") 
	# testdf$APOEe4
	# oneway <- lm(logNFL ~ Lifetime_AD_binary, testdf) 
	# twoway <- lm(logNFL ~ APOEe4 + Lifetime_AD_binary, testdf) 
	# anova(oneway, twoway)

	# #Then check ality in both AD+ and AD- and overall
	# ADposCBSdf1nfl <- subset(CBSdf1nfl, Lifetime_AD_binary=="AD Positive")
	# ADnegCBSdf1nfl <- subset(CBSdf1nfl, Lifetime_AD_binary=="AD Negative")
	# shapiro.test(ADposCBSdf1nfl$logNFL) #al
	# shapiro.test(ADnegCBSdf1nfl$logNFL) #nonal
	# leveneTest(logNFL ~ Lifetime_AD_binary, CBSdf1nfl) #Heterodasticity

	#Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdf1nfl, x = "Age", y = "NFL", add = "reg.line")+ #Not related to age
  	# 	stat_regline_equation(aes()) 
	# ggscatter(CBSdf1nfl, x = "Age", y = "NFL", color = "Lifetime_AD_binary", add = "reg.line")+
  	# 	stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdf1nfl$NFL, CBSdf1nfl$Age) #no corr
	# cor.test(CBSdf1nfl$logNFL, CBSdf1nfl$Age) #no corr

	# #Inclusion of Duration as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdf1nfl, x = "LP2_Disease_Duration", y = "NFL", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(CBSdf1nfl, x = "LP2_Disease_Duration", y = "NFL", color = "Lifetime_AD_binary", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdf1nfl$NFL, CBSdf1nfl$LP2_Disease_Duration) #no corr
	# cor.test(CBSdf1nfl$logNFL, CBSdf1nfl$LP2_Disease_Duration) #no corr

	# #Inclusion of Abeta as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "abeta", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "abeta", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdf1nfl$abeta, CBSdf1nfl$NFL)
	# cor.test(CBSdf1nfl$abeta, CBSdf1nfl$logNFL)
	# cor.test(CBSdf1nfl$logabeta, CBSdf1nfl$NFL)
	# summary(lm(CBSdf1nfl$abeta ~ CBSdf1nfl$logNFL)) #no linear relationship

	# #Inclusion of ptau as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "ptau", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "ptau", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdf1nfl$ptau, CBSdf1nfl$NFL)
	# cor.test(CBSdf1nfl$ptau, CBSdf1nfl$logNFL)
	# summary(lm(CBSdf1nfl$ptau ~ CBSdf1nfl$logNFL)) #no linear relationship

	# #Inclusion of ptau as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "ttau", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdf1nfl, x = "NFL", y = "ttau", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdf1nfl$ttau, CBSdf1nfl$NFL)
	# cor.test(CBSdf1nfl$ttau, CBSdf1nfl$logNFL)
	# summary(lm(CBSdf1nfl$ttau ~ CBSdf1nfl$logNFL)) #no linear relationship

	# # Data are not ally distributed so be careful with residuals
	# #Many options were tried: most were eliminated as model significance way >0.10. Interactions including aforementioned covariates were tried just in case as visualization was ambiguous. 
	# t.test(ADposCBSdf1nfl$logNFL, ADnegCBSdf1nfl$logNFL) #Welch's t-test: does not assume variance equal
	# 	CBSdf1nfl %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(NFL, na.rm=T),3),3), sd=format(round(sd(NFL, na.rm=T),3),3))
	# 	CBSdf1nfl %>% summarize(count=n(), mean=format(round(mean(NFL, na.rm=T),3),3), sd=format(round(sd(NFL, na.rm=T),3),3))
	# CBSdf1nfl %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(median(NFL, na.rm=T),2),2), IQR=format(round(IQR(NFL, na.rm=T),2),2), min=min(NFL, na.rm=T), max=max(NFL, na.rm=T))
	# CBSdf1nfl %>% summarize(count=n(), format(round(median(NFL, na.rm=T),2),2), IQR=IQR(NFL, na.rm=T), min=min(NFL, na.rm=T), max=max(NFL, na.rm=T))




		#TXT: RTQUIC CHISQUARE WITH AD IN CBS ONLY
	# #Compare proportions between AD+ and AD- for RT positivity but in CBS only: Close to significance. 
	# 	table(CBSdf2$AD_binary, CBSdf2$RTQUIC)
	# 	chisq.test(table(CBSdf2$AD_binary, CBSdf2$RTQUIC)) #correct: Yates'continuity correction
	# fisher.test(table(CBSdf2$AD_binary, CBSdf2$RTQUIC)) # Expected count is <5 for one cell
	# 	# cramerV(table(CBSdf2$AD_binary, CBSdf2$RTQUIC))
	# phi(table(CBSdf2$AD_binary, CBSdf2$RTQUIC), digits=6)

	#TXT: RTQUIC CHISQUARE WITH AD IN PSP ONLY
	# #Compare proportions between AD+ and AD- for RT positivity but in CBS only
	# 	chisq.test(table(PSPdf2$AD_binary, PSPdf2$RTQUIC), correct=T) #correct: Yates'continuity correction
	# fisher.test(table(PSPdf2$AD_binary, PSPdf2$RTQUIC)) # Expected count is <5 for one cell
	# 	table(PSPdf2$AD_binary, PSPdf2$RTQUIC)
	# cramerV(table(PSPdf2$AD_binary, PSPdf2$RTQUIC))
