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
# library(multcomp) #multiple comparisons
library(car) #very useful for Anova(), Levene, etc. Please note that base R anova() and car:Anova() differ for:
##contrast option already set as c("contr.sum", "contr.poly") for car::Anova(). Also, car::Anova() does type II SS
##method as default, whereas base R anova() does type I SS method (same results as summary(aov))
library(emmeans) #estimated marginal means of a model
library(modelbased) #extract some parameters and estimates from a model
library(performance) #check_normality
library(rcompanion) #Cramer's V
library(caret)
library(psych)
# library(ggstatsplot) #plot figures
library(ggplot2) #plot figures with stats
library(ggpubr) #ggscatter




###############################################################################################################################
													#SETUP OF THE DATA
###############################################################################################################################


#EXCLUDE SUBJECTS
dfRS <- subset(df, ID!="T102, T192")
dfRSDW <- subset(df, ID!="T102, T192" & ID!="FTLD_010") #for some analysis, need to exclude the MAPT carrier
dfMVRS <- subset(df, ID!="2308_v1" & ID!="T102, T192") #for some analysis, need to exclude the potential false negative too


##CREATE SOME OF THE SUBSETS USED LATER (MOSTLY FOR OUTLIER IDENTIFICATION)
##CREATE A FUNCTION FOR THIS LATER
#df
CBSdf <- subset(df, DX_APD=="CBS") #for the description of AD+/- cohort
PSPdf <- subset(df, DX_APD=="PSP") #for the description of AD+/- cohort
ADposdf <- subset(df, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdf <- subset(df, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdf <- subset(df, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdf <- subset(df, APOEe4=="Negative") #for the description of APOE+/- cohort

#dfRS
RTposdfRS <- subset(dfRS, RTQUIC=="SAA positive") #for the main results on RT+/-
RTnegdfRS <- subset(dfRS, RTQUIC=="SAA negative") #for the main results on RT+/-
CBSdfRS <- subset(dfRS, DX_APD=="CBS") #for the main results on RT+/-
PSPdfRS <- subset(dfRS, DX_APD=="PSP") #for the main results on RT+/-
ADposdfRS <- subset(dfRS, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdfRS <- subset(dfRS, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdfRS <- subset(dfRS, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdfRS <- subset(dfRS, APOEe4=="Negative") #for the description of APOE+/- cohort

#dfMVRS
RTposdfMVRS <- subset(dfMVRS, RTQUIC=="SAA positive") #for the main results on RT+/-
RTnegdfMVRS <- subset(dfMVRS, RTQUIC=="SAA negative") #for the main results on RT+/-
CBSdfMVRS <- subset(dfMVRS, DX_APD=="CBS") #for the main results on RT+/-
PSPdfMVRS <- subset(dfMVRS, DX_APD=="PSP") #for the main results on RT+/-
ADposdfMVRS <- subset(dfMVRS, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegdfMVRS <- subset(dfMVRS, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposdfMVRS <- subset(dfMVRS, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegdfMVRS <- subset(dfMVRS, APOEe4=="Negative") #for the description of APOE+/- cohort

##NORMALIZE THE DATA FOR MLR
#dfRS
procdfRS <- preProcess(as.data.frame(dfRS), method=c("range")) #All numerical data are put in a range of 0-1
normdfRS <- predict(procdfRS, as.data.frame(dfRS))

RTposnormdfRS <- subset(normdfRS, RTQUIC=="SAA positive") #for the main results on RT+/-
RTnegnormdfRS <- subset(normdfRS, RTQUIC=="SAA negative") #for the main results on RT+/-
CBSnormdfRS <- subset(normdfRS, DX_APD=="CBS") #for the main results on RT+/-
PSPnormdfRS <- subset(normdfRS, DX_APD=="PSP") #for the main results on RT+/-
ADposnormdfRS <- subset(normdfRS, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegnormdfRS <- subset(normdfRS, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposnormdfRS <- subset(normdfRS, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegnormdfRS <- subset(normdfRS, APOEe4=="Negative") #for the description of APOE+/- cohortRTposdfRS <- subset(dfRS, RTQUIC=="SAA positive") #for the main results on RT+/-

#dfMVRS
procdfMVRS <- preProcess(as.data.frame(dfMVRS), method=c("range"))
normdfMVRS <- predict(procdfMVRS, as.data.frame(dfMVRS))

RTposnormdfMVRS <- subset(normdfMVRS, RTQUIC=="SAA positive") #for the main results on RT+/-
RTnegnormdfMVRS <- subset(normdfMVRS, RTQUIC=="SAA negative") #for the main results on RT+/-
CBSnormdfMVRS <- subset(normdfMVRS, DX_APD=="CBS") #for the main results on RT+/-
PSPnormdfMVRS <- subset(normdfMVRS, DX_APD=="PSP") #for the main results on RT+/-
ADposnormdfMVRS <- subset(normdfMVRS, Lifetime_AD_binary=="AD Positive") #for the description of AD+/- cohort
ADnegnormdfMVRS <- subset(normdfMVRS, Lifetime_AD_binary=="AD Negative") #for the description of AD+/- cohort
APOEposnormdfMVRS <- subset(normdfMVRS, APOEe4=="Positive") #for the description of APOE+/- cohort
APOEnegnormdfMVRS <- subset(normdfMVRS, APOEe4=="Negative") #for the description of APOE+/- cohortRTposdfRS <- subset(dfRS, RTQUIC=="SAA positive") #for the main results on RT+/-
RTnegnormdfMVRS <- subset(normdfMVRS, RTQUIC=="SAA negative") #for the main results on RT+/-


###############################################################################################################################
														#METHODS
###############################################################################################################################



##########################			DEMOGRAPHICS ON INITIAL DATASET (IE: NO EXCLUSIONS)		  #################################
###############################################################################################################################



#COMPARISONS OF DEMOGRAPHICS FOR DX:
#################################### 

#TXT: TOTAL NUMBER + SEX%
# df %>% count(DX_APD)
# df %>% group_by(DX_APD) %>% count(Sex)

#TXT: ONSET OF DISEASE
	# shapiro.test(CBSdf$Onset) #normal
	# shapiro.test(PSPdf$Onset) #normal
	# var.test(Onset ~ DX_APD, data = df) #homoscedasticity
# t.test(df$Onset ~ df$DX_APD, var.equal=TRUE) 
# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
# df %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))

#TXT: ONSET OF PARKINSONISM
	# shapiro.test(CBSdf$Park_onset) #normal
	# shapiro.test(PSPdf$Park_onset) #normal
	# hist(CBSdf$Park_onset)
	# hist(PSPdf$Park_onset)
	# var.test(Park_onset ~ DX_APD, data = df) #homoscedasticity
# t.test(df$Park_onset ~ df$DX_APD, var.equal=TRUE) 
# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
# df %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

#TXT: AGE AT LP
	# shapiro.test(CBSdf$Age) #normal
	# shapiro.test(PSPdf$Age) #normal
	# var.test(Age ~ DX_APD, data = df) #homoscedasticity
# t.test(df$Age ~ df$DX_APD, var.equal=TRUE)
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


	#TABLE1: APOE
	# df %>% group_by(DX_APD) %>% count(APOEe4)
	# table(df$APOEe4, df$DX_APD)
	# chisq.test(table(df$APOEe4, df$DX_APD), correct=F)


	#TABLE1: AD
	# df %>% group_by(DX_APD) %>% count(Lifetime_AD_binary)
	# table(df$Lifetime_AD_binary, df$DX_APD)
	# chisq.test(table(df$Lifetime_AD, df$DX_APD), correct=F)


	#TABLE1: ABETA42
	##Even if just a descriptive table, for the biomarkers prefer to remove outliers. It only really makes a difference for the NFL though. 
		#Even if it is just a descriptive analysis, remove outlier to better represent the data
	# 	vec1 <- boxplot(abeta ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	# 		stripchart(abeta ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	boxplot(abeta ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	# 		stripchart(abeta ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	dfabeta <- subset(df, abeta!=vec1[1]) #remove outlier

	# 	###Then check normality and heterodasticity
	# 	CBSdfabeta <- subset(dfabeta, DX_APD=="CBS")
	# 	PSPdfabeta <- subset(dfabeta, DX_APD=="PSP")
	# 	shapiro.test(CBSdfabeta$logabeta) #normal
	# 	shapiro.test(PSPdfabeta$logabeta) #normal
	# 	leveneTest(logabeta ~ DX_APD, data = dfabeta) #homoscedasticity

	# 	# t.test(dfabeta$logabeta ~ dfabeta$DX_APD, var.equal=TRUE) 
	# aov <- aov(logabeta ~ Age + DX_APD, dfabeta) 
	# Anova(aov, type="II") #Compare with type III
	# dfabeta %>% summarize(count=n(), format(round(mean(abeta, na.rm=T),2),2), sd=sd(abeta, na.rm=T)) #Rounds up the sd for some reason
	# dfabeta %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(abeta, na.rm=T),2),2), sd=sd(abeta, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdfabeta$abeta)
	# sd(PSPdfabeta$abeta)



	# #TABLE1: PTAU181
	# 	#Even if it is just a descriptive analysis, remove outlier to better represent the data
	# 	vec1 <- boxplot(ptau ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	# 		stripchart(ptau ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	boxplot(ptau ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	# 		stripchart(ptau ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	dfptau <- subset(df, ptau!=vec1[1] & ptau!=vec1[2]) #remove outlier

	# 	###Then check normality and heterodasticity
	# 	CBSdfptau <- subset(dfptau, DX_APD=="CBS")
	# 	PSPdfptau <- subset(dfptau, DX_APD=="PSP")
	# 	shapiro.test(CBSdfptau$logptau) #normal
	# 	shapiro.test(PSPdfptau$logptau) #normal
	# 	leveneTest(logptau ~ DX_APD, data = dfptau) #homoscedasticity

	# 	t.test(dfptau$logptau ~ dfptau$DX_APD, var.equal=TRUE) 
	# aov <- aov(logptau ~ Age + DX_APD, dfptau) 
	# Anova(aov, type="II") #Compare with type III
	# dfptau %>% summarize(count=n(), format(round(mean(ptau, na.rm=T),2),2), sd=sd(ptau, na.rm=T)) #Rounds up the sd for some reason
	# dfptau %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ptau, na.rm=T),2),2), sd=sd(ptau, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdfptau$ptau)
	# sd(PSPdfptau$ptau)



	# #TABLE1: TOTAL TAU
	# #	#Even if it is just a descriptive analysis, remove outlier to better represent the data
	# 	vec1<- boxplot(ttau ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	# 		stripchart(ttau ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	vec2<- boxplot(ttau ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	# 		stripchart(ttau ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	dfttau <- subset(df, ttau!=vec1[1] & ttau!=vec2[1]) #remove outlier

	#	#Then check normality and heterodasticity
	# 	CBSdfttau <- subset(dfttau, DX_APD=="CBS")
	# 	PSPdfttau <- subset(dfttau, DX_APD=="PSP")
	# 	shapiro.test(CBSdfttau$logttau) #normal
	# 	shapiro.test(PSPdfttau$logttau) #normal
	# 	leveneTest(logttau ~ DX_APD, data = dfttau) #homoscedasticity

	# 	t.test(dfttau$logttau ~ dfttau$DX_APD, var.equal=TRUE) 
	# aov <- aov(logttau ~ Age + DX_APD, dfttau) 
	# Anova(aov, type="II") #Compare with type III
	# dfttau %>% summarize(count=n(), format(round(mean(ttau, na.rm=T),2),2), sd=sd(ttau, na.rm=T)) #Rounds up the sd for some reason
	# dfttau %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ttau, na.rm=T),2),2), sd=sd(ttau, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdfttau$ttau)
	# sd(PSPdfttau$ttau)



	# #TABLE1: ATI
	#	#Even if it is just a descriptive analysis, remove outlier to better represent the data
	# 	boxplot(ATI ~ DX_APD, data= CBSdf, col = "white")$out #identify outliers in each diagnosis. First, look at CBS: there is one so attribute its value to vector.   
	# 		stripchart(ATI ~ DX_APD, data = CBSdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# 	boxplot(ATI ~ DX_APD, data= PSPdf, col = "white")$out #Now identify outliers in PSP: since there is none, no need to attribute to a vector. 
	# 		stripchart(ATI ~ DX_APD, data = PSPdf, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	#  	#Then check normality and heterodasticity
	# 	shapiro.test(CBSdf$ATI) #nonnormal
	# 	shapiro.test(PSPdf$ATI) #nonnormal
	# 	shapiro.test(df$ATI) #normal
	# 	hist(df$ATI)
	# 	leveneTest(ATI ~ DX_APD, data = df) #homoscedasticity

	# 	wilcox.test(df$ATI ~ df$DX_APD, paired=F)
	# aov <- aov(ATI ~ Age + DX_APD, df) 
	# Anova(aov, type="II") #Compare with type III
	# 	check_normality(aov) #Check because the ATI was bimodal including within CBS
	# df %>% summarize(count=n(), format(round(mean(ATI, na.rm=T),2),2), sd=sd(ATI, na.rm=T)) #Rounds up the sd for some reason
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(ATI, na.rm=T),2),2), sd=sd(ATI, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdf$ATI)
	# sd(PSPdf$ATI)



	# TABLE1: NFL
	# Even if it is just a descriptive analysis, remove outlier to better represent the data
		# boxplot <- boxplot(NFL ~ DX_APD, data= CBSdf, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
		# # boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
		# threshold <- min(max(CBSdf$NFL,na.rm=T), as.numeric(quantile(CBSdf$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdf$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
		# dfnfl <- subset(df, (NFL<threshold & DX_APD=="CBS") | (DX_APD=="PSP"))

		# boxplot <- boxplot(NFL ~ DX_APD, data= PSPdf, col = "white")
		# # boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
		# threshold <- min(max(PSPdf$NFL,na.rm=T), as.numeric(quantile(PSPdf$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (PSPdf$NFL)*3)))
		# dfnfl <- subset(dfnfl, (NFL<threshold & DX_APD=="PSP") | (DX_APD=="CBS")) 
	
		# Then check normality and heterodasticity
		# CBSdfnfl <- subset(dfnfl, DX_APD=="CBS")
		# PSPdfnfl <- subset(dfnfl, DX_APD=="PSP")
		# shapiro.test(CBSdfnfl$logNFL) #normal
		# shapiro.test(PSPdfnfl$logNFL) #normal
		# shapiro.test(dfnfl$logNFL) #normal
		# leveneTest(logNFL ~ DX_APD, data = dfnfl) #homoscedasticity

		# t.test(dfnfl$logNFL ~ dfnfl$DX_APD, var.equal=TRUE) 
	# aov <- aov(logNFL ~ Age + DX_APD, dfnfl) 
	# Anova(aov, type="II")
	# dfnfl %>% summarize(count=n(), format(round(mean(NFL, na.rm=T),2),2), sd=sd(NFL, na.rm=T)) #Rounds up the sd for some reason
	# dfnfl %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(NFL, na.rm=T),2),2), sd=sd(NFL, na.rm=T)) #Rounds up the sd for some reason
	# sd(CBSdfnfl$NFL)
	# sd(PSPdfnfl$NFL)


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
	# 	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(MOCA_Z, na.rm=T),2),2), IQR=IQR(MOCA_Z, na.rm=T), min=min(MOCA_Z, na.rm=T), max=max(MOCA_Z, na.rm=T))
	# 	# df %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	# df %>% summarize(count=n(), format(round(median(MOCA_Z, na.rm=T),2),2), IQR=IQR(MOCA_Z, na.rm=T), min=min(MOCA_Z, na.rm=T), max=max(MOCA_Z, na.rm=T))


		#Z-SCORE for ALL COGNITIVE SCORES: NO NEED IN TABLE 1 (REDUNDANT)
			# shapiro.test(CBSdf$Cognitive_Z) #not normal
			# shapiro.test(PSPdf$Cognitive_Z) #not normal
			# leveneTest(Cognitive_Z ~ DX_APD, data = df) #heterodasticity
		# wilcox.test(df$Cognitive_Z ~ df$DX_APD, paired=F)
		# df %>% group_by(DX_APD) %>% summarize(count=n(), format(round(median(Cognitive_Z, na.rm=T),2),2), IQR=IQR(Cognitive_Z, na.rm=T), min=min(Cognitive_Z, na.rm=T), max=max(Cognitive_Z, na.rm=T))
		# df %>% summarize(count=n(), format(round(median(Cognitive_Z, na.rm=T),2),2), IQR=IQR(Cognitive_Z, na.rm=T), min=min(Cognitive_Z, na.rm=T), max=max(Cognitive_Z, na.rm=T))





# ###############################################################################################################################
# 														#RESULTS
# ###############################################################################################################################


##################################		ASYN SAA POSITIVITY AND DEMOGRAPHICS		 ###########################################
################################################################################################################################


#COMPARISONS OF DEMOGRAPHICS FOR RTQUIC:
######################################## 

#TXT: TOTAL NUMBER RTQUIC + SEX% WITHIN DIAGNOSIS
	# df %>% group_by(RTQUIC) %>% count(DX_APD)
# dfMVRS %>% group_by(RTQUIC) %>% count(DX_APD)
# sum(dfMVRS$DX_APD=="CBS" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="M")
# sum(dfMVRS$DX_APD=="PSP" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="M")

	#BONFERRONI CALCULATION FOR AGE/ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO RTQUIC
	# 0.05/3 #0.01666667

#TXT: AGE
#No fundamental change by adding APOEe4 (which has interaction with diagnosis)
	# shapiro.test(RTposdfMVRS$Age) #normal
	# shapiro.test(RTnegdfMVRS$Age) #normal
	# var.test(Age ~ RTQUIC, data = dfMVRS) #homoscedasticity
	# t.test(dfMVRS$Age ~ dfMVRS$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Age ~ DX_APD*RTQUIC, dfMVRS)  
	# Anova(aov2, type="III")
# aov <- aov(Age ~ DX_APD + RTQUIC, dfMVRS)
# Anova(aov, type="II")
# dfMVRS %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Age, na.rm=T), sd=sd(Age, na.rm=T))

#TXT: AGE OF ONSET 
	# shapiro.test(RTposdfMVRS$Onset) #normal
	# shapiro.test(RTnegdfMVRS$Onset) #normal
	# var.test(Onset ~ RTQUIC, data = dfMVRS) #homoscedasticity
	# t.test(dfMVRS$Onset ~ dfMVRS$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Onset ~ DX_APD*RTQUIC, dfMVRS)  
	# Anova(aov2, type="III")
# aov <- aov(Onset ~ DX_APD + RTQUIC, dfMVRS)
# Anova(aov, type="II")
# dfMVRS %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))
# dfMVRS %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))


#TXT: AGE AT PARKINSONISM ONSET 
	# shapiro.test(RTposdfMVRS$Park_onset) #normal
	# shapiro.test(RTnegdfMVRS$Park_onset) #normal
	# var.test(Park_onset ~ RTQUIC, data = dfMVRS) #homoscedasticity
	# t.test(dfMVRS$Park_onset ~ dfMVRS$RTQUIC, var.equal=TRUE)
	# aov2 <- aov(Park_onset ~ DX_APD*RTQUIC, dfMVRS)  
	# Anova(aov2, type="III")
# aov <- aov(Park_onset ~ DX_APD + RTQUIC, dfMVRS)
# Anova(aov, type="II")
# dfMVRS %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(Park_onset, na.rm=T), sd=sd(Park_onset, na.rm=T))
# dfMVRS %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(Onset, na.rm=T), sd=sd(Onset, na.rm=T))


#TXT: TYPE OF ONSET
	# fisher.test(table(dfMVRS$Parkinsonian_onset, dfMVRS$RTQUIC))
# table(dfMVRS$Parkinsonian_onset, dfMVRS$RTQUIC)
	# fisher.test(table(dfMVRS$Cognitive_onset, dfMVRS$RTQUIC))
	# fisher.test(table(dfMVRS$Language_onset, dfMVRS$RTQUIC))

#TXT: SEX
# fisher.test(table(dfMVRS$Sex, dfMVRS$RTQUIC)) # Expected count is <5 for one cell





##################################			CHARACTERIZATION OF THE AD+ GROUP		  ###########################################
##################################################################################################################################


# #COMPARISONS OF DEMOGRAPHICS FOR AD:
# ####################################


#TXT: PREVALENCE OF DX BY AD GROUP
# Here remove CBS-HIV as looking at DX #USE LIFETIME_AD
# chisq.test(table(dfRS$Lifetime_AD_binary, dfRS$DX_APD), correct=F)
# table(dfRS$Lifetime_AD_binary, dfRS$DX_APD) 
# sum(dfRS$Lifetime_AD_binary=="AD Positive" & dfRS$DX_APD=="CBS" & dfRS$Sex=="M")
# sum(dfRS$Lifetime_AD_binary=="AD Positive" & dfRS$DX_APD=="PSP" & dfRS$Sex=="M")


#TXT: PREVALENCE OF APOEe4 BY AD GROUP
#Here interested in AD and vars independent of DX so keeping CBS-HIV
# dfRS %>% group_by(Lifetime_AD_binary) %>% count(APOEe4)
	# chisq.test(table(dfRS$Lifetime_AD_binary, dfRS$APOEe4), correct=T)

	
	#BONFERRONI CALCULATION FOR ONSET AGE/PARKINSONISM AGE COMPARISONS ATTRIBUTABLE TO AD 
	# 0.05/2 #0.025

#TXT: DIFFERENCES IN AGE OF ONSET BETWEEN AD+ AND AD-
#Here remove CBS-HIV as looking at DX: dfRS dataset. 
	# shapiro.test(CBSdfRS$Onset) #normal
	# shapiro.test(PSPdfRS$Onset) #normal
	# shapiro.test(ADposdfRS$Onset) #normal
	# shapiro.test(ADnegdfRS$Onset) #normal
	# shapiro.test(APOEposdfRS$Onset) #normal
	# shapiro.test(APOEnegdfRS$Onset) #normal
	# leveneTest(Onset ~ DX_APD*APOEe4*Lifetime_AD_binary, dfRS) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

	# Compare models with Ftest. Cannot include APOE as there is missing data. 
	# test1 <- lm(Onset ~ Lifetime_AD_binary, dfRS)
	# test2 <- lm(Onset ~ DX_APD, dfRS)
	# test3 <- lm(Onset ~ DX_APD + Lifetime_AD_binary, dfRS)
	# anova(test1, test3)
	# anova(test2, test3) #dx_APD by itself is better thna with AD
	# dfRSapoe <- dfRS[!is.na(dfRS$APOEe4), ]
	# test1 <- lm(Onset ~ APOEe4, dfRSapoe)
	# test2 <- lm(Onset ~ DX_APD, dfRSapoe)
	# test3 <- lm(Onset ~ DX_APD +APOEe4, dfRSapoe)
	# anova(test1, test3) 
	# anova(test2, test3) #dx_APD by itself is better thna with APOEe4. Also tested with Sex.

	#For the sake of clarity and because the SME analysis does not indicate APOEe4 to be a very important factor, report lm2 instead. 
	# aov <- aov(Onset ~ DX_APD*APOEe4 + Lifetime_AD_binary, dfRS) #The model reported as we are interested to describe the effect (here non-effect) of AD
	# Anova(aov, type="III") #Since it includes interaction, report III
	# emmeans(aov, ~ DX_APD:APOEe4, opt.digits=T) #opt.digits: if F, base R getOption setting used. However, this may be more precision than justified using SE. 
	# pairs(emmeans(aov,~ DX_APD:APOEe4))

# aov2 <- aov(Onset ~ DX_APD + Lifetime_AD_binary, dfRS) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(aov2)
# Anova(aov2, type="II") 
	# shapiro.test(residuals(aov2)) #Careful as not normal
	# dfRS %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
	# dfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))


#TXT: DIFFERENCES IN AGE OF ONSET AT PARKINSONISM BETWEEN AD+ AND AD-
#Here remove CBS-HIV as looking at DX: dfRS dataset. 
	# shapiro.test(CBSdfRS$Park_onset) #normal
	# shapiro.test(PSPdfRS$Park_onset) #normal
	# shapiro.test(ADposdfRS$Park_onset) #normal
	# shapiro.test(ADnegdfRS$Park_onset) #normal
	# shapiro.test(APOEposdfRS$Park_onset) #normal
	# shapiro.test(APOEnegdfRS$Park_onset) #normal
	# leveneTest(Park_onset ~ DX_APD*APOEe4*Lifetime_AD_binary, dfRS) #Homoscedasticity. Specify saturated model, ie includes the interaction term even if aovmodel does not

# aov2 <- aov(Park_onset ~ DX_APD + Lifetime_AD_binary, dfRS) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(aov2)
# Anova(aov2, type="II") 
	# shapiro.test(residuals(aov2))
	# dfRS %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
	# dfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))



#TXT: DIFFERENCES IN COGNITIVE SCORES BETWEEN AD+ and AD- 
# Here remove CBS-HIV as looking at DX
	# boxplot(MOCA_Z ~ DX_APD, data= CBSnormdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = CBSnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ DX_APD, data= PSPnormdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = PSPnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADposnormdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADposnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADnegnormdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADnegnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEposnormdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEposnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEnegnormdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEnegnormdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	# Then check normality of my group for AD and DX groups, and homoscedasticity
	# shapiro.test(CBSnormdfRS$MOCA_Z) 
	# shapiro.test(PSPnormdfRS$MOCA_Z) #nonnormal
	# shapiro.test(ADposnormdfRS$MOCA_Z)
	# shapiro.test(ADnegnormdfRS$MOCA_Z) #nonnormal
	# shapiro.test(APOEposnormdfRS$MOCA_Z)
	# shapiro.test(APOEnegnormdfRS$MOCA_Z) #nonnormal
	# leveneTest(MOCA_Z ~ DX_APD*APOEe4*Lifetime_AD_binary, normdfRS) #Homoscedasticity. 

		#Observing some bimodality in the z-scores. Explore:
		# hist(normdfRS$MOCA_Z)
		# testdf <- subset(normdfRS, Onset<=65)
		# hist(testdf$MOCA_Z)
		# testdf <- subset(normdfRS, Onset>65)
		# hist(testdf$MOCA_Z) #Cannot tell source of bimodality.
		# Seems it is just some coincidence + artefact of the nature of z-scores (the left skew is normal for a demented cohort)
		# As a result, preference for MLR.

	# #Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# cor.test(normdfRS$MOCA_Z, normdfRS$Age) #correlated
	# summary(lm(normdfRS$MOCA_Z ~ normdfRS$Age)) #linear relationship
	# ggscatter(normdfRS, x = "Age", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# 	stat_regline_equation(aes()) 
	# ggscatter(normdfRS, x = "Age", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# 	stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(normdfRS, x = "Age", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
  	# 	stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(normdfRS, x = "Age", y = "MOCA_Z", color = "Lifetime_AD_binary", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
  	# 	stat_regline_equation(aes(color = Lifetime_AD_binary))
	##NOTE: Interaction of Age with either DX or AD on Cognitive z-score (CBS/AD+ show much more steep differences in cognitive impairment
			#relative to demographics depending on age. IE young AD+ has much worse CI than a young AD-, but an old AD+ and an old AD- are comparable.
			#This points to AD+ subjects being much more affected at young age.  

	#Inclusion of Duration as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# cor.test(normdfRS$MOCA_Z, normdfRS$LP2_Disease_Duration) #no correlation
	# summary(lm(normdfRS$MOCA_Z ~ normdfRS$LP2_Disease_Duration)) #no linear relationship

	##Compare models with Ftest
	# test1 <- lm(MOCA_Z ~ DX_APD, normdfRS)
	# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, normdfRS)
	# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD, normdfRS) #both additions are beneficial so model should have AD and DX
	# anova(test1, test3)
	# anova(test2, test3)

	# #Inclusion of APOE in model:
	# dfRSapoe <- normdfRS[!is.na(normdfRS$APOEe4), ]
	# test1 <- lm(MOCA_Z ~ APOEe4, dfRSapoe)
	# test2 <- lm(MOCA_Z ~ DX_APD, dfRSapoe)
	# test3 <- lm(MOCA_Z ~ DX_APD +APOEe4, dfRSapoe)
	# anova(test1, test3) 
	# anova(test2, test3) #DX by itself is better than with APOEe4

# Do not include interaction of DX_APD with AD because too few PSP-AD subjects. 
#Include the interaction of age with AD though. 
# mlr <- lm(MOCA_Z ~ DX_APD + Lifetime_AD_binary*Age, normdfRS) #The model reported as we are interested to describe the effect (here non-effect) of AD
# summary(mlr)
	# shapiro.test(residuals(mlr))
# emmeans(mlr, ~ Lifetime_AD_binary, opt.digits=T) #opt.digits: if F, base R getOption setting used. However, this may be more precision than justified using SE. 
	
	## ALL COGNITIVE SCORES FOR CONFIRMATION
	# mlr <- lm(Cognitive_Z ~ DX_APD + Lifetime_AD_binary*Age, normdfRS)
	# summary(mlr)



###############################			ASYN SAA POSITIVITY AND ASSOCIATION WITH AD  	######################################
##############################################################################################################################


# DESCRIPTION OF AD/RTQUIC RELATIONSHIP WITH FREQUENCY DATA
###########################################################


#TXT: TOTAL NUMBER RTQUIC + SEX% WITHIN AD
	# df %>% count(RTQUIC)
	# dfMVRS %>% count(RTQUIC)
# table(dfMVRS$RTQUIC, dfMVRS$AD_binary)
# dfMVRS %>% count(AD_binary)
# sum(dfMVRS$AD_binary=="AD Positive" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="M")
# sum(dfMVRS$AD_binary=="AD Positive" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="F")
# sum(dfMVRS$AD_binary=="AD Negative" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="M")
# sum(dfMVRS$AD_binary=="AD Negative" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Sex=="F")

	#BONFERRONI CALCULATION FOR FISHER TESTS USED FOR ASYN/AD RELATIONSHIP IN OVERALL DATASET AND YO DATASET 
	# 0.05/2 #0.025

#TXT: RTQUIC CHISQUARE WITH AD
# Compare proportions between AD+ and AD- for RT positivity
	# table(dfMVRS$AD_binary, dfMVRS$RTQUIC)
	# chisq.test(table(dfMVRS$AD_binary, dfMVRS$RTQUIC), correct=T) #correct: Yates'continuity correction
# fisher.test(table(dfMVRS$AD_binary, dfMVRS$RTQUIC)) # Expected count is <5 for one cell. Could use one-way but kind of weak to do so. (alternative="greater")
	# cramerV(table(dfMVRS$AD_binary, dfMVRS$RTQUIC))
# phi(table(dfMVRS$AD_binary, dfMVRS$RTQUIC), digits=6)


#TXT: RTQUIC CHISQUARE WITH AD IN YOD ONLY
#Analysis of AD*RTQUIC association in the young cohort
# dfMVRSyoung <- subset(dfMVRS, 65>=Onset)
# nrow(dfMVRSyoung)
# fisher.test(table(dfMVRSyoung$AD_binary, dfMVRSyoung$RTQUIC)) # Expected count is <5 for one cell
	# cramerV(table(dfMVRSyoung$AD_binary, dfMVRSyoung$RTQUIC))
# phi(table(dfMVRSyoung$AD_binary, dfMVRSyoung$RTQUIC), digits=6)
	# dfMVRSyoung$Onset

# #TXT: RTQUIC CHISQUARE WITH AD IN LOD ONLY
# dfMVRSold <- subset(dfMVRS, Onset>65)
# nrow(dfMVRSold)
# fisher.test(table(dfMVRSold$AD_binary, dfMVRSold$RTQUIC)) # Expected count is <5 for one cell
	# cramerV(table(dfMVRSold$AD_binary, dfMVRSold$RTQUIC))
# phi(table(dfMVRSold$AD_binary, dfMVRSold$RTQUIC), digits=6)
	# dfMVRSold$Onset


#TXT: PROPORTION OF AD- RTQUIC+ SUBJECTS in YOD and LOD
# sum(dfMVRS$AD_binary=="AD Positive" & dfMVRS$RTQUIC=="SAA positive" & 65>=dfMVRS$Onset)
# sum(dfMVRS$RTQUIC=="SAA positive" & 65>=dfMVRS$Onset)
# sum(dfMVRS$AD_binary=="AD Positive" & dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Onset>65)
# sum(dfMVRS$RTQUIC=="SAA positive" & dfMVRS$Onset>65)



# DESCRIPTION OF AD/RTQUIC RELATIONSHIP WITH ABETA42
#####################################################

#TXT: ABETA COMPARISON BETWEEN ASYN+ VS ASYN- 
#Use preprocessed dataset (normalized so all cont variables fall within range)
	# boxplot(logabeta ~ RTQUIC, data= RTposnormdfMVRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(logabeta ~ RTQUIC, data = RTposnormdfMVRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(logabeta ~ RTQUIC, data= RTnegnormdfMVRS, col = "white")$out #create vector with outlier values (from AD-negative group)
	# 	stripchart(logabeta ~ RTQUIC, data = RTnegnormdfMVRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# normdfMVRSabeta <- normdfMVRS #No positive outlier
	# nrow(normdfMVRSabeta)
	# head(normdfMVRSabeta)

	# RTposnormdfMVRSabeta<- subset(normdfMVRSabeta, RTQUIC=="SAA positive")
	# RTnegnormdfMVRSabeta <- subset(normdfMVRSabeta, RTQUIC=="SAA negative")
	# CBSnormdfMVRSabeta <- subset(normdfMVRSabeta, DX_APD=="CBS")
	# PSPnormdfMVRSabeta <- subset(normdfMVRSabeta, DX_APD=="PSP")

	# shapiro.test(RTposnormdfMVRSabeta$logabeta) #normal
	# shapiro.test(RTnegnormdfMVRSabeta$logabeta)  #normal
	# shapiro.test(CBSnormdfMVRSabeta$logabeta)  #normal
	# shapiro.test(PSPnormdfMVRSabeta$logabeta)  #normal
	# leveneTest(logabeta ~ DX_APD*RTQUIC, normdfMVRSabeta) #homo

	# Compare models with Ftest
	# test1 <- lm(logabeta ~ RTQUIC, normdfMVRSabeta)
	# test2 <- lm(logabeta ~ RTQUIC  + DX_APD, normdfMVRSabeta)
	# test3 <- lm(logabeta ~ RTQUIC + Sex, normdfMVRSabeta) #not adding any value to the model
	# anova(test1, test2) 
	# anova(test1, test3) 

	# Inclusion of Onset as a covariate:
	# Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(normdfMVRSabeta, x = "Onset", y = "logabeta", add = "reg.line")+ 
  	# 	stat_regline_equation(aes()) 
	# ggscatter(normdfMVRSabeta, x = "Onset", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
  	# 	stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(normdfMVRSabeta, x = "Onset", y = "logabeta", color = "RTQUIC", add = "reg.line")+
  	# 	stat_regline_equation(aes(color = RTQUIC))
	# cor.test(normdfMVRSabeta$logabeta, normdfMVRSabeta$Onset) #corr
	# summary(lm(normdfMVRSabeta$logabeta ~ normdfMVRSabeta$Onset)) #linear relationship 

	# Inclusion of NFL as a covariate:
	# Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(normdfMVRSabeta, x = "NFL", y = "logabeta", add = "reg.line")+ 
  	# 	stat_regline_equation(aes()) 
	# ggscatter(normdfMVRSabeta, x = "NFL", y = "logabeta", color = "DX_APD", add = "reg.line")+ 
  	# 	stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(normdfMVRSabeta, x = "NFL", y = "logabeta", color = "RTQUIC", add = "reg.line")+
  	# 	stat_regline_equation(aes(color = RTQUIC))
	# cor.test(normdfMVRSabeta$logabeta, normdfMVRSabeta$NFL) #corr
	# summary(lm(normdfMVRSabeta$logabeta ~ normdfMVRSabeta$NFL)) #linear relationship 

	# Inclusion of other variables as covariate:
	# summary(lm(normdfMVRSabeta$logabeta ~ normdfMVRSabeta$LP2_Disease_Duration)) #no linear relationship
	# summary(lm(normdfMVRSabeta$logabeta ~ normdfMVRSabeta$ttau)) #no linear relationship
	# summary(lm(normdfMVRSabeta$logabeta ~ normdfMVRSabeta$ptau)) #no linear relationship

#Model
# mlr <- lm(logabeta ~ Onset*RTQUIC + DX_APD + NFL, normdfMVRSabeta) 
# summary(mlr) 

	# Diagnostics of the model run
	# check_normality(mlr) #Ok
	# durbinWatsonTest(mlr) #Ok

#Simple slopes for onset: 
# emtrends(mlr,pairwise ~  RTQUIC, var="Onset")


#Main effects: 
# emmeans(mlr, ~ RTQUIC) #adjusted means: cannot be interpreted due to interaction (slopes crossing each other)
# emmeans(mlr, ~ DX_APD) #adjusted means. Ok because no interaction. 




###################################	BINARY LOGISTIC REGRESSION INCLUDING CLINICAL ###########################################
##############################################################################################################################

#TXT: BLR
normdfMVRS <- normdfMVRS %>%
        		mutate(RTQUIC = case_when(RTQUIC == "SAA positive" ~ 1,
                       RTQUIC == "SAA negative" ~ 0)) %>%
        		data.frame() #Convert to dataframe to facilitate operations.

	sum(normdfMVRS$RTQUIC ==1)

library(boot)
library(perm)

#Model is selected based on previous analyses (logabeta*Onset), simple comparisons in Supp material (Gait/RBD_binary).
#Model was run with and without AD/DX_APD status and all values remained very similar. 
blr <- glm(RTQUIC ~ DX_APD + logabeta*Onset + RBD_binary + Gait, data= normdfMVRS, family = "binomial")
summary(blr)
	
	#model diagnostics
	AIC(blr)
	pscl::pR2(blr)["McFadden"]
	caret::varImp(blr)	

	# car::vif(blr)
	# 	blronset <- glm(RTQUIC ~ DX_APD + Onset + RBD_binary + Gait, data= normdfMVRS, family = "binomial")
	# 		car::vif(blronset)
	# 	blrabeta <- glm(RTQUIC ~ DX_APD + logabeta + RBD_binary + Gait, data= normdfMVRS, family = "binomial")
	# 		car::vif(blrabeta)



# # logodds <- blr$linear.predictors
# # boxTidwell(logodds ~ dfMVRS$DX_APD)
# # boxTidwell(logodds ~ dfMVRS$RBD_binary)
# # boxTidwell(logodds ~ dfMVRS$Gait)

# # residuals <- residuals(blr)
# # hist(residuals)
# # qqnorm(residuals)
# # shapiro.test(residuals)
# durbinWatsonTest(blr) #Check the residuals are independent (multiple regression assumption)
# AIC(blr) 
# durbinWatsonTest(blr2) #Check the residuals are independent (multiple regression assumption)
# AIC(blr2) 


# #6 people have RBD and are RTQUIC positive + 1 is excluded from model but is positive on autopsy (MV) + 1 is negative for RTQUIC

# exp(4.88278) #RBD. 81.60169. Odds of someone with RBD being RTquic+ increases by 8000%. IE x80. 
# exp(-2.85663) #Gait. 0.09253855. Odds of someone with gait issues being RTQUIC+ decreases by 10% compared to someone without. 
##### exp(0.15122) #Age: 1.16. #The value indicates that as age increase by one more unit, then the odds of being SAA+ increases by 16%

# ggscatter(dfMVRS, x = "Age", y = "logabeta", color = "RTQUIC", add = "reg.line")+
  # stat_regline_equation(aes(color = RTQUIC))
# summary(lm(logabeta ~ RTQUIC*Age, normdfMVRS))


###############################################	ASYN SAA POSITIVITY AND NFL ##################################################
##############################################################################################################################


#Visualize as we expect very positive skew
# boxplot(NFL ~ RTQUIC, data= dfMVRS, col = "white")$out
# stripchart(NFL ~ RTQUIC,data = df, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE) #There's clearly a difference

# boxplot <- boxplot(NFL ~ RTQUIC, data= RTposdfMVRS, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
# # boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
# threshold <- min(max(RTposdfMVRS$NFL,na.rm=T), as.numeric(quantile(CBSdfMVRS$RTposdfMVRS, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTposdfMVRS$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
# min(max(RTposdfMVRS$NFL,na.rm=T), as.numeric(quantile(RTposdfMVRS$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (RTposdfMVRS$NFL)*3))) 
# dfnfl <- subset(df, (NFL<threshold & AD=="Positive") | (AD=="Negative")) #make sure you don't accidentally erase the AD negative here, or filter their values
# # sort(dfnfl$NFL)

# # Remove outliers for DX
# boxplot <- boxplot(NFL ~ DX_APD, data= CBSdfMVRS, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
# 	stripchart(NFL ~ DX_APD, data = CBSdfMVRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# # boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
# threshold <- min(max(CBSdfMVRS$NFL,na.rm=T), as.numeric(quantile(CBSdfMVRS$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (CBSdfMVRS$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
# # dfRS[!((dfRS$NFL<threshold & dfRS$DX_APD=="CBS")| (dfRS$DX_APD=="PSP")), c("ID", "NFL")] 
# dfMVRSnflDX <- subset(dfMVRS, (NFL<threshold & DX_APD=="CBS") | (DX_APD=="PSP"))

# boxplot <- boxplot(NFL ~ DX_APD, data= PSPdfRS, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
# 	stripchart(NFL ~ DX_APD, data = PSPdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# # boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
# threshold <- min(max(PSPdfRS$NFL,na.rm=T), as.numeric(quantile(PSPdfRS$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (PSPdfRS$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
# # dfRSnflDX[!((dfRSnflDX$NFL<threshold & dfRSnflDX$DX_APD=="PSP")| (dfRSnflDX$DX_APD=="CBS")), c("ID", "NFL")] 
# dfRSnflDX <- subset(dfRSnflDX, (NFL<threshold & DX_APD=="PSP") | (DX_APD=="CBS")) 

# cor.test(dfMVRS$NFL, dfMVRS$abeta) #correlated
# cor.test(dfMVRS$NFL, dfMVRS$ptau)
# cor.test(dfMVRS$NFL, dfMVRS$ttau)

# wilcox.test(dfMVRS$NFL ~ dfMVRS$RTQUIC)
# dfMVRS %>% group_by(RTQUIC) %>% summarise(count=n(), median= format(round(median(NFL, na.rm=T),3),3), mean= format(round(mean(NFL, na.rm=T),3),3), sd= format(round(sd(NFL, na.rm=T),3),3)) #format is used for the decimal
# median= format(round(median(NFL, na.rm=T),3),3)
# # OLS?

# boxplot <- boxplot(NFL ~ AD, data= ADnegdf, col = "white") #Same for AD- group
# boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
# threshold<- min(max(ADnegdf$NFL,na.rm=T), as.numeric(quantile(ADnegdf$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (ADnegdf$NFL)*3))) 
# dfnfl <- subset(dfnfl, (NFL<threshold & AD=="Negative") | (AD=="Positive")) 
# sort(dfnfl$NFL)
# boxplot(NFL ~ AD, data= dfnfl, col = "white") 
# shapiro.test(dfnfl$logNFL) #normal for overall dataset
# ADposdf <- subset(dfnfl, AD=="Positive")
# ADnegdf <- subset(dfnfl, AD=="Negative")
# shapiro.test(ADposdf$logNFL) #normal
# shapiro.test(ADnegdf$logNFL) #normal

# aovmodel <- aov(logNFL ~ AD + DX_APD + Age + LP2_Disease_Duration, dfnfl) #No need for interaction given how few PSP-AD there are.
# summary(aovmodel)
# Anova(aovmodel, type='III') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# dfnfl %>% group_by(DX_APD) %>% summarise(count=n(), median= format(round(median(NFL, na.rm=T),3),3), mean= format(round(mean(NFL, na.rm=T),3),3), sd= format(round(sd(NFL, na.rm=T),3),3)) #format is used for the decimal
# CBSdfnfl <- subset(dfnfl, DX_APD=="CBS")
# aovmodel <- aov(logNFL ~ AD + Age + LP2_Disease_Duration, CBSdf) #No need for interaction given how few PSP-AD there are.
# summary(aovmodel)
# Anova(aovmodel, type='II') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# emmeans(aovmodel, ~ AD) #if interested in variables of the fitted aov model, ie SE instead of SD and marginal means
# #https://stats.stackexchange.com/questions/369532/interpreting-the-standard-error-from-emmeans-r
# #https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
# par(mfrow=c(2,2))
# plot(aovmodel)
# par(mfrow=c(1,1))

# posth = glht(model, linfct=mcp(factor='Tukey'))
# summary(posth)
# emmeans(model, ~ factor:covariate)
#perform comparisons
# TukeyHSD(aovmodel, conf.level=.95, which = 'DX_APD')
# TukeyHSD(aovmodel, conf.level=.95, which = 'AD')

############################	TABLE 2 ASYN POSITIVITY AND MOTOR/AUTONOMIC/PD SYMPTOMS   ####################################
#############################################################################################################################

# dfMVRS %>% count(Sex)
# dfMVRS %>% count(AD_binary)
# dfMVRS %>% count(APOEe4)

# dfMVRS %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
# dfMVRS %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
# dfMVRS %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# DO CBS FIRST: COUNTS + ANALYSES
#################################
# RTposCBSdfMVRS <- subset(CBSdfMVRS, RTQUIC=="SAA positive")
# RTnegCBSdfMVRS <- subset(CBSdfMVRS, RTQUIC=="SAA negative")

# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Sex)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(AD_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(APOEe4)

# CBSdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
# CBSdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
# CBSdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# STATISTICAL COMPARISONS FOR AGES AND SEX
# table(CBSdfMVRS$Sex, CBSdfMVRS$RTQUIC)
# chisq.test(table(CBSdfMVRS$Sex, CBSdfMVRS$RTQUIC), correct=F)

# table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell

# table(CBSdfMVRS$APOEe4, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$APOEe4, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell


# shapiro.test(RTposCBSdfMVRS$Age) #normal
# shapiro.test(RTnegCBSdfMVRS$Age) #normal
# var.test(Age ~ RTQUIC, data = CBSdfMVRS) #homo
# t.test(CBSdfMVRS$Age ~ CBSdfMVRS$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposCBSdfMVRS$Onset) #normal
# shapiro.test(RTnegCBSdfMVRS$Onset) #normal
# var.test(Onset ~ RTQUIC, data = CBSdfMVRS) #homo
# t.test(CBSdfMVRS$Onset ~ CBSdfMVRS$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposCBSdfMVRS$Park_onset) #normal
# shapiro.test(RTnegCBSdfMVRS$Park_onset) #normal
# var.test(Park_onset ~ RTQUIC, data = CBSdfMVRS) #homo
# # # wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
# t.test(CBSdfMVRS$Park_onset ~ CBSdfMVRS$RTQUIC, var.equal=TRUE) 

# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Tremor_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(RestTremor)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(LimbRigidity)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Slowness_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Falls_PI)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Gait)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(RBD_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Anosmia_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Constipation_binary)
### CBSdfMVRS %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Sexual_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Urinary_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Bowel_binary)
# CBSdfMVRS %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

# table(CBSdfMVRS$anyPPA, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$anyPPA, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Tremor_binary, CBSdfMVRS$RTQUIC)
# chisq.test(table(CBSdfMVRS$Tremor_binary, CBSdfMVRS$RTQUIC), correct=F)
# table(CBSdfMVRS$RestTremor, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$RestTremor, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$LimbRigidity, CBSdfMVRS$RTQUIC)
# chisq.test(table(CBSdfMVRS$LimbRigidity, CBSdfMVRS$RTQUIC), correct=F)
# table(CBSdfMVRS$Slowness_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Slowness_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Falls_PI, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Falls_PI, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Gait, CBSdfMVRS$RTQUIC)
# chisq.test(table(CBSdfMVRS$Gait, CBSdfMVRS$RTQUIC), correct=F)
# table(CBSdfMVRS$RBD_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$RBD_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Lifetime_Dopa_responder_true, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Lifetime_Dopa_responder_true, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Lifetime_VisualHallucinations_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Lifetime_VisualHallucinations_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Constipation_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Constipation_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Urinary_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Urinary_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(CBSdfMVRS$Bowel_binary, CBSdfMVRS$RTQUIC)
# fisher.test(table(CBSdfMVRS$Bowel_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell


# DO PSP SECOND: COUNTS + ANALYSES
#################################
# RTposPSPdfMVRS <- subset(PSPdfMVRS, RTQUIC=="SAA positive")
# RTnegPSPdfMVRS <- subset(PSPdfMVRS, RTQUIC=="SAA negative")

# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Sex)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(AD_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(APOEe4)
# PSPdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
# PSPdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
# PSPdfMVRS %>% group_by(RTQUIC) %>% summarize(format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))

# STATISTICAL COMPARISONS FOR AGES AND SEX
# table(PSPdfMVRS$Sex, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Sex, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell

# table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell

# table(PSPdfMVRS$APOEe4, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$APOEe4, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
#
# shapiro.test(RTposPSPdfMVRS$Age) #normal
# shapiro.test(RTnegPSPdfMVRS$Age) #normal
# var.test(Age ~ RTQUIC, data = PSPdfMVRS) #homo
# t.test(PSPdfMVRS$Age ~ PSPdfMVRS$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposPSPdfMVRS$Onset) #normal
# shapiro.test(RTnegPSPdfMVRS$Onset) #normal
# var.test(Onset ~ RTQUIC, data = PSPdfMVRS) #homo
# t.test(PSPdfMVRS$Onset ~ PSPdfMVRS$RTQUIC, var.equal=TRUE) 

# shapiro.test(RTposPSPdfMVRS$Park_onset) #normal
# shapiro.test(RTnegPSPdfMVRS$Park_onset) #normal
# var.test(Park_onset ~ RTQUIC, data = PSPdfMVRS) #homo
# # # wilcox.test(df$Park_onset ~ df$DX_APD, paired=F) #Cannot be computed due to ties.
# t.test(PSPdfMVRS$Park_onset ~ PSPdfMVRS$RTQUIC, var.equal=TRUE) 


# PSPdfMVRS %>% group_by(RTQUIC) %>% count(anyPPA)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Tremor_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(RestTremor)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(LimbRigidity)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Slowness_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Falls_PI)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Gait)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(RBD_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Lifetime_Dopa_responder_true)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Anosmia_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Lifetime_VisualHallucinations_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Constipation_binary)
# ## PSPdfMVRS %>% group_by(RTQUIC) %>% count(Light_binary) #Need to include but maybe not worth it anyway
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Sexual_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Orthostatism_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Urinary_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Bowel_binary)
# PSPdfMVRS %>% group_by(RTQUIC) %>% count(Thermoregulatory_binary)

# table(PSPdfMVRS$Tremor_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Tremor_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$LimbRigidity, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$LimbRigidity, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Slowness_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Slowness_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Falls_PI, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Falls_PI, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Gait, PSPdfMVRS$RTQUIC)
# chisq.test(table(PSPdfMVRS$Gait, PSPdfMVRS$RTQUIC), correct=F)
# table(PSPdfMVRS$RBD_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$RBD_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Lifetime_Dopa_responder_true, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Lifetime_Dopa_responder_true, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Lifetime_VisualHallucinations_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Lifetime_VisualHallucinations_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Constipation_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Constipation_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Orthostatism_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Orthostatism_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Urinary_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Urinary_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
# table(PSPdfMVRS$Bowel_binary, PSPdfMVRS$RTQUIC)
# fisher.test(table(PSPdfMVRS$Bowel_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell


###################################	SURVIVAL ANALYSES ON THE SEVERE MOTOR DISEASE ###########################################
##############################################################################################################################
# library(survival)
# library(lubridate)
# library(ggsurvfit)
# library(gtsummary)
# library(tidycmprsk)
# library(condsurv)

# dfMVRS$Severe_Motor_onset_status_12 <- as.numeric(dfMVRS$Severe_Motor_onset_status_12)
# CBSdfMVRS$Severe_Motor_onset_status_12 <- as.numeric(CBSdfMVRS$Severe_Motor_onset_status_12)
# Surv(dfMVRS$Severe_onset_diff, dfMVRS$Severe_Motor_onset_status_12==2)
# surv1 <- survfit(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ 1, data = dfMVRS)
# str(surv1)

# summary(survfit(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ 1, data = dfMVRS), times = 15) #average surv time for 5 years

# coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ RTQUIC, data = CBSdfMVRS)
# coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ DX_APD, data = dfMVRS)

# surv<-coxph(Surv(Severe_onset_diff, Severe_Motor_onset_status_12) ~ RTQUIC, data = dfMVRS)
# summary(surv)


####################################################	MISCELLANEOUS ########################################################
##############################################################################################################################

# #Exploratory: remove other influential points
# lev <- hat(model.matrix(mlr1)) #The leverage of an obs. measures its ability to move the regression model, ie the amount by which the predicted value would change if the obs. was shifted by one unit
# plot(lev)
# dfMVRStauratio[lev > mean(lev)*3, c("ID", "NFL")]
# cooks <- cooks.distance(mlr1)
# plot(cooks)
# 	abline(h = 4/nrow(dfMVRStauratio), lty = 2, col = "steelblue") # add cutoff line
# dfMVRStauratio[cooks > (4/nrow(dfMVRStauratio)), c("ID", "tauratio", "ptau", "ttau")]
# # vec <- dfMVRStauratio[cooks > (4/nrow(dfRSnflDX)), "ID"]
# as.numeric(names(cooks)[(cooks > (4/nrow(dfRSnfl)))])
# testdf <- subset(dfRSnflDX, ID!=vec[1] & ID!=vec[2] & ID!=vec[3] & ID!=vec[4] & ID!=vec[5]) 

# CBSdfMVRSptau <- subset(dfMVRSptau, DX_APD=="CBS")
# shapiro.test(RTposdfMVRSptau$logptau) #normal
# shapiro.test(RTnegdfMVRSptau$logptau) #normal
# # Run aov model
# # https://statisticsbyjim.com/anova/ancova/
# acovmodel <- aov(logptau ~ AD + DX_APD + RTQUIC + abeta, dfMVRSptau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(acovmodel)
# Anova(acovmodel, type="II") #Since it is an ancova, report III
# Anova(acovmodel, type="III")
# AIC(acovmodel)
# 	dfMVRSptau %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	dfMVRSptau %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	#https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
# 	emmeans(acovmodel, ~ RTQUIC: abeta) #adjusted means
# 	emmeans(acovmodel, ~ DX_APD:abeta) #adjusted means
# acovmodel <- aov(logptau ~ AD + RTQUIC + abeta, CBSdfMVRSptau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(acovmodel)
# Anova(acovmodel, type="II") #Since it is an ancova, report III
# Anova(acovmodel, type="III")
# AIC(acovmodel)
# 	dfMVRSptau %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(ptau, na.rm=T), sd=sd(ptau, na.rm=T))
# 	emmeans(acovmodel, ~ RTQUIC:abeta) #adjusted means

# #Exploratory: remove other influential points
# lev <- hat(model.matrix(aovmodel)) #The leverage of an obs. measures its ability to move the regression model, ie the amount by which the predicted value would change if the obs. was shifted by one unit
# plot(lev)
# dfRSnflDX[lev > mean(lev)*3, c("ID", "NFL")]
# cooks <- cooks.distance(aovmodel)
# plot(cooks)
# 	abline(h = 4/nrow(dfRSnflDX), lty = 2, col = "steelblue") # add cutoff line
# dfRSnflDX[cooks > (4/nrow(dfRSnflDX)), c("ID", "NFL")]
# vec <- dfRSnflDX[cooks > (4/nrow(dfRSnflDX)), "ID"]
# # as.numeric(names(cooks)[(cooks > (4/nrow(dfRSnfl)))])
# testdf <- subset(dfRSnflDX, ID!=vec[1] & ID!=vec[2] & ID!=vec[3] & ID!=vec[4] & ID!=vec[5]) 
# #Identify and remove outliers in Ttau: here doing it directly on logged data to avoid unnecessary removal of data given it is normal enough that log procedure fixes issue. 
# boxplot(RTposdfMVRS$logttau ~ RTposdfMVRS$RTQUIC)$out #no outlier
# 	stripchart(logttau ~ RTQUIC, data = RTposdfMVRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# vec <- boxplot(RTnegdfMVRS$logttau ~ RTnegdfMVRS$RTQUIC)$out #2 highest values are outliers
# 	stripchart(logttau ~ RTQUIC, data = RTnegdfMVRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# dfMVRSttau <- subset(dfMVRS, logttau!=vec[1]) #remove outlier in logged data
# boxplot(dfMVRSptau$logttau ~ dfMVRSttau$RTQUIC)$out
# 	stripchart(logttau ~ RTQUIC, data = dfMVRSttau, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
# RTposdfMVRSttau <- subset(dfMVRSttau, RTQUIC=="SAA positive")
# RTnegdfMVRSttau <- subset(dfMVRSttau, RTQUIC=="SAA negative")
# shapiro.test(RTposdfMVRSttau$logttau) #normal
# shapiro.test(RTnegdfMVRSttau$logttau) #normal
# #Run aov model
# aovmodel <- aov(logttau ~ AD + DX_APD + RTQUIC + Age + abeta + LP2_Disease_Duration, dfMVRSttau) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# summary(aovmodel)
# AIC(aovmodel)
# Anova(aovmodel, type='II') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# dfMVRSttau %>% group_by(RTQUIC) %>% summarize(count=n(), mean= format(round(mean(ttau, na.rm=T),3),3), sd= format(round(sd(ttau, na.rm=T),3),3))
# dfMVRSttau %>% group_by(DX_APD) %>% summarize(count=n(), mean= format(round(mean(ttau, na.rm=T),3),3), sd= format(round(sd(ttau, na.rm=T),3),3))

# vec <- boxplot(dfMVRS$abeta ~ dfMVRS$RTQUIC)$out
# dfMVRSabeta <- subset(dfMVRS, abeta!=vec[1] & abeta!=vec[2]) #remove outlier in raw data
# boxplot(dfMVRSabeta$logabeta ~ dfMVRSabeta$RTQUIC)$out
# 	aovmodel <- aov(logabeta ~ AD + DX_APD + RTQUIC + Age + ptau + ttau, dfMVRSabeta) #No need for interaction given how few PSP-AD there are. ttau not included due to correlation.
# 	summary(aovmodel)
# 	Anova(aovmodel, type='II') #Compare the fitted model to a null model. No need for posthoc given these are 2-level factors only.
# 	dfMVRSabeta %>% group_by(RTQUIC) %>% summarize(count=n(), mean=mean(abeta, na.rm=T), sd=sd(abeta, na.rm=T))
# 	dfMVRSabeta %>% group_by(DX_APD) %>% summarize(count=n(), mean=mean(abeta, na.rm=T), sd=sd(abeta, na.rm=T))
# # wilcox.test(dfMVRS$abeta ~ dfMVRS$RTQUIC)
# #Ordinal logistic regression?
# boxplot(abeta ~ RTQUIC, data = dfMVRSabeta, col = "white")
# stripchart(abeta ~ RTQUIC,
#            data = dfMVRSabeta,
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
	# Data are not normally distributed so be careful with residuals
	# Using lm to compare Rsquare to that of test models below (there is no difference with aov)
	# mlr <- lm(Cognitive_Z ~ Lifetime_AD_binary + DX_APD + Age, dfRS)
	# mlr2 <- lm(Cognitive_Z ~Age*Lifetime_AD_binary + DX_APD, dfRS) #Technically best model but since interaction remains ns I think for sake of simplicity it is okay to select mlr instead. 
	# mlr3 <- lm(Cognitive_Z ~Age*DX_APD + Lifetime_AD_binary, dfRS)
	# summary(mlr) 
	# summary(mlr2) 
	# summary(mlr3) 

	# ##Diagnostics of the model run
	# residuals <- residuals(mlr)
	# hist(residuals)
	# qqnorm(residuals)
	# shapiro.test(residuals)
	# check_normality(mlr) #Very non normal: model not fittint too well the data
	# durbinWatsonTest(mlr) #Check the residuals are independent (multiple regression assumption)

	##Means to report for model mlr: report EMMs + SE 
	## dfRS %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Cognitive_Z, na.rm=T),3),3), sd=format(round(sd(Cognitive_Z, na.rm=T),3),3))
	# emmeans(mlr, ~ AD, opt.digits=T) #adjusted means
	## emmeans(mlr, ~ DX_APD, opt.digits=T) #adjusted means


# MULT LINEAR REGRESSION  TO COMPARE MOCA Z SCORES BETWEEN DX AND AD WITH DURATION AS COVARIATE
# #Here remove CBS-HIV as looking at DX
	# boxplot(MOCA_Z ~ DX_APD, data= CBSdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = CBSdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ DX_APD, data= PSPdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ DX_APD, data = PSPdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADposdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADposdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ Lifetime_AD_binary, data= ADnegdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ Lifetime_AD_binary, data = ADnegdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEposdfRS, col = "white")$out #identify outliers: in AD positive df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEposdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# boxplot(MOCA_Z ~ APOEe4, data= APOEnegdfRS, col = "white")$out #identify outliers: in AD negative df, none. 
	# 	stripchart(MOCA_Z ~ APOEe4, data = APOEnegdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)

	# #Then check normality of my group for AD and DX groups, and homoscedasticity
	# shapiro.test(CBSdfRS$MOCA_Z)
	# shapiro.test(PSPdfRS$MOCA_Z) #not normal
	# shapiro.test(ADposdfRS$MOCA_Z)
	# shapiro.test(ADnegdfRS$MOCA_Z) #not normal
	# shapiro.test(APOEposdfRS$MOCA_Z)
	# shapiro.test(APOEnegdfRS$MOCA_Z) #not normal
	# leveneTest(MOCA_Z ~ DX_APD*APOEe4*Lifetime_AD_binary, dfRS) #Homoscedasticity. 

	#Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(dfRS, x = "Age", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(dfRS, x = "Age", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(dfRS, x = "Age", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
	#   stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(dfRS, x = "Age", y = "MOCA_Z", color = "Lifetime_AD_binary", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
	#   stat_regline_equation(aes(color = Lifetime_AD_binary))
	# cor.test(dfRS$MOCA_Z, dfRS$Age)
	# summary(lm(dfRS$MOCA_Z ~ dfRS$Age)) #linear relationship
	##NOTE: possible interaction of Age with either DX or AD on Cognitive z-score (CBS/AD+ show much more steep differences in cognitive impairment
			#relative to demographics depending on age. IE young AD+ has much worse CI than a young AD-, but an old AD+ and an old AD- are comparable.
			#This points to AD+ subjects being much more affected at young age. Not sure if makes sense to include in this model). 

	#Inclusion of Duration as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(dfRS, x = "LP2_Disease_Duration", y = "MOCA_Z", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(dfRS, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "DX_APD", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes(color = DX_APD)) 
	# ggscatter(dfRS, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "APOEe4", add = "reg.line")+ #Okay
	#   stat_regline_equation(aes(color = APOEe4)) 
	# ggscatter(dfRS, x = "LP2_Disease_Duration", y = "MOCA_Z", color = "AD", add = "reg.line")+#Same as for DX, AD+ seem to be much worse when young
	#   stat_regline_equation(aes(color = AD))
	# cor.test(dfRS$MOCA_Z, dfRS$LP2_Disease_Duration) #No relationship

	##Compare models with Ftest
	# test1 <- lm(MOCA_Z ~ DX_APD, dfRS)
	# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, dfRS)
	# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD, dfRS) #both additions are beneficial so model should have AD and DX
	# anova(test1, test3)
	# anova(test2, test3)

	# Inclusion of APOE in model:
	# dfRSapoe <- dfRS[!is.na(dfRS$APOEe4), ]
	# test1 <- lm(MOCA_Z ~ APOEe4, dfRSapoe)
	# test2 <- lm(MOCA_Z ~ DX_APD, dfRSapoe)
	# test3 <- lm(MOCA_Z ~ DX_APD +APOEe4, dfRSapoe)
	# anova(test1, test3) 
	# anova(test2, test3) #DX by itself is better than with APOEe4

# Using lm to compare Rsquare to that of test models below (there is no difference with aov)
# mlr <- lm(MOCA_Z ~ Lifetime_AD_binary + DX_APD + Age, dfRS)
	# mlr2 <- lm(MOCA_Z ~Age*Lifetime_AD_binary + DX_APD, dfRS) #Technically best model but since interaction remains ns I think for sake of simplicity it is okay to select mlr instead. 
	# mlr3 <- lm(MOCA_Z ~Age*DX_APD + Lifetime_AD_binary, dfRS)
# summary(mlr) #type I (no interaction so order does not matter)
	# summary(mlr2) 
	# summary(mlr3) 

# # ##Diagnostics of the model run
# hist(residuals)
# qqnorm(residuals)
# shapiro.test(residuals(mlr))
# check_normality(mlr) #Very non normal: model not fittint too well the data
# durbinWatsonTest(mlr) #Check the residuals are independent (multiple regression assumption)

##Means to report for model mlr: report EMMs + SE 
## dfRS %>% group_by(AD) %>% summarize(count=n(), mean=format(round(mean(Cognitive_Z, na.rm=T),3),3), sd=format(round(sd(Cognitive_Z, na.rm=T),3),3))
# emmeans(mlr, ~ Lifetime_AD_binary, opt.digits=T) #adjusted means
## emmeans(mlr, ~ DX_APD, opt.digits=T) #adjusted means





	##### TABLES2 #####
	#Remove CBS-HIV patient
	
	# CBSdfRS %>% count(Lifetime_AD_binary)

	# #TABLE1: AGE
	# 	shapiro.test(ADposCBSdfRS$Age) #normal
	# 	shapiro.test(ADnegCBSdfRS$Age) #normal
	# 	var.test(Age ~ Lifetime_AD_binary, data = CBSdfRS) #homoscedasticity
	# t.test(CBSdfRS$Age ~ CBSdfRS$Lifetime_AD_binary, var.equal=TRUE)
	# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))
	# CBSdfRS %>% summarize(count=n(), format(round(mean(Age, na.rm=T),2),2), sd=sd(Age, na.rm=T))

	# #TABLE1:SEX
	# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% count(Sex)
	# table(CBSdfRS$Sex, CBSdfRS$Lifetime_AD_binary)
	# chisq.test(table(CBSdfRS$Sex, CBSdfRS$Lifetime_AD_binary), correct=F)


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			# shapiro.test(ADposCBSdfRS$Onset) #normal
			# shapiro.test(ADnegCBSdfRS$Onset) #normal
			# var.test(Onset ~ Lifetime_AD_binary, data = CBSdfRS) #homoscedasticity
		# t.test(CBSdfRS$Onset ~ CBSdfRS$Lifetime_AD_binary, var.equal=TRUE) 
		# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
		# CBSdfRS %>% summarize(count=n(), format(round(mean(Onset, na.rm=T),2),2), sd=sd(Onset, na.rm=T))
		
		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		# 	shapiro.test(ADposCBSdfRS$Park_onset) #normal
		# 	shapiro.test(ADnegCBSdfRS$Park_onset) #normal
		# 	var.test(Park_onset ~ Lifetime_AD_binary, data = CBSdfRS) #homoscedasticity
		# t.test(CBSdfRS$Park_onset ~ CBSdfRS$Lifetime_AD_binary, var.equal=TRUE) 
		# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))
		# CBSdfRS %>% summarize(count=n(), format(round(mean(Park_onset, na.rm=T),2),2), sd=sd(Park_onset, na.rm=T))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			# shapiro.test(ADposCBSdfRS$Cognitive_Z) #normal
			# shapiro.test(ADnegCBSdfRS$Cognitive_Z) #nonnormal
			# leveneTest(Cognitive_Z ~ Lifetime_AD_binary, data = CBSdfRS) #homoscedasticity
		# wilcox.test(CBSdfRS$Cognitive_Z ~ CBSdfRS$Lifetime_AD_binary, paired=F)
		# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(Cognitive_Z, na.rm=T),2),2), sd=sd(Cognitive_Z, na.rm=T))
		# CBSdfRS %>% summarize(count=n(), format(round(mean(Cognitive_Z, na.rm=T),2),2), sd=sd(Cognitive_Z, na.rm=T))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% count(APOEe4)
		# table(CBSdfRS$APOEe4, CBSdfRS$Lifetime_AD_binary)
		# fisher.test(table(CBSdfRS$APOEe4, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell


			###SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
	###COMPARISONS OF BIOMARKERS FOR AD:
	##################################################

			##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		##Look only in CBS: should be no interaction of APOEe4 but possible effect of AD
		# CBSdfRSonsetapoe <- CBSdfRSonset[!is.na(CBSdfRSonset$APOEe4), ]
			# test1 <- lm(Onset ~ Lifetime_AD_binary, CBSdfRSonsetapoe)
			# test2 <- lm(Onset ~ APOEe4, CBSdfRSonsetapoe)
			# test3 <- lm(Onset ~ APOEe4+ Lifetime_AD_binary, CBSdfRSonsetapoe)
			# anova(test1, test3)
			# anova(test2, test3) 
			# ADposCBSdfRSonset <- subset(CBSdfRSonset, Lifetime_AD_binary=="AD Positive")
			# ADnegCBSdfRSonset <- subset(CBSdfRSonset, Lifetime_AD_binary=="AD Negative")
			# shapiro.test(ADposCBSdfRSonset$Onset) #normal
			# shapiro.test(ADnegCBSdfRSonset$Onset) #normal
			# var.test(Onset ~ Lifetime_AD_binary, data = CBSdfRSonset) #homoscedasticity
		# t.test(CBSdfRSonset$Onset ~ CBSdfRSonset$Lifetime_AD_binary, var.equal=TRUE)
		# CBSdfRSonset %>% group_by(DX_APD) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))
		# CBSdfRSonset %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Onset, na.rm=T),3),3), sd=format(round(sd(Onset, na.rm=T),3),3))


		##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
		##Look only in CBS: should be no interaction of APOEe4 but possible effect of AD
		# CBSdfRSonsetapoe <- CBSdfRSonset[!is.na(CBSdfRSonset$APOEe4), ]
		# test1 <- lm(Park_onset ~ Lifetime_AD_binary, CBSdfRSonsetapoe)
		# test2 <- lm(Park_onset ~ APOEe4, CBSdfRSonsetapoe)
		# test3 <- lm(Park_onset ~ APOEe4+ Lifetime_AD_binary, CBSdfRSonsetapoe)
		# anova(test1, test3)
		# anova(test2, test3) 
		# ADposCBSdfRSonset <- subset(CBSdfRSonset, Lifetime_AD_binary=="AD Positive")
		# ADnegCBSdfRSonset <- subset(CBSdfRSonset, Lifetime_AD_binary=="AD Negative")
		# shapiro.test(ADposCBSdfRSonset$Park_onset) #normal
		# shapiro.test(ADnegCBSdfRSonset$Park_onset) #normal
		# var.test(Park_onset ~ Lifetime_AD_binary, data = CBSdfRSonset) #homoscedasticity
		# t.test(CBSdfRSonset$Park_onset ~ CBSdfRSonset$Lifetime_AD_binary, var.equal=TRUE)
		# CBSdfRSonset %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))
		# CBSdfRSonset %>% summarize(count=n(), mean=format(round(mean(Park_onset, na.rm=T),3),3), sd=format(round(sd(Park_onset, na.rm=T),3),3))


	##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
	# table(CBSdfRS$Parkinsonian_onset, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Parkinsonian_onset, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$Cognitive_onset, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Cognitive_onset, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$Language_onset, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Language_onset, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$anyPPA, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$anyPPA, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$Tremor_binary, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Tremor_binary, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$Myoclonus, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Myoclonus, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell
	# table(CBSdfRS$Language_onset, CBSdfRS$Lifetime_AD_binary)
	# fisher.test(table(CBSdfRS$Language_onset, CBSdfRS$Lifetime_AD_binary)) # Expected count is <5 for one cell



			##SUPPLEMENTARY MATERIAL: TABLE CBS AD+ VS AD-
			#Already got rid of outliers earlier. 
			# CBSdfRSapoe <- CBSdfRS[!is.na(CBSdfRS$APOEe4), ]
			# test1 <- lm(MOCA_Z ~ APOEe4, CBSdfRSapoe)
			# test2 <- lm(MOCA_Z ~ Lifetime_AD_binary, CBSdfRSapoe)
			# test3 <- lm(MOCA_Z ~ Lifetime_AD_binary +APOEe4, CBSdfRSapoe)
			# anova(test1, test3) 
			# anova(test2, test3) #DX by itself is better than with APOEe
			# shapiro.test(ADposCBSdfRS$MOCA_Z) #normal
			# shapiro.test(ADnegCBSdfRS$MOCA_Z) #normal
			# leveneTest(MOCA_Z ~ Lifetime_AD_binary, data = CBSdfRS) #homoscedasticity
			# t.test(CBSdfRS$MOCA_Z ~ CBSdfRS$Lifetime_AD_binary, var.equal=TRUE) 
			# CBSdfRS %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
			# CBSdfRS %>% summarize(count=n(), format(round(mean(MOCA_Z, na.rm=T),2),2), sd=sd(MOCA_Z, na.rm=T))
	

	# #MULT LINEAR REGRESSION TO COMPARE NFL LEVELS
	# #For NFL, another approach is chosen, as the data are expected to be very right-skewed. Threshold for outlier is changed.

	# Remove outliers for AD
	#First identify and remove outlier from the AD+ group only. 
	# boxplot <- boxplot(NFL ~ Lifetime_AD_binary, data= ADposCBSdfRS, col = "white") #For NFL, chose a more tolerant threshold for outliers. Q3+IQR*3 threshold instead of IQR*1.5
	# 	stripchart(NFL ~ Lifetime_AD_binary, data = ADposCBSdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# #boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
	# threshold <- min(max(ADposCBSdfRS$NFL,na.rm=T), as.numeric(quantile(ADposCBSdfRS$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (ADposCBSdfRS$NFL)*3))) #reports the value Q3+ IQR*3 (3 is very tolerant threshold)
	# #dfRS[!((dfRS$NFL<threshold & dfRS$Lifetime_AD_binary=="AD Positive")| (dfRS$Lifetime_AD_binary=="AD Negative")), c("ID", "NFL")] 
	# CBSdfRSnfl <- subset(CBSdfRS, (NFL<threshold & Lifetime_AD_binary=="AD Positive") | (Lifetime_AD_binary=="AD Negative")) #make sure you don't accidentally erase the AD negative here, or filter their values

	# ##Then identify and remove outlier from the AD- group only
	# boxplot <- boxplot(NFL ~ Lifetime_AD_binary, data= ADnegCBSdfRS, col = "white") #Same for AD- group
	# 	stripchart(NFL ~ Lifetime_AD_binary, data = ADnegCBSdfRS, method = "jitter", pch = 19, col = 2:4, vertical = TRUE, add = TRUE)
	# #boxplot$stats #[1,] lower whisker, [3,] median, [5,] upper whisker
	# threshold<- min(max(ADnegCBSdfRS$NFL,na.rm=T), as.numeric(quantile(ADnegCBSdfRS$NFL, 0.75, na.rm=T)) + (IQR(na.rm=T, (ADnegCBSdfRS$NFL)*3))) 
	# #dfRSnflAD[!((dfRSnflAD$NFL<threshold & dfRSnflAD$Lifetime_AD_binary=="AD Negative")| (dfRSnflAD$Lifetime_AD_binary=="AD Positive")), c("ID", "NFL")] 
	# CBSdfRSnfl <- subset(CBSdfRSnfl, (NFL<threshold & Lifetime_AD_binary=="AD Negative") | (Lifetime_AD_binary=="AD Positive")) 


	#Compare models withand without AD after removing the AD outliers with the Ftest
	# test1 <- lm(logNFL ~ Lifetime_AD_binary, CBSdfRSnfl)
	# test2 <- lm(logNFL ~ Sex, CBSdfRSnfl)
	# test3 <- lm(logNFL ~ Lifetime_AD_binary + Sex, CBSdfRSnfl)
	# anova(test1, test3)
	# anova(test2, test3)

	# #Oneway seems better but let's perform F-test to have full picture. 
	# testdf <- CBSdfRSnfl %>% drop_na("APOEe4") 
	# testdf$APOEe4
	# oneway <- lm(logNFL ~ Lifetime_AD_binary, testdf) 
	# twoway <- lm(logNFL ~ APOEe4 + Lifetime_AD_binary, testdf) 
	# anova(oneway, twoway)

	# #Then check normality in both AD+ and AD- and overall
	# ADposCBSdfRSnfl <- subset(CBSdfRSnfl, Lifetime_AD_binary=="AD Positive")
	# ADnegCBSdfRSnfl <- subset(CBSdfRSnfl, Lifetime_AD_binary=="AD Negative")
	# shapiro.test(ADposCBSdfRSnfl$logNFL) #normal
	# shapiro.test(ADnegCBSdfRSnfl$logNFL) #nonnormal
	# leveneTest(logNFL ~ Lifetime_AD_binary, CBSdfRSnfl) #Heterodasticity

	#Inclusion of Age as a covariate:
	#Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdfRSnfl, x = "Age", y = "NFL", add = "reg.line")+ #Not related to age
  	# 	stat_regline_equation(aes()) 
	# ggscatter(CBSdfRSnfl, x = "Age", y = "NFL", color = "Lifetime_AD_binary", add = "reg.line")+
  	# 	stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdfRSnfl$NFL, CBSdfRSnfl$Age) #no corr
	# cor.test(CBSdfRSnfl$logNFL, CBSdfRSnfl$Age) #no corr

	# #Inclusion of Duration as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdfRSnfl, x = "LP2_Disease_Duration", y = "NFL", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
	#   stat_regline_equation(aes()) 
	# ggscatter(CBSdfRSnfl, x = "LP2_Disease_Duration", y = "NFL", color = "Lifetime_AD_binary", add = "reg.line")+ #Relative to their demographics, old CBS tend to fare less bad than young CBS, which is not seen in PSP. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdfRSnfl$NFL, CBSdfRSnfl$LP2_Disease_Duration) #no corr
	# cor.test(CBSdfRSnfl$logNFL, CBSdfRSnfl$LP2_Disease_Duration) #no corr

	# #Inclusion of Abeta as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "abeta", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "abeta", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdfRSnfl$abeta, CBSdfRSnfl$NFL)
	# cor.test(CBSdfRSnfl$abeta, CBSdfRSnfl$logNFL)
	# cor.test(CBSdfRSnfl$logabeta, CBSdfRSnfl$NFL)
	# summary(lm(CBSdfRSnfl$abeta ~ CBSdfRSnfl$logNFL)) #no linear relationship

	# #Inclusion of ptau as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "ptau", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "ptau", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdfRSnfl$ptau, CBSdfRSnfl$NFL)
	# cor.test(CBSdfRSnfl$ptau, CBSdfRSnfl$logNFL)
	# summary(lm(CBSdfRSnfl$ptau ~ CBSdfRSnfl$logNFL)) #no linear relationship

	# #Inclusion of ptau as a covariate:
	# #Check independence of DV and covariate (not expected in observational data) + homogeneity of regression slopes. 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "ttau", add = "reg.line")+ #For overall group, NFL not related to age
  	# stat_regline_equation(aes()) 
	# ggscatter(CBSdfRSnfl, x = "NFL", y = "ttau", color = "Lifetime_AD_binary", add = "reg.line")+ #There seems to be interaction. Possibly PSP have lower levels in older group. CBS slight increase as expected with age. 
  	# stat_regline_equation(aes(color = Lifetime_AD_binary)) 
	# cor.test(CBSdfRSnfl$ttau, CBSdfRSnfl$NFL)
	# cor.test(CBSdfRSnfl$ttau, CBSdfRSnfl$logNFL)
	# summary(lm(CBSdfRSnfl$ttau ~ CBSdfRSnfl$logNFL)) #no linear relationship

	# # Data are not normally distributed so be careful with residuals
	# #Many options were tried: most were eliminated as model significance way >0.10. Interactions including aforementioned covariates were tried just in case as visualization was ambiguous. 
	# t.test(ADposCBSdfRSnfl$logNFL, ADnegCBSdfRSnfl$logNFL) #Welch's t-test: does not assume variance equal
	# 	CBSdfRSnfl %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), mean=format(round(mean(NFL, na.rm=T),3),3), sd=format(round(sd(NFL, na.rm=T),3),3))
	# 	CBSdfRSnfl %>% summarize(count=n(), mean=format(round(mean(NFL, na.rm=T),3),3), sd=format(round(sd(NFL, na.rm=T),3),3))
	# CBSdfRSnfl %>% group_by(Lifetime_AD_binary) %>% summarize(count=n(), format(round(median(NFL, na.rm=T),2),2), IQR=format(round(IQR(NFL, na.rm=T),2),2), min=min(NFL, na.rm=T), max=max(NFL, na.rm=T))
	# CBSdfRSnfl %>% summarize(count=n(), format(round(median(NFL, na.rm=T),2),2), IQR=IQR(NFL, na.rm=T), min=min(NFL, na.rm=T), max=max(NFL, na.rm=T))




		#TXT: RTQUIC CHISQUARE WITH AD IN CBS ONLY
	# #Compare proportions between AD+ and AD- for RT positivity but in CBS only: Close to significance. 
	# 	table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC)
	# 	chisq.test(table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC)) #correct: Yates'continuity correction
	# fisher.test(table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC)) # Expected count is <5 for one cell
	# 	# cramerV(table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC))
	# phi(table(CBSdfMVRS$AD_binary, CBSdfMVRS$RTQUIC), digits=6)

	#TXT: RTQUIC CHISQUARE WITH AD IN PSP ONLY
	# #Compare proportions between AD+ and AD- for RT positivity but in CBS only
	# 	chisq.test(table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC), correct=T) #correct: Yates'continuity correction
	# fisher.test(table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC)) # Expected count is <5 for one cell
	# 	table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC)
	# cramerV(table(PSPdfMVRS$AD_binary, PSPdfMVRS$RTQUIC))

#########################################################################################################
										   #DOCUMENTATION AND REFERENCES
###############################################################################################################################


#https://study.com/skill/learn/how-to-calculate-expected-counts-for-the-chi-square-test-for-goodness-of-fit-explanation.html
#https://stats.stackexchange.com/questions/298/in-linear-regression-when-is-it-appropriate-to-use-the-log-of-an-independent-va
#https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
##https://rdrr.io/cran/car/man/Anova.html
#https://study.com/skill/learn/how-to-calculate-expected-counts-for-the-chi-square-test-for-goodness-of-fit-explanation.html

#ANCOVA ASSUMPTIONS
# https://www.researchgate.net/post/How_can_continue_ANCOVA_when_assumption_of_homogeneity_of_regression_slopes_is_violated
#https://www.theanalysisfactor.com/ancova-assumptions-when-slopes-are-unequal/ #when slopes are unequal
# https://www.theanalysisfactor.com/assumptions-of-ancova/ #On the independence of vars and understanding criticism

##REALLY GOOD EXPLANATION OF WHY AOV == LM USING BINARY FACTORIAL AOV AS EXAMPLEA
#https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Learning_Statistics_with_R_-_A_tutorial_for_Psychology_Students_and_other_Beginners_(Navarro)/16%3A_Factorial_ANOVA/16.06%3A_ANOVA_As_a_Linear_Model
##We use aov() when we would like to fit an ANOVA model but the mode of distribution of SS to different effects is type I, which
##cannot be used by car::Anova(). 
##Base R anova() and car:Anova() differ for: contrast option already set as c("contr.sum", "contr.poly") for car::Anova() +
##car::Anova() does type II SS method as default, whereas base R anova() does type I SS method (same results as summary(aov))
## We also use anova() when we would like to compare the fit of nested regression models to determine if a regression model with a
##certain set of coefficients offers a significantly better fit than a model with only a subset of the coefficients.
#https://stats.stackexchange.com/questions/20452/how-to-interpret-type-i-type-ii-and-type-iii-anova-and-manova/20455#20455
#https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
# https://utstat.utoronto.ca/reid/sta442f/2009/typeSS.pdf# https://bookdown.org/steve_midway/DAR/understanding-anova-in-r.html#introduction-2
#https://md.psych.bio.uni-goettingen.de/mv/unit/lm_cat/lm_cat_unbal_ss_explained.html
#https://seriousstats.wordpress.com/2020/05/13/type-ii-and-type-iii-sums-of-squares-what-should-i-choose/
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.utstat.utoronto.ca/reid/sta442f/2009/typeSS.pdf
#https://towardsdatascience.com/anovas-three-types-of-estimating-sums-of-squares-don-t-make-the-wrong-choice-91107c77a27a
#http://mypage.concordia.ca/faculty/pperesne/BIOL_422_680/lecture-10-types-of-sum-of-squares.html
#http://dwoll.de/rexrepos/posts/ancova.html
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://web.pdx.edu/~newsomj/mvclass/ho_ancova.pdf

#MULTIPLE LINEAR REGRESSION
#https://www.statology.org/durbin-watson-test-r/ #DW test for independence of residuals
# https://www.statology.org/multiple-linear-regression-assumptions/
#https://www.statology.org/how-to-report-regression-results/ #Report mlr
# https://stats.stackexchange.com/questions/3200/is-adjusting-p-values-in-a-multiple-regression-for-multiple-comparisons-a-good-i]
