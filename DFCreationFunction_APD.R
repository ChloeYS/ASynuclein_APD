# FILENAME: DFCreationFunction_APD.R

#Last updated 16 May 2024
##var.func.1 creates the APD dataset that is used in the manuscript. It performs basic QC and formatting of the data. 
## The function is shared for informative purposes regarding format of the data/later analyses.
## The dataframe resulting of the var.func.1() call was saved as a cvs, which is the document that would be shared upon request. 
## The original datasheet cannot be shared due to containing potentially sensitive information. 



###############################################################################################################################
#SOURCE PACKAGES AND FUNCTIONS
###############################################################################################################################

source('Functions.R') #Source the functions from file in same repository (for now)
	#tidyverse() loaded in Functions.R


library(lubridate) #for is.Date()


###############################################################################################################################
#SUBJECT INCLUSION
###############################################################################################################################

var.func.1 <- function(df) {

#Make sure only subjects to be included in the paper are here. 
  df <- subset(df, as.numeric(Include_APDs)==1)


###############################################################################################################################
#CHECK DATES ARE IN PROPER FORMAT
###############################################################################################################################

 #Make sure that date variables are indeed dates. Usually issues arise in DOB values due to format. Reformatting original data in Excel may be required.

date.vars <- df[, c("DOB_dd.mmmm.yy", "RTQUIC_2_Date_dd.mmmm.yy", "First_Visit_Date_dd.mmmm.yy", "Last_Visit_Date_dd.mmmm.yy")]
    if(sum(sapply(date.vars, is.Date)) != 4) {
      df$DOB <- as.Date(as.character(df$DOB_dd.mmmm.yy), format="%Y-%m-%d")
      df$Date <- as.Date(as.character(df$RTQUIC_2_Date_dd.mmmm.yy), format="%Y-%m-%d") #Dates used for clinical tp
      df$First_Visit <- as.Date(as.character(df$First_Visit_Date_dd.mmmm.yy), format="%Y-%m-%d")
      df$Last_Visit <- as.Date(as.character(df$Last_Visit_Date_dd.mmmm.yy), format="%Y-%m-%d")
      date.vars <- df[, c("DOB", "Date", "First_Visit", "Last_Visit")]
      sapply(date.vars, is.Date) 
    }


###############################################################################################################################
#CHECK NUMERIC VALUES ARE IN THE PROPER FORMAT
###############################################################################################################################

# Make sure that numeric variables are indeed numeric. Usually issues arise in NFL values due to format. Reformatting original data in Excel may be required.

df$ID_Age <- as.numeric(df$ID_Age) #Manually change so doesnt need to be included in code below (will use it to compare to ID_Age_TXT after to make sure all is good)
num.vars <- c("Education", "ID_Age_TXT", "Onset_age", "Park_onset", "ptau_2",
				"ttau_2", "abeta_2", "ATI_2", "NFL_2", "LP2_Cognitive_Z.score",
				"LP2_MOCA_Z.score", "LP2_MOCA_total", "Lag_hours", "ThTmax") #Doesnt include Age so i can compare it later. 
df.num <- df[, num.vars] #dataframe with the numerical variables only
df.num <- sapply(df.num, as.numeric) #Apply as.numeric to all the values that should be numeric. This introduces NAs. 

cols <- names(df) %in% num.vars #Find the names of the columns that are numerical. "Cols" is a class logical: it says whether or not the name of each df column is within the num.vars vector. 
df.nonnum <- df[!cols] #dataframe without the numerical variables. 

df <- cbind(df.nonnum, df.num)

cat("as.numeric() call:", "\n",
	"Check potential issues when using cbind:", "\n",
	"1- duplicate column names;",  "\n",
	"2- different row numbers, etc.", "\n")

# SANITY CHECK
print(nrow(df.nonnum))
print(nrow(df.num))
print(nrow(df))

print(names(df))


###############################################################################################################################
#CREATE NEW VARIABLES: AGES AND DURATIONS
###############################################################################################################################

df <- df %>%

#"Start with "Age", for which defensive coding is implemented below
mutate(Age= (as.numeric(Date - DOB))/365) %>%
mutate(Age_Calculation_SanityCheck= Age - ID_Age_TXT) %>%
mutate(Age_cbind_SanityCheck= Age - ID_Age) %>%

mutate(Followup_duration= (as.numeric(Last_Visit- First_Visit))/365) %>%  #Followup duration for tables of demographic

mutate(LP2_Disease_Duration= Age - Onset_age) %>% #Disease duration by the time of LP
mutate(LP2_Park_duration= Park_onset - Onset_age) %>% #Disease duration by the time of LP

#To know how long the subject was followed after onset of disease: can be useful for survival analysis (although rn already accountd for in form of Severe_Motor_Status variable)
mutate(Age_at_last_visit= (as.numeric(Last_Visit- DOB))/365) %>%  #Age at the last visit, in order to determine how long each person was followed in total, for the survival analysis
mutate(Duration_last_visit= (Age_at_last_visit - Onset_age)) %>%  #Duration between disease onset and last age checked up on

data.frame() #Convert to dataframe to facilitate operations

# DEFENSE

##Below is for the comparison of the age calculated from DOB vs age entered in excel spreadsheet after dataframe subsetting + cbind call above
vec<- rep(1.1, nrow(df))
vec2<- rep(-1.1, nrow(df))
test <- df[df$Age_Calculation_SanityCheck>vec | df$Age_Calculation_SanityCheck<vec2,]
test2 <- df[df$Age_cbind_SanityCheck>vec | df$Age_cbind_SanityCheck<vec2,]
if (nrow(test)>0 | nrow(test2)>0) {
	print(test)
	print(test2)
    print("Check ages carefully as there seems to be a discrepancy")
    }


###############################################################################################################################
#CHECK CATEGORICAL VALUES ARE IN THE PROPER FORMAT        
###############################################################################################################################

df$ID.factor <- as.factor(df$ID) #Manually change so doesnt need to be included in code below (will use it to compare to ID after to make sure all is good)
cat.vars <- c("Sex", "PPA", "Onset_type", "DX_APD", "Race_ethnicity",
			  "RTQUIC_lifetime", "AD_lifetime_ATHENA", "LP2_auton_signs", "LP2_RBD_plus_other", "LP2_Anosmia",
			  "LP2_Sexual_dysfunction", "LP2_Constipation", "LP2_Urinary", "LP2_Bowel_Incontinence","APOEe4_alleles",
			  "LP2_gait", "LP2_falls_PI","LP2_retropulsion","LP2_tremor","LP2_slowness",
			  "LP2_oculomotor", "Lifetime_oculomotor", "LP2_rigidity", "LP2_apraxia","Lifetime_apraxia",
			  "Lifetime_Dopa","Lifetime_Dopa_responder_true")

df.cat <- df[, cat.vars] #dataframe with the categorical variables only
df.cat <- sapply(df.cat, as.factor) #Apply as.numeric to all the values that should be numeric. This introduces NAs. 

cols <- names(df) %in% cat.vars #Find the names of the columns that are numerical. "Cols" is a class logical: it says whether or not the name of each df column is within the num.vars vector. 
df.noncat <- df[!cols] #dataframe without the numerical variables. 

df <- cbind(df.noncat, df.cat)

cat("as.factor() call:", "\n",
	"Check potential issues when using cbind:", "\n",
	"1- duplicate column names;",  "\n",
	"2- different row numbers, etc.", "\n")

# SANITY CHECK
print(nrow(df.noncat))
print(nrow(df.cat))
print(nrow(df))

print(names(df))


# DEFENSE
if (sum(df$ID.factor == df$ID) != 67) {
	n <- sum(df$ID.factor != df$ID)
	cat("WARNING: Possible error when re-merging dataframe after as.numeric or as.factor call. There are at least: ",
		n, " potential errors when merging using cbind(). Consider other sources of error too.", "\n")
	}


###############################################################################################################################
# CREATE NEW VARIABLES: BINARY STRINGS
##WARNING: BE VERY CAREFUL WITH NAs OR NON-BINARY VARS. IT WILL BINARIZE THE DATA AND LUMP ANY NA WITH THE NON-REF LEVEL
###############################################################################################################################

      ##Demographics
      df <-string.var(df, "Onset_type", "Parkinsonism", "Parkinsonian_onset")
      df <-string.var(df, "Onset_type", "Cognitive", "Cognitive_onset")
      df <-string.var(df, "Onset_type", "Language", "Language_onset")

      ##Biomarkers
      df <-string.var(df, "AD_lifetime_ATHENA", "Positive", "AD") 

      ##Autonomic or PD-related non motor
      df <-string.var(df, "LP2_RBD_plus_other", "Yes", "RBD_binary")
      df <-string.var(df, "LP2_Anosmia", "Yes", "Anosmia_binary") 
      df <-string.var(df, "LP2_Sexual_dysfunction", "Yes", "Sexual_binary") 
      df <-string.var(df, "LP2_Constipation", "Yes", "Constipation_binary") 
      df <-string.var(df, "LP2_Urinary", "Yes", "Urinary_binary") 
      df <-string.var(df, "LP2_Bowel_Incontinence", "Yes", "Bowel_binary") 

      ##Neurological/psych history
      df <-string.var(df, "Lifetime_Hallucinations", "Yes", "Lifetime_Hallucinations_binary")            
      df <-string.var(df, "Lifetime_Hallucinations", "Visual", "Lifetime_VisualHallucinations_binary")

      ##Medication
      df <-string.var(df, "Lifetime_Dopa", "Yes", "Lifetime_Dopa_binary")

      #Motor
      df <-string.var(df, "LP2_tremor", "Yes", "Tremor_binary")
      df <-string.var(df, "LP2_tremor", "action", "ActionTremor")
      df <-string.var(df, "LP2_tremor", "rest", "RestTremor")

      df <-string.var(df, "LP2_slowness", "Yes", "Slowness_binary")

      df <-string.var(df, "LP2_oculomotor", "Yes", "OM_binary")
      df <-string.var(df, "LP2_oculomotor", "vertical", "VerticalOM")
      df <-string.var(df, "Lifetime_oculomotor", "vertical", "Lifetime_VerticalOM")

      df <-string.var(df, "LP2_rigidity", "Yes", "Rigidity_binary")
      df <-string.var(df, "LP2_rigidity", "axial", "AxialRigidity")
      df <-string.var(df, "LP2_rigidity", "limb", "LimbRigidity")


      df <-string.var(df, "Lifetime_apraxia", "Yes", "Lifetime_apraxia_binary")


# ###############################################################################################################################
# # CREATE NEW VARIABLES: COMBINE VARIABLES BASED ON STRINGS
# ##############################################################################################################################

# df <- df %>%

#           ##Onset_type_simplified: multiple domains instead of "Cognitive, Language" etc
#           mutate(Onset_type_simplified= Onset_type) %>%
#           mutate(Onset_type_simplified= case_when(grepl(",", Onset_type_simplified) ~ "Multiple domains",
#                                                   TRUE ~ as.character(Onset_type_simplified))) %>% #TRUE in this case captures all the values excluded from previous conditions,
#                                                                                                     #in which case case_when applies the new value, which is actually the original value.
          
#           ##PPA: whether a PPA was diagnosed or not
#           mutate(anyPPA= PPA) %>%
#           mutate(anyPPA = case_when(anyPPA == "svPPA" ~ "Yes",
#                                     anyPPA == "PNFA" ~ "Yes",
#                                     anyPPA == "lvPPA" ~ "Yes",
#                                     anyPPA == "PPA" ~ "Yes",
#                                     TRUE ~ as.character(anyPPA))) %>%

#           ##DXRTQUIC: combines RTQUIC and DX_APD strings
#           mutate(DXRTQUIC = as.factor(paste(DX_APD, RTQUIC, sep="_"))) %>%

#           data.frame() #Convert to dataframe to facilitate operations.

# ###############################################################################################################################
# # CREATE NEW VARIABLES: CREATE RATIO AND LOGS
# ##############################################################################################################################
# #http://www.sthda.com/english/wiki/survival-analysis-basics
# df <- df %>%     
#           mutate(tauratio= ptau/ttau) %>%
#           mutate(logNFL= log(NFL)) %>%
#           mutate(logptau= log(ptau)) %>%
#           mutate(logabeta= log(abeta)) %>%
#           mutate(logttau= log(ttau)) %>%
#           mutate(logtauratio= log(tauratio)) %>%

#           data.frame() #Convert to dataframe to facilitate operations.

# ###############################################################################################################################
# # CREATE NEW VARIABLES: CHANGE LEVELS OF CATEGORICAL VARS FOR TABLES
# ##############################################################################################################################
  
# df <- df %>%     

#           ##GeneticFTLD: Instead of TRUE and False have yes and no
#           mutate(GeneticFTLD= GeneticFTLD) %>%
#           mutate(GeneticFTLD = case_when(GeneticFTLD == TRUE ~ "Yes",
#                                           GeneticFTLD == FALSE ~ "No")) %>%

#           ##APOEe4: Instead of TRUE and False have Positive and Negative
#           mutate(APOEe4= APOEe4) %>%
#           mutate(APOEe4 = case_when(APOEe4 == "Yes" ~ "Positive",
#                                     APOEe4 == "No" ~ "Negative")) %>%
#           mutate(APOEe4= as.factor(APOEe4)) %>%

#           ##Parkinsonism_onset: Instead of TRUE and False have Positive and Negative
#           mutate(Parkinsonian_onset= Parkinsonian_onset) %>%
#           mutate(Parkinsonian_onset = case_when(Parkinsonian_onset == TRUE ~ "Yes",
#                                               Parkinsonian_onset == FALSE ~ "No")) %>%

#           ##Cognitive_onset: Instead of TRUE and False have Positive and Negative
#           mutate(Cognitive_onset= Cognitive_onset) %>%
#           mutate(Cognitive_onset = case_when(Cognitive_onset == TRUE ~ "Yes",
#                                               Cognitive_onset == FALSE ~ "No")) %>%

#           ##Language_onset: Instead of TRUE and False have Positive and Negative
#           mutate(Language_onset= Language_onset) %>%
#           mutate(Language_onset = case_when(Language_onset == TRUE ~ "Yes",
#                                               Language_onset == FALSE ~ "No")) %>%

#           ##Dysphagia_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Dysphagia_binary= Dysphagia_binary) %>%
#           mutate(Dysphagia_binary = case_when(Dysphagia_binary == TRUE ~ "Yes",
#                                              Dysphagia_binary == FALSE ~ "No")) %>%


#           ##Thermoregulatory_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Thermoregulatory_binary= Thermoregulatory_binary) %>%
#           mutate(Thermoregulatory_binary = case_when(Thermoregulatory_binary == TRUE ~ "Yes",
#                                                      Thermoregulatory_binary == FALSE ~ "No")) %>%


#           ##Orthostatism_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Orthostatism_binary= Orthostatism_binary) %>%
#           mutate(Orthostatism_binary = case_when(Orthostatism_binary == TRUE ~ "Yes",
#                                                 Orthostatism_binary == FALSE ~ "No")) %>%


#           ##LP2_tremor_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Tremor_binary= Tremor_binary) %>%
#           mutate(Tremor_binary = case_when(Tremor_binary == TRUE ~ "Yes",
#                                           Tremor_binary == FALSE ~ "No")) %>%

#           ##LP2_resttremor_binary: Instead of TRUE and False have Positive and Negative
#           mutate(RestTremor= RestTremor) %>%
#           mutate(RestTremor = case_when(RestTremor == TRUE ~ "Yes",
#                                         RestTremor == FALSE ~ "No")) %>%

#           ##LP2_slowness_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Slowness_binary= Slowness_binary) %>%
#           mutate(Slowness_binary = case_when(Slowness_binary == TRUE ~ "Yes",
#                                             Slowness_binary == FALSE ~ "No")) %>%

#           ##LP2_rigidity_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Rigidity_binary= Rigidity_binary) %>%
#           mutate(Rigidity_binary = case_when(Rigidity_binary == TRUE ~ "Yes",
#                                             Rigidity_binary == FALSE ~ "No")) %>%

#           ##LP2_limbrigidity_binary: Instead of TRUE and False have Positive and Negative
#           mutate(LimbRigidity= LimbRigidity) %>%
#           mutate(LimbRigidity = case_when(LimbRigidity == TRUE ~ "Yes",
#                                           LimbRigidity == FALSE ~ "No")) %>%

#           ##LP2_axialrigidity_binary: Instead of TRUE and False have Positive and Negative
#           mutate(AxialRigidity= AxialRigidity) %>%
#           mutate(AxialRigidity = case_when(AxialRigidity == TRUE ~ "Yes",
#                                           AxialRigidity == FALSE ~ "No")) %>%

#           ##LP2_verticaloculomotor_binary: Instead of TRUE and False have Positive and Negative
#           mutate(VerticalOM= VerticalOM) %>%
#           mutate(VerticalOM = case_when(VerticalOM == TRUE ~ "Yes",
#                                         VerticalOM == FALSE ~ "No")) %>%

#           ##Lifetime_VerticalOM: Instead of TRUE and False have Positive and Negative
#           mutate(Lifetime_VerticalOM= Lifetime_VerticalOM) %>%
#           mutate(Lifetime_VerticalOM = case_when(Lifetime_VerticalOM == TRUE ~ "Yes",
#                                                  Lifetime_VerticalOM == FALSE ~ "No")) %>%


#           ##Dystonia_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Dystonia_binary= Dystonia_binary) %>%
#           mutate(Dystonia_binary = case_when(Dystonia_binary == TRUE ~ "Yes",
#                                             Dystonia_binary == FALSE ~ "No")) %>%

#           ##LimbDystonia: Instead of TRUE and False have Positive and Negative
#           mutate(LimbDystonia= LimbDystonia) %>%
#           mutate(LimbDystonia = case_when(LimbDystonia == TRUE ~ "Yes",
#                                           LimbDystonia == FALSE ~ "No")) %>%

#           ##Lifetime_apraxia_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Lifetime_apraxia_binary= Lifetime_apraxia_binary) %>%
#           mutate(Lifetime_apraxia_binary = case_when(Lifetime_apraxia_binary == TRUE ~ "Yes",
#                                                  Lifetime_apraxia_binary == FALSE ~ "No")) %>%

#           ##Hypomimia_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Hypomimia_binary= Hypomimia_binary) %>%
#           mutate(Hypomimia_binary = case_when(Hypomimia_binary == TRUE ~ "Yes",
#                                               Hypomimia_binary == FALSE ~ "No")) %>%

#           ##Lifetime_RBD_binary: Instead of TRUE and False have Positive and Negative
#           mutate(RBD_binary= RBD_binary) %>%
#           mutate(RBD_binary = case_when(RBD_binary == TRUE ~ "Yes",
#                                         RBD_binary == FALSE ~ "No")) %>%

#           ##RLS_binary: Instead of TRUE and False have Positive and Negative
#           mutate(RLS_binary= RLS_binary) %>%
#           mutate(RLS_binary = case_when(RLS_binary == TRUE ~ "Yes",
#                                         RLS_binary == FALSE ~ "No")) %>%

#           ##AD_binary: Instead of TRUE and False have Positive and Negative
#           mutate(AD_binary= AD_binary) %>%
#           mutate(AD_binary = case_when(AD_binary == TRUE ~ "AD Positive",
#                                         AD_binary == FALSE ~ "AD Negative")) %>%

#           ##AD_Brink_binary: Instead of TRUE and False have Positive and Negative
#           mutate(AD_Brink_binary= AD_Brink_binary) %>%
#           mutate(AD_Brink_binary = case_when(AD_Brink_binary == TRUE ~ "AD Positive",
#                                         AD_Brink_binary == FALSE ~ "AD Negative")) %>%

#           ##AD_Brink_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Lifetime_AD_binary= Lifetime_AD_binary) %>%
#           mutate(Lifetime_AD_binary = case_when(Lifetime_AD_binary == TRUE ~ "AD Positive",
#                                                 Lifetime_AD_binary == FALSE ~ "AD Negative")) %>%

#           ##RTQUIC_lifetime: Instead of TRUE and False have Positive and Negative
#           mutate(RTQUIC= RTQUIC) %>%
#           mutate(RTQUIC = case_when(RTQUIC == "Positive" ~ "aSyn-SAA positive", #paste("\U03B1","Syn-SAA+") fails at plot level. expression(alpha*) is incompatible with case_when
#                                     RTQUIC == "Negative" ~ "aSyn-SAA negative")) %>%
#           mutate(RTQUIC= as.factor(RTQUIC)) %>%

#           ##Lifetime_VisualHallucinations_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Lifetime_VisualHallucinations_binary= Lifetime_VisualHallucinations_binary) %>%
#           mutate(Lifetime_VisualHallucinations_binary = case_when(Lifetime_VisualHallucinations_binary == TRUE ~ "Yes",
#                                                                   Lifetime_VisualHallucinations_binary == FALSE ~ "No")) %>%

#           ##Lifetime_Dopa_binary: Instead of TRUE and False have Positive and Negative
#           mutate(Lifetime_Dopa_binary= Lifetime_Dopa_binary) %>%
#           mutate(Lifetime_Dopa_binary = case_when(Lifetime_Dopa_binary == TRUE ~ "Yes",
#                                               Lifetime_Dopa_binary == FALSE ~ "No")) %>%

#           data.frame() #Convert to dataframe to facilitate operations.


###############################################################################################################################
# SAVE DATAFRAME
##############################################################################################################################
return(df)

} #var.func.1

