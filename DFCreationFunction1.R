##DFCreationFunction1.R
##var.func.1 creates the APD dataset


##NOTES FOR LATER: ADD THE NEW VARIABLES

###############################################################################################################################
#SOURCE PACKAGES AND FUNCTIONS
###############################################################################################################################

source('Functions.R') #Source the functions from file in same repository (for now)

###############################################################################################################################
#PREPARE FOR MODIFICATION OF THE APD DATASET ONLY
###############################################################################################################################

var.func.1 <- function(df) { #Here I'm creating a function specific to this dataset so I can reuse it in another dataset-specific analysis

  #Get rid of anyone who is not CBS or PSP
  df <- subset(df, as.numeric(Include_APDs)==1)

###############################################################################################################################
#CHECK DATES ARE IN PROPER FORMAT
###############################################################################################################################

  #Make sure that date variables are indeed dates. Usually issues arise in DOB values due to format. Reformatting original data in Excel may be required.
  library(lubridate) #for is.Date
  date.vars <- df[, c("DOB_dd.mmmm.yy", "RTQUIC_2_Date_dd.mmmm.yy", "First_Visit_Date_dd.mmmm.yy", "Last_Visit_Date_dd.mmmm.yy", "Severe_Motor_onset_dd.mmmm.yy")]
    if(sum(sapply(date.vars, is.Date)) != 5) {
      df$DOB <- as.Date(as.character(df$DOB_dd.mmmm.yy), format="%Y-%m-%d")
      df$Date <- as.Date(as.character(df$RTQUIC_2_Date_dd.mmmm.yy), format="%Y-%m-%d") #Dates used for clinical tp
      df$First_Visit <- as.Date(as.character(df$First_Visit_Date_dd.mmmm.yy), format="%Y-%m-%d")
      df$Last_Visit <- as.Date(as.character(df$Last_Visit_Date_dd.mmmm.yy), format="%Y-%m-%d")
      df$Severe_Motor_onset <- as.Date(as.character(df$Severe_Motor_onset_dd.mmmm.yy), format="%Y-%m-%d")
      date.vars <- df[, c("DOB", "Date", "First_Visit", "Last_Visit", "Severe_Motor_onset")]
      sapply(date.vars, is.Date) 
    }


###############################################################################################################################
#CHECK NUMERIC VALUES ARE IN THE PROPER FORMAT
###############################################################################################################################

# Make sure that numeric variables are indeed numeric. Usually issues arise in NFL values due to format. Reformatting original data in Excel may be required.
 num.vars <- df[, c("Education", "ID_Age_TXT", "Onset_age", "Park_onset", "Death_age",
                    "ptau_2", "ttau_2", "abeta_2", "ATI_2", "NFL_2", 
                    "YKL40_2","LP2_Cognitive_Z.score", "LP2_MOCA_Z.score", "LP2_MOCA_total", "LP2_MOCA_VS",
                    "LP2_MOCA_Language","LP2_MOCA_delayed","LP2_Orientation", "LP2_Fwd_correct", "LP2_Backwd_correct",
                    "LP2_Fluency_category_animals","LP2_TMTA","LP2_TMTB", "LP2_Craft_immediate_z.score","LP2_Craft_delayed_z.score",
                    "LP2_Benson_recall", "LP2_CDR_SOB","LP2_CDR_Total", "LP2_CVLT_10min_delay_correct", "LP2_IRI_ECS",
                    "LP2_IRI_PT", "LP2_BIS_Total", "LP2_RSMS_Total", "LP2_PSPRS", "LP2_auton_signs_n")] #"LP2_Conf_Naming",
      if(sum(sapply(num.vars, is.numeric)) != 35) {

#COULD BE A FOR-LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #Demographics & clinic
        df$Education <- as.numeric(df$Education)
        df$ID_Age_TXT <- as.numeric(df$ID_Age_TXT)
        df$Onset <- as.numeric(df$Onset_age)
        df$Park_onset <- as.numeric(df$Park_onset)
        df$Death_age <- as.numeric(df$Death_age)

        #Biomarkers
        df$ptau <- as.numeric(df$ptau_2)
        df$ttau <- as.numeric(df$ttau_2)
        df$abeta <- as.numeric(df$abeta_2)
        df$ATI <- as.numeric(df$ATI_2)
        df$NFL <- as.numeric(df$NFL_2)
        df$YKL <- as.numeric(df$YKL40_2)

        #Neuropsych  
        df$Cognitive_Z <- as.numeric(df$LP2_Cognitive_Z.score)
        df$MOCA_Z <- as.numeric(df$LP2_MOCA_Z.score)
        df$MOCA_total <- as.numeric(df$LP2_MOCA_total_corrected)
        df$MOCA_VS <- as.numeric(df$LP2_MOCA_VS)
        df$MOCA_Language <- as.numeric(df$LP2_MOCA_Language)
        df$MOCA_delayed <- as.numeric(df$LP2_MOCA_delayed)
        df$Orientation <- as.numeric(df$LP2_Orientation) #NEEDS TO BE REFORMATTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        df$Fwd <- as.numeric(df$LP2_Fwd_correct)
        df$Backwd <- as.numeric(df$LP2_Backwd_correct)
        df$Fluency <- as.numeric(df$LP2_Fluency_category_animals)
        df$TMTA <- as.numeric(df$LP2_TMTA)
        df$TMTB <- as.numeric(df$LP2_TMTB)
        df$Craft_immediate <- as.numeric(df$LP2_Craft_immediate_z.score)
        df$Craft_delayed <- as.numeric(df$LP2_Craft_delayed_z.score)
        df$Benson_recall <- as.numeric(df$LP2_Benson_recall)
        # df$LP2_Conf_Naming <- as.numeric(df$LP2_Conf_Naming) #NEEDS TO BE CREATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        df$CDR_SOB <- as.numeric(df$LP2_CDR_SOB)
        df$CDR_Total <- as.numeric(df$LP2_CDR_Total)
        df$CVLT <- as.numeric(df$LP2_CVLT_10min_delay_correct)
        df$IRI_ECS <- as.numeric(df$LP2_IRI_ECS)
        df$IRI_PT <- as.numeric(df$LP2_IRI_PT)
        df$BIS <- as.numeric(df$LP2_BIS_Total)
        df$RSMS <- as.numeric(df$LP2_RSMS_Total)

        #Motor tests
        df$PSPRS <- as.numeric(df$LP2_PSPRS)

        #Clinic
        df$Auton_n <- as.numeric(df$LP2_auton_signs_n)

        # sapply(num.vars, is.numeric) #In case of failure of as.numeric, may need to check the specific variable failing to convert properly
      }


###############################################################################################################################
#CREATE NEW VARIABLES: AGES AND DURATIONS
###############################################################################################################################

      df <- df %>%

          #age, for which defensive coding is implemented below
          mutate(Age= (as.numeric(Date - DOB))/365) %>%
          mutate(Age_SanityCheck= Age - ID_Age_TXT) %>%

          mutate(Followup_duration= (as.numeric(Last_Visit- First_Visit))/365) %>%  #Followup duration for tables of demographic

          mutate(LP2_Disease_Duration= Age - Onset) %>% #Disease duration by the time of LP
          mutate(LP2_Park_duration= Park_onset - Onset) %>% #Disease duration by the time of LP

          #Survival analysis variables: First, calculate age at onset of severe motor impairment. Then calculate duration of disease without severe impairment. 
          mutate(Severe_onset_age= (as.numeric(Severe_Motor_onset - DOB))/365) %>%
          mutate(Severe_onset_diff = Severe_onset_age - Onset) %>% #This is necessary as a variable that is later used to determine survival: ie, time without having yet severe disease
          
          #To know how long the subject was followed after onset of disease: can be useful for survival analysis (although rn already accountd for in form of Severe_Motor_Status variable)
          mutate(Age_at_last_visit= (as.numeric(Last_Visit- DOB))/365) %>%  #Age at the last visit, in order to determine how long each person was followed in total, for the survival analysis
          mutate(Duration_last_visit= (Age_at_last_visit - Onset)) %>%  #Duration between disease onset and last age checked up on

          data.frame() #Convert to dataframe to facilitate operations

          #DEFENSIVE CODING TO MAKE SURE NO ONE HAS ABERRANT AGE VALUES
          vec<- rep(1.1, nrow(df))
          vec2<- rep(-1.1, nrow(df))
          test <- df[df$Age_SanityCheck>vec | df$Age_SanityCheck<vec2,]
          if (nrow(test)>0) {
            print(test)
            print("Check ages carefully as there seems to be a discrepancy")
          }


###############################################################################################################################
#CHECK CATEGORICAL VALUES ARE IN THE PROPER FORMAT         LP2_Orthostatism_plus_other LP2_Bowel_Incontinence  LP2_RBD_plus_other
###############################################################################################################################

# Make sure that numeric variables are indeed numeric. Usually issues arise in NFL values due to format. Reformatting original data in Excel may be required.
 cat.vars <- df[, c("Sex", "Onset_type", "RTQUIC_lifetime", "AD_lifetime_ATHENA", "AD_2_ATHENA", 
                    "AD_2_LITERATURE", "AD_2_BRINK", "DX_Jabbari", "DX_Lifetime", "DX_APD", #AD_2_Abeta
                    "DX_FTLD", "DX_pathol","LP2_auton_signs", "LP2_Dysphagia", "LP2_Sexual_dysfunction",
                    "LP2_Constipation", "LP2_Urinary","LP2_Bowel_Incontinence", "LP2_Hyperhidrosis_Thermoregulatory_plus_other","LP2_Orthostatism_plus_other",
                    "LP2_RBD_plus_other", "LP2_Anosmia","APOEe4_alleles", "Stroke_TIA", "HeartAttack_CardiacArrest",
                    "Cardiovascular_other","Smoking","DIAB","HYT", "HYPERCHOL",
                    "TBI", "Seizure_baseline","Arthritis","Medical_other", "LP2_Delusions", 
                    "Lifetime_Hallucinations","LP2_gait", "LP2_falls_PI","LP2_retropulsion","LP2_tremor",
                    "LP2_slowness", "LP2_oculomotor", "Lifetime_oculomotor", "LP2_rigidity","LP2_dystonia", 
                    "LP2_apraxia","Lifetime_apraxia", "LP2_myoclonus", "LP2_alien_limb","LP2_hypomimia",
                    "Medication_1_anx","Medication_2_antidep", "Medication_3_ACheL","Medication_4_memantine","Medication_5_AP",
                    "Medication_6_Mood", "Medication_7_Stim", "Medication_8_Vasc","Medication_9_Sleep","Medication_10_Mov",
                    "Lifetime_Dopa","Lifetime_Dopa_responder_true", "Severe_Motor_onset_status_12")]

      if(sum(sapply(cat.vars, is.factor)) != 63) {

#COULD BE A FOR-LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #Demographics
        df$Sex <- as.factor(df$Sex)
        df$DX_Jabbari <- as.factor(df$DX_Jabbari)
        df$DX_Lifetime <- as.factor(df$DX_Lifetime)
        df$DX_APD <- as.factor(df$DX_APD)
        df$DX_FTLD <- as.factor(df$DX_FTLD)
        df$DX_pathol <- as.factor(df$DX_pathol)
        df$APOEe4 <- as.factor(df$APOEe4_alleles)
        df$Onset_type <- as.factor(df$Onset_type)
        df$Severe_status <- as.factor(df$Severe_Motor_onset_status_12)

        #Biomarkers
        df$RTQUIC <- as.factor(df$RTQUIC_lifetime)
        df$Lifetime_AD <- as.factor(df$AD_lifetime_ATHENA)
        df$AD <- as.factor(df$AD_2_ATHENA)
        df$AD_lit <- as.factor(df$AD_2_LITERATURE)
        df$AD_Brink <- as.factor(df$AD_2_BRINK)

        #Autonomic or PD-related non motor
        df$Auton_binarized <- as.factor(df$LP2_auton_signs)
        df$Dysphagia <- as.factor(df$LP2_Dysphagia) #Cannot be used as is
        df$Sexual <- as.factor(df$LP2_Sexual_dysfunction)
        df$Constipation <- as.factor(df$LP2_Constipation)
        df$Urinary <- as.factor(df$LP2_Urinary)
        df$Bowel <- as.factor(df$LP2_Bowel_Incontinence)
        df$Thermoregulatory <- as.factor(df$LP2_Hyperhidrosis_Thermoregulatory_plus_other) #Cannot be used as is
        df$Orthostatism <- as.factor(df$LP2_Orthostatism_plus_other) #Cannot be used as is
        df$RBD_plus <- as.factor(df$LP2_RBD_plus_other) #Cannot be used as is
        df$Anosmia <- as.factor(df$LP2_Anosmia)
        
        #Vascular RF
        df$Stroke_TIA <- as.factor(df$Stroke_TIA)
        df$HeartAttack_CardiacArrest <- as.factor(df$HeartAttack_CardiacArrest)
        df$Cardiovascular_other <- as.factor(df$Cardiovascular_other)
        df$Smoking <- as.factor(df$Smoking)
        df$DIAB <- as.factor(df$DIAB)
        df$HYT <- as.factor(df$HYT)
        df$HYPERCHOL <- as.factor(df$HYPERCHOL)

        #Neurological/psych history
        df$TBI <- as.factor(df$TBI)
        df$Seizure_baseline <- as.factor(df$Seizure_baseline)
        df$Delusions <- as.factor(df$LP2_Delusions)
        df$Lifetime_Hallucinations <- as.factor(df$Lifetime_Hallucinations)

        #Other history
        df$Arthritis <- as.factor(df$Arthritis)
        df$Medical_other <- as.factor(df$Medical_other)

        #Motor
        df$Gait <- as.factor(df$LP2_gait)
        df$Falls_PI <- as.factor(df$LP2_falls_PI)
        df$Retropulsion <- as.factor(df$LP2_retropulsion)
        df$Tremor <- as.factor(df$LP2_tremor)
        df$Slowness <- as.factor(df$LP2_slowness)
        df$OM <- as.factor(df$LP2_oculomotor)
        df$Lifetime_OM <- as.factor(df$Lifetime_oculomotor)
        df$Rigidity <- as.factor(df$LP2_rigidity)
        df$Dystonia <- as.factor(df$LP2_dystonia)
        df$Apraxia <- as.factor(df$LP2_apraxia)
        df$Lifetime_apraxia <- as.factor(df$Lifetime_apraxia)
        df$Myoclonus <- as.factor(df$LP2_myoclonus)
        df$Alien <- as.factor(df$LP2_alien_limb)
        df$Hypomimia <- as.factor(df$LP2_hypomimia)

        #Medication
        df$Medication_1_anx <- as.factor(df$Medication_1_anx)
        df$Medication_2_antidep <- as.factor(df$Medication_2_antidep)
        df$Medication_3_ACheL <- as.factor(df$Medication_3_ACheL)
        df$Medication_4_memantine <- as.factor(df$Medication_4_memantine)
        df$Medication_5_AP <- as.factor(df$Medication_5_AP)
        df$Medication_6_Mood <- as.factor(df$Medication_6_Mood)
        df$Medication_7_Stim <- as.factor(df$Medication_7_Stim)
        df$Medication_8_Vasc <- as.factor(df$Medication_8_Vasc)
        df$Medication_9_Sleep <- as.factor(df$Medication_9_Sleep)
        df$Medication_10_Mov <- as.factor(df$Medication_10_Mov)
        df$Lifetime_Dopa <- as.factor(df$Lifetime_Dopa)
        df$Lifetime_Dopa_responder_true <- as.factor(df$Lifetime_Dopa_responder_true)

        # sapply(cat.vars, is.factor) #In case of failure of as.factor, may need to check the specific variable failing to convert properly

      }


###############################################################################################################################
# CREATE NEW VARIABLES: BINARY STRINGS
##WARNING: BE VERY CAREFUL WITH NAS OR NON-BINARY VARS. IT WILL BINARIZE THE DATA AND LUMP ANY NA WITH THE NON-REF LEVEL
###############################################################################################################################

      ##Demographics
      df <-string.var(df, "DX_Lifetime", "genetic", "GeneticFTLD")
      df <-string.var(df, "Onset_type", "Parkinsonism", "Parkinsonian_onset")
      df <-string.var(df, "Onset_type", "Cognitive", "Cognitive_onset")
      df <-string.var(df, "Onset_type", "Language", "Language_onset")

      ##Biomarkers
      df <-string.var(df, "Lifetime_AD", "Positive", "Lifetime_AD_binary") 
      df <-string.var(df, "AD_lit", "Positive", "AD_lit_binary") 
      df <-string.var(df, "AD", "Positive", "AD_binary") 
      df <-string.var(df, "AD_Brink", "Positive", "AD_Brink_binary") 

      ##Autonomic or PD-related non motor
      df <-string.var(df, "Dysphagia", "Yes", "Dysphagia_binary")
      df <-string.var(df, "Thermoregulatory", "Yes", "Thermoregulatory_binary")
      df <-string.var(df, "Orthostatism", "Yes", "Orthostatism_binary")
      df <-string.var(df, "RBD_plus", "Yes", "RBD_binary")
      df <-string.var(df, "RBD_plus", "RLS", "RLS_binary")
      df <-string.var(df, "Anosmia", "Yes", "Anosmia_binary") 
      df <-string.var(df, "Sexual", "Yes", "Sexual_binary") 
      df <-string.var(df, "Constipation", "Yes", "Constipation_binary") 
      df <-string.var(df, "Urinary", "Yes", "Urinary_binary") 
      df <-string.var(df, "Bowel", "Yes", "Bowel_binary") 

      ##Neurological/psych history
      df <-string.var(df, "Lifetime_Hallucinations", "Yes", "Lifetime_Hallucinations_binary")            
      df <-string.var(df, "Lifetime_Hallucinations", "Visual", "Lifetime_VisualHallucinations_binary")

      ##Medication
      df <-string.var(df, "Lifetime_Dopa", "Yes", "Lifetime_Dopa_binary")

      #Motor
      df <-string.var(df, "Tremor", "Yes", "Tremor_binary")
      df <-string.var(df, "Tremor", "action", "ActionTremor")
      df <-string.var(df, "Tremor", "rest", "RestTremor")

      df <-string.var(df, "Slowness", "Yes", "Slowness_binary")

      df <-string.var(df, "OM", "Yes", "OM_binary")
      df <-string.var(df, "OM", "vertical", "VerticalOM")
      df <-string.var(df, "Lifetime_OM", "vertical", "Lifetime_VerticalOM")

      df <-string.var(df, "Rigidity", "Yes", "Rigidity_binary")
      df <-string.var(df, "Rigidity", "axial", "AxialRigidity")
      df <-string.var(df, "Rigidity", "limb", "LimbRigidity")

      df <-string.var(df, "Dystonia", "Yes", "Dystonia_binary")
      df <-string.var(df, "Dystonia", "limb", "LimbDystonia")

      df <-string.var(df, "Lifetime_apraxia", "Yes", "Lifetime_apraxia_binary")

      df <-string.var(df, "Hypomimia", "Yes", "Hypomimia_binary") #Counts poker face as no hypomimia


###############################################################################################################################
# CREATE NEW VARIABLES: COMBINE VARIABLES BASED ON STRINGS
##############################################################################################################################

df <- df %>%

          ##Onset_type_simplified: multiple domains instead of "Cognitive, Language" etc
          mutate(Onset_type_simplified= Onset_type) %>%
          mutate(Onset_type_simplified= case_when(grepl(",", Onset_type_simplified) ~ "Multiple domains",
                                                  TRUE ~ as.character(Onset_type_simplified))) %>% #TRUE in this case captures all the values excluded from previous conditions,
                                                                                                    #in which case case_when applies the new value, which is actually the original value.
          
          ##PPA: whether a PPA was diagnosed or not
          mutate(anyPPA= PPA) %>%
          mutate(anyPPA = case_when(anyPPA == "svPPA" ~ "Yes",
                                    anyPPA == "PNFA" ~ "Yes",
                                    anyPPA == "lvPPA" ~ "Yes",
                                    anyPPA == "PPA" ~ "Yes",
                                    TRUE ~ as.character(anyPPA))) %>%

          data.frame() #Convert to dataframe to facilitate operations.

###############################################################################################################################
# CREATE NEW VARIABLES: CREATE RATIO AND LOGS
##############################################################################################################################
#http://www.sthda.com/english/wiki/survival-analysis-basics
df <- df %>%     
          mutate(tauratio= ptau/ttau) %>%
          mutate(logNFL= log(NFL)) %>%
          mutate(logptau= log(ptau)) %>%
          mutate(logabeta= log(abeta)) %>%
          mutate(logttau= log(ttau)) %>%
          mutate(logtauratio= log(tauratio)) %>%

          data.frame() #Convert to dataframe to facilitate operations.

###############################################################################################################################
# CREATE NEW VARIABLES: CHANGE LEVELS OF CATEGORICAL VARS FOR TABLES
##############################################################################################################################
  
df <- df %>%     

          ##GeneticFTLD: Instead of TRUE and False have yes and no
          mutate(GeneticFTLD= GeneticFTLD) %>%
          mutate(GeneticFTLD = case_when(GeneticFTLD == TRUE ~ "Yes",
                                          GeneticFTLD == FALSE ~ "No")) %>%

          ##APOEe4: Instead of TRUE and False have Positive and Negative
          mutate(APOEe4= APOEe4) %>%
          mutate(APOEe4 = case_when(APOEe4 == "Yes" ~ "Positive",
                                    APOEe4 == "No" ~ "Negative")) %>%
          mutate(APOEe4= as.factor(APOEe4)) %>%

          ##Parkinsonism_onset: Instead of TRUE and False have Positive and Negative
          mutate(Parkinsonian_onset= Parkinsonian_onset) %>%
          mutate(Parkinsonian_onset = case_when(Parkinsonian_onset == TRUE ~ "Yes",
                                              Parkinsonian_onset == FALSE ~ "No")) %>%

          ##Cognitive_onset: Instead of TRUE and False have Positive and Negative
          mutate(Cognitive_onset= Cognitive_onset) %>%
          mutate(Cognitive_onset = case_when(Cognitive_onset == TRUE ~ "Yes",
                                              Cognitive_onset == FALSE ~ "No")) %>%

          ##Language_onset: Instead of TRUE and False have Positive and Negative
          mutate(Language_onset= Language_onset) %>%
          mutate(Language_onset = case_when(Language_onset == TRUE ~ "Yes",
                                              Language_onset == FALSE ~ "No")) %>%

          ##Dysphagia_binary: Instead of TRUE and False have Positive and Negative
          mutate(Dysphagia_binary= Dysphagia_binary) %>%
          mutate(Dysphagia_binary = case_when(Dysphagia_binary == TRUE ~ "Yes",
                                             Dysphagia_binary == FALSE ~ "No")) %>%


          ##Thermoregulatory_binary: Instead of TRUE and False have Positive and Negative
          mutate(Thermoregulatory_binary= Thermoregulatory_binary) %>%
          mutate(Thermoregulatory_binary = case_when(Thermoregulatory_binary == TRUE ~ "Yes",
                                                     Thermoregulatory_binary == FALSE ~ "No")) %>%


          ##Orthostatism_binary: Instead of TRUE and False have Positive and Negative
          mutate(Orthostatism_binary= Orthostatism_binary) %>%
          mutate(Orthostatism_binary = case_when(Orthostatism_binary == TRUE ~ "Yes",
                                                Orthostatism_binary == FALSE ~ "No")) %>%


          ##LP2_tremor_binary: Instead of TRUE and False have Positive and Negative
          mutate(Tremor_binary= Tremor_binary) %>%
          mutate(Tremor_binary = case_when(Tremor_binary == TRUE ~ "Yes",
                                          Tremor_binary == FALSE ~ "No")) %>%

          ##LP2_resttremor_binary: Instead of TRUE and False have Positive and Negative
          mutate(RestTremor= RestTremor) %>%
          mutate(RestTremor = case_when(RestTremor == TRUE ~ "Yes",
                                        RestTremor == FALSE ~ "No")) %>%

          ##LP2_slowness_binary: Instead of TRUE and False have Positive and Negative
          mutate(Slowness_binary= Slowness_binary) %>%
          mutate(Slowness_binary = case_when(Slowness_binary == TRUE ~ "Yes",
                                            Slowness_binary == FALSE ~ "No")) %>%

          ##LP2_rigidity_binary: Instead of TRUE and False have Positive and Negative
          mutate(Rigidity_binary= Rigidity_binary) %>%
          mutate(Rigidity_binary = case_when(Rigidity_binary == TRUE ~ "Yes",
                                            Rigidity_binary == FALSE ~ "No")) %>%

          ##LP2_limbrigidity_binary: Instead of TRUE and False have Positive and Negative
          mutate(LimbRigidity= LimbRigidity) %>%
          mutate(LimbRigidity = case_when(LimbRigidity == TRUE ~ "Yes",
                                          LimbRigidity == FALSE ~ "No")) %>%

          ##LP2_axialrigidity_binary: Instead of TRUE and False have Positive and Negative
          mutate(AxialRigidity= AxialRigidity) %>%
          mutate(AxialRigidity = case_when(AxialRigidity == TRUE ~ "Yes",
                                          AxialRigidity == FALSE ~ "No")) %>%

          ##LP2_verticaloculomotor_binary: Instead of TRUE and False have Positive and Negative
          mutate(VerticalOM= VerticalOM) %>%
          mutate(VerticalOM = case_when(VerticalOM == TRUE ~ "Yes",
                                        VerticalOM == FALSE ~ "No")) %>%

          ##Lifetime_VerticalOM: Instead of TRUE and False have Positive and Negative
          mutate(Lifetime_VerticalOM= Lifetime_VerticalOM) %>%
          mutate(Lifetime_VerticalOM = case_when(Lifetime_VerticalOM == TRUE ~ "Yes",
                                                 Lifetime_VerticalOM == FALSE ~ "No")) %>%


          ##Dystonia_binary: Instead of TRUE and False have Positive and Negative
          mutate(Dystonia_binary= Dystonia_binary) %>%
          mutate(Dystonia_binary = case_when(Dystonia_binary == TRUE ~ "Yes",
                                            Dystonia_binary == FALSE ~ "No")) %>%

          ##LimbDystonia: Instead of TRUE and False have Positive and Negative
          mutate(LimbDystonia= LimbDystonia) %>%
          mutate(LimbDystonia = case_when(LimbDystonia == TRUE ~ "Yes",
                                          LimbDystonia == FALSE ~ "No")) %>%

          ##Lifetime_apraxia_binary: Instead of TRUE and False have Positive and Negative
          mutate(Lifetime_apraxia_binary= Lifetime_apraxia_binary) %>%
          mutate(Lifetime_apraxia_binary = case_when(Lifetime_apraxia_binary == TRUE ~ "Yes",
                                                 Lifetime_apraxia_binary == FALSE ~ "No")) %>%

          ##Hypomimia_binary: Instead of TRUE and False have Positive and Negative
          mutate(Hypomimia_binary= Hypomimia_binary) %>%
          mutate(Hypomimia_binary = case_when(Hypomimia_binary == TRUE ~ "Yes",
                                              Hypomimia_binary == FALSE ~ "No")) %>%

          ##Lifetime_RBD_binary: Instead of TRUE and False have Positive and Negative
          mutate(RBD_binary= RBD_binary) %>%
          mutate(RBD_binary = case_when(RBD_binary == TRUE ~ "Yes",
                                        RBD_binary == FALSE ~ "No")) %>%

          ##RLS_binary: Instead of TRUE and False have Positive and Negative
          mutate(RLS_binary= RLS_binary) %>%
          mutate(RLS_binary = case_when(RLS_binary == TRUE ~ "Yes",
                                        RLS_binary == FALSE ~ "No")) %>%

          ##AD_binary: Instead of TRUE and False have Positive and Negative
          mutate(AD_binary= AD_binary) %>%
          mutate(AD_binary = case_when(AD_binary == TRUE ~ "AD Positive",
                                        AD_binary == FALSE ~ "AD Negative")) %>%

          ##AD_Brink_binary: Instead of TRUE and False have Positive and Negative
          mutate(AD_Brink_binary= AD_Brink_binary) %>%
          mutate(AD_Brink_binary = case_when(AD_Brink_binary == TRUE ~ "AD Positive",
                                        AD_Brink_binary == FALSE ~ "AD Negative")) %>%

          ##AD_Brink_binary: Instead of TRUE and False have Positive and Negative
          mutate(Lifetime_AD_binary= Lifetime_AD_binary) %>%
          mutate(Lifetime_AD_binary = case_when(Lifetime_AD_binary == TRUE ~ "AD Positive",
                                                Lifetime_AD_binary == FALSE ~ "AD Negative")) %>%

          ##RTQUIC_lifetime: Instead of TRUE and False have Positive and Negative
          mutate(RTQUIC= RTQUIC) %>%
          mutate(RTQUIC = case_when(RTQUIC == "Positive" ~ "aSyn-SAA positive", #paste("\U03B1","Syn-SAA+") fails at plot level. expression(alpha*) is incompatible with case_when
                                    RTQUIC == "Negative" ~ "aSyn-SAA negative")) %>%
          mutate(RTQUIC= as.factor(RTQUIC)) %>%

          ##Lifetime_VisualHallucinations_binary: Instead of TRUE and False have Positive and Negative
          mutate(Lifetime_VisualHallucinations_binary= Lifetime_VisualHallucinations_binary) %>%
          mutate(Lifetime_VisualHallucinations_binary = case_when(Lifetime_VisualHallucinations_binary == TRUE ~ "Yes",
                                                                  Lifetime_VisualHallucinations_binary == FALSE ~ "No")) %>%

          ##Lifetime_Dopa_binary: Instead of TRUE and False have Positive and Negative
          mutate(Lifetime_Dopa_binary= Lifetime_Dopa_binary) %>%
          mutate(Lifetime_Dopa_binary = case_when(Lifetime_Dopa_binary == TRUE ~ "Yes",
                                              Lifetime_Dopa_binary == FALSE ~ "No")) %>%

          data.frame() #Convert to dataframe to facilitate operations.


###############################################################################################################################
# SAVE DATAFRAME
##############################################################################################################################
return(df)

} #var.func.1


