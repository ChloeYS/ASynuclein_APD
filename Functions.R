#Last updated 16 May 2024
##Only kept functions useful for the APD manuscript


# LIBRARIES #

## DATAFRAME MANIPULATION ##
library(tidyverse) #https://tidyverse.tidyverse.org/
  # library(dplyr) #part of tidyverse
  # library(forcats) # for factors (categorical data) #part of tidyverse
  # library(lubridate)#part of tidyverse

#library(janitor) #https://www.rdocumentation.org/packages/janitor/versions/2.1.0
# library(skimr) #Has function skim() which gives overview of data
# library(DataExplorer) # for plot_str
# library(GGally) #generates plots for the whole dataframe to give to quick overview of relationships in data
# library(lubridate) #for is.Date

## TABLES ##
#library(formattable)
#library(data.table)
# library(tableone) # for descriptive stats table one
# library(broom) #prints results of statistical tests in dataframe

## STATISTICAL TESTS ##
# library(car) #Anova() #Levene's test
# library(lme4)
# library(onewaytests) #Bartlett-Forsythe test
library(emmeans) #allows pairwise comparison for post hoc of ANCOVA

## FITTING MODELS ##
# library(boot)
# library(perm)

## POSTHOC TESTS ##
# library(multcomp) #multiple comparisons
# library(FSA) #Dunn's test with bonferroni correction
# library(pwr) #for power analysis
# library(rstatix) #wrapper for anova/ancova


## FIGURES ##
library(ggpubr) # for ggplot themes
ibrary(ggplot2) #part of tidyverse
# library(grid) #for arranging figures
# library(gridExtra) #for arranging figures


## REPORTS IN R ##
# library(knitr) #https://hbctraining.github.io/Training-modules/Rmarkdown/lesson.html
# library("webshot")
# webshot::install_phantomjs() #install once in order to export the formattable object
# library("htmltools")



# DEFINE REFERENCES LISTS #

#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
#https://www.colorhexa.com/999999
cbPalette <- c("gray46", "orangered1", "purple1", "green4", "blue") 
cbPalette2 <- c("blue", "hotpink", "slateblue1", "darkgoldenrod1", "chocolate4")
cbPalette3 <- c("blue", "slateblue1")



# FUNCTIONS USED IN THE MANUSCRIPT #

## READ.FUNC() ##
##read.func is a function that assigns the data in a given file to a dataframe
#then prints the name of the file & returns the data.
read.func <- function(file, df.name) { #file is the name of a file in the same directory, and df.name is a character string
             if (!file.exists(file)) { #then we check that the input file exists
               stop("The first argument, data file, cannot be found.")
             } else {
                 df <- read.csv(file, header=T, na.strings="") #creates dataframe from file data
                 df <- filter(df, rowSums(is.na(df)) != ncol(df)) # Apply filter function in order to remove completely empty rows. Basically, only keeps rows for which the number of NAs is different from the number of columns
                 cat('The file:', file, 'is read into the dataframe:', df.name, '\n')
                 return(df) #if run independently, read.func will print the dataframe
               } #else
             }

## STRING.VAR() ##
string.var <- function(df, working.col, string, new.var.name) {
  working.col <- df[, working.col]
  string.vec <- string #Create a new vector, contains only one element, "string". This is because names() only works on vectors
  names(string.vec) <- c(new.var.name) #The vector items are named, but there is only one item in vector anyway
  df <- cbind(df, data.frame(lapply(string.vec, grepl, working.col))) #at each row of working col, we search for string.vec value. The stringvec is converted into a Boolean, but keeps the same name (ie new.var). We cbind that new vector to the df.
  return(df)
}


## VARIANCE.FUNC() ##
##We first check that variance is equal between the groups. Based on the results
# from our normality assessment, we choose the Levene test which is robust enough
#to departure from normality but doesnt assume complete departure.
variance.func <- function (DV, IV, test=2) { #test can be Bartlett, Levene, or Brown-Forsythe
  if (test=="1") { #then perform the Bartlett's test
   cat("Performing Bartlett's test for homoscedasticity.\n")
   bart <- bartlett.test(DV ~ IV) #object of type list
   return(bart) #Returns the p-value of the test of homogeneity of variance. p=bart[[3]]
 } else {
    if (test=="2") { #then perform Levene's test
     cat("Performing Levene's test for homoscedasticity.\n")
     library(car) #needed for leveneTest
     lev <- leveneTest(DV ~ IV) #object of type list
     return(lev) #p=lev[[3]]
   } else {
      if (test=="3") {
       library(onewaytests)
       cat("Performing Brown-forsythe's F test for homoscedasticity.\n")
       brown <- bf.test(DV ~ IV) #object of type list
       return(brown) #p=brown[[4]])
      }
   }
  }
}


## TTEST_ASSUMPTIONS() ##

##We first check that variance is equal between the groups. Based on the results
ttest_assumptions <- function(df, col1, col2) {
                     DV <- df[, col1]
                     IV <- df[, col2]
                     norm <- shapiro.test(DV)
                      if (norm[[2]] <= 0.05) { #If Shapiro-Wilk's test is significant
                          var <- variance.func(DV, IV) #Do Levene's test by default, robust to normality deviation
                          if (var$Pr[1] <= 0.05) { #If Levene's test is significant
                            cat("Normality & homoscedasticity of ", col1, " violated. T-tests should not be run. Use Mann-Whitney-Wilcoxon instead. \n")
                            cat("Normality: Shapiro Wilk's test significant with:  \n")
                            print(norm)
                            cat("Homoscedasticity: Levene's test significant with:  \n")
                            print(var)
                            return("MWU")
                          } else cat("Only normality of ", col1, " violated. T-tests should not be run. \n")
                                 cat("Normality: Shapiro Wilk's test significant with with:  \n")
                                 print(norm)
                                 cat("Homoscedasticity: Levene's test not significant with:  \n")
                                 print(var)
                                 return("MWU")
                        } else var <- variance.func(DV, IV, 1) #whole test #Do Bartlett's test if Shapiro-wilks is not significant
                             # variance.func(DV, IV, 3) #need to fix Brown-Forsythe
                             if (var[[3]] <= 0.05) { #If Bartlett's test is significant
                                cat("Only homoscedasticity of ", col1, " violated. Welch's t-test can be run with the Satterthwaite/Welch's method.\n")
                                cat("Normality: Shapiro Wilk's test not significant with:  \n")
                                print(norm)
                                cat("Homoscedasticity: Bartlett's test significant with:  \n")
                                print(var)
                                return("welch")
                              } else cat("All assumptions of ", col1, " respected. Student's t-test can be run without concern using the pooled variance. \n")
                                      cat("Normality: Shapiro Wilk's test not significant with:  \n")
                                      print(norm)
                                      cat("Homoscedasticity: Bartlett's test not significant with:  \n")
                                      print(var)
                                      return("student")
                      }

## TTEST.FUNC() ##
ttest.func <- function(var, vec, df) { #3 arguments: value of IV, vector of DV, df
for (value in vec) { #for-loop that tests each variable
  cat("Performing t-tests or MWU assumption testing and analysis on", value ,"in ", var, "\n")
  test <- ttest_assumptions(df, value, var)
  DV <- df[, value]
  IV <- df[, var]
    if (test=="student") {
        cat("Perform Student's t-test with pooled variance \n")
        print(t.test(DV ~ IV, var.equal=T))
    } else if (test=="welch") {
            cat("Perform Welch's adaptation of the Student's t-test with Welch's/Satterthwaite estimate of the df  \n")
            print(t.test(DV ~ IV, var.equal=F))
        } else cat("Perform MWU \n")
               print(wilcox.test(DV ~ IV, exact=F)) #Exact=F is necessary to deal with the possible ties in ranking of the pairs
} #for-loop
}#function


