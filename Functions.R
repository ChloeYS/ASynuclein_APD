#Functions.R

#################################################################################
###########################         INDEX         ###############################

#######################         BASIC FUNCTIONS         ########################
###1. READ.FUNC
###2. HIGHEST FREQUENCY OF A GIVEN LABEL FOR ONE VARIABLE

#############################      FIT MODELS      #############################
####1. FIT A LINEAR MODEL
####2. FIT A QUADRATIC MODEL
####3.PLOT A LINEAR OR A QUADRATIC MODEL
####4.PLOT A LINEAR AND A QUADRATIC MODEL ON EACH OTHER
####5. ASSESS THE GOODNESS OF FIT OF A LINEAR OR QUADRATIC MODEL


#######################       EXPORT & SAVE FUNCTIONS      #####################
####1. EXPORT TABLES
####2. FORMATTER FUNCTION TO PICK COLOR BASED ON P-VALUE

##############################       FIGURES      ##############################
####1. BOXPLOTS
####2. RADAR PLOTS


#################################################################################
############################      PACKAGES         #############################

##INSTALL NEW PACKAGES
##Traditional installation is with CRAN selection:
#install.packages('', dependencies = TRUE)

##But on my laptop anaconda has been set up, causing issues. In these cases,
#consider running a conda installation.
##https://towardsdev.com/install-r-in-conda-8b2033ec3d4f
# conda install -c r r-tidyverse


##DATAFRAME MANIPULATION
library(tidyverse)
  #https://tidyverse.tidyverse.org/
  # library(dplyr) #part of tidyverse
  # library(forcats) # for factors (categorical data) #part of tidyverse
  # library(lubridate)#part of tidyverse
  
library(janitor)
  #https://www.rdocumentation.org/packages/janitor/versions/2.1.0
# library(skimr) #Has function skim() which gives overview of data\
# library(DataExplorer) # for plot_str
# library(GGally) #generates plots for the whole dataframe to give to quick overview of relationships in data
# library(lubridate) #for is.Date

##TABLES
library(formattable)
library(data.table)
# library(tableone) # for descriptive stats table one
# library(broom) #prints results of statistical tests in dataframe

##STATISTICAL TESTSs
# library(car) #Anova() #Levene's test
# library(lme4)
# library(onewaytests) #Bartlett-Forsythe test

##FITTING MODELS
# library(boot)
# library(perm)

##POSTHOC TESTSs
# library(multcomp) #multiple comparisons
# library(FSA) #Dunn's test with bonferroni correction
# library(pwr) #for power analysis
# library(rstatix) #wrapper for anova/ancova
# library(emmeans) #allows pairwise comparison for post hoc of ANCOVA


##FIGURES
# library(ggpubr) # for ggplot themes
# library(ggplot2) #part of tidyverse
# library(grid) #for arranging figures
# library(gridExtra) #for arranging figures


##REPORTS IN R
# library(knitr)
#https://hbctraining.github.io/Training-modules/Rmarkdown/lesson.html
# library("webshot")
  # webshot::install_phantomjs() #install once in order to export the formattable object
# library("htmltools")

################################################################################
#####################       DEFINE REFERENCES TOOLS     ########################
#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
#https://www.colorhexa.com/999999
cbPalette <- c("gray46", "orangered1", "purple1", "green4", "blue")
cbPalette2 <- c("blue", "hotpink", "slateblue1", "darkgoldenrod1", "chocolate4")
cbPalette3 <- c("blue", "slateblue1")

#################################################################################
#######################         BASIC FUNCTIONS         ########################

###1. READ.FUNC
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


###2. HIGHEST FREQUENCY OF A GIVEN LABEL FOR ONE VARIABLE
##highest.frequency is a function, which, for a given categorical variable, identifies the label with the highest number of a given label of a second variable.
#In other words, for a df where var 1 = neighbourhoods, and var2 = COVID outcomes, we can identify the neighbourhood that has the most COVID fatalities specifically.
#col1: a vector of all values of the first variable. Categorical.
#col2: a vector of all values of the second variable. Categorical.
#label: here, unfortunately, the label of interest in var 2 ("Positive") is hardcoded.
#frequencies: a table presenting all combinations of outcomes by neighborhoods. Contains one FATAL column which is the number of FATAL outcomes per neighborhood.
#fatalities.ranked: vector of the FATAL column in frequencies, ranked from highest to lowest. Its first element is the highest number of fatalities for one neighborhood.
#neighbourhoods.ranked: vector with neighbourhood names ranked based on the number of fatalities. Its first element is the neighborhood with the most fatalities.

highest.frequency <- function(df, working.col1, working.col2) { #requires 1 df, 2 vars
    working.col1 <- df[, working.col1] #subsets values of variable 1
    working.col2 <- df[, working.col2] #subsets values of variable 2
    frequencies <- table(df[, c(working.col1, working.col2)]) #frequency table indicating how many times var2 labels happen for var 1, eg: how many times each outcome happens for each neighborhood.
    ranked <- (sort(frequencies[, 'Positive'], decreasing=T)) #keeps only the label of interest (Positive) and sorts var1 by frequency
    ranked.names <- names(ranked) #keeps the names of var1, in the same order
    cat(ranked[1], 'with', ranked.names[1], 'RT QuIC positive cases.', '\n')
}

###3. CREATE NEW VARIABLE BASED ON STRING IN ONE VARIABLE
string.var <- function(df, working.col, string, new.var.name) {
  working.col <- df[, working.col]
  string.vec <- string #Create a new vector, contains only one element, "string". This is because names() only works on vectors
  names(string.vec) <- c(new.var.name) #The vector items are named, but there is only one item in vector anyway
  df <- cbind(df, data.frame(lapply(string.vec, grepl, working.col))) #at each row of working col, we search for string.vec value. The stringvec is converted into a Boolean, but keeps the same name (ie new.var). We cbind that new vector to the df.
  return(df)
}

################################################################################
#######################         STATS FUNCTIONS         ########################

####1. CHECK HOMOSCEDASTICITY
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

####2. CHECK T-TEST ASSUMPTIONs
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

####3. RUN T-TEST OR MWU
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


####4. CHECK ANCOVA ASSUMPTIONs
##We first check that variance is equal between the groups. Based on the results
#CHeck for outliers
#Check for regression homogeneity


################################################################################
############################      FIT MODELS      ##############################

####1. FIT A LINEAR MODEL
##lin.mod.func returns a linear model for a given set of data.
lin.mod.func <- function (dataframe, col1, col2) {
  DV <- dataframe[, col1]
  IV <- dataframe[, col2]
  lin.mod <- lm(DV ~ IV, dataframe) #Creates linear model where TBIs predicts other injuries
  return(lin.mod) #Returns the model
}

####2. FIT A QUADRATIC MODEL
##quad.mod.func returns a quadratic model for a given set of data.
quad.mod.func <- function (dataframe, col1, col2) {
  DV <- dataframe[, col1]
  IV <- dataframe[, col2]
  square.IV <- IV^2 #Assign the square of the IV to a new variable
  quad.mod <- lm(DV ~ IV + square.IV, dataframe) #Creates quadratic model where TBIs predicts other injuries
  return(quad.mod) #Return the model so it can be used in other contexts
}

####3.PLOT A LINEAR OR A QUADRATIC MODEL
##plot.func takes an already fitted model (linear or quadratic), its variables, and graphs them.
plot.func <- function (model, dataframe, col1, col2) { #give already fitted model and a dataframe
  DV <- dataframe[, col1]
  IV <- dataframe[, col2]
  if (length(coef(model))==2) { #if there are 2 coefficients for the model (ie it's the linear model)
    plot(IV, DV, ylab=col1, xlab=col2, main="Linear model fit", lwd=2,col="black")
      abline(model, col="red") #Same as model
  } else {
      if (length(coef(model))==3) { #if there are 3 coefficients (ie it's the quadratic model)
        plot(IV,DV)
        axis.x <- seq(min(IV), max(IV), len=length(IV)) #Create x axis covering range of IV
        axis.y <- model$coef %*% rbind(1,axis.x,axis.x^2) #Each value of y axis is recomputed from the coefficients and the x values
        lines(axis.x,axis.y, ylab=col1, xlab=col2,  main="Quadratic model fit", lwd=2,col="blue") #Blue line that represents quadratic fit to the data is added to real data
      }
  }
}

####4.PLOT A LINEAR AND A QUADRATIC MODEL ON EACH OTHER
##comp.mod.func compares the linear and the quadratic models to each other.
comp.mod.func <- function (dataframe, col1, col2) {
   lin.mod <- lin.mod.func(dataframe, col1, col2) #Computes the linear model
   quad.mod <- quad.mod.func(dataframe, col1, col2) #Computes the quadratic model
   plot.func(quad.mod, dataframe, col1, col2) #Plots the quadratic curve
   abline(lin.mod, col="red") #Adds the linear relationship
 }

####5. ASSESS THE GOODNESS OF FIT OF A LINEAR OR QUADRATIC MODEL
##AnalyzeModel generates plots that help assess the goodness of fit of a given model.
#It also identifies data points potentially disproportionally influencing the model & returns them in a list for the given tolerance level.
AnalyzeModel <- function (dataframe, model, col1, col2, tolerance=0.25) { #tolerance has to be between 0 and 1

   cat("The set tolerance level is", tolerance, ".\n")
   cat("Plotting leverage, Cook's distance, and deviation from normality for residuals and Studentized residuals.\n")

   #Assign variables to use
   DV <- dataframe[, col1] #For manipulation purposes
   IV <- dataframe[, col2]
   dataframe <- cbind(dataframe, model$fitted) #Enables us to use the model fitted values based on indices

   ##Cook's distance from the model:
   cooks <- cooks.distance(model) #Computes Cook's distance
   high.cook <- dataframe[cooks>tolerance, ] #identifies high Cook's distance data
   high.cook.index <- as.numeric(rownames(high.cook)) #Identifies index of points with high Cook
   plot(cooks, ylab="Cook's distance", main="Cook's distance within the model") #Plots Cook's distance
   points(high.cook.index, cooks[cooks>tolerance], col="red") #Discriminates high Cook's distance points

   ##Leverage in the model:
   leverage <- hat(model.matrix(model)) #Computes leverage
   plot(leverage, ylab="Leverage", main="Leverage within the model") #Plots leverage
   points(high.cook.index, leverage[high.cook.index], col="red") #Discriminates high Cook's distance points

   par(mfrow=c(2,2)) #Organizes visualization of plots into 2 rows, 2 columns per page

   ##Deviation from normality of the residuals:
   plot(IV, model$res, xlab=col2, ylab="Residuals") #residuals vs explanatory variable
    points(high.cook[ ,2], model$res[high.cook.index], col="red") #Discriminates high Cook's distance points
   plot(DV, model$res, xlab=col1, ylab="Residuals") #residuals vs explained variable
    points(high.cook[ ,1], model$res[high.cook.index], col="red") #Discriminates high Cook's distance points
   plot(model$fitted, model$res, xlab="Predicted values", ylab="Residuals") #residuals vs predicted values
    points(high.cook[ ,3], model$res[high.cook.index], col="red") #Discriminates high Cook's distance points
   qqnorm(model$res, ylab="Residuals") #Plot of normality of residuals
    qqline(model$res) #Add line representing normality to the plot of residuals

   ##Deviation from normality of the Studentized residuals:
   plot(IV, rstudent(model), xlab=col2, ylab="Studentized residuals") #studentized residuals vs explanatory variable
    points(high.cook[ ,2], rstudent(model)[high.cook.index], col="red") #Discriminates high Cook's distance points
   plot(DV, rstudent(model), xlab=col1, ylab="Studentized residuals") #studentized residuals vs explained variable
    points(high.cook[ ,1], rstudent(model)[high.cook.index], col="red") #Discriminates high Cook's distance points
   plot(model$fitted, rstudent(model), xlab="Predicted values", ylab="Studentized residuals")  #studentized residuals vs predicted values
    points(high.cook[ ,3], rstudent(model)[high.cook.index], col="red") #Discriminates high Cook's distance points
   qqnorm(rstudent(model), ylab="Studentized residuals")
    qqline(rstudent(model))

   par(mfrow=c(1,1)) #Reset visualization to one figure per page

   ##Histograms of the residuals:
   hist(model$res, xlab="Residuals")
   hist(rstudent(model), xlab="Studentized residuals")

   ##Return list of suspicious points:
   cat("Index of points with Cook's distance above", tolerance, " :\n")
   return(high.cook.index) #Returns a list of points that had Cook,s distance above tolerance threshold
}

################################################################################
#######################      EXPORT & SAVE FUNCTIONS      ######################
####1. EXPORT TABLES
####2. FORMATTER FUNCTION TO PICK COLOR BASED ON P-VALUE

####1. EXPORT TABLES
#https://stackoverflow.com/questions/38833219/command-for-exporting-saving-table-made-with-formattable-package-in-r
#https://tableone.readthedocs.io/_/downloads/en/latest/pdf/
export_formattable <- function(f, file, width = "100%", height = NULL,
                               background = "white", delay = 0.2) {
                      library("webshot")
                      library("htmltools")
                      w <- as.htmlwidget(f, width = width, height = height)
                      path <- html_print(w, background = background, viewer = NULL)
                      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
                      webshot(url, file = file, selector = ".formattable_widget", delay = delay)
                      }

####2. FORMATTER FUNCTIONS TO PICK COLOR BASED ON P-VALUE
#https://clarewest.github.io/blog/post/making-tables-shiny/
p_formatter <- formatter("span", style= x ~ ifelse(as.numeric(x)<=0.05, "color: forestgreen", "color: black"))

# p_formatter <- function() {formatter("span", style= ~ style(color= ifelse(as.numeric(p)<=0.05, "green", "black")))}

p_formatter_0.1 <- formatter("span",
                             style= x ~
                             ifelse(as.numeric(x)<=0.05, "font-size:15px; color: forestgreen", #struggling to get the font-weight to vary based on conditions and to show noticeable variation
                             ifelse(x<=0.1, "font-size:15px; color: orange",
                                            "font-size:15px; color: black")))

color_formatter <- function() {
    formatter("span",
        style = ~ style(color = ifelse(History == "Current", "red",
                                ifelse(History == "Remote", "orange",
                                "black"))), font.weight = "bold")
        }



##############################       FIGURES      ##############################
####1. BOXPLOTS
####2. RADAR PLOTS
#Captions and legends: https://ggplot2.tidyverse.org/reference/labs.html
#
# #   # annotate("segment", x = 1, xend = 2, y = 2575, yend=2575, alpha=0.8) +
# #   # geom_text(aes(label="****"), x = 1.5, y = 2600, linewidth=7) +

####1. BOXPLOTS

box_simple <- function(df, working.col1, working.col2, title="", axis1, axis2, figname="Fig.png", palette=cbPalette) {
     
      working.col1 <- df[, working.col1] #subsets values of variable 1
      working.col2 <- df[, working.col2] #subsets values of variable 2

      Fig <- df %>%

      #Main structure
        ggplot(aes(x = working.col1, y = working.col2, color= working.col1, fill = working.col1)) + 
        scale_fill_manual(values=cbPalette) +
        scale_color_manual(values=cbPalette) +
        geom_boxplot(alpha=0.2) +
        geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=1) +

        #Theme
        theme_classic() +

        #Title
        ggtitle(title) +
        theme(plot.title = element_text(size = rel(2), hjust = 0.5)) +

        #Axes
         labs(x=axis1, y=axis2) +
         theme(axis.title.y = element_text(size = rel(1.5)),
          axis.title.x = element_text(size = rel(1.5)),
          axis.text.x = element_text(hjust = 0.5, size=13,color="black")) + 

         #Legend
         # theme(legend.position="none") + #can also be "bottom", etc
         theme(legend.key.size = unit(1, 'cm'),
               legend.key.height = unit(1, 'cm'), 
               legend.key.width = unit(1, 'cm'), 
               legend.text = element_text(size=10),
               legend.title = element_blank())
              # legend.title = element_text(size=14), 
             # scale_colour_discrete(labels = c('Women', 'Men'))

     ggsave(Fig, filename = figname, bg= "transparent")
     
     return(Fig)

}

box_simple_diffpoints <- function(df, working.col1, working.col2, working.col3, title="", axis1, axis2, figname="Fig.png", palette=cbPalette) {
     
      working.col1 <- df[, working.col1] #subsets values of variable 1
      working.col2 <- df[, working.col2] #subsets values of variable 2
      working.col3 <- df[, working.col3] #subsets values of variable 3


      Fig <- df %>%

      #Main structure
        ggplot(aes(x = working.col1, y = working.col2, color= "black", fill = working.col1)) + 
        scale_fill_manual(values=cbPalette) +
        scale_color_manual(values=cbPalette) +
        geom_boxplot(alpha=0.2) +
        geom_point(x=working.col1,
        # geom_point(position=position_jitterdodge(dodge.width=1), #the flaw with the jitterdodge approach is that some points remain on axis of boxplot and are thus in alpha=0.2 without color
                     alpha=1, aes(color=working.col3)) +

        #Theme
        theme_classic() +

        #Title
        ggtitle(title) +
        theme(plot.title = element_text(size = rel(2), hjust = 0.5)) +

        #Axes
         labs(x=axis1, y=axis2) +
         theme(axis.title.y = element_text(size = rel(1.5)),
          axis.title.x = element_text(size = rel(1.5)),
          axis.text.x = element_text(hjust = 0.5, size=13,color="black")) + 

         #Legend
         # theme(legend.position="none") + #can also be "bottom", etc
         theme(legend.key.size = unit(1, 'cm'),
               legend.key.height = unit(1, 'cm'), 
               legend.key.width = unit(1, 'cm'), 
               legend.text = element_text(size=10),
               legend.title = element_blank())
              # legend.title = element_text(size=14), 
             # scale_colour_discrete(labels = c('Women', 'Men'))

     ggsave(Fig, filename = figname, bg= "transparent")
     
     return(Fig)

}

box_groupped <- function(df, working.col1, working.col2, working.col3, title="", axis1, axis2, figname="Fig.png", palette=cbPalette) {
     
      working.col1 <- df[, working.col1] #subsets values of variable 1
      working.col2 <- df[, working.col2] #subsets values of variable 2
      working.col3 <- df[, working.col3] #subsets values of variable 3


      Fig <- df %>%

      #Main structure
        ggplot(aes(x = working.col1, y = working.col2, color= working.col3, fill = working.col3)) + 
        scale_fill_manual(values=cbPalette) +
        scale_color_manual(values=cbPalette) +
        geom_boxplot(alpha=0.2) +
        geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=1) +

        #Theme
        theme_classic() +

        #Title
        ggtitle(title) +
        theme(plot.title = element_text(size = rel(2), hjust = 0.5)) +

        #Axes
         labs(x=axis1, y=axis2) +
         theme(axis.title.y = element_text(size = rel(1.5)),
          axis.title.x = element_text(size = rel(1.5)),
          axis.text.x = element_text(hjust = 0.5, size=13,color="black")) + 

         #Legend
         # theme(legend.position="none") + #can also be "bottom", etc
         theme(legend.key.size = unit(1, 'cm'),
               legend.key.height = unit(1, 'cm'), 
               legend.key.width = unit(1, 'cm'), 
               legend.text = element_text(size=10),
               legend.title = element_blank())
              # legend.title = element_text(size=14), 
             # scale_colour_discrete(labels = c('Women', 'Men'))

     ggsave(Fig, filename = figname, bg= "transparent")
     
     return(Fig)


}

####2. RADAR PLOTS

radar <- function(title, df, color) {

  #Plot the radar chart 
  radar <- radarchart(df,

          #Axis with %
              axistype=1, #by default is empty. 1 means you end up with % on each line of the web. 
              axislabcol='grey', #color of the legend on the axis of %
          # caxislabels=c(), #Vector with the values to put on the axis (instead of 25% etc),

          #Title and labels
          title=title, vlcex=0.8, # vlcex: Group labels size

          #General layout: web
          cglty=1, #type of web line: 1 is continuous, 2 is disturbing, 3 is pointillé
          cglcol='grey', #web color
          cglwd=1, #net width

          #General layout: "stain"
          pcol=color, #line color for outline
          pfcol = scales::alpha(color, 0.5),
          plwd=2 #line weight
          )

          ggplotly(radar)

}