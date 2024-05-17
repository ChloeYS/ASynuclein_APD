# FILENAME: Functions.R

#Last updated 16 May 2024
##Only kept functions useful for the APD manuscript


# USEFUL LIBRARIES #
library(tidyverse) #https://tidyverse.tidyverse.org/


# DEFINE REFERENCES LISTS #

#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
#https://www.colorhexa.com/999999
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #colorblind-friendly palette. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette_RTQUIC <- c("#E69F00", "#999999") #colorblind-friendly palette. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette_DX_APD <- c("#56B4E9", "#CC79A7") #colorblind-friendly palette. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette



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

                

