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

## CHOOSEX2.FUNC() ##
chooseX2.func <- function(table) {
  total <- table[1] + table[2] + table[3] + table[4]
  cell1 <- (table[1] + table[2]) * (table[1] + table[3])
  cell2 <- (table[1] + table[2]) * (table[2] + table[4])
  cell3 <- (table[1] + table[3]) * (table[3] + table[4])
  cell4 <- (table[3] + table[4]) * (table[2] + table[4])
  if ((cell1/total) <5 | (cell2/total) <5 | (cell3/total) <5 | (cell4/total) <5 ) {
  return("fisher") }
  	else {
  	return("chisquare")
  	}
  }           


## CREATE_BEAUTIFUL_RADARCHART() ##
create_beautiful_radarchart <- function(data, color= "red", plty=1, vlabels = colnames(data), vlcex= 0.7, caxislabels= NULL, title= NULL, ...){
  radarchart(data,
  axistype = 1,
  pcol= color, #Color of the lines of the radar plot
  pfcol = scales::alpha(color, 0.3), #Fill color, based on color
  plwd= 1, #line weight
  plty= plty, #line type for the radar
  cglcol = "#000000", #line color for the grid
  cglty = 1,#line type for the grid
  cglwd = 0.8,#line width for the grid
  axislabcol = "grey",  #color of axis label and numbers. Default is “blue”.
  vlcex = vlcex, #Label size: add in function call.
  vlabels = vlabels, #Column names for variable names. No need to change
  caxislabels= caxislabels, #Character vector to be used as labels on the center axis
  title= title, #title of the plot
  ...
  )
}

