# ASynuclein_APD repository

## Description
The code in this repository allowed us to obtain the values published in "Cerebrospinal fluid α-synuclein seeding activity in patients with atypical parkinsonian disorders". It is meant mostly for others to be able to review the steps followed in the analysis (eg: testing of assumptions, evaluation of model diagnostics, model selection) as well as the exact R functions used in the code. If necessary, it can be used as a way to share further information about the analysis that was not included in the original manuscript. 

## Content
The repository contains the following files in addition to the README:

1. Functions.R: stereotypic calls for basic functions that get repeated throughout the code. Also defines arbitrary lists that need to stay consistent such as palettes, etc. 
2. DataQC_PROCESSING.R: Function that creates the dataframe to be used for analysis. Creates new variables, standardizes the format, etc. 
3. Analysis.R: Actual analysis file, including any step toward data exploration (data distribution, etc). Contains information that could not be included in the manuscript, but also can be updated if necessary upon feedback in the future. Creates the Analysis_Terminal_Output.pdf file when run in bash using Rscript function. 
4. Analysis_Rmarkdown.Rmd: Rmarkdown of Analysis.R which creates the Analysis_Rmarkdown.pdf file.
5. Other files: figures and Rplots created with Analysis.R. 

## Contact information:
If further information is required, or for advice regarding the statistics, please reach out to the corresponding author Maria Carmela Tartaglia (carmela.lastname@uhn.ca) and Github repository owner Chloe Anastassiadis (firstname.lastname@mail.utoronto.ca).

