
These files accompany the paper: "Viral load dynamics of SARS-CoV-2 Delta and Omicron variants
following multiple vaccine doses and previous infection".
This files contains code that produces the main figures (1-3) and main regression results.



System requirements:
This code ran on R version 4.0.3 (2020-10-10) using Rstudio 1.4.1103 and mingw32 operating system.

R packages requirement: 
mltools, ggpubr, ggplot2, gridExtra, grid, jtools, ggstance, lme4, quantreg, forcats,
tidyquant, glmnet, VGAM, scales, dplyr, mvtnorm, lubridate ,DescTools, MASS.
Please type   source('libraries.R')   to upload them (after installation). 

Data availability:
Since the dataset is individual-level, it cannot be shared, due to privacy and Israel Ministry of Health's restrictions.

Code Availability:
Code for producing main text figures and linear regression results is openly available. 
Full code's version will be shared upon request. 


Running time:
Running main.R on the full dataset took approximately 1 min.



Files include:

libraries.R: load R packages
main.R: main file that runs all functions
prepData.R: preliminary data preperation for further analysis
prepFigTable1.R: process table for Fig.1 
prepFigTable2.R: process table for Fig.2 
prepFigTable3.R: process table for Fig.3
code_output_.xlsx: table output for figures 1-3, each tab represents a different figure.
ctTable.RData: a data table example (with 100 rows). Due to sparsity of the data, the code will not run properly on that table.


Extended tables - SI: a folder that contains all SI tables in xlsx format.
Figures: a folder that contains all main text and SI figures in JPEG format. 
 