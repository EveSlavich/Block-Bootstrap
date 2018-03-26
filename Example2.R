## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----results = "hide", echo=TRUE, message=FALSE--------------------------
load("PlayData_for_Examples.RData") ### load the data for the examples
source("LoadFunctions.R") #### Source the functions required to run the block bootstrap


#Install packages if required
if(!require(mvabund)) { install.packages("mvabund", repos = "http://cran.us.r-project.org"); 
require(mvabund) }
library(mvabund)



## ----pressure, echo=TRUE-------------------------------------------------
set.seed(42)

lookuptables.folderpathname = "LookupTables/" 

######## Get a bootID matrix (takes a minute or two)
BootID.example2 = BlockBootID(x = x ,
                              y = y,
                              block_Ls = 0.1,
                              NBoot = 500,
                              Grid_space = 0.01,
                              lookuptables.folderpath =  lookuptables.folderpathname)

## ----variogram, echo=TRUE------------------------------------------------
responseMultiSpecies=mvabund(multispecies_dat[,1:20]) #20 species multivariate reponse
mod.1  = manyglm(responseMultiSpecies~temperature, data = multispecies_dat,family="binomial")
mod.2  = manyglm(responseMultiSpecies~temperature*treatment, data = multispecies_dat,family="binomial")
anova.results_naive = anova(mod.1, mod.2)
anova.results = anova(mod.1, mod.2, bootID=BootID.example2, resamp="case")
?anova.manyglm


