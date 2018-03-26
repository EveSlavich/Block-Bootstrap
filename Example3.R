## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,tidy.opts=list(width.cutoff=40))


## ----results = "hide", echo=TRUE, message=FALSE--------------------------
load("PlayData_for_Examples.RData") ### load the data for the examples
source("LoadFunctions.R") #### Source the functions required to run the block bootstrap


#Install packages if required
if(!require(mvabund)) { install.packages("mvabund", repos =
"http://cran.us.r-project.org"); require(mvabund) }

library(mvabund)



## ----pressure, echo=TRUE-------------------------------------------------
set.seed(42)
treatment =rbinom(500,1,0.5)
multispecies_dat = data.frame(multispecies_dat,treatment)

## ----funct, echo=TRUE----------------------------------------------------

MultSpeciesTestStat = function(dat){
responseMultiSpecies = mvabund(dat[,1:20])

null.model = manyglm (responseMultiSpecies~1, data=dat, family="binomial")
alt.model = manyglm (responseMultiSpecies~temperature*treatment, data=dat, family="binomial")
Test.stat = sum(alt.model$two.loglike) - sum(null.model$two.loglike)
;
Test.stat
}

MultSpeciesTestStat (multispecies_dat)


## ----  echo=TRUE---------------------------------------------------------

lookuptables.folderpathname = "LookupTables/" 


## ----  echo=TRUE---------------------------------------------------------


Results.Example3 = BlockBootApply (x = x ,y = y ,block_Ls = c(0,0.05,0.1,0.2), Grid_space = 0.01 ,
dat = multispecies_dat , Stat.function = MultSpeciesTestStat,  tuning_block_length = 2, NBoot = 2, 
method.block.length.select = "Lahiri",  type="SE", lookuptables.folderpath=lookuptables.folderpathname)


L.example3  = Results.Example3$Lahiri_block_size
#The chosen block size parameter
L.example3


## ---- echo=TRUE----------------------------------------------------------
BootID.example3 = BlockBootID(x = x ,
y = y,
block_Ls = L.example3,
NBoot = 500,
Grid_space = 0.01,
lookuptables.folderpath =  lookuptables.folderpathname, shape="disc")

#######

responseMultiSpecies = mvabund(multispecies_dat[,1:20])

mod.full  = manyglm(responseMultiSpecies~temperature*treatment, data = multispecies_dat,
family="binomial")

anova(mod.full, bootID = BootID.example3, resamp="case")

##########

