## ----setup, include=FALSE------------------------------------------------

knitr::opts_chunk$set(echo = TRUE,tidy.opts=list(width.cutoff=60))


## ----results = "hide", echo=TRUE, message=FALSE--------------------------
load("PlayData_for_Examples.RData") ### load the data for the examples
source("LoadFunctions.R") #### Source the functions required to run the block bootstrap


#Install packages if required
if(!require(pROC)) { install.packages("pROC", repos = "http://cran.us.r-project.org"); require(pROC) }

library(pROC)


## ----pressure, echo=TRUE-------------------------------------------------
set.seed(42)

fit1 = glm (PresAbsresponse ~ temperature, data = dat, family="binomial")

auc(response=dat$PresAbsresponse, predictor= predict(fit1))



## ----funct, echo=TRUE----------------------------------------------------

GetAUC = function(dat){
fit1 = glm ("PresAbsresponse~temperature", data=dat, family="binomial")
AUC = auc(response=dat$PresAbsresponse, predictor= predict(fit1))
;
AUC
}


## ----variogram, echo=TRUE------------------------------------------------


ini.vals <- expand.grid(seq(0,10,l=100), seq(0,1,l=100))


variogram_block_length = select_b_length_off_variogram_envelope(x,y,resids = fit1$resid,
max.dist = 1.2, breaks=c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.5),
ini=ini.vals)
variogram_block_length


## ----  echo=TRUE---------------------------------------------------------

lookuptables.folderpathname = "LookupTables/" 


## ----  echo=TRUE, warning=FALSE, message=FALSE, eval=T-------------------

Results.example4 = BlockBootApply (x = x ,y = y ,block_Ls = c(0,0.05,0.1,0.2), Grid_space = 0.01 ,
dat = dat , Stat.function = GetAUC,  tuning_block_length = 2, NBoot = 500,
method.block.length.select = "Lahiri", type="SE", lookuptables.folderpath=lookuptables.folderpathname)


## ---- echo=TRUE----------------------------------------------------------

Results.example4$Empirical.MSE

## ---- echo=TRUE----------------------------------------------------------
Results.example4$Lahiri_block_size

## ---- echo=TRUE----------------------------------------------------------
Results.example4$SE.estimate

