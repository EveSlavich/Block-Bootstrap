## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----results = "hide", echo=TRUE, message=FALSE--------------------------
load("PlayData_for_Examples.RData") ### load the data for the examples
source("LoadFunctions.R") #### Source the functions required to run the block bootstrap


#Install packages if required
if(!require(pROC)) { install.packages("pROC", repos = "http://cran.us.r-project.org"); require(pROC) }

library(pROC)


## ----pressure, echo=TRUE-------------------------------------------------
set.seed(42)

fit1 = glm (PresAbsresponse ~ temperature, data = dat, family="binomial")


## ----funct, echo=TRUE----------------------------------------------------
GetAUC = function(dat){
fit1 = glm ("PresAbsresponse~temperature", data=dat, family="binomial")
AUC = auc(response=dat$PresAbsresponse, predictor= predict(fit1))
;
AUC
}


## ----  echo=TRUE---------------------------------------------------------

lookuptables.folderpathname = "LookupTables/" 


## ----  echo=TRUE, warning=FALSE, message=FALSE, eval = TRUE--------------


batchID = 1
assign (paste0("Results.example4.batch", batchID ), 
        BlockBootApply (x = x ,y = y ,block_Ls = c(0,0.05,0.1,0.2), Grid_space = 0.01 , 
                        dat = dat, Stat.function = GetAUC, tuning_block_length = 2, 
                        NBoot = 250, method.block.length.select = "Lahiri",type="SE",
                        lookuptables.folderpath=lookuptables.folderpathname))



## ----  echo=TRUE, warning=FALSE, message=FALSE, eval = TRUE--------------
batchID = 2
assign (paste0("Results.example4.batch", batchID ), 
        BlockBootApply (x = x ,y = y ,block_Ls = c(0,0.05,0.1,0.2), Grid_space = 0.01 , 
                        dat = dat, Stat.function = GetAUC, tuning_block_length = 2, 
                        NBoot = 250, method.block.length.select = "Lahiri",type="SE",
                        lookuptables.folderpath=lookuptables.folderpathname))

## ----  echo=TRUE, warning=FALSE, message=FALSE, eval = TRUE--------------

batchIDs=1:2

list_of_boot_batches = lapply(paste0("Results.example4.batch", batchIDs), get)

Results.example4.all = CombineBatchesofBootstraps(list_of_boot_batches=
list_of_boot_batches )


## ---- echo=TRUE, eval = TRUE---------------------------------------------

Results.example4.all$Empirical.MSE

## ---- echo=TRUE, eval = TRUE---------------------------------------------
Results.example4.all$Lahiri_block_size

## ---- echo=TRUE, eval = TRUE---------------------------------------------
Results.example4.all$SE.estimate

