## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,tidy.opts=list(width.cutoff=60))
setwd("\\\\myfiles.unsw.edu.au@SSL/DavWWWRoot/Staff078/z3354192/MovingBlockBootstrap/For Publication/Code for users")

## ----results = "hide", echo=TRUE, message=FALSE--------------------------
load("PlayData_for_Examples.RData") ### load the data for the examples

#Let's explore the data visually 
plot(response~temperature, data=dat)

#add the line of best fit
fit1 = lm ("response~temperature", data=dat)
abline(fit1,col="red")
summary(fit1)

#view the spatial pattern of response, covariate and residuals
library(ggplot2)
qplot(x, y, colour = dat$response)
qplot(x, y, colour = dat$temperature)
qplot(x, y, colour = fit1$residuals) #seems to be some spatial autocorrelation


source("LoadFunctions.R") #### Source the functions required to run the block bootstrap




## ----pressure, echo=TRUE-------------------------------------------------
set.seed(42)


fit1 = lm (response ~ temperature, data = dat)
summary(fit1)
abline(fit1,col="red")

## ----variogram, echo=TRUE------------------------------------------------

ini.vals <- expand.grid(seq(0,10,l=100), seq(0,1,l=100)) #inital values for variofit
variogram_block_length = select_b_length_off_variogram_envelope(x,y,resids = fit1$resid, 
max.dist = 0.4, breaks=c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.5),ini=ini.vals)


variogram_block_length


## ----funct, echo=TRUE----------------------------------------------------



BetaCoeff = function(dat){
  fit1 = lm ("response~temperature", data=dat)
  Beta = fit1$coefficients["temperature"]
  Beta
}


## ----  echo=TRUE---------------------------------------------------------

lookuptables.folderpathname = "LookupTables/" 


## ----  echo=TRUE, warning=FALSE, message=FALSE---------------------------
Results = BlockBootApply (x = x ,y = y ,block_Ls = c(0,0.05,0.1, 0.2), Grid_space = 0.01 ,                           dat = dat ,
                          Stat.function = BetaCoeff,  tuning_block_length = 4, 
                          NBoot = 500, method.block.length.select = "Lahiri",
                          type="SE", lookuptables.folderpath=lookuptables.folderpathname)



## ---- echo=TRUE----------------------------------------------------------

Results$Empirical.MSE

## ---- echo=TRUE----------------------------------------------------------
Results$Lahiri_block_size

## ---- echo=TRUE----------------------------------------------------------
Results$SE.estimate

