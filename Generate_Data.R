#Simulate multispecies_dat
load(file="/home/z3354192/hdrive/z3354192/MovingBlockBootStrap/Resubmission/New Code/Motivating Example/Sig.RData")
library(RandomFields)
library(spatstat)
library(mvabund)
Thomas1  = rThomas(0.03, 2, 3, win = owin(c(0, 100), c(0, 100)))[1:500, ]

x =Thomas1$x
y =Thomas1$y
#species responses and covariates   

nsites=500
nspec=22 #number of species
range_autocor1 = c(20) #range parameter for spatial covariate x1
range_autocor2 = c(10)  #range parameter for spatial random effect x2 to induce spatial autocorrelation
range_autocor3 = c(10)  #range parameter for treatment


beta1.assemblage = rnorm(nspec,1,1) #coeff of x1
beta2.assemblage = rnorm(nspec,1,1) #coeff of x2
range2.assemblage = sapply(rnorm(nspec,range_autocor2,1), function(x){max(x, 0.000001)})
sigma_e =1

temperature = RandomFields::RFsimulate(
  model = RandomFields::RMexp(var = 1, scale = range_autocor1),
  x = x,
  y = y,
  grid = F,
  n = 1
)$variable1

#simaulte presence absence data
response = matrix(nrow=nsites, ncol=nspec)



library(MASS)
err_spec = mvrnorm(n=1, mu=rep(0,nspec), Sigma=Sig[1:nspec,1:nspec])


for (speci in 1:nspec){
  x2 = RandomFields::RFsimulate(
    model = RandomFields::RMexp(var = 1, scale = range2.assemblage[speci]),
    x =x,
    y = y,
    grid = F,
    n = 1
  )$variable1
  
  #uncorrelated species
  err = rnorm(nsites, 0,1)
  beta0.assemblage =-2
  linear.response = beta0.assemblage + beta1.assemblage[speci]  %*% x1 + beta2.assemblage[speci]  %*%    x2 +sigma_e * err +err_spec[speci]
  
  count.response = rnbinom(n = nsites, size=1, mu=exp(linear.response) )
  count.response = round(exp(linear.response) )
  response[,speci] = count.response
}

response = mvabund(response) 

treatment  = factor(cut_number(
RandomFields::RFsimulate(
    model = RandomFields::RMexp(var = 1, scale = range_autocor3),
    x = x,
    y = y,
    grid = F,
    n = 1
  )$variable1,2),labels=c("control","treatment"))

multispecies_dat = cbind(response, treatment, temperature)
save(multispecies_dat, file="/home/z3354192/hdrive/z3354192/MovingBlockBootStrap/Resubmission/Code for Users/multispecies_dat.RData")
