dat_poa_pluslayers$x = dat_poa_pluslayers$MGAEasting.y
dat_poa_pluslayers$y = dat_poa_pluslayers$MGANorthin

Generate_dat = function(dat,subsample_, formula1, family1, spec, seed){
set.seed(seed)
x = dat$x
y = dat$y
sites = cbind(x,y)
names(dat)[which(names(dat)==spec)]="binary.response"
ft0 = glm(formula1, family=family1, data=dat)
beta0 = ft0$coeff

subsample_sites=sample(nrow(sites),subsample_,replace=F)
sites=sites[subsample_sites,]
mint5 = dat$mint5[subsample_sites]
dat.run = data.frame(sites,mint5)
names(dat.run) = c("long","lat","mint5")
list(dat.run,beta0)
}

dat.run0 = Generate_dat(
dat = dat_poa_pluslayers,
subsample_=nsites,
seed=nsites+1,
formula1 =as.formula("binary.response~1"),
family1 = "binomial",
spec = "Cymbopogon.refractus")

dat.run = dat.run0[[1]]
xB=rep(dat.run0[[2]], nrow(dat.run))