
Generate_dat_orig_sites=function(dat_poa_all,effect,subsample_, range_autocor1=20000,auto_cor_in_cov=TRUE,seed, ...){
  #if(hasArg(seed)){set.seed(seed)}
  set.seed(seed)
  x = dat_poa_all$x
  y = dat_poa_all$y
  #scale1=30000
  sites = cbind(x,y)
  
  
  spec_intercept = runif(1)
  climate_slope = runif(1)
 
  if(hasArg(subsample_)==FALSE){subsample_=nrow(sites)}
  subsample_sites=sample(nrow(sites),subsample_,replace=F)
  sites=sites[subsample_sites,]
  nsite=nrow(sites)
 
  if(auto_cor_in_cov==TRUE){
  mint5 = RandomFields::RFsimulate(model = RandomFields::RMexp(var=1, scale=range_autocor1), x=sites[,"x"], y=sites[,"y"],grid=F, n=1)
  } else {
  mint5 = runif(nrow(sites))  
  }
  X.alt = data.frame(mint5)
  mint5  = as.vector(data.frame(mint5)[[1]])
  xB=matrix(nrow=nsite, ncol=nspec)
  xB_true=xB
  xB_false=xB

    for (i in 1:nsite){
      for (j in 1:nspec){
        xB_true[i,j] = spec_intercept[j] +  climate_slope*mint5[i]
		xB_false[i,j] = spec_intercept[j] 
      }
    }
	  if(effect==TRUE){
	  xB=xB_true
  }
	  if(effect==FALSE){
	   xB=xB_false
  }
  
  

  #xB=as.vector(xB)[subsample_sites]
 
  dat_short=data.frame(mint5,sites)
  dat_long=data.frame(mint5,sites)
  dat.run = data.frame(mint5)
  
  
  dat.run=data.frame(dat.run,lat=sites[,"y"], long=sites[,"x"])

  names(dat.run)[1]="mint5"
  
  return(list(xB=xB, X.alt=X.alt , dat_short=dat_short , dat_long=dat_long, dat.run=dat.run, mint5=mint5, sites=sites, xB_false=xB_false,xB_true=xB_true))
}

auto_cor_in_cov = TRUE
effect=FALSE
nspec=1
dat_poa_pluslayers$x = dat_poa_pluslayers$MGAEasting.y
dat_poa_pluslayers$y = dat_poa_pluslayers$MGANorthin
Gen.dat = Generate_dat_orig_sites(dat_poa_pluslayers, effect=effect, scale1=100,subsample_=nsites,auto_cor_in_cov=auto_cor_in_cov, seed = nsites+1)
 xB = Gen.dat$xB
 xB_true = Gen.dat$xB_true
  xB_false = Gen.dat$xB_false
 X.alt = Gen.dat$X.alt
 dat_short = Gen.dat$dat_short
 dat_long = Gen.dat$dat_long
 dat.run  = Gen.dat$dat.run
 mint5  = Gen.dat$mint5
 sites  = data.frame(Gen.dat$sites)
 #save(xB,X.alt,dat_short,dat_long,mint5,sites,xB_true,xB_false, dat.run, file=paste0("Gen.dat_orig.SampleSize",sample_size,auto_cor_in_cov_name,".RData"))