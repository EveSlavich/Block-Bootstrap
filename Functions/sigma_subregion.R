sigma_subregion = function(x,y,m_x, m_y, dat,block_L, scale1,function_to_repeat, function_on_reps,NBoot, ...){
  start.time <- Sys.time()
  # break region into subregions
  # bootstrap on subregions
  # calculate parameter
  # calculate RMSE of parameters
  #bins = break_into_subregions (x=x, y=y, m_y=m_y ,m_x=m_x)
  #bins1 = unique(bins[,3])
  sigma_subregion =c()
  x_subrange = diff(range(x))/m_x
  y_subrange =  diff(range(y))/m_y
  reg_size=diff(range(x))*diff(range(y))
  subreg_size = reg_size - x_subrange*y_subrange
  n_blocks = round(subreg_size/(block_L^2))
  print(paste0("nblocks ",n_blocks))
  if (block_L>0){

    if(file.exists(paste0("DataMod/lookup_subregions_mx_",m_x,"my_",m_y,"block_L",block_L,".RData"))){
      load(paste0("DataMod/lookup_subregions_mx_",m_x,"my_",m_y,"block_L",block_L,".RData"))
    }else{
    lookup = create_moving_block_lookup_subregions(x=dat$x, y=dat$y,block_L=block_L, scale1=scale1, m_x=m_x, m_y=m_y)
    save(lookup, file=paste0("DataMod/lookup_subregions_mx_",m_x,"my_",m_y,"block_L",block_L,".RData"))
    }
    lookup_bins = lookup$bin.ind
    lookup_tab = lookup$ind
    site_bins = lookup$site_bins
    nref= lookup$nref
  }
  end.time <- Sys.time()
  time.taken.tonow <- start.time - end.time
  print("Lookuptable created")
  print(time.taken.tonow)
  
  nbins = length(levels(unlist(lookup_bins)))
  print(paste("nbins", nbins))
  for (i in 1:nbins){ 
    if(i ==1){
      start.time <- Sys.time()
    }
    #subset dat
    print(i)
    bin_i = levels(unlist(lookup_bins))[i]
    dat_i = dat[which(site_bins != bin_i),]
    nref_i = sum(nref[which(nref$ref.bin!=bin_i),"nref"])
    #create lookuptable
    if (block_L>0){
      #resample blocks
      samples = resample_blocks_subregion( NBoot = NBoot, lookup_table = lookup_tab, block_L=block_L,  bin_i=bin_i, site_bins=site_bins, nref = nref_i, n_blocks=n_blocks)
      
    }
    if(block_L==0){
      Nrow=dim(dat_i)[1]
      samples = matrix(sample(c(1:Nrow),size=Nrow*NBoot, replace=TRUE),nrow= Nrow, ncol= NBoot)
    }
    print(paste("created samples subregion",i))
    
    
    #samples
    
    results = bootstrap_wrapper2(dat, function_to_repeat =function_to_repeat, new_sample = samples, NBoot=NBoot, function_on_reps=function_on_reps,...)
    sigma_subregion[i] = results[[1]]
    if(i ==1){
      end.time <- Sys.time()
      time.taken.tonow <- start.time - end.time
      print(paste("the 1st subregion time was",time.taken.tonow ))
    }
  }
  ;sigma_subregion
  
}

MSE_sigma_subregions = function (sigma_subregion,sigma_guess){
  mean((sigma_subregion-sigma_guess)^2)
}
VAR_sigma_subregions = function (sigma_subregion,sigma_guess){
  var((sigma_subregion))
}
MEAN_sigma_subregions = function (sigma_subregion,sigma_guess){
  mean(sigma_subregion)
}
BIAS_sigma_subregions = function (sigma_subregion,sigma_guess){
  mean(sigma_subregion)-sigma_guess
}
min_MSE = function(x){
which(x==min(x, na.rm=T))
}
