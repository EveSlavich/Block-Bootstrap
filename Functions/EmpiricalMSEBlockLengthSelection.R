which.min=function(x){
  ind=NA
  if(length(na.omit(x))>0){
  minx=min(x, na.rm=T)
  min.index = which(x==minx)
  ind = max(min.index)
  }
  ind
}

EmpiricalMSEBlockLengthSelection=function(x,y,m_x, m_y,  dat, Grid_space, shape, NBoot, sampling_type, Stat.function,lookuptables.folderpath,block_Ls ,boot.reps.of.Stat.function,tuning_block_length , type,subregion.division,...){
  
  sigma_stat_subregion = BlockBoot_apply_subregions(
    x = x ,
    y = y ,
    m_x = m_x ,
    m_y = m_y ,
    dat = dat,
    block_Ls = block_Ls,
    Grid_space = Grid_space,
    shape = shape,
    NBoot = NBoot,
    sampling_type = sampling_type,
    Stat.function = Stat.function,
    lookuptables.folderpath = lookuptables.folderpath,
   type=type,
   subregion.division=subregion.division,
    ...
  )

    Empirical.BIAS2 = matrix(nrow = length(block_Ls), ncol = length(block_Ls))
      Empirical.VAR = matrix(nrow = length(block_Ls), ncol = length(block_Ls))
	  Empirical.BIAS2.unscaled = matrix(nrow = length(block_Ls), ncol = length(block_Ls))
	  Empirical.VAR.unscaled = matrix(nrow = length(block_Ls), ncol = length(block_Ls))
  n_in_subregions=c()
  Lahiri_chosen_block_length = c()
  Lahiri_chosen_block_length.unscaled= c()
  Lahiri_chosen_block_length_iter= c()
  Lahiri_chosen_block_length_iter.unscaled= c()
  SE.estimate = c()
  SE.estimate.iter= c()
  SE.estimate.unscaled.iter= c()
    SE.estimate.unscaled= c()
  for (subregion in 1:c(m_x*m_y)){
    bins = break_into_subregions(x , y ,  m_y = m_y , m_x = m_x)
	if (subregion.division=="mutually.exclusive"){
    n_in_subregions [subregion] =  table(bins$bin)[levels(bins$bin)[subregion]]
	}else{
	n_in_subregions [subregion] =  length(x) - table(bins$bin)[levels(bins$bin)[subregion]]
	}
        if( is.na(n_in_subregions [subregion])){ n_in_subregions [subregion] = 0}
  }
  
  sigma_stat_hat  = lapply(boot.reps.of.Stat.function, sd.unlist)
  for (L in 1:length(block_Ls)) {
    for (plug_in in 1:length(block_Ls)) {
                      Empirical.BIAS2 [L, plug_in]       = ((mean(sqrt(n_in_subregions/ sum(n_in_subregions))*sigma_stat_subregion[L,], na.rm = T) -
                                    sigma_stat_hat[[plug_in]]) ^ 2   )
                       Empirical.VAR [L, plug_in]       =      var(sqrt(n_in_subregions/sum(n_in_subregions))*sigma_stat_subregion[L,], na.rm=T)
					   
					   Empirical.BIAS2.unscaled [L, plug_in]       = ((mean(sigma_stat_subregion[L,], na.rm = T) -
                                    sigma_stat_hat[[plug_in]]) ^ 2   )
                       Empirical.VAR.unscaled [L, plug_in]       =      var(sigma_stat_subregion[L,], na.rm=T)

    }
  }
   Empirical.MSE =Empirical.BIAS2+Empirical.VAR
    Empirical.MSE.unscaled =Empirical.BIAS2+Empirical.VAR
  if(!hasArg("tuning_block_length")){print("Error, needs tuning block length")}
  if (hasArg("tuning_block_length")) {
    for (tuning.ind in 1:length(tuning_block_length)){
      if(tuning_block_length[tuning.ind] %in% block_Ls == FALSE){
        print ("invalid tuning block length provided")
        tuning_block_length[tuning.ind]  = NA
        SE.estimate[tuning.ind] = NA
        Lahiri_chosen_block_length[tuning.ind] = NA
      }else{

            
      tuning_block_length_ = which(block_Ls == tuning_block_length[tuning.ind])
        SE.estimate[tuning.ind] = NA
        Lahiri_chosen_block_length[tuning.ind] = block_Ls[which.min(Empirical.MSE[, tuning_block_length_])]
		 Lahiri_chosen_block_length.unscaled[tuning.ind] = block_Ls[which.min(Empirical.MSE.unscaled[, tuning_block_length_])]
        iterated_Empirical_MSE_block_length = iteratively_plug_in (Empirical.MSE=Empirical.MSE, tuning_block_length_=tuning_block_length_, block_Ls=block_Ls)
       Lahiri_chosen_block_length_iter[tuning.ind] = iterated_Empirical_MSE_block_length
         Lahiri_chosen_block_length_iter.unscaled[tuning.ind] =   iteratively_plug_in (Empirical.MSE=Empirical.MSE.unscaled, tuning_block_length_=tuning_block_length_, block_Ls=block_Ls)                                                 
        try({
          SE.estimate[tuning.ind] = 
            sd.unlist(boot.reps.of.Stat.function[[paste0("blength", Lahiri_chosen_block_length[tuning.ind])]])
			SE.estimate.iter[tuning.ind] = 
            sd.unlist(boot.reps.of.Stat.function[[paste0("blength", Lahiri_chosen_block_length_iter[tuning.ind])]])
			SE.estimate.unscaled.iter[tuning.ind] = 
            sd.unlist(boot.reps.of.Stat.function[[paste0("blength", Lahiri_chosen_block_length_iter.unscaled[tuning.ind])]])
			SE.estimate.unscaled[tuning.ind] = 
            sd.unlist(boot.reps.of.Stat.function[[paste0("blength", Lahiri_chosen_block_length.unscaled[tuning.ind])]])
        })
      }
      }
    }  else{
        tuning_block_length = NA
        Lahiri_chosen_block_length = NA
        SE.estimate[tuning.ind] = NA
        print(
          "Warning, invalide tuning block length provided for Lahiri subregion method of block length selection"
        )
      }
    rownames(Empirical.MSE) = paste0("blocklength", block_Ls)
    colnames(Empirical.MSE) = paste0("tuninglength", block_Ls)
    list(
      EmpiricalMSE_block_size = Lahiri_chosen_block_length,
	  EmpiricalMSE_block_size.iter = Lahiri_chosen_block_length_iter,
	  EmpiricalMSE_block_size.iter.unscaled= Lahiri_chosen_block_length_iter.unscaled,
	  EmpiricalMSE_block_size.unscaled =	 Lahiri_chosen_block_length.unscaled,
      Empirical.MSE = Empirical.MSE,
Empirical.BIAS2 = Empirical.BIAS2,
Empirical.VAR = Empirical.VAR,
      SE.estimate.EmpiricalMSE = SE.estimate,
	  SE.estimate.EmpiricalMSE.iter = SE.estimate.iter,
  	  SE.estimate.EmpiricalMSE.unscaled.iter = SE.estimate.unscaled.iter,
  	  SE.estimate.EmpiricalMSE.unscaled =  SE.estimate.unscaled,
      sigma_stat_subregion = sigma_stat_subregion,
      n_in_subregions=n_in_subregions,
      sigma_stat_hat =sigma_stat_hat,
      block_Ls=block_Ls,
      tuning_block_length=tuning_block_length
    )

}


      iteratively_plug_in = function(Empirical.MSE, tuning_block_length_,block_Ls){
mins = apply(Empirical.MSE,2,which.min)
block_length_index = c()
for (i in 1:5){
block_length_index[i] = mins[tuning_block_length_]
tuning_block_length_= mins[tuning_block_length_]
}
block_length_index = max(block_length_index[3:5])
block_Ls[block_length_index]
}