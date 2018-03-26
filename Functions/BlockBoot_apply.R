
BlockBootApply = function (x ,
                           y,
                           block_Ls,
                           Stat.function,
                           dat,
                           NBoot = 500,
                           NBoot.subregions = 200,
                           Grid_space,
                           type = "SE",
                           method.block.length.select = NA,
                           lookuptables.folderpath = NA,
                           Stat.function.optimising = NA,
                           m_x = 3 ,
                           m_y = 3 ,
                           shape = "disc",
                           sampling_type = "area",
                           long_format_required = FALSE,
                           tuning_block_length = c(),
                           breaks,
                           ini.variogram.pars, 
                           max.dist,
                           plug_in_tuning_blocksize1,
                           plug_in_tuning_blocksize2,
                           V,
                           subregion.division = "mutually.exclusive",
                           ...) {
 
  boot.reps.of.Stat.function = list() #initialise place to store results
  boot.reps.of.Stat.function.optimising = list() 
  results = list()
  if(is.na(Stat.function.optimising)){Stat.function.optimising = Stat.function}
  #set the spacing between sampling points
  if (missing(Grid_space)) {
 Grid_space = min(subset(block_Ls, block_Ls>0)) / 3
  }
  
  
  for (L in 1:length(block_Ls)) {
    block_L = block_Ls[L]
    if (block_L > 0) {
      if (is.na(lookuptables.folderpath) == FALSE) {
      print("has foldername, checking for lookuptable")
        #check if a lookup table has been created to speed things up, if it has load it, and check the lookup table is for data with same x,y, coordinates
        load_file_if_exists(
          paste0(
            lookuptables.folderpath,
            "lookup_table",
            "_L",
            block_L,
            "_grid_space_",
            Grid_space,
            "_sampling_type_",
            "sites",
            "_",
            shape,
            ".RData"
          )
        )
        if (exists("lookup.coords")) {
          if (identical (lookup.coords$lookup.x , x) == FALSE |
              identical (lookup.coords$lookup.y , y) == FALSE) {
            rm(lookup_table, lookup.coords, envir = .GlobalEnv)
            print("loaded lookup, lookup coords don't match, creating new")
			      new_sample = resample_blocks_by_area(
        x = x,
        y = y,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = sampling_type,
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
		#  ...
      )###will create a lookup table
          }else{
                      print("using lookup")
		        new_sample = resample_blocks_by_area(
        x = x,
        y = y,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = sampling_type,
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
		lookup_table=lookup_table,
      #  ...
      )#will use existing lookuptable}
        }
      } else{print("did not load")
	  new_sample = resample_blocks_by_area(
        x = x,
        y = y,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = sampling_type,
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
		#  ...
      )###will create a lookup table
	  }
      
      #Create the new sample (BlockBootID)

      
    }

      if (is.na(lookuptables.folderpath)) {	
print("code under developement")
	}
	}
    #If block_L =0, do an iid bootstrap
    if (block_L == 0) {
      new_sample = list()
      for (i in 1:NBoot) {
        new_sample [[i]] = sample(1:length(x), size = length(x), replace = T)
      }
    }
    if (long_format_required == TRUE) {
      new_sample =lapply(1:length(new_sample),function(i){BlockResampleData_longformat(new_sample[[i]], nspec=nspec)})

    }

    boot.reps.of.Stat.function.optimising[[L]] = bootstrap_wrapper(
      dat = dat,
      function_to_repeat = Stat.function.optimising,
      new_sample = new_sample,
      NBoot = NBoot,
	  type=type,
      ...
      #DELETE THIS LINE
      # family = family,
      #  formula1 = formula1,
      # formula2 = formula2
    )
     if (identical (Stat.function, Stat.function.optimising) == FALSE){
    boot.reps.of.Stat.function[[L]] = bootstrap_wrapper(
      dat = dat,
      function_to_repeat = Stat.function,
      new_sample = new_sample,
      NBoot = NBoot,
	  type=type,
       ...
    #DELETE THIS LINE
     # family = family,
    #  formula1 = formula1,
     # formula2 = formula2
    )
    
    }else{
     boot.reps.of.Stat.function[[L]] =  boot.reps.of.Stat.function.optimising[[L]]
    }
  }
  names(boot.reps.of.Stat.function.optimising) = paste0("blength", block_Ls)
  if(identical(Stat.function.optimising, Stat.function)==FALSE){
  names(boot.reps.of.Stat.function) = paste0("blength", block_Ls)
  results$boot.reps.of.Stat.function = boot.reps.of.Stat.function
  }else{
    boot.reps.of.Stat.function = boot.reps.of.Stat.function.optimising
  }
  Stat = Stat.function(dat, ...)
  k = length(Stat)
  Variogram.block.size = NA
  results$boot.reps.of.Stat.function.optimising = boot.reps.of.Stat.function.optimising
  results$Stat = Stat

  
  #######RUN BLOCK LENGTH SELECTION
  if ("Variogram" %in% method.block.length.select) {
    resids = variogram.function(dat,...)
    Variogram.block.size = select_b_length_off_variogram_envelope(
      x = x,
      y = y,
      resids = resids,
      breaks = breaks,
      ini = ini.variogram.pars,
      plot = FALSE,
      max.dist=max.dist,
      ...
    )

    results$Variogram.block.size =round_to_vector(Variogram.block.size[[1]],block_Ls)
    results$variogram_range_parameter = Variogram.block.size[[2]]
  }
  
  Plug.In.block.size = NA
  results$Plug.In.block.size = NA
  if ("PlugIn" %in% method.block.length.select) {
    ####CODE THIS SECTION
    Plug.In.block.size = Select.Plug.In.block.size(plug_in_tuning_blocksize1 = plug_in_tuning_blocksize1,
                                                   plug_in_tuning_blocksize2 = plug_in_tuning_blocksize2, 
                                                   boot.reps.of.Stat.function = boot.reps.of.Stat.function,
                                                   V=V)
    results$Plug.In.block.size =round_to_vector(Plug.In.block.size,block_Ls)
  }
 
  if ("EmpiricalMSE" %in% method.block.length.select) {
    ######RUN LAHIRI BL SELECTION
    if(is.na(Variogram.block.size)==FALSE){
      if(is.na(Plug.In.block.size)==FALSE){
        tuning_block_length1 = c(
          tuning_block_length,
          results$Variogram.block.size,
          results$Plug.In.block.size
        )
        
      }else{
        tuning_block_length1 = c(
          tuning_block_length,
          results$Variogram.block.size
    )
      }}else{
      
        if(is.na(Plug.In.block.size)==FALSE){
          tuning_block_length1 = c(
            tuning_block_length,
            results$Plug.In.block.size
          )
          
        }else{
          tuning_block_length1 = c(
            tuning_block_length
          )
        }
      }
      
    

    results.1 = EmpiricalMSEBlockLengthSelection(
      x = x,
      y = y,
      m_x = m_x,
      m_y = m_y,
      dat = dat,
      Grid_space=  Grid_space,
      shape = shape,
      NBoot =  NBoot.subregions,
      sampling_type = sampling_type,
      Stat.function = Stat.function.optimising,
      lookuptables.folderpath = lookuptables.folderpath,
      block_Ls = block_Ls,
      boot.reps.of.Stat.function = boot.reps.of.Stat.function.optimising,
      tuning_block_length = tuning_block_length1,
	  type=type,
	   subregion.division = subregion.division,
      ...
    )
    if( "Variogram.block.size" %in% method.block.length.select & "Plug.In.block.size" %in% method.block.length.select){
    names(results.1$EmpiricalMSE_block_size) = 
      paste0("tuning = ",c(tuning_block_length,
      "Variogram.block.size",
      "Plug.In.block.size"
    ))
    names(results.1$SE.estimate.EmpiricalMSE) = 
      paste0("tuning = ",c(tuning_block_length,
                           "Variogram.block.size",
                           " Plug.In.block.size"
      ))
      }else if("Variogram.block.size" %in% method.block.length.select & !("Plug.In.block.size" %in% method.block.length.select)){
      names(results.1$EmpiricalMSE_block_size) = 
      paste0("tuning = ",c(tuning_block_length,
      "Variogram.block.size",
         ))
    names(results.1$SE.estimate.EmpiricalMSE) = 
      paste0("tuning = ",c(tuning_block_length,
                           "Variogram.block.size",
                           
      ))
      }else if(!("Variogram.block.size" %in% method.block.length.select) & "Plug.In.block.size" %in% method.block.length.select){
      names(results.1$EmpiricalMSE_block_size) = 
      paste0("tuning = ",c(tuning_block_length,
      " Plug.In.block.size",
         ))
    names(results.1$SE.estimate.EmpiricalMSE) = 
      paste0("tuning = ",c(tuning_block_length,
                            "Plug.In.block.size",
                           
      ))
      }else {
      names(results.1$EmpiricalMSE_block_size) = 
      paste0("tuning = ",c(tuning_block_length
         ))
    names(results.1$SE.estimate.EmpiricalMSE) = 
      paste0("tuning = ",c(tuning_block_length
                          
                           
      ))
      }
       

    results=c(results.1, results)
  }
  
  
  if (type == "Pval") {
    pval_block_sizes = c()
    for (block_size in 1:length(block_Ls)) {
      distribution_of_T = unlist(boot.reps.of.Stat.function.optimising[[block_size]])
      pval_block_sizes[block_size] = (length(which(distribution_of_T > Stat)) +
                                        1) / (length(distribution_of_T) + 1)
    }
    results$pval_block_sizes = pval_block_sizes
    names(pval_block_sizes) = paste("block size", block_Ls)
    if ("EmpiricalMSE" %in% method.block.length.select) {
      ind = c()
      pval_EmpiricalMSE = c()
	   pval_EmpiricalMSE.iter = c()
	    pval_EmpiricalMSE.unscaled = c()
		 pval_EmpiricalMSE.unscaled.iter = c()
      for (i in 1:length( results$EmpiricalMSE_block_size)){
        if(is.na( results$EmpiricalMSE_block_size[i])){pval_EmpiricalMSE[i] = NA}else{
      ind[i] = which(names(pval_block_sizes) == paste("block size",
                  results$EmpiricalMSE_block_size[i]))
      pval_EmpiricalMSE[i] = pval_block_sizes[ind[i]]
        }
		        if(is.na( results$EmpiricalMSE_block_size.iter.unscaled[i])){pval_EmpiricalMSE.unscaled.iter[i] = NA}else{
      ind[i] = which(names(pval_block_sizes) == paste("block size",
                  results$EmpiricalMSE_block_size.iter.unscaled[i]))
      pval_EmpiricalMSE.unscaled.iter[i] = pval_block_sizes[ind[i]]
        }
				        if(is.na( results$EmpiricalMSE_block_size.iter[i])){pval_EmpiricalMSE.iter[i] = NA}else{
      ind[i] = which(names(pval_block_sizes) == paste("block size",
                  results$EmpiricalMSE_block_size.iter[i]))
      pval_EmpiricalMSE.iter[i] = pval_block_sizes[ind[i]]
        }
				        if(is.na( results$EmpiricalMSE_block_size.unscaled[i])){pval_EmpiricalMSE.unscaled[i] = NA}else{
      ind[i] = which(names(pval_block_sizes) == paste("block size",
                  results$EmpiricalMSE_block_size.unscaled[i]))
      pval_EmpiricalMSE.unscaled[i] = pval_block_sizes[ind[i]]
        }
		        
		}
      names(pval_EmpiricalMSE) = paste0("tuning = ",c(tuning_block_length,
                                                      "Variogram.block.size",
                                                      " Plug.In.block.size"
      ))
     
      results$pval_EmpiricalMSE =  pval_EmpiricalMSE
	  results$pval_EmpiricalMSE.iter =  pval_EmpiricalMSE.iter
	  results$pval_EmpiricalMSE.unscaled =  pval_EmpiricalMSE.unscaled
	  results$pval_EmpiricalMSE.unscaled.iter =  pval_EmpiricalMSE.unscaled.iter
    }

    
    if ("Variogram" %in% method.block.length.select){
      if(is.na( results$Variogram.block.size)){pval_Variogram= NA}else{
        pval_Variogram = pval_block_sizes[ which(names(pval_block_sizes) == paste("block size",results$Variogram.block.size))]
      }
      results$pval_Variogram =   pval_Variogram
    }
    
    if ("PlugIn" %in% method.block.length.select){
      if(is.na( results$Plug.In.block.size)){pval_PlugIn= NA}else{
        pval_PlugIn = pval_block_sizes[ which(names(pval_block_sizes) == paste("block size",results$Plug.In.block.size))]
      }
      results$pval_PlugIn =   pval_PlugIn
    }
    
  }
  if (type == "SE") {
    se_block_sizes = c()
    for (block_size in 1:length(block_Ls)) {
    distribution_of_T = unlist(boot.reps.of.Stat.function[[block_size]])
    se_block_sizes[block_size] = sd(distribution_of_T)
    }
    results$se_block_sizes = se_block_sizes
    names(se_block_sizes) = paste("block size", block_Ls)
    
    if ("PlugIn" %in% method.block.length.select){
      boot.reps = boot.reps.of.Stat.function[[paste0("blength", results$Plug.In.block.size)]][[1]]
      if(k>1){ 
        results$SE.estimate.PlugIn= sapply(1:k, function(i){sd.unlist(sapply(boot.reps, function(x){x[i]}))})
      }else{      results$SE.estimate.PlugIn = sd.unlist(boot.reps)}    }
    
    if ("Variogram" %in% method.block.length.select){
      boot.reps = boot.reps.of.Stat.function[[paste0("blength", results$Variogram.block.size)]][[1]]
      if(k>1){ 
        results$SE.estimate.Variogram= sapply(1:k, function(i){sd.unlist(sapply(boot.reps, function(x){x[i]}))})
      }else{      results$SE.estimate.Variogram = sd.unlist(boot.reps)}
    }
    
    if ("EmpiricalMSE" %in% method.block.length.select){
      results$SE.estimate.EmpiricalMSE_block_size = list()
      for (ind in 1: length (results$EmpiricalMSE_block_size)){
      boot.reps = boot.reps.of.Stat.function[[paste0("blength", results$EmpiricalMSE_block_size[ind])]][[1]]
      if(k>1){ 
        results$SE.estimate.EmpiricalMSE_block_size[[ind]]= sapply(1:k, function(i){sd.unlist(sapply(boot.reps, function(x){x[i]}))})
      }else{      results$SE.estimate.EmpiricalMSE_block_size[[ind]] = sd.unlist(boot.reps)}
      }
    }
  }
  
results
}
