

BlockBoot_apply_subregions = function (x,y,m_x,m_y, dat, block_Ls, Grid_space,shape,NBoot,Stat.function,sampling_type,  lookuptables.folderpath,type,subregion.division,...){
print(subregion.division) 
 if(subregion.division == "mutually.exclusive"){
  bins = break_into_subregions(x , y ,  m_y = m_y , m_x = m_x)
  bins.levels = levels(bins$bin)
  sigma_stat_subregion = matrix(nrow=length(block_Ls), ncol=(m_x*m_y))
  rownames(sigma_stat_subregion) = block_Ls
  
  for(L in 1:length(block_Ls)){
    for (subregion in 1:(m_x*m_y)){
    print(paste0("subregion",subregion))
        print(paste0("mxy",m_x,m_y))
      block_L = block_Ls[L]
      dat_subregion  = dat [which(bins$bin == bins.levels[subregion]),]
      x_subregion =  x [which(bins$bin == bins.levels[subregion])]
      y_subregion =  y [which(bins$bin == bins.levels[subregion])]
      if(nrow(dat_subregion)>0){
      if (block_L>0){
	            rm(lookup_table,lookup.coords,envir=.GlobalEnv)
        if (is.na(lookuptables.folderpath) ==FALSE){ #check if a lookup table has been created to speed things up, if it has load it, and check the lookup table is for data with same x,y, coordinates
          load_file_if_exists(paste0(lookuptables.folderpath,"lookup_table_subregion",subregion,"_L",block_L,"_grid_space_",Grid_space,"_sampling_type_","sites","_",shape,".RData"))
          if(exists("lookup.coords")){
          print("using existing coords...")
          if( identical ( lookup.coords$lookup.x , x_subregion ) == FALSE |  identical ( lookup.coords$lookup.y , y_subregion ) == FALSE ){
                    print("wrong sites... creating new lookup")
           
			new_sample_subregion = resample_blocks_by_area(NBoot = NBoot, x=x_subregion , y=y_subregion , block_L=block_Ls[L],Grid_space = Grid_space, area_or_sites =sampling_type,shape=shape,lookup_tablename=paste0("lookup_table_subregion",subregion),  lookuptables.folderpath = lookuptables.folderpath, ...)
          }else{
                    
		  new_sample_subregion = resample_blocks_by_area(NBoot = NBoot,lookup_table= lookup_table,  x=x_subregion , y=y_subregion , block_L=block_Ls[L],Grid_space = Grid_space, area_or_sites =sampling_type,shape=shape, ...)
		  }
		    rm(lookup_table,lookup.coords,envir=.GlobalEnv)
 
        }else{
         print("creating new lookup")
	  new_sample_subregion = resample_blocks_by_area(
         x = x_subregion,
        y = y_subregion,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = sampling_type,
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
                lookup_tablename=paste0("lookup_table_subregion",subregion)
		#  ...
      )###will create a lookup table
	  }}else{print("code under developement")}
           
      }
      #If block_L =0, do an iid bootstrap
      if (block_L==0){
        new_sample_subregion = list()
        for(i in 1:NBoot){
          new_sample_subregion [[i]] = sample(1:length(x_subregion), size=length(x_subregion), replace =T) 
        }
      }
      
      boot.reps.of.Stat.function_subregion = bootstrap_wrapper(dat = dat_subregion, function_to_repeat = Stat.function, new_sample = new_sample_subregion, NBoot=NBoot,type=type,...)
     sigma_stat_subregion[L,subregion] = sd(unlist(boot.reps.of.Stat.function_subregion))
    }else{
         sigma_stat_subregion[L,subregion] = 0
    }
    }
  }
  }else{
  ### subregion division not mutually exclusive
    bins = break_into_subregions(x , y ,  m_y = m_y , m_x = m_x)
  bins.levels = levels(bins$bin)
  sigma_stat_subregion = matrix(nrow=length(block_Ls), ncol=(m_x*m_y))
  rownames(sigma_stat_subregion) = block_Ls
  
  for(L in 1:length(block_Ls)){
    for (subregion in 1:(m_x*m_y)){
       block_L = block_Ls[L]
      dat_subregion  = dat [which(bins$bin != bins.levels[subregion]),]
      x_subregion =  x [which(bins$bin != bins.levels[subregion])]
      y_subregion =  y [which(bins$bin != bins.levels[subregion])]
      if(nrow(dat_subregion)>0){
      if (block_L>0){
	            rm(lookup_table,lookup.coords,envir=.GlobalEnv)
        if (is.na(lookuptables.folderpath) ==FALSE){ #check if a lookup table has been created to speed things up, if it has load it, and check the lookup table is for data with same x,y, coordinates
          load_file_if_exists(paste0(lookuptables.folderpath,"lookup_table_subregion_overlap",subregion,"_L",block_L,"_grid_space_",Grid_space,"_sampling_type_","sites","_",shape,".RData"))
          if(exists("lookup.coords")){
          print("using existing coords...")
          if( identical ( lookup.coords$lookup.x , x_subregion ) == FALSE |  identical ( lookup.coords$lookup.y , y_subregion ) == FALSE ){
                    print("wrong sites... creating new lookup")
           
			new_sample_subregion = resample_blocks_by_area(NBoot = NBoot, x=x_subregion , y=y_subregion , block_L=block_Ls[L],Grid_space = Grid_space, area_or_sites =sampling_type,shape=shape,lookup_tablename=paste0("lookup_table_subregion_overlap",subregion),  lookuptables.folderpath = lookuptables.folderpath, ...)
          }else{
                    
		  new_sample_subregion = resample_blocks_by_area(NBoot = NBoot,lookup_table= lookup_table,  x=x_subregion , y=y_subregion , block_L=block_Ls[L],Grid_space = Grid_space, area_or_sites =sampling_type,shape=shape, lookup_tablename=paste0("lookup_table_subregion_overlap",subregion),...)
		  }
		    rm(lookup_table,lookup.coords,envir=.GlobalEnv)
 
        }else{
         print("creating new lookup")
	  new_sample_subregion = resample_blocks_by_area(
         x = x_subregion,
        y = y_subregion,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = sampling_type,
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
                lookup_tablename=paste0("lookup_table_subregion_overlap",subregion)
		#  ...
      )###will create a lookup table
	  }}else{print("code under developement")}
           
      }
      #If block_L =0, do an iid bootstrap
      if (block_L==0){
        new_sample_subregion = list()
        for(i in 1:NBoot){
          new_sample_subregion [[i]] = sample(1:length(x_subregion), size=length(x_subregion), replace =T) 
        }
      }
      
      boot.reps.of.Stat.function_subregion = bootstrap_wrapper(dat = dat_subregion, function_to_repeat = Stat.function, new_sample = new_sample_subregion, NBoot=NBoot,type=type,...)
     sigma_stat_subregion[L,subregion] = sd(unlist(boot.reps.of.Stat.function_subregion))
    }else{
         sigma_stat_subregion[L,subregion] = 0
    }
    }
  }
  }
  
  ;
  sigma_stat_subregion
}

