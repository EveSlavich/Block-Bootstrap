


BlockBootID = function (x ,
                           y,
                           block_Ls,
                           NBoot = 500,
                           Grid_space,
                           lookuptables.folderpath = NA,
                           shape = "disc",

                           
                           ...) {
  
  #set the spacing between sampling points
  if(missing(Grid_space)) {
    Grid_space = min(block_Ls)/3
  }
    block_L = block_Ls[1]
  
    if (block_L > 0) {
      if (is.na(lookuptables.folderpath) ==FALSE){ #check if a lookup table has been created to speed things up, if it has load it, and check the lookup table is for data with same x,y, coordinates
        load_file_if_exists(paste0(lookuptables.folderpath,"lookup_table","_L",block_L,"_grid_space_",Grid_space,"_",shape,".RData"))
        
        if(exists("lookup.coords")){
          if( identical ( lookup.coords$lookup.x , x ) == FALSE |  identical ( lookup.coords$lookup.y , y ) == FALSE ){
            rm(lookup_table,lookup.coords,envir=.GlobalEnv)
          }
        }
      }
      
      #Create the new sample (BlockBootID)
      new_sample = resample_blocks_by_area(
        x = x,
        y = y,
        NBoot = NBoot,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = "sites",
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
        ...
      )
      
    }
    #If block_L =0, do an iid bootstrap
    if (block_L == 0) {
      new_sample = list()
      for (i in 1:NBoot) {
        new_sample [[i]] = sample(1:length(x), size = length(x), replace = T)
      }
    }
  new_sample = t(as.data.frame(new_sample))
  row.names(new_sample) = paste0("boot.",1:NBoot)
  ;
  new_sample
}
