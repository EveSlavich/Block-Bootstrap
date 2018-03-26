resample_blocks_subregion = function ( NBoot, lookup_table, block_L, nref, bin_i, site_bins, n_blocks) {
  new.samples.index=list()
  
  
  for (i in 1:NBoot){
    sites = c()
    for ( j in 1:n_blocks){
      random_block = sample ( 1:nref , 1 ,replace=TRUE)
      if (random_block <= length(lookup_table)){
        sites_in_block = unlist (lookup_table [[ random_block ]] )
        sites_in_block_bins = site_bins [sites_in_block ]
        sites_in_block_and_subregion = sites_in_block[which(sites_in_block_bins !=bin_i)]
        sites = c( sites, sites_in_block_and_subregion )
      }
    }
    new.samples.index [[i]] = sites
  }
  new.samples.index
}
