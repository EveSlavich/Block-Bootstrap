create_ref_subregions = function(scale1, x, y, m_x, m_y){
  #scale 1 is the spacing of the grid points
  #x,y are coordinates of the sites
  
  #order x,y by x then y
  o=order(x,partial=y)
  
  x_range = max(x)-min(x)
  y_range = max(y)-min(y)
  nx = floor ( x_range/scale1 )
  ny = floor ( y_range/scale1 )
  minx = min ( x )
  miny = min ( y )
  x0 = minx + ( x_range - nx * scale1 )/2
  y0 = miny + ( y_range - ny * scale1 )/2
  
  xM = x0 + nx*scale1
  yM = y0 + ny*scale1
  
  seqx = c ( 0, seq( 1:nx ))
  seqy = c ( 0, seq( 1:ny ))
  #the x and y grid coordinates
  finegrid_x = x0 + seqx * scale1
  finegrid_y = y0 + seqy * scale1
  
  
  
  #the x and y grid coordinates to resample
  
  ref = merge ( finegrid_x, finegrid_y, all.x = TRUE, all.y = TRUE )
  x_subrange = diff(range(x))/m_x
  y_subrange =  diff(range(y))/m_y
  x_breaks = min(x) +  c(0:m_x)*x_subrange
  y_breaks = min(y) +  c(0:m_y)*y_subrange
  x_bin = cut(finegrid_x, x_breaks, labels=c(1:m_x), include.lowest=T)
  y_bin = cut(finegrid_y, y_breaks, labels=c(1:m_y), include.lowest=T)
  #x_bin = cbind(x_bin, finegrid_x)
  #y_bin = cbind(y_bin, finegrid_y)
  ref.bins = merge ( x_bin, y_bin, all.x = TRUE, all.y = TRUE )
  bin = paste0(ref.bins[,"x"], ref.bins[,"y"])
  ref = cbind(ref,bin)
  x_bin_orig_sites = cut(x, x_breaks, labels=c(1:m_x), include.lowest=T)
  y_bin_orig_sites = cut(y, y_breaks, labels=c(1:m_y), include.lowest=T)
  site_bins = paste0(x_bin_orig_sites,y_bin_orig_sites)
  #bin = paste0(x_bin, y_bin)
  return(list(ref=ref, site_bins=site_bins))
  
}


