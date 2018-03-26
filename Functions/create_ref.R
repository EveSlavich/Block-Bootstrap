create_ref = function(scale1, x, y,block_L){
  #scale 1 is the spacing of the grid points
  #x,y are coordinates of the sites
  print(length(x))
  print(length(y))
  
  #order x,y by x then y
  o=order(x,partial=y)
  
  x_range = max(x)-min(x)+2*block_L
  y_range = max(y)-min(y)+2*block_L
  nx = floor ( x_range/scale1 )
  ny = floor ( y_range/scale1 )
  minx = min ( x )-block_L
  miny = min ( y )-block_L
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
  return(ref)
}