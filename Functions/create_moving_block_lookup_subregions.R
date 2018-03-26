create_moving_block_lookup_subregions=function(x, y,block_L, scale1, ref, m_x, m_y){
  if(missing(scale1)) {
    scale1 = block_L/3
  }
  if(missing(ref)) {
    ref_ = create_ref_subregions(scale1, x, y, m_x, m_y)
    ref=ref_$ref
    site_bins=ref_$site_bins
  }
  ref.bin = unique(ref[,3])
  nref = c()
  for(r in 1:length(ref.bin)){
    nref [r] = length(which(ref[,3]==ref.bin[r]))
  }
  nref = data.frame(ref.bin, nref)
  
  #  scale1=5000
  #  block_L=20000
  ind = list()
  bin.ind = list()
  l = 1
  for ( i in 1:dim(ref)[1]){
    x1 = ref [ i , "x" ] - block_L
    x2 = ref [ i , "x" ] + block_L
    y1 = ref [ i , "y" ] - block_L
    y2 = ref [ i , "y" ] + block_L
    #sites_in_proximity.ind = which(sapply(x, Within,x1,x2)==TRUE & sapply(y, Within,y1,y2)==TRUE)
    sites_in_proximity.ind  = which(x>x1 & x<x2 &y>y1 &y<y2)
    sites_in_proximity = list()
    sites_in_proximity$x = x[sites_in_proximity.ind]
    sites_in_proximity$y = y[sites_in_proximity.ind]
    sites_in_block=c()
    if (length(sites_in_proximity.ind)>0){
      for (j in 1:length(sites_in_proximity.ind)){
        site_dist = sqrt ( ( sites_in_proximity$x[j] - ref$x[i])^2  + (sites_in_proximity$y[j] - ref$y[i])^2 ) #distance from site to grid cell
        if(site_dist<block_L){
          sites_in_block = c(sites_in_block,sites_in_proximity.ind[j])
        }
      }
    }
    #update ind if sites_in_block is non empty
    if (length (sites_in_block) > 0){
      ind[[l]]=sites_in_block
      bin.ind[[l]]=ref[i, "bin"]
      l=l+1
    }
  }
  return(list(ind=ind,  bin.ind=bin.ind, site_bins=site_bins, nref=nref))
}

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

