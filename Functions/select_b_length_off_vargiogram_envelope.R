try(library(geoR))
select_b_length_off_variogram_envelope = function(x,y, resids, breaks, ini, max.dist, plot=FALSE,...){
  v1 = variog(coords = data.frame(x,y), data = resids,breaks=breaks,...)
  values=c(NA,NA)
  try({v2 = variofit(v1,ini=ini,max.dist=max.dist, weights="cressie")
  if(plot==TRUE){  
    plot(v1)
    abline(v=v2$practical)
    lines( v2 )
  }
  values = c(v2$practical, v2$cov.pars[2])
  })
  ;
 values
}
