bootstrap_wrapper = function(dat, function_to_repeat, new_sample, NBoot, function_on_reps=return_reps,type="SE",...){
   start.time <- Sys.time()
  if(missing(NBoot)){
    NBoot = length(new_sample)
  }else{
    NBoot = min( NBoot, length(new_sample))
  }
  #print (paste("Running Block Bootstrap with", NBoot, "resamples"))
  function_resample = as.list(rep(NA, NBoot))
  for (i in 1:NBoot){
   # if((i%%10) == 1){
      #print every tenth i
    #  print(paste0( "Boot iter = ",i))
    #}
    try({
      dat1 = dat [ new_sample[[i]], ]
      ind1 = new_sample[[i]]
	  if(type=="Pval"){ function_resample[[i]] = function_to_repeat (dat1 , ind1, ...)}
      if(type=="SE"){ function_resample[[i]] = function_to_repeat (dat1 ,  ...)}
      })
  }
  end.time <- Sys.time()
  time.taken <- start.time - end.time
  #print(time.taken)
  #print(value)
  list(function_resample)
}
