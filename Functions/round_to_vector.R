round_to_vector = function(x,v){
  x1=c()
  for (i in 1:length(x)){
    x1[i] = v[which(order(abs(x[i]-v))==1)]
  }
  ;
  x1
}

