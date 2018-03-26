Select.Plug.In.block.size = function( plug_in_tuning_blocksize1,
                                       plug_in_tuning_blocksize2,
                                     boot.reps.of.Stat.function,
                                     V,r=2){
  
  Sigma2_b1 =sd.unlist(boot.reps.of.Stat.function[[paste0("blength",plug_in_tuning_blocksize1)]][[1]])
  Sigma2_b2 = sd.unlist(boot.reps.of.Stat.function[[paste0("blength",plug_in_tuning_blocksize2)]][[1]])
  Sigma2_2b2 =sd.unlist(boot.reps.of.Stat.function[[paste0("blength",2*plug_in_tuning_blocksize2)]][[1]])
  B0 = 2* plug_in_tuning_blocksize2 *(Sigma2_b2 - Sigma2_2b2)
  l0 = (1.5*r/(r+2)) * (B0^2*V/r*Sigma2_b1^2)^(1/(r+2))
  ;
  l0
}  