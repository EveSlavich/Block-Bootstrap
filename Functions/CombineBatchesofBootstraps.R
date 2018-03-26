

CombineBatchesofBootstraps =function( list_of_boot_batches,...){
  N.batches = length(list_of_boot_batches)
  n_in_subregions = list_of_boot_batches[[1]]$n_in_subregions
  block_Ls= list_of_boot_batches[[1]]$block_Ls
  n_subregions=ncol(list_of_boot_batches[[1]]$sigma_stat_subregion)
  sigma2_stat_subregion_all = array(dim= c(length(block_Ls),n_subregions,N.batches))
  sigma2_stat_hat_all = matrix(nrow=N.batches, ncol=length(block_Ls))
  tuning_block_length=list_of_boot_batches[[1]]$tuning_block_length
  for (j in 1:N.batches){
    sigma_stat_subregion = list_of_boot_batches[[j]]$sigma_stat_subregion
    sigma2_stat_subregion_all[,,j] =sigma_stat_subregion^2 #block sizes x nsubregion sized matrix
    sigma2_stat_hat_all[j,]= unlist( list_of_boot_batches[[j]]$sigma_stat_hat)^2
}
  
  sigma_stat_subregion_all = sqrt(apply(sigma2_stat_subregion_all,c(1,2),mean))
  sigma_stat_hat_all = sqrt(apply( sigma2_stat_hat_all,2,mean))
  Empirical.MSE = matrix(nrow = length(block_Ls), ncol = length(block_Ls))
  for (L in 1:length(block_Ls)) {
    for (plug_in in 1:length(block_Ls)) {
      Empirical.MSE[L, plug_in] = (mean(sqrt(n_in_subregions)*sigma_stat_subregion_all[L,], na.rm = T) -
                                     sqrt(sum(n_in_subregions))* sigma_stat_hat_all[plug_in]) ^ 2  +var(sqrt(n_in_subregions)*sigma_stat_subregion_all[L,])
    }
  }

#if (is.na(tuning_block_length)==FALSE) {
  SE.estimate = NA
  Lahiri_chosen_block_length = block_Ls[which(rank(Empirical.MSE[, tuning_block_length]) ==
                                                1)]
  try({
    SE.estimate = sd.unlist(lapply(1:N.batches, function(j){list_of_boot_batches[[j]]$boot.reps.of.Stat.function[[paste0("blength", Lahiri_chosen_block_length)]]}))
  })
#}  else{
#  Lahiri_chosen_block_length = NA
 # SE.estimate = NA
 # print(
#    "Warning, no tuning block length provided for Lahiri subregion method of block length selection"
#  )
#}
rownames(Empirical.MSE) = paste0("blocklength", block_Ls)
colnames(Empirical.MSE) = paste0("tuninglength", block_Ls)
results = list(
#boot.reps.of.Stat.function = boot.reps.of.Stat.function,
Lahiri_block_size = Lahiri_chosen_block_length,
Empirical.MSE = Empirical.MSE,
SE.estimate = SE.estimate,
sigma_stat_subregion = sigma_stat_subregion_all,
n_in_subregions=n_in_subregions,
sigma_stat_hat =sigma_stat_hat_all,
block_Ls=block_Ls
)

results
}