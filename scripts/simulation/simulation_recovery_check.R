## Multivariate chains
runs <- read.csv("~/net/home/ferret/inputs/run_tracker.csv",stringsAsFactors=FALSE)
mcmcPars <- c("adaptive_period"=1000000,"iterations"=7000000,"opt_freq"=1000,"thin"=1000,"save_block"=1000,"popt"=0.234)
thin <- 1
burnin <- mcmcPars["adaptive_period"]
top_wd <- getwd()
all_res <- NULL
for(i in 1:nrow(runs)){
  runName <- runs$runName[i]
  print(runName)
  parTab_file <- paste0("~/net/home/ferret/inputs/parTab_new/",runs$parTab_file[i],".csv")
  parTab <- read.csv(parTab_file,stringsAsFactors=FALSE)
  
  setwd(top_wd)
  chain_wd <- paste0("~/net/home/ferret/outputs_sim_fixed/",runName)
  setwd(chain_wd)
  useMulti <- TRUE
  if(runName %in% c("A","S","V","D")) useMulti <- FALSE
  chains <- load_mcmc_chains(chain_wd,parTab,asList=FALSE,convertMCMC=FALSE,TRUE,1,mcmcPars["adaptive_period"],useMulti)[[2]]
  chains <- chains[,1:(ncol(chains)-1)]
  results <- apply(chains, 2, function(x) quantile(x,c(0,0.025,0.5,0.975,1)))
  true_par <- parTab[which(parTab$fixed == 0),"values"]
  results <- as.data.frame(t(rbind(results,true_par)))
  results$correct_quant <- (results$true_par > results$`2.5%` & results$true_par < results$`97.5%`)
  results$correct_limits <- (results$true_par > results$`0%` & results$true_par < results$`100%`)
  results$runName <- runName
  results$par <- parTab[which(parTab$fixed == 0),"names"]
  results$type <- parTab[which(parTab$fixed == 0),"type"]
  results$strain <- parTab[which(parTab$fixed == 0),"strain"]
  results$exposure <- parTab[which(parTab$fixed == 0),"exposure"]
  all_res <- rbind(all_res, results)
}

sim_results_fixed <- plyr::ddply(all_res,"runName", function(x) sum(x$correct_quant)/nrow(x))
