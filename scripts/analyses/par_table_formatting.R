## Modifies the all_estimates table to get nicer formatting
runs <- read.csv("~/Documents/Ferret_Model/inputs/run_tracker_formatted.csv",stringsAsFactors=FALSE)
par_estimates <- read.csv("~/Documents/new_mcmc1/par_estimates_combined.csv",stringsAsFactors=FALSE)
par_estimates <- par_estimates[,-2]
runs[runs$form=="C","form"] <- "Competitive"
runs[runs$form=="I","form"] <- "Isolated"

runs[runs$antigenic_seniority=="N","antigenic_seniority"] <- "Absent"
runs[runs$antigenic_seniority=="Y","antigenic_seniority"] <- "Present"

runs[runs$cr=="A","cr"] <- "Universal"
runs[runs$cr=="T","cr"] <- "Type specific"

runs[runs$priming=="Y","priming"] <- "Present"
runs[runs$priming=="N","priming"] <- "Absent"

runs[runs$types==3,"types"] <- "3 types"
runs[runs$types==6,"types"] <- "6 types"

runs[runs$wane=="M","wane"] <- "Monophasic"
runs[runs$wane=="B","wane"] <- "Biphasic"

runs[runs$y0_mod=="Y","y0_mod"] <- "Present"
runs[runs$y0_mod=="N","y0_mod"] <- "Absent"

runs <- runs[,c(1,4:10)]
colnames(runs) <- c("runID","Model Form","Antigenic Seniority","Cross Reactivity",
                    "Priming","Typed exposures","Waning","Titre dependent boosting")
write.table(runs,"~/Documents/Ferret_Model/inputs/run_key.csv",sep=",",row.names=FALSE)


type_names <- c("all"="All","vacc"="TIV","adj"="TIV + adjuvant","infection"="Infection",
                "vacc1"="TIV 1","vacc2"="TIV 2","adj1"="TIV 1 + adjuvant","adj2"="TIV 2 + adjuvant",
                "infection1"="Infection 1","infection2"="Infection 2","Priming"="Priming","priming"="Priming")
ordered_names <- c("All","TIV","TIV + adjuvant","Infection",
                   "TIV 1","TIV 2","TIV 1 + adjuvant","TIV 2 + adjuvant",
                   "Infection 1","Infection 2","Priming")
par_names <- c("mod1"="ρ1", "mod2"="ρ2","mod3"="ρ3","mod4"="ρ4",
               "mod"="ρ1", "mod.1"="ρ2","mod.2"="ρ3","mod.3"="ρ4",
               "adjmu"="μ(1-dp)","adjmu.1"="μ(1-dp)", "adjmu.2"="μ(1-dp)", 
               "adjmu.3"="μ(1-dp)", "adjmu.4"="μ(1-dp)", "adjmu.5"="μ(1-dp)",
               "mu"="μ","sigma"="σ","S"="sd","beta"="β",
               "c"="c","lnlike"="Log likelihood","dp"="dp","ts"="ts","m"="m","y0_mod"="γ", "boost_limit"="y_limit")

par_names_2 <- c("mu"="μ","mu.1"="μ1","mu.2"="μ2","mu.3"="μ3","mu.4"="μ4","mu.5"="μ5",
                  "m"="m","m.1"="m1","m.2"="m2","m.3"="m3","m.4"="m4","m.5"="m5",  
                 "ts"="ts","ts.1"="ts1","ts.2"="ts2","ts.3"="ts3","ts.4"="ts4","ts.5"="ts5",  
                 "sigma"="σ","sigma.1"="σ1","sigma.2"="σ2","sigma.3"="σ3","sigma.4"="σ4","sigma.5"="σ5",  
                 "dp"="dp","dp.1"="dp1","dp.2"="dp2","dp.3"="dp3","dp.4"="dp4","dp.5"="dp5",  
                 "y0_mod"="γ","boost_limit"="y_limit","S"="sd","beta"="β","mod1"="ρ1", "mod2"="ρ2","mod3"="ρ3","mod4"="ρ4",
                 "mod"="ρ1", "mod.1"="ρ2","mod.2"="ρ3","mod.3"="ρ4",
                 "γ"="y0_mod","y_limit"="boost_limit","lnlike"="Log likelihood","c"="c"
                 )


par_estimates$par_name <- par_names[par_estimates$par_name]
par_estimates$type <- type_names[par_estimates$type]
colnames(par_estimates) <- c("runID","Parameter name","Exposure Type","Mean","Median","Mode","2.5% CI","97.5% CI")

par_estimates <- merge(runs,par_estimates, by="runID")
write.table(par_estimates,"~/Documents/new_mcmc1/parameter_estimates.csv",sep=",",row.names=FALSE)


########
## Format convergence/WAIC table
########

convergence_table <- read.csv("~/Documents/new_mcmc1/convergence_check.csv",stringsAsFactors=FALSE)
convergence_table$ess_names <- par_names_2[convergence_table$ess_names]
convergence_table$psrf_names <- par_names_2[convergence_table$psrf_names]
#convergence_table$ess_names_univ <- par_names_2[convergence_table$ess_names_univ]
#convergence_table$psrf_names_univ <- par_names_2[convergence_table$psrf_names_univ]

min_ess <- NULL
min_ess_par <- NULL
max_gelman_par <- NULL
max_gelman <- NULL
gelman_mpsrf <- NULL

for(i in 1:nrow(convergence_table)){
  if((convergence_table[i,"multi_chain_fine"] & convergence_table[i,"manual"] != "univ") |
     convergence_table[i,"manual"] == "multi"){
      min_ess[i] <- convergence_table[i,"ess_vector_multi"]
      min_ess_par[i] <- convergence_table[i,"ess_names_multi"]
      max_gelman_par[i] <- convergence_table[i,"psrf_names_multi"]
      max_gelman[i] <- convergence_table[i,"gelman.psrf_multi"]
      gelman_mpsrf[i] <- convergence_table[i,"gelman.mpsrf_multi"]
    
  }
  else if(convergence_table[i,"manual"] == "univ" | 
          (convergence_table[i,"multi_chain_fine"] != TRUE & convergence_table[i,"univ_chain_fine"] == TRUE)){
    min_ess[i] <- convergence_table[i,"ess_vector_univ"]
    min_ess_par[i] <- convergence_table[i,"ess_names_univ"]
    max_gelman_par[i] <- convergence_table[i,"psrf_names_univ"]
    max_gelman[i] <- convergence_table[i,"gelman.psrf_univ"]
    gelman_mpsrf[i] <- convergence_table[i,"gelman.mpsrf_univ"]
  }
  else{
    min_ess[i] <- convergence_table[i,"ess_vector_multi"]
    min_ess_par[i] <- convergence_table[i,"ess_names_multi"]
    max_gelman_par[i] <- convergence_table[i,"psrf_names_multi"]
    max_gelman[i] <- convergence_table[i,"gelman.psrf_multi"]
    gelman_mpsrf[i] <- convergence_table[i,"gelman.mpsrf_multi"]
  }
}

convergence_table <- convergence_table[,c(1,20,21)]
convergence_table <- cbind(convergence_table,min_ess_par,min_ess,max_gelman_par,max_gelman,gelman_mpsrf)

convergence_table <- convergence_table[,c(1,11,12,3:7)]


colnames(convergence_table) <- c("runID","BIC","WAIC",
                                  "Minimum ESS Parameter","Minimum ESS Value",
                                 "Max Gelman PSRF Parameter","Max Gelman PSRF","Gelman MPSRF")
convergence_table$δWAIC <- convergence_table$WAIC - min(convergence_table$WAIC)
convergence_table$δBIC <- convergence_table$BIC - min(convergence_table$BIC)

convergence_table <- merge(runs,convergence_table,by="runID")
write.table(convergence_table,"~/Documents/new_mcmc1/waic_table.csv",sep=",",row.names=FALSE)
