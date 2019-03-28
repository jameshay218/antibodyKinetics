## Modifies the all_estimates table to get nicer formatting
runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker_all.csv",stringsAsFactors=FALSE) ## File location of run key table
par_estimates <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/parameter_estimates.csv",stringsAsFactors=FALSE) ## File location of parameter estimate results

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
#write.table(runs,"~/Documents/Ferret_Model/inputs/run_key.csv",sep=",",row.names=FALSE)


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
               "c"="c","lnlike"="Log likelihood","dp"="dp","ts"="ts","m"="m","y0_mod"="γ", "boost_limit"="y_limit",
               "tau"="τ")

par_names_2 <- c("mu"="μ1","mu.1"="μ2","mu.2"="μ3","mu.3"="μ4","mu.4"="μ5","mu.5"="μ6",
                  "m"="m1","m.1"="m2","m.2"="m3","m.3"="m4","m.4"="m5","m.5"="m6",  
                 "ts"="ts1","ts.1"="ts2","ts.2"="ts3","ts.3"="ts4","ts.4"="ts5","ts.5"="ts6",  
                 "sigma"="σ1","sigma.1"="σ2","sigma.2"="σ3","sigma.3"="σ4","sigma.4"="σ5","sigma.5"="σ6",  
                 "dp"="dp1","dp.1"="dp2","dp.2"="dp3","dp.3"="dp4","dp.4"="dp5","dp.5"="dp6",  
                 "y0_mod"="γ","boost_limit"="y_limit","S"="sd","beta"="β","mod1"="ρ1", "mod2"="ρ2","mod3"="ρ3","mod4"="ρ4",
                 "mod"="ρ1", "mod.1"="ρ2","mod.2"="ρ3","mod.3"="ρ4",
                 "γ"="y0_mod","y_limit"="boost_limit","lnlike"="Log likelihood","c"="c",
                 "tau"="τ"
                 )


par_estimates$par_name <- par_names[par_estimates$par_name]
par_estimates$type <- type_names[par_estimates$type]
colnames(par_estimates) <- c("runID","runName","Parameter name","Exposure Type","Mean","Median","Mode","2.5% CI","97.5% CI","ESS","Rhat","Rhat upper 95% CI")

par_estimates <- merge(runs,par_estimates, by="runID")
write.table(par_estimates,"~/Drive/Influenza/Ferret/PLOS Comp Bio/combined_results/main_par_estimates.csv",sep=",",row.names=FALSE)

########
## Format convergence/WAIC table
########

convergence_table <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/main_results/convergence_check.csv",stringsAsFactors=FALSE)
convergence_table <- convergence_table[convergence_table$runID %in% runs$runID,]
convergence_table$ess_names <- par_names_2[convergence_table$ess_names]
convergence_table$psrf_names <- par_names_2[convergence_table$psrf_names]
#convergence_table$ess_names_univ <- par_names_2[convergence_table$ess_names_univ]
#convergence_table$psrf_names_univ <- par_names_2[convergence_table$psrf_names_univ]

min_ess <- NULL
min_ess_par <- NULL
max_gelman_par <- NULL
max_gelman <- NULL
gelman_mpsrf <- NULL

convergence_table <- convergence_table[,c(1,2,13,14,15,16,17,3,4,5,6,7)]
colnames(convergence_table) <- c("runID","Run name","WAIC","ELPD LOO","ELPD LOO SE","P LOO","P LOO SE",
                                  "Minimum ESS Parameter","Minimum ESS Value",
                                 "Max Gelman PSRF Parameter","Max Gelman PSRF","Gelman MPSRF")
convergence_table$δWAIC <- convergence_table$WAIC - min(convergence_table$WAIC)
convergence_table$`δELPD LOO` <- max(convergence_table$`ELPD LOO`) - convergence_table$`ELPD LOO`
convergence_table <- merge(runs,convergence_table,by="runID")
convergence_table <- convergence_table[order(convergence_table$`δELPD LOO`),]
write.table(convergence_table,"~/Drive/Influenza/Ferret/PLOS Comp Bio/combined_results/model_comparison_table.csv",sep=",",row.names=FALSE)
