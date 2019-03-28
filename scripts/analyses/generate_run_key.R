

all_res <- expand.grid("Model form"=c("Competitive"),
            "Antigenic Seniority"=c("Absent","Present"),
            "Cross Reactivity"=c("Universal","Type specific"),
            "Priming"=c("Absent","Present"),
            "Typed exposures"=c("6 types","3 types"),
            "Waning"=c("Monophasic","Biphasic"),
            "Titre dependent boosting"=c("Present","Absent"))

all_res <- all_res[order(all_res$`Antigenic Seniority`, all_res$Priming,
              all_res$`Titre dependent boosting`, all_res$Waning, 
              all_res$`Typed exposures`,all_res$`Cross Reactivity`
             ),]
all_res$runID <- 1:nrow(all_res)
runs <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/run_tracker_all.csv",stringsAsFactors = FALSE)
all_res <- merge(all_res, runs[,c("runID","runName")])
write.table(all_res,"~/Documents/Ferret_Model/antibodyKinetics/inputs/run_key.csv",sep=",",row.names=FALSE)
