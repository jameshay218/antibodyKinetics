## Calculating model weights
## Read in model comparison results
res <- read.csv("~/Drive/Influenza/Ferret/PLOS Comp Bio/combined_results/model_comparison_table.csv",stringsAsFactors=TRUE)

## Read in ELPD estimates post cross validation
load("~/Drive/Influenza/Ferret/PLOS Comp Bio/cross_validation_results/all_loo_estimates_cv.RData")
all_loo_estimates <- old_loo_estimates

res$Waning <- factor(res$Waning, levels=c("Monophasic","Biphasic"))
res1 <- res[res$δELPD.LOO < 10,]
res2 <- res[res$δELPD.LOO > 100,]
apply(res1[,4:9], 2, table)
apply(res2[,4:9], 2, table)

## Crude GLM of how well each model mechanism predicts ELPD
fit <- glm(data=res, ELPD.LOO~Antigenic.Seniority + Cross.Reactivity + Priming + Typed.exposures + 
             Waning + Titre.dependent.boosting)
summary(fit)
confint(fit)
## Calculate ELPD weights manually (not fully checked)
res$elpd_weight <- (exp(res$ELPD.LOO - 0.5*res$ELPD.LOO.SE))/sum(exp(res$ELPD.LOO - 0.5*res$ELPD.LOO.SE)) ## with correction
res$elpd_weight1 <- (exp(res$ELPD.LOO))/sum(exp(res$ELPD.LOO))

#############
## CALCULATE PSEUDO-BMA+ WEIGHTS USING LOO PACKAGE
library(loo)
res <- res[order(res$runID),]
lpd1 <- lapply(all_loo_estimates[1:64],function(x) x$pointwise[,1])
lpd1 <- do.call("cbind",lpd1)
res$pseudo_bma <- pseudobma_weights(lpd1)
res$stacking_weights <- stacking_weights(lpd1)
res$pbma <- pseudobma_weights(lpd1,BB=FALSE)

## Mechanism weights based on calculated ELPD
sum(res[res$Priming == "Present","elpd_weight"])
sum(res[res$Titre.dependent.boosting == "Present","elpd_weight"])
sum(res[res$Typed.exposures == "6 types","elpd_weight"])
sum(res[res$Waning == "Biphasic","elpd_weight"])
sum(res[res$Cross.Reactivity == "Type specific","elpd_weight"])
sum(res[res$Antigenic.Seniority == "Present","elpd_weight"])

## From Pseudo-BMA+
sum(res[res$Priming == "Present","pseudo_bma"])
sum(res[res$Titre.dependent.boosting == "Present","pseudo_bma"])
sum(res[res$Typed.exposures == "6 types","pseudo_bma"])
sum(res[res$Waning == "Biphasic","pseudo_bma"])
sum(res[res$Cross.Reactivity == "Type specific","pseudo_bma"])
sum(res[res$Antigenic.Seniority == "Present","pseudo_bma"])

## From WAIC
sum(res[res$Priming == "Present","waic_weight"])
sum(res[res$Titre.dependent.boosting == "Present","waic_weight"])
sum(res[res$Typed.exposures == "6 types","waic_weight"])
sum(res[res$Waning == "Biphasic","waic_weight"])
sum(res[res$Cross.Reactivity == "Type specific","waic_weight"])
sum(res[res$Antigenic.Seniority == "Present","waic_weight"])

## Create a plot to look at mechanism distribution by model order
res_binary <- res[,c("runID","Antigenic.Seniority","Cross.Reactivity","Priming","Typed.exposures","Waning","Titre.dependent.boosting")]
res_binary$Antigenic.Seniority <- res_binary$Antigenic.Seniority == "Present"
res_binary$Cross.Reactivity <- res_binary$Cross.Reactivity == "Type specific"
res_binary$Priming <- res_binary$Priming == "Present"
res_binary$Typed.exposures <- res_binary$Typed.exposures == "6 types"
res_binary$Waning <- res_binary$Waning == "Biphasic"
res_binary$Titre.dependent.boosting <- res_binary$Titre.dependent.boosting == "Present"
res_binary$order <- 1:nrow(res_binary)
res_bin_1 <- reshape2::melt(res_binary,id.vars=c("runID","order"))
res_bin_1$value <- as.numeric(res_bin_1$value)
library(ggplot2)

ggplot(res_bin_1) + geom_raster(aes(x=order,y=variable,fill=value))
