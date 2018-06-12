## Install the necessary packages from Rich's github
#devtools::install_github(c("gaborcsardi/progress","dide-tools/context","richfitz/queuer","dide-tools/didewin"))

## Change location to network drive
setwd("~/net/home/ferret")

## Submit credentials to didewin cluster
options(didehpc.credentials = "~/.smbcredentials",
        #didehpc.cluster = "fi--didemrchnb")
        didehpc.cluster = "fi--dideclusthn")

src <- provisionr::package_sources(local = c("~/Documents/Ferret_Model/antibodyKinetics"))
sources <- c("~/net/home/ferret/scripts/full_run.R")

## Setup contexts
context::context_log_start()
root <- "contexts"
packages <- c("antibodyKinetics","coda","ggplot2","plyr","data.table")

ctx <- context::context_save(packages=packages,path=root, sources=sources,package_sources=src)

## Submit setup to cluster
obj1 <- didehpc::queue_didehpc(ctx)
