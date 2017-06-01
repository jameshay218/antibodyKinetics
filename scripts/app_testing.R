setwd("~/Documents/Ferret_Model/antibodyKinetics")
library(devtools)
library(shiny)
load_all()
paramViewer()


df.1 <- data.frame("id"=c(1,2,3,4,5),"values"=c(21,53215,32,34,5),"strain"=c("A","B","C","D","E"),stringsAsFactors=FALSE)
df.2 <- data.frame("id"=c(1,2,3,6,7,8),strain=c("A","B","C","A","B","C"),wow=seq(1,60,by=1),stringsAsFactors=FALSE)
