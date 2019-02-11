## Load the package
library(antibodyKinetics)
devtools::load_all("~/Documents/Ferret_Model/antibodyKinetics")
library(ggplot2)

## Vector of names input parameters to solve the model -
## note that for the cpp implementation, the order of inputs matters
## run `parameter_descriptions()` for further details
parameter_descriptions()

pars <- c("lower_bound"=0,"S"=1,"EA"=0,"MAX_TITRE"=13,
          "mu"=8,"tp"=12,"dp"=0.5,"ts"=10,"m"=0.003,"beta"=0.6, "c"=4,
          "sigma"=1,"y0_mod"=-10000,"boost_limit"=0,"tau"=0.05,
          "order"=1,"primed"=0,"mod"=1,
          "x"=0,"t_i"=10,"y0"=0,"eff_y0"=0)

## Vector of times to solve the model over
times <- seq(0,100,by=1)

## Solve the model using both the R and Cpp implementations 
titre_trajectory_cpp <- model_trajectory_cpp(pars,times)
titre_trajectory_R <- model_trajectory(pars, times)

## Compare two the solvers, identical results
plot(titre_trajectory_cpp, type='l', col="blue")
lines(titre_trajectory_R,col="red")


runName <- "CYTY6BN"
options <- convert_runName_to_options(runName)

exampleParTab <- read.csv("~/Documents/Ferret_Model/antibodyKinetics/inputs/parTabs/G_parTab.csv",stringsAsFactors=FALSE)
data(exampleParTab)
parTab <- parTab_modification(exampleParTab,options)
data(exampleExposureTab)

times <- seq(0,100,by=1)

## Cpp implementation
## Create a function pointer that only needs the unlabelled vector 
## of model parameters and the observation time vector
f <- create_model_group_func_cpp(parTab,exampleExposureTab,version="model",form=options$form,typing = TRUE,cross_reactivity = options$cr)
y <- f(parTab$values, times)
## Cpp implementation retains less labelling for speed, so 
## corresponding groups/strains must be added back to the data frame
cpp_labels <- as.data.frame(expand.grid("strain"=LETTERS[1:5],"group"=1:5))
y <- cbind(cpp_labels, as.data.frame(y))
y <- reshape2::melt(y, id.vars=c("strain","group"))
y$variable <- as.numeric(y$variable)
p <- ggplot(y) + geom_line(aes(x=variable,y=value,col=strain)) + facet_wrap(~group) + theme_bw()


## Attach titre data
data(ferret_titres)

## Only use the columns containing titres
dat <- as.matrix(ferret_titres[,4:ncol(ferret_titres)])

## The first row of the data should give the sampling times
dat <- rbind(c(0,21,37,49,70),dat)
rownames(dat) <- NULL

## Create a function pointer to solve the posterior probability
posterior_func <- create_model_group_func_cpp(exampleParTab,exampleExposureTab, dat=dat, version="posterior",
                                              form=options$form,typing=TRUE,cross_reactivity = TRUE)
posterior_func(exampleParTab$values)
