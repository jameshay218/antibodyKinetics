## Load the package
library(antibodyKinetics)
library(ggplot2)

## Vector of names input parameters to solve the model -
## note that for the cpp implementation, the order of inputs matters
## run `parameter_descriptions()` for further details
parameter_descriptions()

pars <- c("lower_bound"=0,"S"=1,"EA"=0,"MAX_TITRE"=13,
          "mu"=4,"tp"=12,"dp"=0.5,"ts"=10,"m"=0.003,"beta"=0.02, "c"=4,
          "sigma"=0.01,"y0_mod"=-10000,"boost_limit"=0,
          "primed"=0,"mod"=1,
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

data(exampleParTab)
parTab <- parTab_modification(exampleParTab,options)
data(exampleExposureTab)

times <- seq(0,100,by=1)

## R implementation
## Create a function pointer that only needs the unlabelled vector 
## of model parameters and the observation time vector
f1 <- create_model_group_func(parTab,exampleExposureTab,form=options$form,typing = TRUE,cross_reactivity = options$cr)
y1 <- f1(parTab$values, times)
y1 <- reshape2::melt(y1, id.vars=c("times","group"))    
p1 <- ggplot(y1) + geom_line(aes(x=times,y=value,col=variable)) + facet_wrap(~group)+theme_bw()

## Cpp implementation
## Create a function pointer that only needs the unlabelled vector 
## of model parameters and the observation time vector
f2 <- create_model_group_func_cpp(parTab,exampleExposureTab,version="model",form=options$form,typing = TRUE,cross_reactivity = options$cr)
y2 <- f2(parTab$values, times)
## Cpp implementation retains less labelling for speed, so 
## corresponding groups/strains must be added back to the data frame
cpp_labels <- as.data.frame(expand.grid("strain"=LETTERS[1:5],"group"=1:5))
y2 <- cbind(cpp_labels, as.data.frame(y2))
y2 <- reshape2::melt(y2, id.vars=c("strain","group"))
y2$variable <- as.numeric(y2$variable)
p2 <- ggplot(y2) + geom_line(aes(x=variable,y=value,col=strain)) + facet_wrap(~group) + theme_bw()


## Attach titre data
data(ferret_titres)

## Only use the columns containing titres
dat <- as.matrix(ferret_titres[,4:ncol(ferret_titres)])

## The first row of the data should give the sampling times
dat <- rbind(c(0,21,36,49,70),dat)
rownames(dat) <- NULL

## Create a function pointer to solve the posterior probability
posterior_func <- create_model_group_func_cpp(exampleParTab,exampleExposureTab, dat=dat, version="posterior",
                                              form=options$form,typing=TRUE,cross_reactivity = TRUE)
posterior_func(exampleParTab$values)