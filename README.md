# antibodyKinetics
> Fitting antibody dynamics models to HI titre data

## Summary
This package includes all of the code, data and instructions required to reproduce the analysis in the "antibody kinetics model in ferrets" paper, which you have probably been linked from. The current draft of the paper can be found on [biorxiv](https://www.biorxiv.org/content/early/2018/09/07/411751), where the full model methodology and results are described. The purpose of this project is to generate models of antibody kinetics that describe antibody kinetics following varied exposures; to fit these models to observed haemagglutination inhibition (HI) titre data; and to compare the ability of these models in fitting the observed data.

This README covers the following topics, with links to more detailed vignettes where appropriate:

1. Installation and requirements, including a link to the shiny app vignette
2. Running the basic model
3. Work flow, including changing the model structure and generating input tables
4. Link to model fitting vignette
5. Links to scripts for reproducing analyses

**NOTE:** all of the code in this README can be found in a single R script in [`scripts/readme_code.R`](https://github.com/jameshay218/antibodyKinetics/blob/master/scripts/readme_code.R).

## 1. Installation & requirements
Installation of the package itself is straightforward:
```r
devtools::install_github("jameshay218/antibodyKinetics")
library(antibodyKinetics)
```
However, additional packages are required to use all of the features.

#### 1.1 MCMC code
Whilst the base package comes with its own MCMC functionality, the model is well suited for use with the [`lazymcmc`](https://github.com/jameshay218/lazymcmc) package. In particular, the branch with parallel tempering [`https://github.com/jameshay218/lazymcmc/tree/parallel_tempering`](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering) thanks to [Ada Yan](https://github.com/ada-w-yan). Given the nature and amount of data used here, the posterior distributions are often multi-modal. Parallel tempering samples from such distributions more efficiently, and is therefore the recommended MCMC algorithm for this system. 

Please visit this [site](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering) and clone or download the full package. Although the documentation on that site applies to MCMC code without parallel tempering, the interface for the user is pretty much the same.

#### 1.2 Shiny app
There is a shiny app in this package which can be used to generate expected titre trajectories with user specified parameters and exposure schedules. This app can then be used to download .csv files with parameters and exposure timings to be fed directly into the model. It is unlikely that anyone will need to use this app, but a vignette explaining its use can be found here:

[https://jameshay218.github.io/antibodyKinetics/inst/doc/paramViewer.html](https://jameshay218.github.io/antibodyKinetics/inst/doc/paramViewer.html)

If you'd like to run this app (eg. see what titre trajectories would be generated given your own model structure and exposure schedule), the `shiny`, `shinyBS` and `rhandsontable` packages are required. The app can then be run with:

```r
paramViewer()
```

## 2. Running the basic model
Below is a very basic example of generating a single antibody titre trajectory from a given parameter set. This particular trajectory has a single exposure 10 days after the start of the experiment.
```r
## Load the package
library(antibodyKinetics)

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

## Compare the 2 solvers with identical results
plot(titre_trajectory_cpp, type='l', col="blue")
lines(titre_trajectory_R,col="red")
```

## 3. Work flow
The full work flow for the analysis is quite complex to replicate from scratch, but all of the inputs and scripts used are included. The examples below demonstrate how to run the entire model fitting process for a single model variant. 

A table describing the inputs and options used for each model variant can be found in `inputs/run_tracker_all.csv`. The examples below replicate the inputs for runID 62. `inputs/run_key.csv` describes which model mechanisms are included in which runs.

**NOTE** `stringsAsFactors` should be `FALSE` when reading in these parameter tables.

```r
## Load package
library(antibodyKinetics)
```

#### 3.1 Generation of input parameters and exposure tables
The purpose of this project is to explore different assumptions and immunological mechanisms in generating observed antibody titres. The format of input exposures and parameters therefore depends on the assumptions that the user wishes to make and the exposure schedule used. Making the correct input parameter and exposure tables is not completely trivial, and should be done using the `paramViewer()`, a vignette explaining its use can be found [here](https://jameshay218.github.io/antibodyKinetics/inst/doc/paramViewer.html). However, example parameter tables and exposure tables can be found as attached R objects and in the `inputs` folder (including all those used in this analysis). Users should refer to `inputs/run_tracker_all.csv` for a description of which parameter and expsure tables in `inputs` correspond to which model.

The `runName` identifier and corresponding options can be generated as follows:
```r
runName <- "CYTY6BN"
options <- convert_runName_to_options(runName)
```

##### 3.1.1 Example parameter input table
Refer to the documentation for this example table with `?exampleParTab`.
```r
data(exampleParTab)
```
The parameter table can be split into 4 key "sections":

 * The top 4 parameters correspond to boundary conditions and observation error parameters, as described in `describe_parameters()`
 * The next *n* parameters correspond to antibody dynamics parameters specific to either each exposure type, each exposure ID (if not applying to "all") or each exposure-strain specific set of interaction parameters
 * After these type-specific parameters, entries correspond to additional boosting elicited by priming, parameters describing the cross reactivity, titre-dependent boosting and antigenic seniority parameters
 * Parameters named "x" correspond to antigenic distances between each pair of strains in "exposure" and "strain", assuming symmetric antigenic distance (A -> B == B -> A)
 * OUTDATED: the final set of entries, `mod` correspond to the ρ parameters modifying each boost depending on its order. ie. the proportion of the 2nd boost experienced is μ*ρ1. Note that this was used as the previous implementation of antigenic seniority, and has since been superceded by the τ parameter.

Some minor tweaks are needed to fix, remove or modify some parameters depending on the options used. The parameter table can be correctly formatted using the following function:

```r
parTab <- parTab_modification(exampleParTab,options)
```
##### 3.1.2 Example exposure input table
Refer to the documentation for this example table with `?exampleExposureTab`. It is recommended that users use the shiny app to understand how this table is generated, if exposure schedules beyond those provided in the `inputs/exposureTabs` folder are needed.
```r
data(exampleExposureTab)
```
#### 3.2 Solving the model
```r
library(ggplot2)
times <- seq(0,100,by=1)

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
```

#### 3.3 Posterior function
Solving the posterior probability for a set of model parameters given the ferret titre data is straightforward. The only difference relative to the model solving function is to pass the titre data as an argument to `create_model_group_func_cpp`, and to specify "posterior" in the `version` argument.
```r
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
```

#### 3.4 Changing the model structure options
Changing the model structure is not as simple as switching some input flags on or off, because the parameters used and interactions of these parameters changes depending on the assumed model. For ease of implementation, different input tables are used for different models, with input flags being used to inform the model code which parameters to expect. 

Users are referred to:

1. [`inputs/run_tracker_all.csv`](https://github.com/jameshay218/antibodyKinetics/blob/master/inputs/run_tracker_all.csv) to find which parameter and exposure tables are used for which model variants, indexed by `runName`
2. [`inputs/run_key.csv`](https://github.com/jameshay218/antibodyKinetics/blob/master/inputs/run_key.csv) for a key showing which mechanisms are included in each `runName`
3. The shiny app (above), which allows users to generate their own parameter and exposure tables for a given set of assumptions

To test a particular model, use the appropriate `runName` in `inputs/run_key.csv` and use the parameter and exposure table files described in `run_tracker.csv` in place of `exampleParTab` and `exampleExposureTab` above.

## 4. Links to further vignettes
The above examples should be sufficient to give a feel for what the package is doing and how it is doing it. The below vignettes give more detail regarding how to reproduce the analysis in the paper:

1. [Model fitting vignette](https://jameshay218.github.io/antibodyKinetics/inst/doc/model_fitting.html)
2. [Model comparison scripts](https://github.com/jameshay218/antibodyKinetics/tree/master/scripts/analyses). Running the code in `all_analyses.R` will generate and save files with:
     
      a) Convergence diagnostics;
      b) Whether chains need to be visually inspected for convergence or rerun (eg. if bimodal posterior distributions detected);
      c) WAIC, elpd loo and p loo (from the `loo` R package);
      d) A table of posterior mean, median, mode, 95% credible intervals, effective sample size and Rhat values for each estimated parameter for each model (very large!);
      e) A table of residuals (model predicted values - observed values) for each model;
      f) R objects with all `loo` estimates and pareto k estimates (from the [loo](https://cran.r-project.org/web/packages/loo/index.html) R package)
3. [Figure generation scripts](https://github.com/jameshay218/antibodyKinetics/blob/master/scripts/figures): each script corresponds to a figure in the manuscript. However, running these relies on previously generating all of the MCMC chains as described above. Note that some minor additional formatting and labelling has been done to the figures in the manuscript. Note also that some figures (methods) were created using Inkscape, and therefore do not have accompanying scripts.
4. [Parameter table formatting script](https://github.com/jameshay218/antibodyKinetics/blob/master/scripts/analyses/par_table_formatting.R): this script was used to improve the formatting/symbols used in the output files from the model comparison scripts.


## License

GPL-3 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).
