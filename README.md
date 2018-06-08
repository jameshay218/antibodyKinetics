# antibodyKinetics
> Fitting antibody dynamics models to HI titre data

## Summary
This package includes all of the code, data and instructions required to reproduce the analysis in the "antibody kinetics model in ferrets" paper, which you have probably been linked from. The current draft of the paper can be found on [biorxiv](), where the full model methodology and results are described. The purpose of this project is to generate models of antibody kinetics that describe antibody kinetics following varied exposures; to fit these models to observed HI titre data; to compare the ability of these models in fitting the observed data.

## Installation
Installation of the package itself is straightforward:
```r
devtools::install_github("jameshay218/antibodyKinetics")
library(antibodyKinetics)
```
However, additional packages are required to use all of the features.

### MCMC
Whilst the base package comes with its own MCMC functionality, the model is well suited for use with the [`lazymcmc`](https://github.com/jameshay218/lazymcmc) package. In particular, the branch with parallel tempering [`https://github.com/jameshay218/lazymcmc/tree/parallel_tempering`](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering) thanks to [Ada Yan](https://github.com/ada-w-yan). This is because, given the nature and amount of data used here, the posterior distribution is multi modal. Parallel tempering samples from such distributions more efficiently. 

Please visit this [site](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering) and clone or download the full package. Although the documentation on that site applies to MCMC code without parallel tempering, the interface for the user is pretty much the same.

### Shiny app
There is a shiny app in this package which can be used to generate expected titre trajectories with user specified parameters and exposure schedules. This app can then be used to download .csv files with parameters and exposure timings to be fed directly into the model. It is unlikely that anyone will need to use this app, but a vignette explaining its use can be found here:

(vignette link)[placeholder]

If you'd like to run this app (eg. see what titre trajectories would be generated given your own model structure and exposure schedule), the `shiny`, `shinyBS` and `rhandsontable` packages are required. The app can then be run with:

```r
paramViewer()
```

## Work flow
The full work flow for the analysis is quite complex. It can be broken down to the following steps:

### 1. Generation of input parameters and exposure tables
The purpose of this project is to explore different assumptions and immunological mechanisms in generating observed antibody titres. The format of input exposures and parameters therefore depends on the assumptions that the user wishes to make and the exposure schedule used. 

### 2. 


## License

GPL-3 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).
