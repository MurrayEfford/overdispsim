# Simulations of Overdispersion in SECR

An R package to perform and document simulations of spatially explicit 
capture--recapture in which activity centres are overdispersed 
(Efford and Fletcher 2024). Several wrapper functions are provided that call 
functions from **secr** and **secrdesign** to do the real work.

The R code to run simulations is in the Rmarkdown file overdispsim-vignette.rmd 
in the vignettes folder. 

Simulation results are archived on Zenodo (Efford 2025).

The package may be installed in R using
```
devtools::install_github("MurrayEfford/overdispsim",  build_vignettes = TRUE)
```

(this takes some time to run, as it rebuilds the vignette)

or
 
```
install.packages("overdispsim", repos = "https://MurrayEfford.r-universe.dev")
```

(faster as binaries are pre-built).

The vignette may be viewed by installing the package and typing 

```
vignette('overdispsim-vignette', 'overdispsim')
```

or by downloading from Zenodo.

### References

Efford, M. G. (2025) Simulations of overdispersed activity centres in spatially explicit capture-recapture [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14873669

Efford MG, Fletcher D 2024. Effect of spatial overdispersion on confidence intervals for population density estimated by spatial capture-recapture. bioRxiv https://doi.org/10.1101/2024.03.12.584742 
