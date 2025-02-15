# Simulations of Overdispersion in SECR

An R package to perform and document simulations of spatially explicit 
capture--recapture in which activity centres are overdispersed 
(Efford and Fletcher 2024). Several wrapper functions are provided that call 
functions from **secr** and **secrdesign** to do the real work.

Simulation results are archived on Zenodo (Efford 2025).

The package may be installed in R using
```
devtools::install_github("MurrayEfford/overdispsim")
```

The R code to run simulations is in the Rmarkdown file overdispsim-vignette.rmd 
in the vignettes folder. The best way to view the script is to install the 
package and type 
```
vignette('overdispsim-vignette', 'overdispsim')
```

### References

Efford, M. G. (2025) Simulations of overdispersed activity centres in spatially explicit capture-recapture [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14873670

Efford MG, Fletcher D 2024. Effect of spatial overdispersion on confidence intervals for population density estimated by spatial capture-recapture. bioRxiv https://doi.org/10.1101/2024.03.12.584742 
