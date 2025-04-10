###############################################################################
## package 'overdispsim'
## onAttach.R
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('overdispsim'), "")
    packageStartupMessage("This is overdispsim ", version,
                           ". For overview type vignette('overdispsim-vignette', 'overdispsim')")
    setparameters()
}