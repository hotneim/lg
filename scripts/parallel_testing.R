#### Parallelization test

rm(list=ls())
install_github("hotneim/lg")
library(lg)
source("scripts/profiling_setup.R")

tt <- proc.time()
test=dlg(lg_object,grid)
proc.time()-tt

tt <- proc.time()
test=dlg_par(lg_object,grid,num.cores = 10)
proc.time()-tt


detach("package:lg", unload=TRUE)
library(devtools)
install_github("hotneim/lg")
library(lg)
tt <- proc.time()
test=dlg(lg_object,grid)
proc.time()-tt

