#### Parallelization test on Windows laptop with 4 cores.

rm(list=ls())
install_github("martinju/lg")
library(lg)
library(foreach)
source("scripts/profiling_setup.R")

tt <- proc.time()
test1=dlg(lg_object,grid)
proc.time()-tt

tt <- proc.time()
test2=dlg_par(lg_object,grid,num.cores = 4)
proc.time()-tt

detach("package:lg", unload=TRUE)
library(devtools)
install_github("hotneim/lg")
library(lg)
tt <- proc.time()
test3=dlg(lg_object,grid)
proc.time()-tt
