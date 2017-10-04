#### Parallelization test on Windows laptop with 4 cores.

rm(list=ls())
library(devtools)
install_github("martinju/lg")
library(lg)
source("scripts/profiling_setup.R")

tt <- proc.time()
test1.1=dlg(lg_object,grid)
proc.time()-tt
#user  system elapsed
#15.01    0.00   15.04
test1.2=dlg(lg_object,grid+0.5)
proc.time()-tt
#user  system elapsed
#30.47    0.00   30.52

tt <- proc.time()
test2.1=dlg_par(lg_object,grid,num.cores = 4)
proc.time()-tt
#user  system elapsed
#1.00    0.01    6.30
test2.2=dlg_par(lg_object,grid+0.5,num.cores=4)
proc.time()-tt
#user  system elapsed
#2.21    0.04   13.01

detach("package:lg", unload=TRUE)
library(devtools)
install_github("hotneim/lg")
library(lg)
tt <- proc.time()
test3.1=dlg(lg_object,grid)
proc.time()-tt
#user  system elapsed
#15.85    0.02   15.86
test3.2=dlg(lg_object,grid+0.5)
proc.time()-tt
#user  system elapsed
#31.18    0.02   31.23


# Checking that all versions are identical
all.equal(test1.1,test2.1)
all.equal(test1.2,test2.2)

all.equal(test2.1,test3.1)
all.equal(test2.2,test3.2)

#> all.equal(test1.1,test2.1)
#[1] TRUE
#> all.equal(test1.2,test2.2)
#[1] TRUE
#>
#  > all.equal(test2.1,test3.1)
#[1] TRUE
#> all.equal(test2.2,test3.2)
#[1] TRUE


### CONCLUSION: The non-parallelized new version maybe slightly faster than the original (or perhaps inaccurate measuring)
### The parallelized new version is 2.4 times faster than the original using 4 cores.

