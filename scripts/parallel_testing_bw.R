#### Parallelization test on Windows laptop with 4 cores.

rm(list=ls())
library(devtools)
install_github("martinju/lg")
library(lg)

#### Test setup ####
set.seed(123)
n <- 500
n.grid <- 1000
dim.data <- 3

x <- NULL
grid <- NULL
for (i in 1:dim.data){
  x <- cbind(x,rnorm(n))
  grid <- cbind(grid,rnorm(n))
}

### Running the tests

aa=proc.time()
lg_object.par <- lg(x,bw_method = "cv",num_cores_bw_cv = 4)
proc.time()-aa
#user  system elapsed
#0.21    0.01   12.06

aa=proc.time()
lg_object <- lg(x,bw_method ="cv")
proc.time()-aa
#user  system elapsed
#28.30    0.35   29.19

detach("package:lg", unload=TRUE)
library(devtools)
install_github("hotneim/lg")
library(lg)
tt <- proc.time()
lg_object.old.package <- lg(x,bw_method ="cv")
proc.time()-tt
#user  system elapsed
#27.28    0.08   27.50

all.equal(lg_object.par,lg_object)
# TRUE
all.equal(lg_object.par,lg_object.old.package)
# TRUE


###############


### CONCLUSION: The paralelized version is about 2.4 times faster than the non-parallelized version in a test using a 4 cores.
# The results are identical.
