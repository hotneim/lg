#### Parallelization test on Windows laptop with 4 cores.

rm(list=ls())
library(devtools)
#install_github("martinju/lg")
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
lg_object.par <- lg(x,bw_method = "cv",parallelize = "cv",num.cores = 4)
proc.time()-aa
#user  system elapsed
#0.59    0.03   22.83


aa=proc.time()
lg_object <- lg(x,bw_method ="cv",parallelize = NULL)
proc.time()-aa
#user  system elapsed
#58.64    0.22   59.40

all.equal(lg_object.par,lg_object)
# TRUE

###############


### CONCLUSION: The paralelized version is about 2.6 times faster than the non-parallelized version in a test using a 4 cores.
# The results are identical.
