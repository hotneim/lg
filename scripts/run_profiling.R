#library(devtools)
#devtools::install_github("hadley/lineprof")
library(lineprof)
source("scripts/profiling_setup.R")
source("scripts/help.functions.R")

# lineprof
#lg.lineprof <- lineprof(lg(x))
#shine(lg.lineprof)

density_estimate.lineprof <- lineprof(dlg(lg_object,grid=grid))
shine(density_estimate.lineprof)




# To attach hidding functions in the lg package
attach(asNamespace("lg"))

shine(density_estimate.lineprof)

dlg_bivariate.lineprof <- lineprof(dlg_bivariate(x <- x[,1:2]))
shine(dlg_bivariate.lineprof)

tt<-Sys.time()
aa=dlg(lg_object,grid=grid)
Sys.time()-tt
#Time difference of 3.146 secs

tt<-Sys.time()
aa=dlg_bivariate(x <- x[,1:2])
Sys.time()-tt
#Time difference of 0.4160001 secs

for
aa=dlg_bivariate(x <- x[,1:2])



### Trying to optimize lik.1par by Rcpp Armadillo
grid_point = x[1,]

bw <- c(1,1)

x1 <- x[,1]
x2 <- x[,2]

h1 <- bw[1]
h2 <- bw[2]

x1_0 <- grid_point[1]
x2_0 <- grid_point[2]

W <- dnorm(x1, mean = x1_0, sd = h1)*dnorm(x2, mean = x2_0, sd = h2)

m1 <- mean(W)
m2 <- mean(W*x1^2)
m3 <- mean(W*x2^2)
m4 <- mean(W*x1*x2)
rho = 0.5

tt<-Sys.time()
for (i in 1:10^6){
  lik.1par(rho = rho,
           m1 = m1,
           m2 = m2,
           m3 = m3,
           m4 = m4,
           x1_0 = x1_0,
           x2_0 = x2_0,
           h1 = h1,
           h2 = h2)
}
Sys.time()-tt
#Time difference of 7.085 secs

library(Rcpp)

cppFunction('double add(double x, double y, double z) {
  double sum = x + y + z;
  return sum;
}')
# add works like a regular R function

add(1.2,2,3)



sourceCpp("src/RcppFuncs.cpp")

tt<-Sys.time()
for (i in 1:10^6){
  lik1parRcpp(rho = rho,
           m1 = m1,
           m2 = m2,
           m3 = m3,
           m4 = m4,
           x1_0 = x1_0,
           x2_0 = x2_0,
           h1 = h1,
           h2 = h2)
}
Sys.time()-tt
#Time difference of 6.722031 secs








################ OLD
# Rprof-test
tmp <- tempfile()
Rprof(tmp, interval = 0.1)

#profvis({
lg_object <- lg(x)  # Put all the running parameters in here.
density_estimate <- dlg(lg_object,grid=grid)
#})

Rprof(NULL)
summaryRprof(tmp)


# Profvis
profvis({
lg_object <- lg(x)  # Put all the running parameters in here.
density_estimate <- dlg(lg_object,grid=grid)
})


