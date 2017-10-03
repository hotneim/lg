library(lg)

set.seed(123)
n <- 500
n.grid <- 1000
dim.data <- 10

x <- NULL
grid <- NULL
for (i in 1:dim.data){
  x <- cbind(x,rnorm(n))
  grid <- cbind(grid,rnorm(n))
}

lg_object <- lg(x)
