
<!-- README.md is generated from README.Rmd. Please edit that file -->
The lg package for calculating local Gaussian correlations in multivariate applications
=======================================================================================

Otneim and Tjøstheim (2017a) describes a new method for estimating multivariate density functions using the concept of local Gaussian correlations. Otneim and Tjøstheim (2017b) expands the idea to the estimation of *conditional* density functions. This package, written for the R programming language, provides a simple interface for implementing these methods in practical problems.

Construction of the lg-object
=============================

Let us illustrate the use of this package by looking at the built-in data set of daily closing prices of 4 major European stock indices in the period 1991-1998. We load the data and transform them to daily returns:

``` r
data(EuStockMarkets)
x <- apply(EuStockMarkets, 2, function(x) diff(log(x)))

# Remove the days where at least one index did not move at all
x <- x[!apply(x, 1, function(x) any(x == 0)),]
```

When using this package, the first task is always to create an lg-object using the lg\_main()-function. This object contains all the estimation parameters that will be used in the estimation step, including the bandwidths. There are three main parameters that one can tune:

-   **The bandwidth selection method**: The bandwidth selection method is controlled through the argument *bw\_method*. It can take one of two values:
    -   Use *bw\_method = "cv"* to use the cross-validation routine described by Otneim and Tjøstheim (2017a). Depending on your system, this method is fairly slow and may take several minutes even for moderately sized data sets.
    -   Use *bw\_method = "plugin"* for plugin bandwidths. This is the **default** method and very quick. It simply sets all bandwidths equal to 1.75\*n^(-1/6), which is derived from the asymptotic convergence rates for the local Gaussian correlations. Both numbers (1.75 and -1/6) can be set manually (see documentation of the lg\_main()-function).
-   **The estimation method:** The method of estimation is controlled through the argument *est\_method*. It can take one of three values:
    -   Otneim and Tjøstheim (2017a) uses a simplified method for multivariate density estimation. The density estimate is a locally Gaussian distribution, with correlations being estimated *locally* and *pairwise*. The data is transformed for marginal standard normality (see next point) and as a consequence, we fix the means and standard deviations to 0 and 1 respectively. To use this estimation method, write *est\_method = "1par"*. This is the **default** method.
    -   Set *est\_method = "5par\_marginals\_fixed"* to estimate local means and local standard deviations marginally, as well as the pairwise local correlations. This is a more flexible method, but its theoretical properties are not (yet) fully understood. This configuration allows for the estimation of multivariate density functions without having to transform the data.
    -   The option *est\_method = "5par"* is reserved to bivariate problems, and is a fully nonparametric estimation method as laid out by Tjøstheim & Hufthammer (2013). This will simply invoke the *localgauss* package (Berentsen et. al., ).
-   **Transformation of the marginals:** This is controlled by the logical argument *transform\_to\_marginal\_normality*. If true, the marginals are transformed to marginal standard normality according to Otneim and Tjøstheim (2017a). This is the **default** method.

See the documentation of the *lg\_main()*-function for further details. We can now construct the lg-object using the default configuration by running

``` r
library(lg)
lg_object <- lg_main(x)
```

Estimation of density functions
===============================

We can then specify a set of grid points and estimate the probability density function of *x* using the *dlg()*-function. We choose a set of grid ponts that go diagonally through R^4, estimate, and plot the result as follows:

``` r
grid <- matrix(rep(seq(-.03, .03, length.out = 100), 4), ncol = 4, byrow = FALSE)
density_estimate <- dlg(lg_object = lg_object, grid = grid)
# plot(grid[,1], density_estimate$f_est, type = "l",
#     xlab = "Diagonal grid point", ylab = "Estimated density")
```

Estimation of conditional densities
===================================

If we want to calculate conditional density functions, we must take care to notice the *order* of the columns in our data set. This is because the estimation routine, implemented in the *clg()*-function, will always assume that the independent variables come first. Looking at the top of our data set:

``` r
head(x)
#>               DAX          SMI          CAC         FTSE
#> [1,] -0.009326550  0.006178360 -0.012658756  0.006770286
#> [2,] -0.004422175 -0.005880448 -0.018740638 -0.004889587
#> [3,]  0.009003794  0.003271184 -0.005779182  0.009027020
#> [4,] -0.001778217  0.001483372  0.008743353  0.005771847
#> [5,] -0.004676712 -0.008933417 -0.005120160 -0.007230164
#> [6,]  0.012427042  0.006737244  0.011714353  0.008517217
```

we see that DAX comes first. Say that we want to estimate the conditional density of DAX, given that SMI = CAC = FTSE = 0. We do that by running

``` r
grid <- matrix(seq(-.03, .03, length.out = 100), ncol = 1)   # The grid must be a matrix
condition <- c(0, 0, 0)                                      # Value of dependent variables
cond_dens_est <- clg(lg_object = lg_object, 
                     grid = grid,
                     condition = condition)
# plot(grid, cond_dens_est$f_est, type = "l",
#     xlab = "DAX", ylab = "Estimated conditional density")
```

If we want to estimate the conditional density of CAC and FTSE given DAX and SMI, for example, we must first shuffle the data so that CAC and FTSE come first, and supply the conditional value for DAX and SMI through the vector *condition*, now having two elements.

References
==========

Berentsen, Geir Drage, Tore Selland Kleppe, and Dag Tjøstheim. "Introducing localgauss, an R package for estimating and visualizing local Gaussian correlation." Journal of Statistical Software 56.1 (2014): 1-18.

Otneim, Håkon, and Dag Tjøstheim. "The locally Gaussian density estimator for multivariate data." Statistics and Computing 27, no. 6 (2017a): 1595-1616.

Otneim, Håkon, and Dag Tjøstheim. "Conditional density estimation using the local Gaussian correlation" Statistics and Computing (2017b): 1-19.

Tjøstheim, Dag, & Hufthammer, Karl Ove (2013). Local Gaussian correlation: a new measure of dependence. Journal of Econometrics, 172(1), 33-48.
