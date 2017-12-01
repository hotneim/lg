
<!-- README.md is generated from README.Rmd. Please edit that file -->
The lg package for calculating local Gaussian correlations in multivariate applications
=======================================================================================

Otneim and Tjøstheim (2017a) describes a new method for estimating multivariate density functions using the concept of local Gaussian correlations. Otneim and Tjøstheim (2017b) expands the idea to the estimation of *conditional* density functions. This package, written for the R programming language, provides a simple interface for implementing these methods in practical problems.

Example
=======

Let us illustrate the use of this package by looking at the built-in data set of daily closing prices of 4 major European stock indices in the period 1991-1998. We load the data and transform them to daily returns:

``` r
data(EuStockMarkets)
x <- apply(EuStockMarkets, 2, function(x) diff(log(x)))
```

Installation instructions
=========================

Overview
========
