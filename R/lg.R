#' \code{lg}: A package for calculating the local Gaussian correlation in
#' multivariate applications.
#'
#' The \code{lg} package provides implementations for the multivariate density
#' estimation and the conditional density estimation methods using local
#' Gaussian correlation as presented in Otneim & Tjøstheim (2017a) and Otneim &
#' Tjøstheim (2017b).
#'
#' The main function is called \code{lg_main}, and takes as argument a data set
#' (represented by a matrix or data frame) as well as various (optional)
#' configurations that is described in detail in the articles mentioned above,
#' and in the documentation of this package. In particular, this function
#' will calculate the bandwidths used for estimation, using either a plugin
#' estimate (default), or a cross validation estimate. If \code{x} is the data
#' set, then the following line of code will create an \code{lg} object using
#' the default configuration, that can be used for density estimation
#' afterwards:
#'
#' \code{lg_object <- lg_main(x)}
#'
#' You can change estimation method, bandwidth selection method and other
#' parameters by using the arguments of the \code{lg_main} function.
#'
#' You can evaluate the multivariate density estimate on a \code{grid} as
#' described in Otneim & Tjøstheim (2017a) using the \code{dlg}-function as
#' follows:
#'
#' \code{dens_est <- dlg(lg_object, grid = grid).}
#'
#' Assuming that the data set has \strong{p} variables, you can evaluate the \emph{conditional} density of the \strong{p - q} first variables (counting from column 1), \emph{given} the remaining \strong{q} variables being equal to \code{condition = c(v1, ..., vq)}, on a \code{grid}, by running
#'
#' \code{conditional_dens_est <- clg(lg_object, grid = grid, condition = condition)}.
#'
#' @references
#'
#'   Otneim, Håkon, and Dag Tjøstheim. "The locally gaussian density estimator for
#'   multivariate data." Statistics and Computing 27, no. 6 (2017a): 1595-1616.
#'
#'   Otneim, Håkon, and Dag Tjøstheim. "Conditional density estimation using
#'   the local Gaussian correlation" Statistics and Computing (2017b): 1-19.
#'
#' @importFrom stats dnorm optim optimise qnorm quantile
#' @importFrom utils combn
#' @docType package
#' @name lg
NULL
