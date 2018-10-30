
# Functions for bandwidth selection
# ---------------------------------

#' Cross-validation for univariate distributions
#'
#' Uses cross-validation to find the optimal bandwidth for a univariate locally
#' Gaussian fit
#'
#' This function provides the univariate version of the Cross Validation
#' algorithm for bandwidth selection described in Otneim & Tjøstheim (2017),
#' Section 4. Let \eqn{\hat{f}_h(x)} be the univariate locally Gaussian density
#' estimate obtained using the bandwidth \eqn{h}, then this function returns the
#' bandwidth that maximizes \deqn{CV(h) = n^{-1} \sum_{i=1}^n \log
#' \hat{f}_h^{(-i)}(x_i),} where \eqn{\hat{f}_h^{(-i)}} is the density estimate
#' calculated without observation \eqn{x_i}.
#'
#' @param x The vector of data points.
#' @param tol The absolute tolerance in the optimization, passed to the
#'   \code{optim}-function using the BFGS-method.
#'
#' @return The function returns a list with two elements: \code{bw} is the
#'   selected bandwidth, and \code{convergence} is the convergence flag returned
#'   by the \code{optim}-function.
#'
#' @examples
#'   x <- rnorm(100)
#'   bw <- bw_select_cv_univariate(x)
#'
#' @references
#'
#' Otneim, Håkon, and Dag Tjøstheim. "The locally gaussian density estimator for
#' multivariate data." Statistics and Computing 27, no. 6 (2017): 1595-1616.
#'
#' @export
bw_select_cv_univariate <- function(x, tol = 10^(-3)) {

    # The estimated (negative) KL-error for a given bandwidth, calculated by
    # cross-validation
    KL <- function(bw) {

        # Leave one out
        objective <- function(j) {
            log(dlg_marginal(x = x[-j],
                             bw = bw,
                             eval_points = x[j])$f_est)
        }
        -mean(do.call(rbind, lapply(X = as.list(1:length(x)),
                                   FUN = objective)),
             na.rm = TRUE)
    }

    # Minimize negative KL
    result <- optim(1, KL, method = "BFGS", control = list(abstol = tol))

    return(list(bw = result$par,
                convergence = result$convergence))
}

#' Cross-validation for bivariate distributions
#'
#' Uses cross-validation to find the optimal bandwidth for a bivariate locally
#' Gaussian fit
#'
#' This function provides an implementation for the Cross Validation algorithm
#' for bandwidth selection described in Otneim & Tjøstheim (2017), Section 4.
#' Let \eqn{\hat{f}_h(x)} be the bivariate locally Gaussian density estimate
#' obtained using the bandwidth \eqn{h}, then this function returns the
#' bandwidth that maximizes \deqn{CV(h) = n^{-1} \sum_{i=1}^n \log
#' \hat{f}_h^{(-i)}(x_i),} where \eqn{\hat{f}_h^{(-i)}} is the density estimate
#' calculated without observation \eqn{x_i}.
#'
#' The recommended use of this function is through the \code{lg_main} wrapper
#' function.
#'
#' @param x The matrix of data points.
#' @param tol The absolute tolerance in the optimization, used by the
#'   \code{optim}-function.
#' @param est_method The estimation method for the bivariate fit. If estimation
#'   method is \code{5par_marginals_fixed}, the marginal bandwidths must be
#'   supplied as well through the argument \code{bw_marginal}. This is
#'   automatically handled by the \code{lg_main} wrapper function.
#' @param bw_marginal The bandwidths for estimation of the marginals if method
#'   \code{5par_fixed_marginals} is used.
#'
#' @return The function returns a list with two elements: \code{bw} is the
#'   selected bandwidths, and \code{convergence} is the convergence flag returned
#'   by the \code{optim}-function.
#'
#' @examples
#'   x <- cbind(rnorm(100), rnorm(100))
#'   bw <- bw_select_cv_univariate(x)
#'
#' @references
#'
#' Otneim, Håkon, and Dag Tjøstheim. "The locally gaussian density estimator for
#' multivariate data." Statistics and Computing 27, no. 6 (2017): 1595-1616.
#'
#' @export
bw_select_cv_bivariate <- function(x,
                                   tol = 10^(-3),
                                   est_method = "1par",
                                   bw_marginal = NULL) {

    # Perform some checks
    x <- check_data(x, dim_check = 2, type = "data")
    check_est_method(est_method)
    if((est_method == "5par_marginals_fixed") & is.null(bw_marginal))
        stop("If estimation method is '5par_marginals_fixed', then bw_marginal must be supplied.")

    # The estimated (negative) KL-error for a given bandwidth, calculated by
    # cross-validation
    KL <- function(bw) {

        # Leave one out
        objective <- function(j) {
            log(dlg_bivariate(x = x[-j,],
                             bw = bw,
                             eval_points = matrix(x[j,], ncol = 2),
                             est_method = est_method,
                             run_checks = FALSE,
                             bw_marginal = bw_marginal)$f_est)
        }

        objective_values <- do.call(rbind, lapply(X = as.list(1:nrow(x)),
                                    FUN = objective))

        -mean(objective_values[abs(objective_values) != Inf],
              na.rm = TRUE)
    }

    # Minimize negative KL
    result <- optim(c(1, 1), KL, control = list(abstol = tol))

    return(list(bw = result$par,
                convergence = result$convergence))

}

#' Plugin bandwidth selection for univariate data
#'
#' Returns a plugin bandwidth for data vectors for use with univariate locally
#' Gaussian density estimation
#'
#' This function takes in a data vector of length \code{n}, and returns a the
#' real number \code{c*n^a}, which is a quick and dirty way of selecting a
#' bandwidth for univariate locally Gaussian density estimation. The number
#' \code{c} is by default set to \code{1.75}, and \code{c = -1/5} is the usual
#' exponent that stems from the asymptotic convergence rate of the density
#' estimate. Recommended use of this function is through the \code{lg_main} wrapper
#' function.
#' @param x The data vector.
#' @param n The number of data points. Can provide only this if we do not want
#'   to supply the entire data vector.
#' @param a A constant, se details.
#' @param c A constant, se details.
#'
#' @return A number, the selected bandwidth.
#'
#' @examples
#'   x <- rnorm(100)
#'   bw <- bw_select_plugin_univariate(x = x)
#'   bw <- bw_select_plugin_univariate(n = 100)
#'
#' @export
bw_select_plugin_univariate <- function(x = NULL,
                                        n = length(x),
                                        c = 1.75,
                                        a = -1/5) {
    c*n^a
}

#' Plugin bandwidth selection for multivariate data
#'
#' Returns a plugin bandwidth for multivariate data matrices for the estimation
#' of local Gaussian correlations
#'
#' This function takes in a data matrix with \code{n} rows, and returns a the
#' real number \code{c*n^a}, which is a quick and dirty way of selecting a
#' bandwidth for locally Gaussian density estimation. The number  \code{c} is by
#' default set to \code{1.75}, and \code{c = -1/6} is the usual exponent, that
#' stems from the asymptotic convergence rate of the density estimate. This
#' function is usually called from the \code{lg_main} wrapper function.
#' @param x The data matrix.
#' @param n The number of data points. Can provide only this if we do not want
#'   to supply the entire data vector.
#' @param a A constant, se details.
#' @param c A constant, se details.
#'
#' @return A number, the selected bandwidth.
#'
#' @examples
#'   x <- cbind(rnorm(100), rnorm(100))
#'   bw <- bw_select_plugin_multivariate(x = x)
#'   bw <- bw_select_plugin_multivariate(n = 100)
#'
#' @export
bw_select_plugin_multivariate <- function(x = NULL,
                                           n = nrow(x),
                                           c = 1.75,
                                           a = -1/6) {
    c*n^a
}

#' Bandwidth selection for local Gaussian correlation.
#'
#' Takes a matrix of data points and returns the bandwidths used for estimating
#' the local Gaussian correlations.
#'
#' This is the main bandwidth selection function within the framework of locally
#' Gaussian distributions as described in Otneim and Tjøstheim (2017). This
#' function takes in a data set of arbitrary dimension, and calculates the
#' bandwidths needed to find the pairwise local Gaussian correlations, and
#' is mainly used by the main \code{lg_main} wrapper function.
#' @param x A matrix or data frame with data, one column per variable, one row
#'   per observation.
#' @param bw_method The method used for bandwidth selection. Must be either
#'   \code{"cv"} (cross-validation, slow, but accurate) or \code{"plugin"}
#'   (fast, but crude).
#' @param est_method The estimation method, must be either "1par", "5par" or
#'   "5par_marginals_fixed", see \code{\link{lg_main}}.
#' @param plugin_constant_marginal The constant \code{c} in \code{cn^a} used for
#'   finding the plugin bandwidth for locally Gaussian marginal density
#'   estimates, which we need if estimation method is "5par_marginals_fixed".
#' @param plugin_exponent_marginal The constant \code{a} in \code{cn^a} used for
#'   finding the plugin bandwidth for locally Gaussian marginal density
#'   estimates, which we need if estimation method is "5par_marginals_fixed".
#' @param plugin_constant_joint The constant \code{c} in \code{cn^a} used for
#'   finding the plugin bandwidth for estimating the pairwise local Gaussian
#'   correlation between two variables.
#' @param plugin_exponent_joint The constant \code{a} in \code{cn^a} used for
#'   finding the plugin bandwidth for estimating the pairwise local Gaussian
#'   correlation between two variables.
#' @param tol_marginal The absolute tolerance in the optimization for finding
#'   the marginal bandwidths when using cross validation.
#' @param tol_joint The absolute tolerance in the optimization for finding the
#'   joint bandwidths when using cross-validation.
#'
#' @return A list with three elements, \code{marginal}  contains the bandwidths
#'   used for the marginal locally Gaussian estimation,
#'   \code{marginal_convergence} contains the convergence flags for the marginal
#'   bandwidths, as returned by the \code{optim} function, and \code{joint}
#'   contains the pairwise bandwidths and convergence flags.
#'
#' @examples
#'   x <- cbind(rnorm(100), rnorm(100), rnorm(100))
#'   bw <- bw_select(x)
#'
#' @references
#'
#' Otneim, Håkon, and Dag Tjøstheim. "The locally gaussian density estimator for
#' multivariate data." Statistics and Computing 27, no. 6 (2017): 1595-1616.
#'
#' @export
bw_select <- function(x,
                      bw_method = "plugin",
                      est_method = "1par",
                      plugin_constant_marginal = 1.75,
                      plugin_exponent_marginal = -1/5,
                      plugin_constant_joint = 1.75,
                      plugin_exponent_joint = -1/6,
                      tol_marginal = 10^(-3),
                      tol_joint = 10^(-3)) {

    # Do some sanity checks of the inputs
    x <- check_data(x, type = "data")
    check_bw_method(bw_method)
    check_est_method(est_method)

    # Dimension and sample size
    d <- ncol(x)
    n <- nrow(x)

    # Initialize the vectors and matrices
    marginal_bandwidths <- rep(NA, d)
    marginal_convergence <- rep(NA, d)

    # If the estimation method is "5par_marginals_fixed" we must first find the
    # marginal bandwidths
    if(est_method == "5par_marginals_fixed") {

        for(i in 1:d) {
            if(bw_method == "cv") {
                result <- bw_select_cv_univariate(x[, i], tol = tol_marginal)
                marginal_bandwidths[i] <- result$bw
                marginal_convergence[i] <- result$convergence

                # Print a warning if the convergence is not ok
                if(result$convergence != 0)
                    warning(paste("Cross validation for marginal bandwidth",
                                  as.character(i),
                                  "did not converge properly"))

            } else if(bw_method == "plugin") {
                marginal_bandwidths[i] <-
                    bw_select_plugin_univariate(n = n,
                                                c = plugin_constant_marginal,
                                                a = plugin_exponent_marginal)

            }
        }
    }

    # Find the joint bandwidths
    joint_bandwidths <- data.frame(t(combn(c(1:d), 2)), 0*t(combn(c(1:d), 2)))
    joint_bandwidths <- data.frame(joint_bandwidths, joint_bandwidths[,4])
    colnames(joint_bandwidths) <- c('x1', 'x2', 'bw1', 'bw2', 'convergence')

    # Iterate over all the pairs
    for(i in 1:nrow(joint_bandwidths)) {

        variables <- c(joint_bandwidths$x1[i], joint_bandwidths$x2[i])

        # Extract the pairs of variables
        bivariate_data <- x[, variables]

        if(bw_method == "cv") {

            result <- bw_select_cv_bivariate(x = bivariate_data,
                                             tol = tol_joint,
                                             est_method = est_method,
                                             bw_marginal = marginal_bandwidths[variables])

            if(result$convergence != 0)
                    warning(paste("Cross valdidation for joint bandwidths",
                                  as.character(variables[1]), "and",  as.character(variables[2]),
                                  "did not converge properly"))

            joint_bandwidths$bw1[i] <- result$bw[1]
            joint_bandwidths$bw2[i] <- result$bw[2]
            joint_bandwidths$convergence[i] <- result$convergence

        } else if(bw_method == "plugin") {

            bw <- bw_select_plugin_multivariate(n = n,
                                             c = plugin_constant_joint,
                                             a = plugin_exponent_joint)
            joint_bandwidths$bw1 <- bw
            joint_bandwidths$bw2 <- bw
            joint_bandwidths$convergence <- NA

        }

    }

    ret <- list(marginal = marginal_bandwidths,
                marginal_convergence = marginal_convergence,
                joint = joint_bandwidths)

    return(ret)

}

#' Create simple bandwidth object
#'
#' Create a simple bandwidths object for local Gaussian correlations
#'
#' This function provides a quick way of producing a bandwidth object that may
#' be used in the \code{lg_main()}-function. The user must specify a bandwidth
#' \code{joint} that is used for all joint bandwidths, and the user may specify
#' \code{marg}, a marginal bandwidth that will be used for all marginal
#' bandwidths. This is needed if the subsequent analyses use
#' \code{est_method = "5par_marginals_fixed"}.
#'
#' The function must know the dimension of the problem, which is achieved by
#' either supplying the data set \code{x} or the number of variables \code{dim}.
#'
#' @param joint Joint bandwidth
#' @param marg Marginal bandwidths
#' @param x The data set
#' @param dim The number of variables
#'
#' @examples
#'
#'   bw_object <- bw_simple(joint = 1, marg = 1, dim = 3)
#'
#' @export

bw_simple <- function(joint = 1, marg = NA, x = NULL, dim = NULL) {

  # Check that we have wither the data or a dimension
  if(is.null(x) & is.null(dim)) {
    stop("Specify either x or dim")
  }

  # Create a dummy data det if only the dimension is specified
  if(is.null(x)) {
    x <- matrix(1, nrow = 2, ncol = dim)
  }

  bw_select(x,
            bw_method = "plugin",
            est_method = "5par_marginals_fixed",
            plugin_constant_marginal = marg,
            plugin_exponent_marginal = 0,
            plugin_constant_joint = joint,
            plugin_exponent_joint = 0)

}
