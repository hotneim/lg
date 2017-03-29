
# Functions for bandwidth selection
# ---------------------------------

#' Cross-validation for univariate distributions
#'
#' Uses cross-valdiation to find the optimal bandwidth for a univariate locally
#' Gaussian fit
#'
#' This function looks at a vector of observations and estimates the bandwidth that
#' minimizes the Kullback-Leibler distance between the true and the estimated densities.
#'
#' @param x The vector of data points.
#' @param tol The absolute tolerance in the optimization
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

#' Cross-validation for a pair of variables
#' 
bw_select_cv_bivariate <- function() {
}

bw_select <- function() {
}
