
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

#' Cross-validation for bivariate distributions
#'
#' Uses cross-valdiation to find the optimal bandwidth for a bivariate locally
#' Gaussian fit
#' 
#' This function looks at a 2-column matrix of observations and estimates the bandwidth that
#' minimizes the Kullback-Leibler distance between the true and the estimated densities.
#'
#' @param x The matrix of data points.
#' @param tol The absolute tolerance in the optimization
#' @param est_method The estimation method for the bivariate fit. If estimation method is
#'   \code{5par_marginals_fixed}, the marginal bandwidths must be supplied a s well through
#'   \code{bw_marginal}
#' @param bw_marginal The bandwidths for estimation of the marginals if method
#'   \code{5par_fixed_marginals} is used 
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
