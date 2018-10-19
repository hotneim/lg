
# Functions for checking the arguments
# ------------------------------------

#' Check the data and grid
#'
#' Checks that the data or grid provided is of the correct form. This function
#' is an auxiliary function that can quickly check that a supplied data set or
#' grid is a matrix or a data frame, and that it has the correct dimension, as
#' defined by the \code{dim_check} parameter. The \code{type} argument is simply
#' a character vector "data" or "grid" that is used for printing error messages.
#'
#' @param x Data or grid
#' @param dim_check How many columns do we expect?
#' @param type Is it the "grid" or "data" for use in error messages.
#'
check_data <- function(x, dim_check = NA, type) {

    # The data an the evaluation points should be a matrix, and if a data frame we
    # convert it to a matrix
    if (is.matrix(x) ) {
    } else if (is.data.frame(x)) {
        x <- as.matrix(x)
    } else {
        stop(paste("The ",
                   type,
                   " must be a matrix or a data frame", sep = ""))
    }

    # Check that ret now has the correct number of columns
    if(!is.na(dim_check)) {
        if(ncol(x) != dim_check) {
            stop(paste("The ",
                       type,
                       " can only have ",
                       as.character(dim_check),
                       " variables", sep = ""))
        }
    }

    # Remove NAs, but issue message if any is found
    x <- na.omit(x)
    if(class(attributes(x)$na.action) == "omit") {
      message(paste("Number of NA observations removed from the data: ",
                    length((attributes(x)$na.action)),
                    ".",
                    sep = ""))
    }

    return(x)
}

#' Check bandwidth vector
#'
#' Checks that the bandwidth vector supplied to the bivariate density function
#' is a numeric vector of length 2.
#'
#' @param bw The bandwidth vector to be checked
check_bw_bivariate <- function(bw) {
     # The bandwidths can only be a numerical vector of two elements
    if(!is.vector(bw)) {
        stop("bw must be a vector")
    } else if(length(bw) != 2) {
        stop("bw must have length 2")
    } else if(!is.numeric(bw)) {
        stop("bw must be numeric")
    }
}

#' Check estimation method
#'
#' Checks that the estimation method is one of the allowed values, currently
#' "1par", "5par" and "5par_marginals_fixed".
#'
#' @param est_method Check if equal to "1par" or "5par"
check_est_method <- function(est_method) {
    if(!(est_method %in% c("1par", "5par", "5par_marginals_fixed")))
        stop("Estimation method must be either '1par', '5par' or '5par_marginals_fixed'")
}

#' Check bw method
#'
#' Checks that the bandwidth method is one of the allowed values, currently "cv"
#' or "plugin".
#'
#' @param bw_method Check if equal to "cv" or "plugin"
check_bw_method <- function(bw_method) {
    if(!(bw_method %in% c("cv", "plugin")))
        stop("Estimation method must be either 'cv' or 'plugin'")
}

#' Check the arguments for the \code{dmvnorm_wrapper} function
#'
#' Checks that the arguments provided to the \code{dmvnorm_wrapper}-function are
#' numerical vectors, all having the same lengths.
#'
#' @inheritParams dmvnorm_wrapper
check_dmvnorm_arguments <- function(eval_points, mu_1, mu_2, sig_1, sig_2, rho) {
    if(!is.numeric(eval_points) |
       !is.numeric(mu_1) |
       !is.numeric(mu_2) |
       !is.numeric(sig_1) |
       !is.numeric(sig_2) |
       !is.numeric(rho))
        stop("All arguments must be numeric")

    if(!(ncol(eval_points) == 2))
        stop("eval_points must be a matrix with 2 columns")

    ref <- nrow(eval_points)
    if(!all(length(mu_1) == ref,
           length(mu_2) == ref,
           length(sig_1) == ref,
           length(sig_2) == ref,
           length(rho) == ref))
        stop("Mismatch in dimensions in dmvnorm_wrapper")
}

#' Check that an object has class "lg"
#'
#' Checks that the provided object has class \code{lg}.
#'
#' @param check_object The object to be checked
check_lg <- function(check_object) {
    if(class(check_object) != "lg")
        stop("Object must be of class 'lg'")
}

