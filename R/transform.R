
# Functions used to transform the data to marginal standard normality
# -------------------------------------------------------------------

#' Transform the marginals of a multivariate data set to standard normality
#' based on the logspline density estimator (Kooperberg and Stone, 1991). See
#' Otneim and Tjøstheim (2017) for details.
#' @param x The data matrix, one row per observation.
#' @return A list containing the transformed data ($transformed_data), and a
#'   function ($trans_new)  that can be used to transform grid points and obtain
#'   normalizing constants for use in density estimation functions
#'
#' @references
#'
#' Kooperberg, Charles, and Charles J. Stone. "A study of logspline density
#'   estimation." Computational Statistics & Data Analysis 12.3 (1991): 327-347.
#'
#' Otneim, Håkon, and Dag Tjøstheim. "The locally gaussian density estimator for
#'   multivariate data." Statistics and Computing 27, no. 6 (2017): 1595-1616.
trans_normal<- function(x) {

    # The sample size and the dimension
    n <- nrow(x)
    d <- ncol(x)

    # Estimate the marginals using the logspline
    estimate_marginal <- function(i) logspline::logspline(x[,i])
    marginal_estimates <- lapply(X = as.list(1:d), FUN = estimate_marginal)

    # Transform each observation vector
    transform_data <- function(i) qnorm(logspline::plogspline(x[,i],
                                                              marginal_estimates[[i]]))
    transformed_data <- matrix(unlist(lapply(X = as.list(1:d),
                                             FUN = transform_data)), ncol = d)

    # The list to be returned
    ret = list(transformed_data = transformed_data)

    # We crate a function that can be used to transform new data and grid points, and
    # to calculate the normalizing constants
    trans_new <- function(new, return_normalizing_constants = TRUE) {

        transform_grid <- function(i)
            qnorm(logspline::plogspline(new[,i], marginal_estimates[[i]]))

        calculate_normalizing_constants <- function(i) {
            logspline::dlogspline(new[,i], marginal_estimates[[i]])/
                dnorm(qnorm(logspline::plogspline(new[,i], marginal_estimates[[i]])))
        }

        normalizing_constants <- matrix(unlist(lapply(X = as.list(1:d),
                                                      FUN = calculate_normalizing_constants)),
                                        ncol = d)


        transformed_grid <- matrix(unlist(lapply(X = as.list(1:d),
                                                 FUN = transform_grid)), ncol = d)

        return(list(trans = transformed_grid,
                    normalizing_constants = normalizing_constants))
    }

    ret$trans_new <- trans_new
    return(ret)
}
