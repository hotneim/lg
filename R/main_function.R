
# The main function that creates the 'lg'-object
# ----------------------------------------------

#' Create an \code{lg} object
#'
#' Create an \code{lg}-object, that can be used to estimate local Gussian correlations,
#' unconditional and conditional densities, local partial correlation and for testing
#' purposes.
#'
#' This is the main function in the package. It lets the user supply a data set and set
#' a number of options, which is then used to prepare an \code{lg} object that can be
#' supplied to other functions in the package, such as \code{dlg} (density estimation),
#' \code{clg} (conditional density estimation), or other tasks in the locally Gaussian
#' universe.
#' @param x A matrix or data frame with data, on column per variable, one row
#'   per observation
#' @param bw_method The method used for bandwidth selection. Must be either
#'   \code{"cv"} (cross-validation, slow, but accurate) or \code{"plugin"} (fast,
#'   but crude)
#' @param est_method The estimation method, must be either "1par", "5par" or
#'   "5par_marginals_fixed"
#' @param transform_to_marginal_normality Logical true if we want to transform our
#'   data to marginal standard normality. This is assumed by method "1par", but can
#'   of course be skipped using this argument if it has been done already.
#' @param bw Bandwidth object if it has already been calculated.
#' @param plugin_constant_marginal The constant \code{c} in \code{cn^a} used for
#'   finding the plugin bandwidth for locally Gaussian marginal density estimates,
#'   which we need if estimation method is "5par_marginals_fixed".
#' @param plugin_exponent_marginal The constant \code{a} in \code{cn^a} used for
#'   finding the plugin bandwidth for locally Gaussian marginal density estimates,
#'   which we need if estimation method is "5par_marginals_fixed".
#' @param plugin_constant_joint The constant \code{c} in \code{cn^a} used for
#'   finding the plugin bandwidth for estimating the pairwise local Gaussian
#'   correlation between two variables.
#' @param plugin_exponent_joint The constant \code{a} in \code{cn^a} used for
#'   finding the plugin bandwidth for estimating the pairwise local Gaussian
#'   correlation between two variables.
#' @param tol_marginal The absolute tolerance in the optimization for finding the
#'   marginal bandwidths
#' @param tol_joint The absolute tolerance in the optimization for finding the
#'   joint bandwidths
#' @export
lg <- function(x,
               bw_method = "plugin",
               est_method = "1par",
               transform_to_marginal_normality = TRUE,
               bw = NULL,
               plugin_constant_marginal = 1.75,
               plugin_constant_joint = 1.75,
               plugin_exponent_marginal = -1/5,
               plugin_exponent_joint = -1/6,
               tol_marginal = 10^(-3),
               tol_joint = 10^(-3)) {

    # Sanity checks
    x <- check_data(x, type = "data")
    check_est_method(est_method)
    if((est_method == "5par") & (ncol(x) != 2)) {
        stop("Data must be bivariate if estimation method is '5par'")
    }
    if((est_method == "1par") & (transform_to_marginal_normality == FALSE)) {
        warning("Estimation method '1par' assumes marginal standard normality.")
    }

    # Return a list
    ret <- list()
    ret$x <- x
    ret$bw_method <- bw_method
    ret$est_method <- est_method

    # Transformation
    if(transform_to_marginal_normality) {
        transformed <- trans_normal(x = x)
        ret$transformed_data <- transformed$transformed_data
        ret$trans_new <- transformed$trans_new
    } else {
        ret$transformed_data = x
        ret$trans_new = NA
    }
    
    # Bandwidth selection
    if(is.null(bw)) {
        bw <- bw_select(ret$transformed_data,
                        bw_method = bw_method,
                        est_method = est_method,
                        plugin_constant_marginal = plugin_constant_marginal,
                        plugin_exponent_marginal =  plugin_exponent_marginal,
                        plugin_constant_joint = plugin_constant_joint,
                        plugin_exponent_joint = plugin_exponent_joint,
                        tol_marginal = tol_marginal,
                        tol_joint = tol_joint)
    }
    ret$bw <- bw
    
    class(ret) <- "lg"

    return(ret)
}
