
# Functions for calculating the local Gaussian partial correlation
# ----------------------------------------------------------------

#' Calculate the local Gaussian partial correlation
#'
#' A function that calculates the local Gaussian partial correlation for a pair
#' of variables, given the values of some conditionaing variables.
#'
#' Explanation....
#'
#'
partial_cor <- function(lg_object, grid = NULL, condition = NULL) {

  # The local partial correlation is (currently) only defined for *pairs* of
  # variables, given a *set* of variables. As in the conditional density
  # estimation routine, we always assume that the two variables for which we
  # want to estimate the local partial correlation are the two first variables
  # (columns) in the data, and that the remaining variables with index > 2 are
  # the conditioning variables. We therefore run an initial check to ensure that
  # the dimensions in the supplied arguments make sense.
  if(is.null(grid)) {
    stop("You must supply a grid where the local partial correlation is to be evaluated.")
  } else {
    if(ncol(grid) != 2) {
      stop("The grid must have exactly two columns.")
    }
  }

  if(is.null(condition)) {
    stop("You must supply the values of the conditioning variables.")
  } else {
    if(length(condition) != (ncol(lg_object$x) - 2))
      stop("The number of conditions must be exactly two less than the number of variables in the data.")
  }

  # We use the function for conditional density estimation to extract the
  # conditional (or partial) covariance matrices.
  clg_object <- clg(lg_object = lg_object, grid = grid, condition = condition)

  # The element c_cov in the clg_object is now a list of the local covariance
  # matrices, which we translate into correlation matrices and pick out the
  # off-diagonal element.
  partial_correlations <- unlist(lapply(clg_object$c_cov, function(S) stats::cov2cor(S)[1,2]))

}
