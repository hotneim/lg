
# Functions for calculating the local Gaussian partial correlation
# ----------------------------------------------------------------

#'Calculate the local Gaussian partial correlation
#'
#'A function that calculates the local Gaussian partial correlation for a pair
#'of variables, given the values of some conditioning variables.
#'
#'This function is a wrapper for the \code{clg}-function (for conditional
#'density estimation) that returns the local conditional, or partial,
#'correlations described by Otneim & Tjøstheim (2018). The function takes as
#'arguments an \code{lg}-object as produced by the main \code{lg_main}-
#'function, a grid of points where the density estimate should be estimated, and
#'a set of conditions.
#'
#'The variables must be sorted before they are supplied to this function. It
#'will always assume that the free variables come before the conditioning
#'variables, see \code{?clg} for details.
#'
#'Assume that X is a stochastic vector with scalar components X1 and X2, and a
#'possibly d-dimensional component X3. This function will thus compute the local
#'*partial* correlation between X1 and X2 given X3 = x3.
#'
#'@param lg_object An object of type \code{lg}, as produced by the
#'  \code{lg_main}-function.
#'@param grid A matrix of grid points, where we want to evaluate the density
#'  estimate. Number of columns *must* be equal to 2.
#'@param condition A vector with conditions for the variables that we condition
#'  upon. Length of this vector *must* be the same as the number of variables in
#'  X3. The function will throw an error of there is any discrepancy in the
#'  dimensions of the \code{grid}, \code{condition} and data set.
#'@param level Specify a level if asymptotic standard deviations and confidence
#'   intervals should be returned. If not, set to \code{NULL}.
#'
#'@return A list containing the local partial Gaussian correlations as well as all the
#'   running parameters that has been used. The elements are:
#'
#'   \itemize{
#'     \item \code{grid} The grid where the estimation was performed, on the
#'           original scale.
#'     \item \code{partial_correlations} The estimated local partial Gaussian
#'           correlations.
#'     \item \code{cond_density} The estimated conditional density of X1 and X2 given
#'           X3, as described by Otneim & Tjøstheim (2018).
#'     \item \code{transformed_grid}: The grid where the estimation was
#'           performed, on the marginal standard normal scale.
#'     \item \code{bw}: The bandwidth object.
#'     \item \code{partial_correlations_sd} Estimated standard deviations of the local
#'           partial Gaussian correlations, as described in a forthcoming paper.
#'     \item \code{partial_correlations_lower} Lower confidence limit based on the
#'           asymptotic standard deviation.
#'     \item \code{partial_correlations_upper} Upper confidence limit based on the
#'           asymptotic standard deviation.
#'   }
#'
#' @examples
#'   # A 3 variate example
#'   x <- cbind(rnorm(100), rnorm(100), rnorm(100))
#'
#'   # Generate the lg-object with default settings
#'   lg_object <- lg_main(x)
#'
#'   # Estimate the local partial Gaussian correlation between X1 and X2 given X3 = 1 on
#'   # a small grid
#'   partial_correlations <- partial_cor(lg_object,
#'                                      grid = cbind(-4:4, -4:4),
#'                                      condition = 1)
#'
#' @references
#'
#'   Otneim, Håkon, and Dag Tjøstheim. "Conditional density estimation using
#'   the local Gaussian correlation" Statistics and Computing 28, no. 2 (2018):
#'   303-321.
#'
#'@export
partial_cor <- function(lg_object, grid = NULL, condition = NULL, level = 0.95) {

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
  partial_correlations <-
    unlist(lapply(clg_object$c_cov, function(S) stats::cov2cor(S)[1,2]))

  # When calculating the gradient, we need the estimated local partial
  # covariance matrices, and they are returned by the clg-function.
  partial_cov <- clg_object$c_cov

  # Return
  ret = list(grid = clg_object$grid,
             partial_correlations = partial_correlations,
             cond_density = clg_object$f_est,
             transformed_grid = clg_object$transformed_grid,
             x = lg_object$x,
             transformed_data = lg_object$transformed_data,
             bw = lg_object$bw)

  # Calculate the asymptotic standard deviation
  # If the level-argument is provided, we also calculate the standard
  # deviations based on asymptotic formulas with corresponding confidence
  # intervals for the local partial correlations.
  if(!is.null(level)) {

    # Build a list of asymptotic covariance matrices for the *ordinary* local
    # correlations. One element per grid point.
    omega <- lapply(as.list(data.frame(t(clg_object$density_object$loc_cor_sd^2))),
                    diag)

    # In the derivation of the gradient of the transformation to partial
    # correlation, we saw that there are three categories of elements. We need
    # to identify them.
    pairs <- clg_object$density_object$pairs
    pair_category <- rep(NA, nrow(clg_object$density_object$pairs))

    # Category 1, the pair (1, 2)
    pair_category[apply((pairs == 1) | (pairs == 2 ), 1, all)] <- 1

    # Category 2, *one* element in the pair is either 1 or 2
    pair_category[apply((pairs == 1) | (pairs == 2 ), 1, any) &
                    !apply((pairs == 1) | (pairs == 2 ), 1, all)] <- 2

    # Category 3, no element is either 1 or 2
    pair_category[!apply((pairs == 1) | (pairs == 2 ), 1, any)] <- 3

    # Initialize the gradient vector. One row per grid point, one column per pair.
    g_grad <- matrix(NA, nrow = nrow(clg_object$grid), ncol = nrow(pairs))

    # Calculate the C-matrices in each grid point.
    C_list <- lapply(as.list(data.frame(t(clg_object$density_object$loc_cor))),
                     make_C,
                     pairs = pairs,
                     p = ncol(lg_object$x))

    # Replace the category 1 entries with the value of the gradient, which is just 1:
    g_grad[, pair_category == 1] <- 1

    # Replace the category 2 entries with the value of the gradient there. We loop over the pairs:
    for(i in (1:nrow(pairs))[pair_category == 2]) {

      # We need the dR12-matrix, defined in the partial correlation paper. It is
      # zero everywhere, except in the position of the correlation in question.
      dR <- matrix(0, nrow = ncol(lg_object$x), ncol = ncol(lg_object$x))
      dR[pairs[i, 1], pairs[i, 2]] <- dR[pairs[i, 2], pairs[i, 1]] <- 1
      dR12 <- dR[1:2, 3:ncol(dR)]
      dR12_list <- lapply(1:nrow(g_grad), function(a) dR12)

      # We then calculate the derivative of the partial covariance matrix with
      # respect to the local Gaussian correlation in question, in the paper
      # denoted by Sigma^(k).
      sigma_k_list <- Map('+', Map('%*%', dR12_list, C_list),
                          Map('%*%', lapply(C_list, t), lapply(dR12_list, t)))

      # Finally, we calculate the gradient from the list of partial covariances
      # and its derivative, using the formula in the paper.
      g_grad[, i] <- unlist(Map(gradient, partial_cov, sigma_k_list))

    }

    # Replace the category 3 enstries with the value of the gradient there. We loop over the pairs:
    for(i in (1:nrow(pairs))[pair_category == 3]) {

      # We need the dR3-matrix, defined in the partial correlation paper. It is
      # zero everywhere, except in the position of the correlation in question.
      dR <- matrix(0, nrow = ncol(lg_object$x), ncol = ncol(lg_object$x))
      dR[pairs[i, 1], pairs[i, 2]] <- dR[pairs[i, 2], pairs[i, 1]] <- 1
      dR3 <- dR[3:ncol(dR), 3:ncol(dR)]
      dR3_list <- lapply(1:nrow(g_grad), function(a) dR3)

      # We then calculate the derivative of the partial covariance matrix with
      # respect to the local Gaussian correlation in question, in the paper
      # denoted by Sigma^(k).
      sigma_k_list <- lapply(Map('%*%', Map('%*%', lapply(C_list, t), dR3_list), C_list), '-')

      # Finally, we calculate the gradient from the list of partial covariances
      # and its derivative, using the formula in the paper.
      g_grad[, i] <- unlist(Map(gradient, partial_cov, sigma_k_list))

    }

    # We now have the value of the gradient of g, and we can calculate the
    # asymptotic variance of the estimated local partial correlation.
    g_grad_list <- lapply(unname(as.list(data.frame(t(g_grad)))), matrix, ncol = 1)

    partial_correlations_sd <-
      sqrt(unlist(Map('%*%', Map('%*%', lapply(g_grad_list, t), omega), g_grad_list)))

    # Calculate the confidence limits based on the asymptotic normality:
    partial_correlations_lower <- partial_correlations + qnorm((1 - level)/2) *
      partial_correlations_sd
    partial_correlations_upper <- partial_correlations - qnorm((1 - level)/2) *
      partial_correlations_sd

    ret$partial_correlations_sd <- partial_correlations_sd
    ret$partial_correlations_lower <- partial_correlations_lower
    ret$partial_correlations_upper <- partial_correlations_upper

    }

  class(ret) <- "partial"
  ret
}
