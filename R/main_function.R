
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

lg <- function() {

    ret <- list()
    class(ret) <- "lg"
}
