
# Functions for plotting the local correlations
# ---------------------------------------------

#' Plot unconditional local correlation maps
#'
#' Plot the estimated local correlation map for a pair of variables
#'
#' This function plots a map of estimated local Gaussian correlations of a
#' specified pair (defaults to the first pair) of variables as produced by the
#' dlg-function. This plot is heavily inspired by the local correlation plots
#' produced by the 'localgauss'-package by Berentsen et. al (2014), but it is
#' here more easily customized and specially adapted to the ecosystem within the
#' \code{lg}-package. The plotting is carried out using the ggplot2-package
#' (Wickham, 2009).
#'
#' @param dlg_object The density estimation object produced by the dlg-function
#' @param pair Integer indicating which pair of variables you want to plot. The
#'   function looks up the corresponding variables in the bandwidth object used
#'   to calculate the dlg object, and you can inspect this in
#'   \code{dlg_object$bw$joint}. Defaults to 1 (the first pair, usually variable
#'   1 against variable 2).
#' @param gaussian_scale Logical, if \code{TRUE} the plot is produced on the
#'   marginal standard Gaussian scale.
#' @param plot_colormap Logical, if \code{TRUE} the plot includes a colormap to
#'   visualize the value of the local correlation.
#' @param plot_obs Logical, if \code{TRUE} the observations are plotted.
#' @param plot_labels Logical, if \code{TRUE} character labels with local
#'   correlation values are plotted.
#' @param plot_legend Logical, if \code{TRUE} a color legend is plotted.
#' @param plot_thres A number between 0 and 1 indicating the threshold value to
#'   be used for not plotting the estimated local correlation in areas with no
#'   data. Uses a quick bivariate kernel density estimate a criterion, and skips
#'   plotting in areas with kernel density estimate less than the fraction
#'   plot_thres of the maximum density estimate. If 0 (default), everything is
#'   plotted, if 1 nothing is plotted. Typical values may be in the
#'   0.001-0.01-range.
#' @param alpha_tile The alpha-value indicating the transparency of the color
#'   tiles. Number between 0 (transparent) and 1 (not transparent).
#' @param alpha_point he alpha-value indicating the transparency of the
#'   observations. Number between 0 (transparent) and 1 (not transparent).
#' @param low_color The color corresponding to correlation equal to -1 (default:
#'   blue).
#' @param high_color The color corresponding to correlation equal to 1 (default:
#'   red).
#' @param break_int Break interval in the color gradient.
#' @param label_size Size of text labels, if plotted.
#' @param font_family Font family used for text labels, if plotted.
#' @param point_size Size of points used for plotting the observations.
#' @param xlim x-limits
#' @param ylim y-limits
#' @param xlab x-label
#' @param ylab y-label
#' @param rholab Label for the legend, if plotted
#' @param main Title of plot
#' @param subtitle Subtitle of plot
#'
#' @references
#'
#' Berentsen, G. D., Kleppe, T. S., & Tjøstheim, D. (2014). Introducing
#' localgauss, an R package for estimating and visualizing local Gaussian
#' correlation. Journal of Statistical Software, 56(1), 1-18.
#'
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
#' York, 2009.
#'
#' @export
corplot <- function(dlg_object,
                    pair = 1,
                    gaussian_scale = FALSE,
                    plot_colormap = TRUE,
                    plot_obs = FALSE,
                    plot_labels = TRUE,
                    plot_legend = FALSE,
                    plot_thres = 0,
                    alpha_tile = .8,
                    alpha_point = .8,
                    low_color = "blue",
                    high_color = "red",
                    break_int = .2,
                    label_size = 3,
                    font_family = "sans",
                    point_size = NULL,
                    xlim = NULL,
                    ylim = NULL,
                    xlab = NULL,
                    ylab = NULL,
                    rholab = NULL,
                    main = NULL,
                    subtitle = NULL) {

    # First, we chack that the supplied dlg_object actually comes from the
    # dlg-function:
    if(class(dlg_object) != "dlg") {
        stop("dlg_object needs to have class 'dlg'")
    }

    # We also check that the pair number is actually a pair in this model
    if(!(pair %in% 1:nrow(dlg_object$bw$joint))) {
      stop(paste("'pair' must be an integer between 1 and ",
                 nrow(dlg_object$bw$joint),
                 " (the number of pairs in the model"))
    }

    # Next, we construct the data frame that contains the local correlations.
    # To do that, we extract the grid for this pair, and remove the duplicates,
    # because the pairwise local correlation will only depend on the pairwise
    # grid points.
    if(gaussian_scale) {
        full_grid <-
          dlg_object$transformed_grid[, unlist(dlg_object$bw$joint[pair, c(1, 2)])]
    } else {
        full_grid <-
          dlg_object$grid[, unlist(dlg_object$bw$joint[pair, c(1, 2)])]
    }

    # Where do we find unique values in the pairwise grid?
    distinct_grid <- !duplicated(full_grid, MARG = 1)

    # The data frame used for plotting the dependence map
    x = full_grid[distinct_grid, 1]
    y = full_grid[distinct_grid, 2]
    rho <- dlg_object$loc_cor[distinct_grid, pair]

    # Create the label for use in plot
    label = ifelse(!is.na(rho),
                   formatC(rho, digits = 2, format = "f", flag = "+"), "")

    # Create another data frame containing the observations
    if(gaussian_scale) {
        obs <-
          data.frame(x = dlg_object$transformed_data[, unlist(dlg_object$bw$joint[pair, 1])],
                     y = dlg_object$transformed_data[, unlist(dlg_object$bw$joint[pair, 2])])
    } else {
        obs <-
          data.frame(x = dlg_object$x[, unlist(dlg_object$bw$joint[pair, 1])],
                     y = dlg_object$x[, unlist(dlg_object$bw$joint[pair, 2])])
    }

    # Set the alpha equal to zero if we do not want the color map. If we do,
    # set selected grid points to no plot, because there is no data nearby.
    if(!plot_colormap) {
        alpha <- rep(0, length(x))
    } else if(plot_thres == 0) {
        alpha <- rep(alpha_tile, length(x))
    } else {
        kernel_density <- ks::kde(obs, eval.points = cbind(x, y))$estimate
        alpha <- rep(alpha_tile, length(x))
        alpha[kernel_density/max(kernel_density) < plot_thres] <- 0
    }

    # Create the data frame containing this information
    plot_data <- data.frame(x = x,
                            y = y,
                            rho = rho,
                            label = label,
                            alpha = alpha, stringsAsFactors = FALSE)


    # Initialize the plot.
    g = ggplot2::ggplot() +
      ggplot2::geom_tile(data = plot_data,
                         mapping = ggplot2:: aes(x = x, y = y, fill = rho), alpha = alpha)

    # Add the color gradient
    gr = ggplot2::scale_fill_gradient2(midpoint = 0,
                                       low = low_color,
                                       high = high_color,
                                       space = "Lab",
                                       limits = c(-1, 1),
                                       breaks = seq(-1, 1, by = break_int))
    g <- g + gr

    # Set axes limits if provided
    if(!is.null(xlim)) {
      g <- g + ggplot2::xlim(xlim)
    }

    if(!is.null(ylim)) {
      g <- g + ggplot2::ylim(ylim)
    }

    # Set axis labels if they are not provided
    if(is.null(xlab)) {
        if(is.null(colnames(dlg_object$x))) {
            x_label <- "X"
        } else {
            x_label <- colnames(dlg_object$x)[1]
        }
    } else {
        x_label <- xlab
    }

    if(is.null(ylab)) {
        if(is.null(colnames(dlg_object$x))) {
            y_label <- "Y"
        } else {
            y_label <- colnames(dlg_object$x)[2]
        }
    } else {
        y_label <- ylab
    }

    if(is.null(main)) {
        main_label <- paste("Local correlations for pair ",
                             pair,
                             " (variables ",
                             dlg_object$bw$joint[pair, 1],
                             " and ",
                             dlg_object$bw$joint[pair, 2],
                             ")", sep = "")
    } else {
      main_label <- main
    }

    if(is.null(rholab)) {
        rho_label <- "Rho"
    } else {
        rho_label <- rholab
    }

    g <- g + ggplot2::labs(fill = rho_label, x = x_label, y = y_label) +
      ggplot2::ggtitle(label = main_label, subtitle = subtitle)

    # Add the observations, if requested
    if(plot_obs) {
      if(is.null(point_size)) {
        g = g + ggplot2::geom_point(data = obs,
                                    mapping = ggplot2::aes(x = x,
                                                           y = y),
                                    alpha = alpha_point)
      } else {
        g = g + ggplot2::geom_point(data = obs,
                                    mapping = ggplot2::aes(x = x,
                                                           y = y),
                                    alpha = alpha_point,
                                    size = point_size)
      }
    }


    # Add the labels, if requested
    if(plot_labels) {
      g = g + ggplot2::geom_text(data = plot_data,
                                 mapping = ggplot2::aes(x = x, y = y, label = label),
                                 size = label_size,
                                 family = font_family,
                                 alpha = ifelse(alpha == 0, 0, 1))
    }

    # Remove the legend, if requested
    if(!plot_legend) {
        g = g + ggplot2::theme(legend.position = 'none')
    }

    return(g)
}
