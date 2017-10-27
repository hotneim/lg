
# Functions for plotting the local correlations
# ---------------------------------------------

#' Plot unconditional local correlations
#'
#' blabla
#'
#' blablajkaklajkl
#'

corplot <- function(dlg_object,
                    pair,
                    plot_map = TRUE,
                    gaussian_scale = FALSE,
                    plot_colormap = TRUE,
                    plot_obs = FALSE,
                    plot_labels = TRUE,
                    plot_legend = FALSE,
                    alpha_tile = .5,
                    alpha_point = .8,
                    low_color = "blue",
                    high_color = "red",
                    break_int = .2,
                    label_size = 3,
                    font_family = "sans",
                    point_size = NULL,
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

    # The data fram used for plotting the dependence map
    x = full_grid[distinct_grid, 1]
    y = full_grid[distinct_grid, 2]

    # Extract the local correlations, or let them have an out-of-range
    # value if we do not want to plot the color map.
    if(plot_colormap) {
        rho = dlg_object$loc_cor[distinct_grid, pair]
    } else {
        rho = rep(99, length(x))
    }

    # Create the label for use in plot
    label = ifelse(!is.na(rho),
                   formatC(rho, digits = 2, format = "f", flag = "+"), "")

    # Create the data frame containing this information
    plot_data <- data.frame(x = x,
                            y = y,
                            rho = rho,
                            label = label, stringsAsFactors = FALSE)

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

    # Initialize the plot.
    g = ggplot2::ggplot() +
      ggplot2::geom_tile(data = plot_data,
                         mapping = ggplot2:: aes(x = x, y = y, fill = rho),
                         alpha = alpha_tile)

    # Add the color gradient
    gr = ggplot2::scale_fill_gradient2(midpoint = 0,
                                       low = low_color,
                                       high = high_color,
                                       space = "Lab",
                                       limits = c(-1, 1),
                                       breaks = seq(-1, 1, by = break_int))
    g <- g + gr

    # Set axis labels id they are not provided
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
                                 family = font_family)
    }

    # Remove the legend, if requested
    if(!plot_legend) {
        g = g + ggplot2::theme(legend.position = 'none')
    }

    return(g)
}

#' Plot all pairs and collect in grid

arrange_corplots <- function(dlg_objects) {

}

