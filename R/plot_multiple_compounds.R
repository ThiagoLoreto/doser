#' Plot Multiple Dose-Response Curves Using ggplot2
#'
#' Creates composite plots showing multiple dose-response curves with professional
#' styling, flexible legend placement, and intelligent defaults. Optimized for
#' scientific publications with black-and-white or color output options.
#'
#' @param results List object returned by \code{\link{fit_dose_response}} containing
#'   dose-response analysis results.
#' @param compound_indices Numeric vector specifying which compounds to plot.
#'   If NULL, plots all available compounds (default: NULL).
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 150)).
#' @param point_shapes Numeric vector of point shapes for each compound. If NULL,
#'   generates optimal shapes automatically based on number of compounds (default: NULL).
#' @param colors Character vector of colors for each compound. If NULL, uses black
#'   for all compounds (default: NULL). Provide custom colors for colored output.
#' @param legend_position Character specifying legend position: "right", "left",
#'   "top", "bottom", or "none" (default: "right").
#' @param show_grid Logical indicating whether to show background grid (default: FALSE).
#' @param show_legend Logical indicating whether to show legend (default: TRUE).
#' @param save_plot Either a file path for saving the plot, or TRUE for automatic naming
#'   (default: NULL, no saving).
#' @param plot_width Plot width in inches for saved plots (default: 10).
#' @param plot_height Plot height in inches for saved plots (default: 8).
#' @param plot_dpi Resolution for saved raster images (default: 600).
#' @param axis_label_size Font size for axis labels (default: 14).
#' @param axis_text_size Font size for axis numbers (default: 14).
#' @param show_error_bars Logical indicating whether to show error bars (default: TRUE).
#' @param error_bar_width Width of error bar ends (default: 0.05).
#' @param plot_title Character string for plot title (default: "Multiple Dose-Response Curves").
#' @param legend_text_size Font size for legend text. If NULL, automatically adjusts
#'   based on number of compounds (default: NULL).
#' @param legend_title_size Font size for legend title (default: 11).
#' @param legend_key_height Height of legend key in points. If NULL, automatically
#'   adjusts based on number of compounds (default: NULL).
#' @param legend_ncol Number of columns in legend. If NULL, automatically determines
#'   optimal number (default: NULL).
#' @param legend_label_wrap Maximum characters per line in legend labels before
#'   wrapping (default: 25).
#' @param plot_ratio Numeric ratio of plot area allocated to main plot vs legend
#'   area (default: 0.7).
#' @param legend_title Character string for legend title. If NULL, no title is shown
#'   (default: NULL).
#'
#' @return Returns a ggplot object with comprehensive metadata stored as attributes:
#' \itemize{
#'   \item \code{compound_names}: Names of plotted compounds
#'   \item \code{compound_indices}: Indices of plotted compounds
#'   \item \code{n_compounds}: Number of compounds plotted
#'   \item \code{point_shapes}: Point shapes used for each compound
#'   \item \code{colors}: Colors used for each compound
#'   \item \code{point_size}: Point size used in plot
#'   \item \code{legend_position}: Legend position used
#'   \item \code{legend_ncol}: Number of columns in legend
#'   \item \code{legend_text_size}: Legend text size used
#'   \item \code{legend_key_height}: Legend key height used
#'   \item \code{legend_title}: Legend title used
#'   \item \code{x_limits}: X-axis limits used
#'   \item \code{plot_dimensions}: Dimensions of the plot (width, height, dpi)
#'   \item \code{file_saved}: Path to saved file if plot was saved
#' }
#'
#' @details
#' This function creates professional multi-curve dose-response plots using ggplot2
#' with intelligent defaults optimized for scientific publications. The default
#' styling uses black lines with different point shapes for optimal differentiation
#' in black-and-white publications, while supporting custom colors for presentations.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Black-and-white optimized}: Default uses black lines with distinct point shapes
#'   \item \strong{Smart point selection}: Automatically chooses optimal point shapes based on compound count
#'   \item \strong{Adaptive sizing}: Point sizes and legend elements adjust based on number of compounds
#'   \item \strong{Intelligent text wrapping}: Automatically wraps long compound names in legend
#'   \item \strong{Professional styling}: Clean, publication-ready appearance with customizable elements
#'   \item \strong{Self-contained}: No external package loading required
#' }
#'
#' \strong{Automatic Adjustments:}
#' \itemize{
#'   \item \strong{Point shapes}: Uses most distinguishable shapes first, recycles intelligently
#'   \item \strong{Point size}: Larger points for few compounds (3.5), smaller for many (2.5)
#'   \item \strong{Legend text}: Smaller text for many compounds (9pt), larger for few (11pt)
#'   \item \strong{Legend columns}: Single column for ≤10 compounds, two columns for >10 compounds
#'   \item \strong{X-axis limits}: Automatically calculated from data with 5% margin
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Default black-and-white with automatic point shapes
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#' 
#' # Plot compounds with default settings (black lines, different point shapes)
#' p <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = c(1, 3, 5, 7),
#'   plot_title = "Selected Compound Comparison"
#' )
#' print(p)
#' 
#' # Example 2: Custom colors for presentation
#' p <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:4,
#'   colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), # ColorBrewer palette
#'   plot_title = "Colored Dose-Response Curves",
#'   legend_title = "Test Compounds"
#' )
#' print(p)
#' 
#' # Example 3: Publication-ready with custom point shapes
#' p <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:6,
#'   point_shapes = c(15, 16, 17, 18, 0, 1), # Specific shapes
#'   plot_title = "Custom Point Shapes",
#'   legend_position = "bottom",
#'   legend_ncol = 3
#' )
#' print(p)
#' 
#' # Example 4: Many compounds with optimized layout
#' p <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:15,
#'   plot_title = "High-Throughput Screening",
#'   legend_text_size = 8,
#'   plot_width = 12,
#'   save_plot = "screening_results.png"
#' )
#' print(p)
#' 
#' # Example 5: Access and use plot metadata
#' p <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = c(2, 4, 6)
#' )
#' 
#' # Extract metadata for reproducibility
#' meta <- attr(p, "metadata")
#' cat("Plotted", meta$n_compounds, "compounds\n")
#' cat("Point shapes:", meta$point_shapes, "\n")
#' cat("X-axis range:", round(meta$x_limits, 2), "\n")
#' 
#' # Save styling information
#' styling_info <- data.frame(
#'   Compound = meta$compound_names,
#'   PointShape = meta$point_shapes,
#'   Color = meta$colors
#' )
#' write.csv(styling_info, "plot_styling.csv", row.names = FALSE)
#' }
#'
#' @seealso
#' \code{\link{fit_dose_response}} for generating input data
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @export






plot_multiple_compounds <- function(results, compound_indices = NULL, 
                                    y_limits = c(0, 150),
                                    point_shapes = NULL,
                                    colors = NULL,
                                    legend_position = "right",
                                    show_grid = FALSE, show_legend = TRUE,
                                    save_plot = NULL, plot_width = 10, 
                                    plot_height = 8, plot_dpi = 600,
                                    axis_label_size = 14,
                                    axis_text_size = 14,
                                    show_error_bars = TRUE,
                                    error_bar_width = 0.05,
                                    plot_title = "Multiple Dose-Response Curves",
                                    legend_text_size = NULL,
                                    legend_title_size = 11,
                                    legend_key_height = NULL,
                                    legend_ncol = NULL,
                                    legend_label_wrap = 25,
                                    plot_ratio = 0.7,
                                    legend_title = NULL) {
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it with: install.packages('ggplot2')")
  }
  
  # Input validation
  if (is.null(compound_indices)) {
    compound_indices <- seq_along(results$detailed_results)
  }
  
  n_compounds <- length(compound_indices)
  
  if (n_compounds == 0) {
    stop("No compounds selected for plotting")
  }
  
  if (max(compound_indices) > length(results$detailed_results)) {
    stop("Compound indices exceed available compounds")
  }
  
  # Validate parameters
  if (plot_ratio <= 0 || plot_ratio >= 1) {
    stop("plot_ratio must be between 0 and 1 (exclusive)")
  }
  
  valid_positions <- c("right", "left", "top", "bottom", "none")
  if (!legend_position %in% valid_positions) {
    stop("legend_position must be one of: '", paste(valid_positions, collapse = "', '"), "'")
  }
  
  # Smart defaults for legend based on number of compounds
  if (is.null(legend_text_size)) {
    legend_text_size <- if (n_compounds > 15) 9 else
      if (n_compounds > 8) 10 else 11
  }
  
  if (is.null(legend_key_height)) {
    legend_key_height <- if (n_compounds > 15) 12 else 15
  }
  
  # Extract compound data and prepare for plotting
  curve_data <- list()
  point_data <- list()
  
  for (i in seq_along(compound_indices)) {
    idx <- compound_indices[i]
    result <- results$detailed_results[[idx]]
    
    if (is.null(result$model) || !isTRUE(result$success)) {
      warning("Compound ", idx, ": no successful model fit - skipping")
      next
    }
    
    compound_name <- strsplit(result$compound, " \\| ")[[1]][1]
    
    # Prepare curve data
    data <- result$data
    valid_data <- data[is.finite(data$log_inhibitor) & is.finite(data$response), ]
    
    if (nrow(valid_data) < 2) {
      warning("Compound ", compound_name, ": insufficient data - skipping")
      next
    }
    
    # Generate fitted curve
    x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)
    pred_df <- data.frame(
      log_inhibitor = x_seq,
      response = predict(result$model, newdata = data.frame(log_inhibitor = x_seq)),
      compound = compound_name,
      compound_index = i
    )
    
    curve_data[[i]] <- pred_df
    
    # Prepare point data with means and SDs
    conc_levels <- unique(valid_data$log_inhibitor)
    
    point_stats <- data.frame(
      log_inhibitor = numeric(0),
      mean_response = numeric(0),
      sd_response = numeric(0),
      compound = character(0),
      compound_index = numeric(0)
    )
    
    for (conc in conc_levels) {
      conc_responses <- valid_data$response[valid_data$log_inhibitor == conc]
      if (length(conc_responses) > 0) {
        point_stats <- rbind(point_stats, data.frame(
          log_inhibitor = conc,
          mean_response = mean(conc_responses, na.rm = TRUE),
          sd_response = sd(conc_responses, na.rm = TRUE),
          compound = compound_name,
          compound_index = i
        ))
      }
    }
    
    point_data[[i]] <- point_stats
  }
  
  # Combine all data
  all_curve_data <- do.call(rbind, curve_data)
  all_point_data <- do.call(rbind, point_data)
  
  if (is.null(all_curve_data) || nrow(all_curve_data) == 0) {
    stop("No valid data to plot")
  }
  
  # Smart point shape selection
  if (is.null(point_shapes)) {
    optimal_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14)
    
    if (n_compounds <= 6) {
      point_shapes <- optimal_shapes[1:n_compounds]
    } else {
      point_shapes <- rep(optimal_shapes, length.out = n_compounds)
    }
  } else {
    point_shapes <- rep(point_shapes, length.out = n_compounds)
  }
  
  # Color logic: default is black, but allows custom colors
  use_colors <- if (!is.null(colors)) {
    if (length(colors) < n_compounds) {
      warning("Number of colors provided (", length(colors), ") is less than number of compounds (", n_compounds, "). Recycling colors.")
      rep(colors, length.out = n_compounds)
    } else {
      colors[1:n_compounds]
    }
  } else {
    rep("black", n_compounds)
  }
  
  # Smart text wrapping for long compound names
  smart_label_wrap <- function(labels, width = legend_label_wrap) {
    sapply(labels, function(label) {
      if (nchar(label) <= width) return(label)
      
      # Try to break at hyphens or underscores first
      if (grepl("[-_]", label)) {
        parts <- strsplit(label, "[-_]")[[1]]
        if (length(parts) > 1 && any(nchar(parts) <= width)) {
          best_break <- which(cumsum(nchar(parts) + 1) <= width)
          if (length(best_break) > 0) {
            break_point <- max(best_break)
            line1 <- paste(parts[1:break_point], collapse = "-")
            line2 <- paste(parts[(break_point + 1):length(parts)], collapse = "-")
            return(paste(line1, line2, sep = "\n"))
          }
        }
      }
      
      # Break by words
      words <- strsplit(label, " ")[[1]]
      if (length(words) > 1) {
        lines <- character(0)
        current_line <- ""
        
        for (word in words) {
          test_line <- if (current_line == "") word else paste(current_line, word)
          if (nchar(test_line) <= width) {
            current_line <- test_line
          } else {
            if (current_line != "") lines <- c(lines, current_line)
            current_line <- word
          }
        }
        if (current_line != "") lines <- c(lines, current_line)
        
        if (length(lines) <= 3) {
          return(paste(lines, collapse = "\n"))
        }
      }
      
      # Simple break for very long single words
      if (nchar(label) > width * 2) {
        part1 <- substr(label, 1, width)
        part2 <- substr(label, width + 1, width * 2)
        part3 <- substr(label, (width * 2) + 1, nchar(label))
        return(paste(part1, part2, part3, sep = "\n"))
      } else {
        mid <- ceiling(nchar(label) / 2)
        space_pos <- gregexpr(" ", substr(label, mid - 5, mid + 5))[[1]]
        if (any(space_pos > 0)) {
          break_point <- mid - 6 + min(space_pos[space_pos > 0])
        } else {
          break_point <- mid
        }
        return(paste(substr(label, 1, break_point), 
                     substr(label, break_point + 1, nchar(label)), sep = "\n"))
      }
    })
  }
  
  # Apply text wrapping
  compound_labels <- unique(all_curve_data$compound)
  wrapped_labels <- smart_label_wrap(compound_labels, legend_label_wrap)
  
  # Calculate optimal X-axis limits
  calculate_x_limits <- function(curve_data, point_data) {
    all_x <- c(curve_data$log_inhibitor, point_data$log_inhibitor)
    valid_x <- all_x[is.finite(all_x)]
    
    if (length(valid_x) == 0) return(c(-10, -2))
    
    x_range <- range(valid_x)
    x_margin <- diff(x_range) * 0.05
    return(c(x_range[1] - x_margin, x_range[2] + x_margin))
  }
  
  x_limits <- calculate_x_limits(all_curve_data, all_point_data)
  
  # Adaptive point size based on number of compounds
  point_size <- if (n_compounds > 15) 2.5 else 
    if (n_compounds > 8) 3 else 3.5
  
  # Final legend title (NULL = no title)
  legend_title_final <- if (!is.null(legend_title)) legend_title else NULL
  
  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::labs(
      x = expression(paste("Log"[10], " Concentration [M]")),
      y = ifelse(results$normalized, "Normalized BRET ratio [%]", "BRET ratio"),
      title = plot_title,
      shape = legend_title_final
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = axis_label_size, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      plot.title = ggplot2::element_text(size = axis_label_size + 2, face = "bold", hjust = 0.5, color = "black"),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = legend_text_size, lineheight = 0.8, color = "black"),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold", color = "black"),
      legend.key.height = ggplot2::unit(legend_key_height, "points"),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = ifelse(show_grid, "grey90", "white")),
      panel.grid.minor = ggplot2::element_line(color = ifelse(show_grid, "grey95", "white")),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Add plot elements with color control
  if (!is.null(colors)) {
    # COLOR MODE: colored lines and points
    p <- p +
      ggplot2::geom_line(data = all_curve_data, 
                         aes(x = log_inhibitor, y = response, group = compound, color = compound),
                         linewidth = 1, alpha = 0.7) +
      ggplot2::geom_point(data = all_point_data,
                          aes(x = log_inhibitor, y = mean_response, shape = compound, color = compound),
                          size = point_size) +
      ggplot2::scale_color_manual(
        values = use_colors,
        labels = wrapped_labels,
        guide = "none"
      ) +
      ggplot2::scale_shape_manual(
        values = point_shapes[1:n_compounds],
        labels = wrapped_labels
      )
    
    if (show_error_bars && nrow(all_point_data) > 0) {
      p <- p +
        ggplot2::geom_errorbar(data = all_point_data,
                               aes(x = log_inhibitor, 
                                   ymin = mean_response - sd_response, 
                                   ymax = mean_response + sd_response,
                                   group = compound, color = compound),
                               width = error_bar_width, 
                               linewidth = 0.5)
    }
  } else {
    # DEFAULT MODE: all black
    p <- p +
      ggplot2::geom_line(data = all_curve_data, 
                         aes(x = log_inhibitor, y = response, group = compound),
                         color = "black", linewidth = 1, alpha = 0.7) +
      ggplot2::geom_point(data = all_point_data,
                          aes(x = log_inhibitor, y = mean_response, shape = compound),
                          color = "black", size = point_size) +
      ggplot2::scale_shape_manual(
        values = point_shapes[1:n_compounds],
        labels = wrapped_labels
      )
    
    if (show_error_bars && nrow(all_point_data) > 0) {
      p <- p +
        ggplot2::geom_errorbar(data = all_point_data,
                               aes(x = log_inhibitor, 
                                   ymin = mean_response - sd_response, 
                                   ymax = mean_response + sd_response,
                                   group = compound),
                               color = "black",
                               width = error_bar_width, 
                               linewidth = 0.5)
    }
  }
  
  # Configure legend columns
  if (!is.null(legend_ncol)) {
    ncol_final <- legend_ncol
  } else {
    ncol_final <- if (n_compounds > 10) 2 else 1
  }
  
  # Configure legend appearance
  if (!is.null(colors)) {
    p <- p + ggplot2::guides(
      shape = ggplot2::guide_legend(
        ncol = ncol_final,
        override.aes = list(
          color = use_colors,
          size = point_size
        )
      )
    )
  } else {
    p <- p + ggplot2::guides(
      shape = ggplot2::guide_legend(
        ncol = ncol_final,
        override.aes = list(
          color = "black",
          size = point_size
        )
      )
    )
  }
  
  # Hide legend if requested
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  # Save plot if requested
  if (!is.null(save_plot)) {
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      indices_str <- paste(compound_indices, collapse = "_")
      filename <- paste0("multiple_compounds_", indices_str, ".png")
    } else {
      stop("save_plot must be either a file path or TRUE for auto-naming")
    }
    
    plot_dir <- dirname(filename)
    if (plot_dir != "." && !dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    ggplot2::ggsave(filename, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi)
    message("Plot successfully saved as: ", normalizePath(filename))
  }
  
  # Return plot and metadata
  metadata <- list(
    compound_names = unique(all_curve_data$compound),
    compound_indices = compound_indices,
    n_compounds = n_compounds,
    point_shapes = point_shapes[1:n_compounds],
    colors = use_colors,
    point_size = point_size,
    legend_position = legend_position,
    legend_ncol = ncol_final,
    legend_text_size = legend_text_size,
    legend_key_height = legend_key_height,
    legend_title = legend_title_final,
    x_limits = x_limits,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    file_saved = if (!is.null(save_plot)) filename else NULL
  )
  
  attr(p, "metadata") <- metadata
  return(p)
}
