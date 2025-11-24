#' Plot Multiple Dose-Response Curves
#'
#' The `plot_multiple_compounds()` function generates a consolidated plot of 
#' fitted dose-response curves for multiple compounds, allowing visual 
#' comparison between different responses. It provides extensive customization 
#' options for colors, shapes, titles, legends, gridlines, and file export.
#'
#' @param results A list containing results from dose-response analysis. 
#'   Expected to include a sublist `detailed_results`, where each 
#'   element contains at least `model`, `data`, `success`, and `compound`.
#' @param compound_indices Numeric vector specifying which compounds to include in the plot
#'   If `NULL` (default), all available compounds in `results$detailed_results` 
#'   are included.
#' @param y_limits Numeric vector of length 2 specifying the y-axis limits
#' @param point_shapes Numeric vector of point shapes for different compounds
#' @param colors Character vector of colors for different compound curves
#'   - `FALSE` (default): all curves and points are black;
#'   - `TRUE`: automatically assigns distinct colors (via `scales::hue_pal()`);
#'   - Character vector of custom colors (recycled if shorter than the number of compounds).
#' @param legend_position Position of the legend: one of `"right"`, `"left"`, 
#'   `"top"`, `"bottom"`, or `"none"`. Default: `"right"`.
#' @param show_grid Logical indicating whether to show background grid lines
#' @param show_legend Logical indicating whether to display the legend
#' @param save_plot Defines whether to save the plot: \code{NULL} (do not save, default), 
#'   \code{TRUE} (automatically saves as PNG with default name), or a file path with extension 
#'   (\code{.png}, \code{.pdf}, \code{.jpeg}, \code{.tiff}, \code{.svg}, \code{.eps}) to save in a specific format.
#' @param plot_width,plot_height Plot dimensions (in inches) when saving.
#' @param plot_dpi Resolution (in DPI) for saved plots. Default: `600`.
#' @param axis_label_size Numeric value for axis title font size
#' @param axis_text_size Numeric value for axis tick label font size
#' @param show_error_bars Logical indicating whether to display error bars around data points
#' @param error_bar_width Numeric value controlling the width of error bars
#' @param plot_title Custom plot title. If `NULL`, a smart title is generated 
#'   automatically based on the compound names.
#' @param legend_text_size Numeric value for legend text font size
#' @param legend_title_size Numeric value for legend title font size
#' @param legend_key_height Height (in points) of legend keys. Automatically adjusted.
#' @param legend_ncol Numeric value specifying number of columns in legend
#' @param legend_label_wrap Maximum character width before legend labels 
#'   automatically wrap to new lines. Default: `25`.
#' @param plot_ratio Numeric value for plot aspect ratio (width/height)
#' @param legend_title Title for the legend (displayed above symbols).
#' @param legend_labels Character vector of custom labels for legend (overrides target/compound names)
#'
#'
#'@importFrom ggplot2 aes
#'
#'
#' @return ggplot object with additional metadata stored as attribute containing:
#'   - compound_names: Original compound names
#'   - compound_indices: Indices used for plotting
#'   - n_compounds: Number of valid compounds plotted
#'   - point_shapes: Point shapes used
#'   - colors: Colors used
#'   - legend_position: Legend position
#'   - legend_ncol: Number of legend columns
#'   - legend_text_size: Legend text size
#'   - legend_key_height: Legend key height
#'   - legend_title: Legend title used
#'   - legend_labels_used: Legend labels actually used
#'   - wrapped_labels: Labels after text wrapping applied
#'   - x_limits: X-axis limits used
#'   - plot_dimensions: Plot dimensions and resolution
#'   - file_saved: Filename if plot was saved
#'   - smart_title_used: Title actually displayed on plot
#'
#' @details
#' This function overlays fitted dose-response curves (based on nonlinear models)
#' together with empirical mean ± SD values for each concentration, allowing
#' direct visual comparison across multiple compounds or experimental conditions.
#'
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
#'   \item \strong{Legend columns}: Single column for <=10 compounds, two columns for >10 compounds
#'   \item \strong{X-axis limits}: Automatically calculated from data with 5% margin
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic plot for all compounds
#' p1 <- plot_multiple_compounds(results)
#' print(p1)
#'
#' # Example 2: Select specific compounds by index
#' p2 <- plot_multiple_compounds(results, compound_indices = c(1, 3, 5))
#'
#' # Example 3: Enable automatic coloring and save the plot
#' p3 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:4,
#'   colors = TRUE,
#'   save_plot = "plots/multi_colored_curves.png"
#' )
#'
#' # Example 4: Customize legend labels and point shapes
#' p4 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:3,
#'   legend_labels = c("Compound A", "Compound B", "Compound C"),
#'   point_shapes = c(15, 17, 19),
#'   legend_title = "Treatments"
#' )
#'
#' # Example 5: Disable error bars and use black-and-white mode
#' p5 <- plot_multiple_compounds(
#'   results,
#'   colors = FALSE,
#'   show_error_bars = FALSE,
#'   show_grid = TRUE,
#'   plot_title = "Curves without error bars"
#' )
#'
#' # Example 6: Place legend below the plot in multiple columns
#' p6 <- plot_multiple_compounds(
#'   results,
#'   legend_position = "bottom",
#'   legend_ncol = 3,
#'   colors = TRUE
#' )
#'
#' # Example 7: Plot without any legend
#' p7 <- plot_multiple_compounds(results, show_legend = FALSE)
#'
#' # Example 8: Fine-tune fonts, limits, and title
#' p8 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:2,
#'   y_limits = c(0, 120),
#'   axis_label_size = 16,
#'   axis_text_size = 13,
#'   plot_title = "Comparison of Two Compounds"
#' )
#' 
#' 
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
#' \code{\link{fit_drc_3pl}} for generating input data
#' \code{\link{fit_drc_4pl}} for generating input data
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @export


plot_multiple_compounds <- function(results, compound_indices = NULL, 
                                    y_limits = c(0, 150),
                                    point_shapes = NULL,
                                    colors = FALSE,
                                    legend_position = "right",
                                    show_grid = FALSE, show_legend = TRUE,
                                    save_plot = NULL, plot_width = 10, 
                                    plot_height = 8, plot_dpi = 600,
                                    axis_label_size = 14,
                                    axis_text_size = 14,
                                    show_error_bars = TRUE,
                                    error_bar_width = 0.05,
                                    plot_title = NULL,
                                    legend_text_size = NULL,
                                    legend_title_size = 11,
                                    legend_key_height = NULL,
                                    legend_ncol = NULL,
                                    legend_label_wrap = 25,
                                    plot_ratio = 0.7,
                                    legend_title = NULL,
                                    legend_labels = NULL) {
  
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
  
  if (!is.null(legend_labels)) {
    if (length(legend_labels) != n_compounds) {
      stop("legend_labels must have the same length as compound_indices (", n_compounds, ")")
    }
  }
  
  if (plot_ratio <= 0 || plot_ratio >= 1) {
    stop("plot_ratio must be between 0 and 1 (exclusive)")
  }
  
  valid_positions <- c("right", "left", "top", "bottom", "none")
  if (!legend_position %in% valid_positions) {
    stop("legend_position must be one of: '", paste(valid_positions, collapse = "', '"), "'")
  }
  
  # Smart defaults for legend
  if (is.null(legend_text_size)) {
    legend_text_size <- if (n_compounds > 15) 9 else if (n_compounds > 8) 10 else 11
  }
  
  if (is.null(legend_key_height)) {
    legend_key_height <- if (n_compounds > 15) 12 else 15
  }
  
  # Extract compound data
  curve_data <- list()
  point_data <- list()
  compound_names_original <- character(n_compounds)
  valid_compound_indices <- numeric(0)
  
  for (i in seq_along(compound_indices)) {
    idx <- compound_indices[i]
    result <- results$detailed_results[[idx]]
    
    if (is.null(result$model) || !isTRUE(result$success)) {
      warning("Compound ", idx, ": no successful model fit - skipping")
      next
    }
    
    compound_name <- strsplit(result$compound, " \\| ")[[1]][1]
    compound_names_original[i] <- compound_name
    valid_compound_indices <- c(valid_compound_indices, i)
    
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
  
  n_valid_compounds <- length(valid_compound_indices)
  
  # Intelligent title generation
  generate_smart_title <- function(compound_names) {
    if (length(compound_names) == 0) {
      return("Multiple Dose-Response Curves")
    }
    
    extract_parts <- function(full_name) {
      parts <- strsplit(full_name, ":")[[1]]
      if (length(parts) == 2) {
        return(list(target = parts[1], compound = parts[2]))
      } else {
        return(list(target = full_name, compound = full_name))
      }
    }
    
    parts_list <- lapply(compound_names, extract_parts)
    targets <- sapply(parts_list, function(x) x$target)
    compounds <- sapply(parts_list, function(x) x$compound)
    
    # Remove _2, _3 suffixes for comparison
    compounds_base <- gsub("_\\d+$", "", compounds)
    
    # Check if all have same compound (ignoring suffixes)
    unique_compounds_base <- unique(compounds_base)
    if (length(unique_compounds_base) == 1) {
      return(unique_compounds_base[1])
    }
    
    # Check if all have same target
    unique_targets <- unique(targets)
    if (length(unique_targets) == 1) {
      return(unique_targets[1])
    }
    
    return("Multiple Dose-Response Curves")
  }
  
  # Apply smart title logic
  if (is.null(plot_title)) {
    plot_title_final <- generate_smart_title(unique(all_curve_data$compound))
  } else {
    plot_title_final <- plot_title
  }
  
  # Point shape selection
  if (is.null(point_shapes)) {
    optimal_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14)
    
    if (n_valid_compounds <= 6) {
      point_shapes <- optimal_shapes[1:n_valid_compounds]
    } else {
      point_shapes <- rep(optimal_shapes, length.out = n_valid_compounds)
    }
  } else {
    point_shapes <- rep(point_shapes, length.out = n_valid_compounds)
  }
  
  # Color logic - handles 3 cases: FALSE (black), TRUE (auto colors), custom colors
  use_colors <- if (is.logical(colors)) {
    if (isTRUE(colors)) {
      # Automatic colorful mode
      scales::hue_pal()(n_valid_compounds)
    } else {
      # Black mode (default when colors = FALSE)
      rep("black", n_valid_compounds)
    }
  } else if (is.character(colors)) {
    # Custom colors mode
    if (length(colors) < n_valid_compounds) {
      warning("Number of colors provided (", length(colors), ") is less than number of valid compounds (", n_valid_compounds, "). Recycling colors.")
      rep(colors, length.out = n_valid_compounds)
    } else {
      colors[1:n_valid_compounds]
    }
  } else {
    # Fallback: black
    rep("black", n_valid_compounds)
  }
  
  # Legend labels
  if (!is.null(legend_labels)) {
    compound_labels <- legend_labels[valid_compound_indices]
  } else {
    compound_labels <- unique(all_curve_data$compound)
  }
  
  # Text wrapping for long labels
  smart_label_wrap <- function(labels, width = legend_label_wrap) {
    sapply(labels, function(label) {
      if (is.na(label) || nchar(label) <= width) return(label)
      
      # Try to break at hyphens or underscores first
      if (grepl("[-_]", label)) {
        parts <- strsplit(label, "[-_]")[[1]]
        if (length(parts) > 1) {
          current_length <- 0
          break_points <- numeric(0)
          for (j in seq_along(parts)) {
            current_length <- current_length + nchar(parts[j]) + 1
            if (current_length <= width + 1) {
              break_points <- c(break_points, j)
            }
          }
          
          if (length(break_points) > 0) {
            break_point <- max(break_points)
            if (break_point < length(parts)) {
              line1 <- paste(parts[1:break_point], collapse = "-")
              line2 <- paste(parts[(break_point + 1):length(parts)], collapse = "-")
              return(paste(line1, line2, sep = "\n"))
            }
          }
        }
      }
      
      # Break by words
      words <- strsplit(label, " ")[[1]]
      if (length(words) > 1) {
        lines <- character(0)
        current_line <- words[1]
        
        for (j in 2:length(words)) {
          test_line <- paste(current_line, words[j])
          if (nchar(test_line) <= width) {
            current_line <- test_line
          } else {
            lines <- c(lines, current_line)
            current_line <- words[j]
          }
        }
        lines <- c(lines, current_line)
        
        if (length(lines) <= 3) {
          return(paste(lines, collapse = "\n"))
        } else {
          mid <- ceiling(length(words) / 2)
          line1 <- paste(words[1:mid], collapse = " ")
          line2 <- paste(words[(mid + 1):length(words)], collapse = " ")
          return(paste(line1, line2, sep = "\n"))
        }
      }
      
      # Force break for very long single words
      if (nchar(label) > width) {
        break_point <- ceiling(nchar(label) / 2)
        part1 <- substr(label, 1, break_point)
        part2 <- substr(label, break_point + 1, nchar(label))
        return(paste(part1, part2, sep = "\n"))
      }
      
      return(label)
    }, USE.NAMES = FALSE)
  }
  
  wrapped_labels <- smart_label_wrap(compound_labels, legend_label_wrap)
  
  # Factor ordering for consistent legend
  all_curve_data$compound <- factor(all_curve_data$compound, 
                                    levels = unique(all_curve_data$compound))
  
  # Map custom legend labels if provided
  if (!is.null(legend_labels)) {
    name_mapping <- setNames(compound_labels, unique(all_curve_data$compound))
    
    all_curve_data$compound_display <- name_mapping[as.character(all_curve_data$compound)]
    all_point_data$compound_display <- name_mapping[as.character(all_point_data$compound)]
    
    all_curve_data$compound_display <- factor(all_curve_data$compound_display, 
                                              levels = compound_labels)
    all_point_data$compound_display <- factor(all_point_data$compound_display, 
                                              levels = compound_labels)
  } else {
    all_curve_data$compound_display <- all_curve_data$compound
    all_point_data$compound_display <- all_point_data$compound
  }
  
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
  
  # Adaptive point size
  point_size <- if (n_valid_compounds > 15) 2.5 else if (n_valid_compounds > 8) 3 else 3.5
  
  # Final legend title
  legend_title_final <- if (!is.null(legend_title)) legend_title else NULL
  
  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::labs(
      x = expression(paste("Log"[10], " Concentration [M]")),
      y = ifelse(results$normalized, "Normalized BRET ratio [%]", "BRET ratio"),
      title = plot_title_final,
      color = legend_title_final
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
  
  # Add plot elements 
  p <- p +
    ggplot2::geom_line(data = all_curve_data, 
                       ggplot2::aes(x = log_inhibitor, y = response, 
                                    group = compound_display, 
                                    color = compound_display), 
                       linewidth = 1, alpha = 0.7) +
    ggplot2::geom_point(data = all_point_data,
                        ggplot2::aes(x = log_inhibitor, y = mean_response, 
                                     shape = compound_display, 
                                     color = compound_display),
                        size = point_size) +
    ggplot2::scale_color_manual(
      values = use_colors,
      labels = wrapped_labels
    ) +
    ggplot2::scale_shape_manual(
      values = point_shapes[1:n_valid_compounds],
      guide = "none" 
    )
  
  if (show_error_bars && nrow(all_point_data) > 0) {
    p <- p +
      ggplot2::geom_errorbar(data = all_point_data,
                             ggplot2::aes(x = log_inhibitor, 
                                          ymin = mean_response - sd_response, 
                                          ymax = mean_response + sd_response,
                                          group = compound_display, 
                                          color = compound_display), 
                             width = error_bar_width, 
                             linewidth = 0.5)
  }
  
  # Configure legend columns
  if (!is.null(legend_ncol)) {
    ncol_final <- legend_ncol
  } else {
    ncol_final <- if (n_valid_compounds > 10) 2 else 1
  }
  
  # Configure legend appearance 
  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(
      ncol = ncol_final,
      override.aes = list(
        shape = point_shapes[1:n_valid_compounds],
        size = point_size,                      
        linetype = 1,                          
        linewidth = 0.8,                       
        fill = NA                                 
      ),
      title = legend_title_final
    )
  )
  
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
    n_compounds = n_valid_compounds,
    point_shapes = point_shapes[1:n_valid_compounds],
    colors = use_colors,
    point_size = point_size,
    legend_position = legend_position,
    legend_ncol = ncol_final,
    legend_text_size = legend_text_size,
    legend_key_height = legend_key_height,
    legend_title = legend_title_final,
    legend_labels_used = compound_labels,
    wrapped_labels = wrapped_labels,
    x_limits = x_limits,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    file_saved = if (!is.null(save_plot)) filename else NULL,
    smart_title_used = plot_title_final,
    color_mode = if(isTRUE(colors)) "colorful" else "black"
  )
  
  attr(p, "metadata") <- metadata
  return(p)
}