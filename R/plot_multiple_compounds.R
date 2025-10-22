#' Plot Multiple Dose-Response Curves in a Single Graph
#'
#' Creates composite plots showing multiple dose-response curves with advanced
#' differentiation methods, professional styling, and flexible legend placement.
#' Ideal for comparing multiple compounds or conditions in pharmacological studies.
#'
#' @param results List object returned by \code{\link{fit_dose_response}} containing
#'   dose-response analysis results.
#' @param compound_indices Numeric vector specifying which compounds to plot.
#'   If NULL, plots all available compounds (default: NULL).
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 150)).
#' @param colors Character vector of colors for each compound. If NULL, generates
#'   appropriate colors automatically (default: NULL).
#' @param line_types Numeric vector of line types for each compound (default: NULL).
#' @param point_shapes Numeric vector of point shapes for each compound (default: NULL).
#' @param differentiation_method Character specifying how to differentiate curves:
#'   "color", "linetype", "pointshape", or "combined" (default: "pointshape").
#' @param legend_position Character specifying legend position: "bottomright", "bottom",
#'   "bottomleft", "left", "topleft", "top", "topright", "right", or "outside"
#'   (default: "outside").
#' @param show_grid Logical indicating whether to show background grid (default: FALSE).
#' @param show_legend Logical indicating whether to show legend (default: TRUE).
#' @param save_plot Either a file path for saving the plot, or TRUE for automatic naming
#'   (default: NULL, no saving).
#' @param plot_width Plot width in inches for saved plots (default: 10).
#' @param plot_height Plot height in inches for saved plots (default: 8).
#' @param plot_dpi Resolution for saved raster images (default: 600).
#' @param axis_label_cex Character expansion factor for axis labels (default: 1.4).
#' @param axis_number_cex Character expansion factor for axis numbers (default: 1.4).
#' @param auto_combine_threshold Numeric threshold for automatically switching to
#'   "combined" differentiation method (default: 12).
#' @param show_error_bars Logical indicating whether to show error bars (default: TRUE).
#' @param error_bar_width Width of error bar ends (default: 0.03).
#' @param error_bar_lwd Line width for error bars (default: 1).
#' @param plot_title Character string for plot title (default: "Multiple Dose-Response Curves").
#' @param legend_cex Character expansion factor for legend text (default: 0.8).
#' @param legend_area_ratio Numeric ratio of plot area allocated to external legend
#'   when legend_position = "outside" (default: 0.25).
#' @param legend_point_cex Point size multiplier in legend (default: 1.0).
#'
#' @return Invisibly returns a list containing comprehensive plot metadata:
#' \itemize{
#'   \item \code{compound_names}: Names of plotted compounds
#'   \item \code{compound_indices}: Indices of plotted compounds
#'   \item \code{n_compounds}: Number of compounds plotted
#'   \item \code{plot_limits}: List with x and y axis limits used
#'   \item \code{differentiation_method}: Method used for curve differentiation
#'   \item \code{styling}: List with colors, line types, and point shapes used
#'   \item \code{error_bars}: Error bar configuration
#'   \item \code{legend_position}: Legend position used
#'   \item \code{legend_settings}: Legend configuration parameters
#'   \item \code{plot_title}: Plot title used
#'   \item \code{file_saved}: Path to saved file if plot was saved
#'   \item \code{file_format}: Format of saved file
#'   \item \code{plot_dimensions}: Dimensions of the plot (width, height, dpi)
#'   \item \code{timestamp}: Time when plot was generated
#' }
#'
#' @details
#' This function creates sophisticated multi-curve dose-response plots with
#' intelligent automatic styling and professional presentation. It automatically
#' handles curve differentiation, legend placement, and output formatting for
#' publication-quality figures.
#'
#' \strong{Automatic Features:}
#' \itemize{
#'   \item \strong{Smart Differentiation}: Auto-switches to combined methods for many compounds
#'   \item \strong{Color Management}: Generates distinct color palettes based on compound count
#'   \item \strong{Limit Calculation}: Automatically determines optimal axis limits from data
#'   \item \strong{Error Bar Handling}: Only shows error bars when meaningful data exists
#'   \item \strong{Layout Optimization}: Adjusts plot layout for external legends
#' }
#'
#' \strong{Differentiation Methods:}
#' \itemize{
#'   \item \strong{color}: Uses distinct colors (optimal for 2-8 compounds)
#'   \item \strong{linetype}: Uses line types 1-6 (repeats after 6 compounds)
#'   \item \strong{pointshape}: Uses point shapes (optimal for 2-15 compounds)
#'   \item \strong{combined}: Uses colors + line types + point shapes (best for 8+ compounds)
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Publication-ready multi-panel comparison
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#' 
#' # Create comparison of different compound classes
#' control_compounds <- c(1, 2, 3)   # Reference compounds
#' test_compounds <- c(4, 5, 6, 7)   # Experimental compounds
#' 
#' # Plot controls
#' plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = control_compounds,
#'   differentiation_method = "color",
#'   colors = c("gray40", "gray60", "gray80"),
#'   plot_title = "Control Compounds",
#'   legend_position = "bottomright"
#' )
#' 
#' # Plot test compounds with bright colors
#' plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = test_compounds,
#'   differentiation_method = "color",
#'   colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), # ColorBrewer colors
#'   plot_title = "Experimental Compounds",
#'   legend_position = "bottomright"
#' )
#' 
#' # Example 2: Black and white publication figure
#' plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:6,
#'   differentiation_method = "combined",
#'   colors = rep("black", 6),  # All black lines
#'   line_types = 1:6,          # Different line types
#'   point_shapes = 15:20,      # Different point shapes
#'   plot_title = "Dose-Response Curves (Black & White)",
#'   save_plot = "bw_figure.tiff",
#'   plot_dpi = 1200
#' )
#' 
#' # Example 3: Large dataset with external legend
#' plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:20,
#'   differentiation_method = "combined",
#'   legend_position = "outside",
#'   legend_area_ratio = 0.4,    # More space for legend
#'   legend_cex = 0.6,           # Smaller legend text
#'   plot_title = "High-Throughput Screening Results",
#'   plot_width = 14,            # Wider plot for many compounds
#'   plot_height = 8
#' )
#' 
#' # Example 4: Custom styling for specific comparisons
#' custom_colors <- c(
#'   "Compound_A" = "#1B9E77",
#'   "Compound_B" = "#D95F02", 
#'   "Compound_C" = "#7570B3",
#'   "Compound_D" = "#E7298A"
#' )
#' 
#' plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = c(2, 5, 8, 11),
#'   colors = custom_colors,
#'   line_types = c(1, 2, 1, 2),      # Alternating line types
#'   point_shapes = c(16, 17, 15, 18), # Distinct point shapes
#'   plot_title = "Selected Compound Comparison",
#'   show_error_bars = TRUE
#' )
#' 
#' # Example 5: Access and use detailed metadata
#' plot_meta <- plot_multiple_compounds(
#'   analysis_results,
#'   compound_indices = 1:8,
#'   differentiation_method = "combined"
#' )
#' 
#' # Create a plot summary report
#' cat("Multi-Curve Plot Summary:\n")
#' cat("Compounds plotted:", plot_meta$n_compounds, "\n")
#' cat("Differentiation method:", plot_meta$differentiation_method, "\n")
#' cat("X-axis range:", round(plot_meta$plot_limits$x_limits, 2), "\n")
#' cat("Y-axis range:", plot_meta$plot_limits$y_limits, "\n")
#' cat("Colors used:", length(unique(plot_meta$styling$colors)), "unique colors\n")
#' cat("Generated:", format(plot_meta$timestamp, "%Y-%m-%d %H:%M"), "\n")
#' 
#' # Save styling information for reproducibility
#' write.csv(
#'   data.frame(
#'     Compound = plot_meta$compound_names,
#'     Color = plot_meta$styling$colors,
#'     LineType = plot_meta$styling$line_types,
#'     PointShape = plot_meta$styling$point_shapes
#'   ),
#'   "plot_styling_reference.csv",
#'   row.names = FALSE
#' )
#' }






plot_multiple_compounds <- function(results, compound_indices = NULL, 
                                    y_limits = c(0, 150), colors = NULL,
                                    line_types = NULL, point_shapes = NULL,
                                    differentiation_method = "pointshape",
                                    legend_position = "outside",
                                    show_grid = FALSE, show_legend = TRUE,
                                    save_plot = NULL, plot_width = 10, 
                                    plot_height = 8, plot_dpi = 600,
                                    axis_label_cex = 1.4, axis_number_cex = 1.4,
                                    auto_combine_threshold = 12,
                                    show_error_bars = TRUE,
                                    error_bar_width = 0.03,
                                    error_bar_lwd = 1,
                                    plot_title = "Multiple Dose-Response Curves",
                                    legend_cex = 0.8,                    
                                    legend_area_ratio = 0.25,         
                                    legend_point_cex = 1.0) {         
  
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
  
  if (!is.character(plot_title) || length(plot_title) != 1) {
    stop("plot_title must be a single character string")
  }
  
  # Validate legend position
  valid_positions <- c("bottomright", "bottom", "bottomleft", "left", 
                       "topleft", "top", "topright", "right", "outside")
  if (!legend_position %in% valid_positions) {
    stop("legend_position must be one of: '", paste(valid_positions, collapse = "', '"), "'")
  }
  
  # Validate legend_area_ratio
  if (legend_area_ratio <= 0 || legend_area_ratio >= 1) {
    stop("legend_area_ratio must be between 0 and 1 (exclusive)")
  }
  
  # Validate differentiation method
  valid_methods <- c("color", "linetype", "pointshape", "combined")
  if (!differentiation_method %in% valid_methods) {
    stop("differentiation_method must be one of: '", paste(valid_methods, collapse = "', '"), "'")
  }
  
  original_method <- differentiation_method
  
  # Auto-switch to combined method for many compounds
  if (n_compounds > auto_combine_threshold && differentiation_method != "combined") {
    if (differentiation_method == "linetype" && n_compounds > 6) {
      differentiation_method <- "combined"
      message("Auto-switching to 'combined' method: ", n_compounds, 
              " compounds exceed the 6 available line types")
    } else if (differentiation_method == "pointshape" && n_compounds > 15) {
      differentiation_method <- "combined"
      message("Auto-switching to 'combined' method: ", n_compounds, 
              " compounds exceed the 15 available point shapes")
    }
  }
  
  # Generate default colors
  if (is.null(colors)) {
    if (differentiation_method %in% c("color", "combined")) {
      colors <- if (n_compounds <= 8) {
        c("blue", "red", "green", "purple", "orange", "brown", "pink", "gray")
      } else if (n_compounds <= 15) {
        c("blue", "red", "green", "purple", "orange", "brown", "pink", "gray",
          "darkblue", "darkred", "darkgreen", "cyan", "magenta", "yellow", "black")
      } else {
        grDevices::rainbow(n_compounds)
      }
    } else {
      colors <- rep("black", n_compounds)
    }
  }
  
  # Ensure sufficient colors
  if (length(colors) < n_compounds) {
    if (differentiation_method %in% c("color", "combined")) {
      warning("Insufficient number of colors, generating additional colors automatically")
      colors <- c(colors, grDevices::rainbow(n_compounds - length(colors)))
    } else {
      colors <- rep(colors, length.out = n_compounds)
    }
  }
  
  # Generate default line types
  if (is.null(line_types)) {
    if (differentiation_method == "linetype" && n_compounds > 6) {
      warning("Only 6 line types available in R. Using repeated line types for ", n_compounds, " compounds")
    }
    line_types <- rep(1:6, length.out = n_compounds)
  } else {
    line_types <- rep(line_types, length.out = n_compounds)
  }
  
  # Generate default point shapes
  if (is.null(point_shapes)) {
    base_point_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14)
    point_shapes <- rep(base_point_shapes, length.out = n_compounds)
  } else {
    point_shapes <- rep(point_shapes, length.out = n_compounds)
  }
  
  # User guidance for many compounds
  if (n_compounds > 6 && original_method == "linetype") {
    message("For ", n_compounds, " compounds with line types:")
    message("  - Only 6 unique line types available in R")
    message("  - Line types will repeat after the 6th compound")
    message("  - Consider using differentiation_method = 'combined' for better distinction")
  }
  
  if (n_compounds > 15 && differentiation_method != "combined") {
    message("Recommendation: For ", n_compounds, " compounds, use differentiation_method = 'combined'")
    message("  This combines colors + line types + point shapes for maximum differentiation")
  }
  
  # Extract compound names
  compound_names_clean <- sapply(compound_indices, function(i) {
    strsplit(results$detailed_results[[i]]$compound, " \\| ")[[1]][1]
  })
  
  # Calculate X-axis limits from all compounds
  all_log_inhibitors <- numeric(0)
  for (idx in compound_indices) {
    result <- results$detailed_results[[idx]]
    if (!is.null(result$model) && isTRUE(result$success)) {
      data <- result$data
      n_rows <- nrow(data) / 2
      control_vals <- data$log_inhibitor[c(1, n_rows)]
      plot_data <- data[!data$log_inhibitor %in% control_vals, ]
      
      if (nrow(plot_data) < 2) plot_data <- data
      
      valid_log_inhibitors <- plot_data$log_inhibitor[is.finite(plot_data$log_inhibitor)]
      all_log_inhibitors <- c(all_log_inhibitors, valid_log_inhibitors)
    }
  }
  
  # Determine x_limits with fallback
  if (length(all_log_inhibitors) == 0) {
    x_limits <- c(-10, -2)
    warning("No valid log_inhibitor values found, using default limits")
  } else {
    x_limits <- range(all_log_inhibitors, na.rm = TRUE)
  }
  
  # Plot saving setup 
  plot_saved <- FALSE
  original_dev <- grDevices::dev.cur()
  filename <- NULL
  
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
    
    file_ext <- tolower(tools::file_ext(filename))
    supported_formats <- c("png", "jpg", "jpeg", "tiff", "pdf", "svg")
    
    if (!file_ext %in% supported_formats) {
      warning("Unsupported format '", file_ext, "'. Using PNG instead.")
      filename <- sub(paste0("\\.", file_ext, "$"), ".png", filename, ignore.case = TRUE)
      file_ext <- "png"
    }
    
    if (grDevices::dev.cur() != 1) {
      grDevices::dev.off()
    }
    
    switch(file_ext,
           png = grDevices::png(filename, width = plot_width, height = plot_height, 
                                units = "in", res = plot_dpi, bg = "white"),
           jpg = grDevices::jpeg(filename, width = plot_width, height = plot_height, 
                                 units = "in", res = plot_dpi, quality = 90, bg = "white"),
           jpeg = grDevices::jpeg(filename, width = plot_width, height = plot_height, 
                                  units = "in", res = plot_dpi, quality = 90, bg = "white"),
           tiff = grDevices::tiff(filename, width = plot_width, height = plot_height, 
                                  units = "in", res = plot_dpi, compression = "lzw", bg = "white"),
           pdf = grDevices::pdf(filename, width = plot_width, height = plot_height, 
                                bg = "white", pointsize = 12),
           svg = grDevices::svg(filename, width = plot_width, height = plot_height, 
                                bg = "white")
    )
    
    plot_saved <- TRUE
    on.exit({
      if (grDevices::dev.cur() != 1) {
        grDevices::dev.off()
      }
      # Restaurar dispositivo original se necessário
      if (original_dev > 1 && grDevices::dev.cur() == 1) {
        grDevices::dev.set(original_dev)
      }
    })
  }
  
  # Adjust layout for external legend
  if (legend_position == "outside") {
    original_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(original_par), add = TRUE)
    
    # Calculate layout based on legend_area_ratio
    plot_ratio <- (1 - legend_area_ratio) / legend_area_ratio
    layout_widths <- c(plot_ratio, 1)
    
    # Create layout with main plot and legend area
    layout_matrix <- matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE)
    graphics::layout(layout_matrix, widths = layout_widths)
    
    # Set margins for main plot (left plot)
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  }
  
  # Create base plot
  graphics::plot(NA, xlim = x_limits, ylim = y_limits,
                 xlab = expression(paste("Log"[10], " Concentration [M]")), 
                 ylab = ifelse(results$normalized, 
                               "Normalized BRET ratio [%]", 
                               "BRET ratio"),
                 main = plot_title,
                 cex.lab = axis_label_cex,
                 cex.axis = axis_number_cex)
  
  if (show_grid) {
    graphics::grid()
  }
  
  # Styling configuration
  styling_config <- list(
    color = list(lty = 1, pch = NA, show_points = FALSE),
    linetype = list(lty = line_types, pch = NA, show_points = FALSE),
    pointshape = list(lty = 1, pch = point_shapes, show_points = TRUE),
    combined = list(lty = line_types, pch = point_shapes, show_points = TRUE)
  )
  
  current_style <- styling_config[[differentiation_method]]
  
  # Error bars function
  add_error_bars <- function(x, y_mean, y_sd, color, bar_width, lwd) {
    if (length(x) == 0 || any(is.na(y_sd))) return()
    
    x_left <- x - bar_width/2
    x_right <- x + bar_width/2
    
    graphics::segments(x_left, y_mean, x_right, y_mean, col = color, lwd = lwd)
    
    for (i in seq_along(x)) {
      if (!is.na(y_sd[i]) && y_sd[i] > 0) {
        y_top <- y_mean[i] + y_sd[i]
        y_bottom <- y_mean[i] - y_sd[i]
        
        graphics::segments(x[i], y_bottom, x[i], y_top, col = color, lwd = lwd)
        graphics::segments(x_left[i], y_top, x_right[i], y_top, col = color, lwd = lwd)
        graphics::segments(x_left[i], y_bottom, x_right[i], y_bottom, col = color, lwd = lwd)
      }
    }
  }
  
  # Plot each compound
  for (i in seq_along(compound_indices)) {
    idx <- compound_indices[i]
    result <- results$detailed_results[[idx]]
    
    if (is.null(result$model) || !isTRUE(result$success)) {
      warning("Compound ", compound_names_clean[i], ": no successful model fit - skipping")
      next
    }
    
    data <- result$data
    n_rows <- nrow(data) / 2
    control_vals <- data$log_inhibitor[c(1, n_rows)]
    plot_data <- data[!data$log_inhibitor %in% control_vals, ]
    
    if (nrow(plot_data) < 2) plot_data <- data
    
    valid_data <- plot_data[is.finite(plot_data$log_inhibitor), ]
    
    if (nrow(valid_data) < 2) {
      warning("Compound ", compound_names_clean[i], ": insufficient data to generate curve")
      next
    }
    
    # Generate and plot fitted curve
    x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)
    pred_df <- data.frame(log_inhibitor = x_seq)
    pred_df$response <- predict(result$model, newdata = pred_df)
    
    graphics::lines(pred_df$log_inhibitor, pred_df$response, 
                    col = colors[i], 
                    lty = ifelse(length(current_style$lty) == 1, 
                                 current_style$lty, current_style$lty[i]),
                    lwd = 2)
    
    # Add points if using point-based differentiation
    if (current_style$show_points) {
      conc_levels <- unique(valid_data$log_inhibitor)
      
      mean_responses <- vapply(conc_levels, function(conc) {
        mean(valid_data$response[valid_data$log_inhibitor == conc], na.rm = TRUE)
      }, numeric(1))
      
      sd_responses <- vapply(conc_levels, function(conc) {
        sd(valid_data$response[valid_data$log_inhibitor == conc], na.rm = TRUE)
      }, numeric(1))
      
      graphics::points(conc_levels, mean_responses, 
                       pch = current_style$pch[i], 
                       col = colors[i], 
                       cex = 1.2)
      
      if (show_error_bars) {
        add_error_bars(conc_levels, mean_responses, sd_responses, 
                       colors[i], error_bar_width, error_bar_lwd)
      }
    }
  }
  
  # Add legend
  if (show_legend) {
    legend_params <- list(
      legend = compound_names_clean,
      bty = "n",
      col = colors,
      cex = legend_cex
    )
    
    if (legend_position == "outside") {
      # Switch to legend area (right panel) with adequate margins
      graphics::par(mar = c(5.1, 0.001, 4.1, 0.001))
      graphics::plot.new()
      
      legend_params$x <- "center"
    } else {
      legend_params$x <- legend_position
    }
    
    # Apply styling based on differentiation method
    if (differentiation_method == "color") {
      legend_params$lwd <- 2
      legend_params$lty <- 1
      legend_params$pch <- NA
    } else if (differentiation_method == "linetype") {
      legend_params$lwd <- 2
      legend_params$lty <- line_types
      legend_params$pch <- NA
    } else if (differentiation_method == "pointshape") {
      legend_params$pch <- point_shapes
      legend_params$pt.cex <- legend_point_cex
      legend_params$pt.lwd <- 1.5
      legend_params$lwd <- 0
      legend_params$lty <- 0
      legend_params$seg.len <- 0
    } else { # combined
      legend_params$lwd <- 2
      legend_params$lty <- line_types
      legend_params$pch <- point_shapes
      legend_params$pt.cex <- legend_point_cex
      legend_params$pt.lwd <- 1.5
    }
    
    do.call(graphics::legend, legend_params)
  }
  
  if (plot_saved) {
    message("Plot successfully saved as: ", normalizePath(filename))
  }
  
  # Return metadata
  invisible(list(
    compound_names = compound_names_clean,
    compound_indices = compound_indices,
    n_compounds = n_compounds,
    plot_limits = list(x_limits = x_limits, y_limits = y_limits),
    differentiation_method = differentiation_method,
    styling = list(colors = colors, line_types = line_types, point_shapes = point_shapes),
    error_bars = list(
      enabled = show_error_bars,
      width = error_bar_width,
      line_width = error_bar_lwd
    ),
    legend_position = legend_position,
    legend_settings = list(
      cex = legend_cex,
      area_ratio = legend_area_ratio,
      point_cex = legend_point_cex
    ),
    plot_title = plot_title,
    file_saved = if (plot_saved) filename else NULL,
    file_format = if (plot_saved) tools::file_ext(filename) else NULL,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    timestamp = Sys.time()
  ))
}
