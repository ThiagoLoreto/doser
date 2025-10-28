#' Plot Dose-Response Curves with Professional Styling
#'
#' Generates publication-quality dose-response plots from analysis results.
#' Supports error bars, fitted curves, IC50 lines, and multiple export formats.
#' @param results List object returned by \code{\link{fit_dose_response}} containing
#'   dose-response analysis results.
#' @param compound_index Numeric index specifying which compound to plot (default: 1).
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 150)).
#' @param point_color Color for data points (default: "black").
#' @param line_color Color for fitted curve (default: "red").
#' @param ic50_line_color Color for IC50 vertical line (default: "gray").
#' @param point_size Size multiplier for data points (default: 1).
#' @param line_width Line width for fitted curve (default: 2).
#' @param error_bar_width Width of error bar ends (default: 0.01).
#' @param show_ic50_line Logical indicating whether to show vertical IC50 line (default: TRUE).
#' @param show_legend Logical indicating whether to show parameter legend (default: TRUE).
#' @param show_grid Logical indicating whether to show background grid (default: FALSE).
#' @param save_plot Either a file path for saving the plot, or TRUE for automatic naming
#'   (default: NULL, no saving).
#' @param plot_width Plot width in inches for saved plots (default: 10).
#' @param plot_height Plot height in inches for saved plots (default: 8).
#' @param plot_dpi Resolution for saved raster images (default: 600).
#' @param axis_label_cex Character expansion factor for axis labels (default: 1.4).
#' @param axis_number_cex Character expansion factor for axis numbers (default: 1.4).
#' @param x_axis_title Custom x-axis title. If NULL, uses default expression .
#' @param y_axis_title Custom y-axis title. If NULL, uses default based on 
#'   normalization status
#'
#' @return Invisibly returns a list containing plot metadata:
#' \itemize{
#'   \item \code{compound_name}: Name of the plotted compound
#'   \item \code{compound_index}: Index of the plotted compound
#'   \item \code{model_success}: Whether model fitting was successful
#'   \item \code{summary_data}: Data frame with summarized plotting data
#'   \item \code{plot_config}: Configuration settings used for plotting
#'   \item \code{y_limits_used}: Y-axis limits actually used
#'   \item \code{data_points}: Number of data points plotted
#'   \item \code{concentration_levels}: Number of concentration levels
#'   \item \code{file_saved}: Path to saved file if plot was saved
#'   \item \code{plot_dimensions}: Dimensions of the plot (width, height, dpi)
#'   \item \code{timestamp}: Time when plot was generated
#' }
#'
#' @details
#' This function creates professional-quality dose-response plots suitable for
#' publications and presentations. Key features include:
#'
#' \strong{Plot Elements:}
#' \itemize{
#'   \item Data points with error bars (standard deviation)
#'   \item Fitted dose-response curve (when model converged)
#'   \item Vertical IC50 line with dashed style
#'   \item Parameter legend with IC50 and R2 values
#'   \item Professional axis labels and formatting
#'   \item Optional background grid
#' }
#'
#' \strong{Supported Export Formats:}
#' \itemize{
#'   \item PNG (high-resolution, recommended for publications)
#'   \item JPEG (good for presentations)
#'   \item TIFF (lossless compression)
#'   \item PDF (vector format, scalable)
#'   \item SVG (vector format, editable)
#' }
#'
#' @examples
#' \dontrun{
#' # Perform dose-response analysis first
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#'
#' # Basic plot for first compound
#' plot_dose_response(analysis_results)
#'
#' # Customized plot with specific styling
#' plot_dose_response(
#'   results = analysis_results,
#'   compound_index = 2,
#'   point_color = "blue",
#'   line_color = "darkred",
#'   show_grid = TRUE,
#'   y_limits = c(0, 200)
#' )
#'
#' # Save plot automatically with compound name
#' plot_dose_response(
#'   results = analysis_results,
#'   compound_index = 1,
#'   save_plot = TRUE  # Saves as "dose_response_CompoundName.png"
#' )
#'
#' # Save plot to specific file with custom dimensions
#' plot_dose_response(
#'   results = analysis_results,
#'   save_plot = "my_plot.pdf",
#'   plot_width = 8,
#'   plot_height = 6
#' )
#'
#' # Access plot metadata
#' plot_info <- plot_dose_response(analysis_results, compound_index = 3)
#' print(plot_info$compound_name)
#' print(plot_info$data_points)
#' }
#'
#' @section Plot Customization:
#' The function provides extensive customization options:
#' \itemize{
#'   \item \strong{Colors}: Customize points, lines, and IC50 marker
#'   \item \strong{Sizes}: Adjust point sizes, line widths, and error bars
#'   \item \strong{Axes}: Control limits, labels, and text sizes
#'   \item \strong{Elements}: Toggle legend, grid, and IC50 line
#'   \item \strong{Export}: Multiple formats with quality control
#' }
#'
#' @section Automatic Features:
#' \itemize{
#'   \item Automatic directory creation for saved plots
#'   \item Smart error bar handling (only shown when meaningful)
#'   \item Graceful handling of failed model fits
#'   \item Professional axis formatting and labeling
#'   \item Compound name extraction and display
#' }
#'
#' @seealso
#' \code{\link{fit_dose_response}} for generating analysis results
#' \code{\link[graphics]{plot}} for base plotting functions
#' \code{\link[grDevices]{png}} for plot export options
#'
#' @export
#'
#' @references
#' For visualization best practices:
#' \itemize{
#'   \item Nature Scientific Figures Guidelines
#'   \item R Graphics Cookbook (O'Reilly)
#' }





plot_dose_response <- function(results, compound_index = 1, y_limits = c(0, 150), 
                               point_color = "black", line_color = "red", 
                               ic50_line_color = "gray", point_size = 1, 
                               line_width = 2, error_bar_width = 0.01,
                               show_ic50_line = TRUE, show_legend = TRUE,
                               show_grid = FALSE, save_plot = NULL, 
                               plot_width = 10, plot_height = 8, plot_dpi = 600,
                               axis_label_cex = 1.4, axis_number_cex = 1.4,
                               x_axis_title = NULL, y_axis_title = NULL,
                               enforce_bottom_threshold = NULL, bottom_threshold = 60) {
  
  # Input validation
  validate_inputs <- function(results, compound_index) {
    if (missing(results)) {
      stop("Argument 'results' is required")
    }
    
    if (!is.list(results) || !"detailed_results" %in% names(results)) {
      stop("Invalid 'results' object. Must be output from fit function")
    }
    
    if (length(results$detailed_results) == 0) {
      stop("No compounds found in results object")
    }
    
    if (compound_index < 1 || compound_index > length(results$detailed_results)) {
      stop("Compound index ", compound_index, " out of range. Must be between 1 and ", 
           length(results$detailed_results))
    }
  }
  
  validate_inputs(results, compound_index)
  
  # Extract compound data
  result <- results$detailed_results[[compound_index]]
  
  if (!is.list(result) || !"data" %in% names(result)) {
    stop("Invalid result structure for compound ", compound_index)
  }
  
  # Use threshold settings from results if not explicitly provided
  if (is.null(enforce_bottom_threshold)) {
    enforce_bottom_threshold <- if (!is.null(results$threshold_settings)) {
      results$threshold_settings$enforce_bottom_threshold
    } else {
      FALSE
    }
  }
  
  # Clean and validate data
  clean_data <- stats::na.omit(result$data)
  required_cols <- c("log_inhibitor", "response")
  if (nrow(clean_data) == 0 || !all(required_cols %in% names(clean_data))) {
    stop("No valid data points available or missing required columns")
  }
  
  # Calculate summary statistics (mean ± SD per concentration)
  calculate_summary_stats <- function(data) {
    summary_data <- do.call(rbind, lapply(split(data, data$log_inhibitor), function(sub_df) {
      data.frame(
        log_inhibitor = unique(sub_df$log_inhibitor),
        mean_response = mean(sub_df$response, na.rm = TRUE),
        sd_response = sd(sub_df$response, na.rm = TRUE),
        n_replicates = nrow(sub_df)
      )
    }))
    
    stats::na.omit(summary_data)
  }
  
  summary_data <- calculate_summary_stats(clean_data)
  
  if (nrow(summary_data) == 0) {
    stop("No valid summary data available for plotting")
  }
  
  # Extract compound name (remove plate info if present)
  compound_name_display <- strsplit(result$compound, " \\| ")[[1]][1]
  
  # Setup plot configuration with flexible axis titles
  setup_plot_config <- function() {
    x_lab <- if (!is.null(x_axis_title)) {
      x_axis_title
    } else {
      expression(paste("Log"[10], " Concentration [M]"))
    }
    
    y_lab <- if (!is.null(y_axis_title)) {
      y_axis_title
    } else {
      ifelse(results$normalized, "Normalized BRET ratio [%]", "BRET ratio")
    }
    
    list(
      x_lab = x_lab,
      y_lab = y_lab,
      point_color = point_color,
      line_color = line_color,
      point_size = point_size,
      line_width = line_width,
      error_bar_width = error_bar_width,
      axis_label_cex = axis_label_cex,
      axis_number_cex = axis_number_cex
    )
  }
  
  plot_config <- setup_plot_config()
  
  # File saving setup with multiple format support
  setup_plot_device <- function() {
    if (is.null(save_plot)) return(NULL)
    
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      safe_name <- gsub("[^a-zA-Z0-9._-]", "_", compound_name_display)
      filename <- paste0("dose_response_", safe_name, ".png")
    } else {
      stop("save_plot must be either a file path or TRUE for auto-naming")
    }
    
    # Create directory if needed
    plot_dir <- dirname(filename)
    if (plot_dir != "." && !dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Determine file format with fallback to PNG
    file_ext <- tolower(tools::file_ext(filename))
    supported_formats <- c("png", "jpg", "jpeg", "tiff", "pdf", "svg")
    
    if (!file_ext %in% supported_formats) {
      warning("Unsupported format '", file_ext, "'. Using PNG instead.")
      filename <- sub(paste0("\\.", file_ext, "$"), ".png", filename, ignore.case = TRUE)
      file_ext <- "png"
    }
    
    # Open appropriate graphics device
    device_func <- switch(file_ext,
                          png = function() grDevices::png(filename, width = plot_width, height = plot_height, 
                                                          units = "in", res = plot_dpi, bg = "white"),
                          jpg =, jpeg = function() grDevices::jpeg(filename, width = plot_width, height = plot_height, 
                                                                   units = "in", res = plot_dpi, quality = 90, bg = "white"),
                          tiff = function() grDevices::tiff(filename, width = plot_width, height = plot_height, 
                                                            units = "in", res = plot_dpi, compression = "lzw", bg = "white"),
                          pdf = function() grDevices::pdf(filename, width = plot_width, height = plot_height, 
                                                          bg = "white", pointsize = 12),
                          svg = function() grDevices::svg(filename, width = plot_width, height = plot_height, 
                                                          bg = "white")
    )
    
    device_func()
    return(list(device = grDevices::dev.cur(), filename = filename))
  }
  
  # Setup plot device and auto-close on exit
  original_dev <- grDevices::dev.cur()
  plot_device_info <- setup_plot_device()
  
  if (!is.null(plot_device_info)) {
    on.exit({
      if (grDevices::dev.cur() == plot_device_info$device) {
        grDevices::dev.off()
        if (grDevices::dev.cur() == 1 && original_dev > 1) {
          grDevices::dev.set(original_dev)
        }
        cat("Plot saved as:", normalizePath(plot_device_info$filename), "\n")
      }
    })
  }
  
  # Core plotting functions
  create_base_plot <- function(title_suffix = "") {
    main_title <- paste(compound_name_display, title_suffix)
    
    # Professional plot settings
    old_par <- graphics::par(
      mar = c(4.5, 5, 3.5, 1.5),  # margins: bottom, left, top, right
      mgp = c(2.8, 0.8, 0),       # axis title, labels, line
      las = 1                      # horizontal labels
    )
    on.exit(graphics::par(old_par))
    
    # Create base scatter plot
    graphics::plot(summary_data$log_inhibitor, summary_data$mean_response, 
                   xlab = plot_config$x_lab, 
                   ylab = plot_config$y_lab, 
                   main = trimws(main_title),
                   pch = 21,
                   bg = plot_config$point_color,
                   col = "black",
                   cex = plot_config$point_size,
                   ylim = y_limits, 
                   xaxt = "n",
                   cex.lab = plot_config$axis_label_cex,
                   cex.axis = plot_config$axis_number_cex,
                   cex.main = plot_config$axis_label_cex * 1.1,
                   font.main = 2,
                   panel.first = if (show_grid) {
                     graphics::grid(col = "grey90", lty = "solid", lwd = 0.7)
                   },
                   bty = "l",
                   tcl = -0.3)
    
    # Custom x-axis with proper formatting
    x_ticks <- pretty(summary_data$log_inhibitor)
    graphics::axis(1, at = x_ticks, 
                   labels = format(x_ticks, digits = 2, nsmall = 1),
                   cex.axis = plot_config$axis_number_cex)
  }
  
  # Add error bars (SD) for points with valid standard deviation
  add_error_bars <- function() {
    valid_mask <- !is.na(summary_data$sd_response) & 
      is.finite(summary_data$sd_response) & 
      summary_data$n_replicates > 1 &
      summary_data$sd_response > 1e-10
    
    if (any(valid_mask)) {
      valid_data <- summary_data[valid_mask, ]
      graphics::arrows(
        x0 = valid_data$log_inhibitor,
        y0 = valid_data$mean_response - valid_data$sd_response,
        x1 = valid_data$log_inhibitor,
        y1 = valid_data$mean_response + valid_data$sd_response,
        angle = 90, 
        code = 3, 
        length = plot_config$error_bar_width, 
        col = plot_config$point_color,
        lwd = 1.2
      )
    }
  }
  
  # Generate smooth fitted curve for plotting - NO EXTRAPOLATION
  generate_fitted_curve <- function(model) {
    x_range <- range(summary_data$log_inhibitor, na.rm = TRUE)
    if (!all(is.finite(x_range))) return(NULL)
    
    # Create sequence ONLY within observed data range
    x_seq <- seq(x_range[1], x_range[2], length.out = 300)
    
    predictions <- tryCatch({
      predict(model, newdata = data.frame(log_inhibitor = x_seq))
    }, error = function(e) NULL)
    
    if (!is.null(predictions)) {
      data.frame(log_inhibitor = x_seq, response = predictions)
    }
  }
  
  add_fitted_curve <- function(model) {
    curve_data <- generate_fitted_curve(model)
    if (!is.null(curve_data)) {
      graphics::lines(curve_data$log_inhibitor, curve_data$response, 
                      col = plot_config$line_color, 
                      lwd = plot_config$line_width,
                      lty = "solid")
    }
  }
  
  # Add vertical line at IC50 position
  add_ic50_line <- function(model) {
    if (!show_ic50_line) return(NULL)
    
    log_ic50 <- tryCatch({
      coefs <- stats::coef(model)
      if ("LogIC50" %in% names(coefs)) coefs["LogIC50"] else NA
    }, error = function(e) NA)
    
    if (is.finite(log_ic50)) {
      graphics::abline(v = log_ic50, lty = 2, col = ic50_line_color, lwd = 1.5)
    }
    return(log_ic50)
  }
  
  # Create legend content with model parameters
  create_legend_content <- function(model = NULL) {
    if (!show_legend) return(NULL)
    
    if (!is.null(model) && isTRUE(result$success)) {
      # Extract model parameters
      log_ic50 <- tryCatch(stats::coef(model)["LogIC50"], error = function(e) NA)
      ic50_value <- if (is.finite(log_ic50)) 10^log_ic50 else NA
      r_squared <- round(result$goodness_of_fit$R_squared, 3)
      
      # Check if IC50 was excluded due to threshold
      ic50_excluded <- FALSE
      if (enforce_bottom_threshold && !is.na(result$parameters$Value[1])) {
        bottom_value <- result$parameters$Value[1]  # Bottom parameter
        ic50_excluded <- bottom_value >= bottom_threshold
      }
      
      legend_text <- c()
      
      if (ic50_excluded) {
        legend_text <- c(legend_text, "LogIC50 = NA")
      } else if (is.finite(log_ic50)) {
        legend_text <- c(legend_text, paste("LogIC50 =", round(log_ic50, 3)))
      } else {
        legend_text <- c(legend_text, "LogIC50 = NA")
      }
      
      if (ic50_excluded) {
        legend_text <- c(legend_text, "IC50 = NA")
      } else if (is.finite(ic50_value)) {
        legend_text <- c(legend_text, paste("IC50 =", sprintf("%.2e", ic50_value)))
      } else {
        legend_text <- c(legend_text, "IC50 = NA")
      }
      
      # Always show R²
      legend_text <- c(legend_text, paste("R² =", r_squared))
      
      return(legend_text)
      
    } else {
      return("Model did not converge")
    }
  }
  
  add_legend <- function(legend_content) {
    if (!is.null(legend_content)) {
      graphics::legend("bottomright", 
                       legend = legend_content,
                       bty = "n",
                       cex = 0.8 * plot_config$axis_label_cex)
    }
  }
  
  # Main plotting logic
  model_success <- !is.null(result$model) && isTRUE(result$success)
  
  if (model_success) {
    create_base_plot()
    add_error_bars()
    add_fitted_curve(result$model)
    add_ic50_line(result$model)
    legend_content <- create_legend_content(result$model)
  } else {
    create_base_plot("(Model failed)")
    add_error_bars()
    legend_content <- create_legend_content()
  }
  
  add_legend(legend_content)
  
  # Return comprehensive metadata
  invisible(list(
    compound_name = compound_name_display,
    compound_index = compound_index,
    model_success = model_success,
    summary_data = summary_data,
    plot_config = plot_config,
    y_limits_used = y_limits,
    data_points = nrow(clean_data),
    concentration_levels = nrow(summary_data),
    file_saved = if (!is.null(plot_device_info)) plot_device_info$filename else NULL,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    timestamp = Sys.time()
  ))
}