#' Plot Dose-Response Curves with Professional Styling
#'
#' Generates publication-quality dose-response plots from analysis results.
#' Supports error bars, fitted curves, IC50 lines, and multiple export formats.
#'
#' @param results List object returned by \code{\link{fit_dose_response}} containing
#'   dose-response analysis results.
#' @param compound_index Numeric index specifying which compound to plot (default: 1).
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 150)).
#' @param point_color Color for data points (default: "black").
#' @param line_color Color for fitted curve (default: "black").
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
                               point_color = "black", line_color = "black", 
                               ic50_line_color = "gray", point_size = 1, 
                               line_width = 2, error_bar_width = 0.01,
                               show_ic50_line = TRUE, show_legend = TRUE,
                               show_grid = FALSE, save_plot = NULL, 
                               plot_width = 10, plot_height = 8, plot_dpi = 600,
                               axis_label_cex = 1.4, axis_number_cex = 1.4) {
  
  # Enhanced input validation
  if (missing(results)) {
    stop("Argument 'results' is required")
  }
  
  if (!is.list(results) || !"detailed_results" %in% names(results)) {
    stop("Invalid 'results' object. Must be output from analyze_dose_response_complete()")
  }
  
  if (length(results$detailed_results) == 0) {
    stop("No compounds found in results object")
  }
  
  if (compound_index < 1 || compound_index > length(results$detailed_results)) {
    stop("Compound index ", compound_index, " out of range. Must be between 1 and ", 
         length(results$detailed_results))
  }
  
  # Extract and validate data
  result <- results$detailed_results[[compound_index]]
  if (!is.list(result) || !"data" %in% names(result)) {
    stop("Invalid result structure for compound ", compound_index)
  }
  
  clean_data <- stats::na.omit(result$data)
  if (nrow(clean_data) == 0) {
    stop("No valid data points available for plotting after removing NAs")
  }
  
  # Enhanced data validation
  if (!all(c("log_inhibitor", "response") %in% names(clean_data))) {
    stop("Data must contain 'log_inhibitor' and 'response' columns")
  }
  
  # Calculate summary statistics using efficient methods
  summary_data <- do.call(rbind, lapply(split(clean_data, clean_data$log_inhibitor), function(sub_df) {
    data.frame(
      log_inhibitor = unique(sub_df$log_inhibitor),
      mean_response = mean(sub_df$response, na.rm = TRUE),
      sd_response = sd(sub_df$response, na.rm = TRUE),
      n_replicates = nrow(sub_df)
    )
  }))
  
  summary_data <- stats::na.omit(summary_data)
  
  if (nrow(summary_data) == 0) {
    stop("No valid summary data available for plotting")
  }
  
  # Process compound name
  compound_name_display <- strsplit(result$compound, " \\| ")[[1]][1]
  
  # Enhanced plot saving with multiple format support
  plot_device <- NULL
  original_dev <- grDevices::dev.cur()
  
  if (!is.null(save_plot)) {
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      safe_name <- gsub("[^a-zA-Z0-9._-]", "_", compound_name_display)
      filename <- paste0("dose_response_", safe_name, ".png")
    } else {
      stop("save_plot must be either a file path or TRUE for auto-naming")
    }
    
    # Create directory if it doesn't exist
    plot_dir <- dirname(filename)
    if (plot_dir != "." && !dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Determine format from extension with fallback
    file_ext <- tolower(tools::file_ext(filename))
    supported_formats <- c("png", "jpg", "jpeg", "tiff", "pdf", "svg")
    
    if (!file_ext %in% supported_formats) {
      warning("Unsupported format '", file_ext, "'. Using PNG instead.")
      filename <- sub(paste0("\\.", file_ext, "$"), ".png", filename, ignore.case = TRUE)
      file_ext <- "png"
    }
    
    # Open appropriate device with professional settings
    switch(file_ext,
           png = grDevices::png(filename, width = plot_width, height = plot_height, 
                                units = "in", res = plot_dpi, bg = "white"),
           jpg =,
           jpeg = grDevices::jpeg(filename, width = plot_width, height = plot_height, 
                                  units = "in", res = plot_dpi, quality = 90, bg = "white"),
           tiff = grDevices::tiff(filename, width = plot_width, height = plot_height, 
                                  units = "in", res = plot_dpi, compression = "lzw", bg = "white"),
           pdf = grDevices::pdf(filename, width = plot_width, height = plot_height, 
                                bg = "white", pointsize = 12),
           svg = grDevices::svg(filename, width = plot_width, height = plot_height, 
                                bg = "white")
    )
    
    plot_device <- grDevices::dev.cur()
    on.exit({
      if (grDevices::dev.cur() == plot_device) {
        grDevices::dev.off()
        # Restore original device if it was different
        if (grDevices::dev.cur() == 1 && original_dev > 1) {
          grDevices::dev.set(original_dev)
        }
        cat("Plot saved as:", normalizePath(filename), "\n")
      }
    })
  }
  
  plot_config <- list(
    x_lab = expression(paste("Log"[10], " Concentration [M]")),
    y_lab = ifelse(results$normalized, 
                   "Normalized BRET ratio [%]", 
                   "BRET ratio"),
    point_color = point_color,
    line_color = line_color,
    point_size = point_size,
    line_width = line_width,
    error_bar_width = error_bar_width,
    axis_label_cex = axis_label_cex,
    axis_number_cex = axis_number_cex
  )
  
  # Enhanced core plotting functions
  create_base_plot <- function(title_suffix = "") {
    main_title <- paste(compound_name_display, title_suffix)
    
    # Professional margins and settings
    old_par <- graphics::par(
      mar = c(4.5, 5, 3.5, 1.5),  # Bottom, Left, Top, Right margins
      mgp = c(2.8, 0.8, 0),       # Axis title, labels, line
      las = 1                      # Horizontal labels
    )
    on.exit(graphics::par(old_par))
    
    # Create base plot with professional styling
    graphics::plot(summary_data$log_inhibitor, summary_data$mean_response, 
                   xlab = plot_config$x_lab, 
                   ylab = plot_config$y_lab, 
                   main = trimws(main_title),
                   pch = 21,                       # Filled circles with border
                   bg = plot_config$point_color,    # Fill color
                   col = "black",                  # Border color
                   cex = plot_config$point_size,
                   ylim = y_limits, 
                   xaxt = "n",
                   cex.lab = plot_config$axis_label_cex,
                   cex.axis = plot_config$axis_number_cex,
                   cex.main = plot_config$axis_label_cex * 1.1,
                   font.main = 2,                  # Bold title
                   panel.first = {
                     if (show_grid) {
                       # Professional grid
                       graphics::grid(col = "grey90", lty = "solid", lwd = 0.7)
                     }
                   },
                   bty = "l",                      # L-shaped box
                   tcl = -0.3)                     # Tick length
    
    # Professional x-axis formatting
    x_ticks <- pretty(summary_data$log_inhibitor)
    graphics::axis(1, at = x_ticks, 
                   labels = format(x_ticks, digits = 2, nsmall = 1),
                   cex.axis = plot_config$axis_number_cex)
  }
  
  # CORRECAO: Funcao add_error_bars atualizada para evitar o erro
  add_error_bars <- function() {
    # Filtra pontos com SD valido, positivo e significativo
    valid_mask <- !is.na(summary_data$sd_response) & 
      is.finite(summary_data$sd_response) & 
      summary_data$n_replicates > 1 &
      summary_data$sd_response > 1e-10  # ? CORRECAO: exclui SD zero ou muito pequeno
    
    if (any(valid_mask)) {
      valid_data <- summary_data[valid_mask, ]
      
      # Verificacao adicional para garantir que as setas tenham comprimento > 0
      arrow_lengths <- 2 * valid_data$sd_response
      non_zero_arrows <- arrow_lengths > 1e-10
      
      if (any(non_zero_arrows)) {
        graphics::arrows(
          x0 = valid_data$log_inhibitor[non_zero_arrows],
          y0 = valid_data$mean_response[non_zero_arrows] - valid_data$sd_response[non_zero_arrows],
          x1 = valid_data$log_inhibitor[non_zero_arrows],
          y1 = valid_data$mean_response[non_zero_arrows] + valid_data$sd_response[non_zero_arrows],
          angle = 90, 
          code = 3, 
          length = plot_config$error_bar_width, 
          col = plot_config$point_color,
          lwd = 1.2
        )
      }
    }
  }
  
  # Enhanced fitted curve generation
  generate_fitted_curve <- function(model) {
    x_range <- range(summary_data$log_inhibitor, na.rm = TRUE)
    if (!all(is.finite(x_range))) return(NULL)
    
    # Extend range slightly for smoother curve edges
    x_padding <- diff(x_range) * 0.08
    x_seq <- seq(x_range[1] - x_padding, x_range[2] + x_padding, length.out = 300)
    
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
  
  create_legend_content <- function(model = NULL) {
    if (!show_legend) return(NULL)
    
    if (!is.null(model) && isTRUE(result$success)) {
      log_ic50 <- tryCatch(stats::coef(model)["LogIC50"], error = function(e) NA)
      ic50_value <- if (is.finite(log_ic50)) 10^log_ic50 else NA
      r_squared <- round(result$goodness_of_fit$R_squared, 3)
      
      c(
        if (is.finite(log_ic50)) paste("LogIC50 =", round(log_ic50, 3)) else "LogIC50 = NA",
        if (is.finite(ic50_value)) paste("IC50 =", sprintf("%.2e", ic50_value)) else "IC50 = NA",
        paste("R2 =", r_squared)
      )
    } else {
      "Model did not converge"
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
  
  # Main plotting execution
  model_success <- !is.null(result$model) && isTRUE(result$success)
  
  if (model_success) {
    # Successful model fit
    create_base_plot()
    add_error_bars()
    add_fitted_curve(result$model)
    add_ic50_line(result$model)
    legend_content <- create_legend_content(result$model)
  } else {
    # Model failed
    create_base_plot("(Model failed)")
    add_error_bars()
    legend_content <- create_legend_content()
  }
  
  add_legend(legend_content)
  
  # Enhanced return of rich metadata
  invisible(list(
    compound_name = compound_name_display,
    compound_index = compound_index,
    model_success = model_success,
    summary_data = summary_data,
    plot_config = plot_config,
    y_limits_used = y_limits,
    data_points = nrow(clean_data),
    concentration_levels = nrow(summary_data),
    file_saved = if (!is.null(plot_device)) filename else NULL,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    timestamp = Sys.time()
  ))
}
