#' Save all dose-response curves from batch analysis
#'
#' @description
#' Generates and saves dose-response curves for all compounds in a batch analysis
#' result object. Plots are saved in the specified format and can be organized
#' by plate, by compound, or in a flat structure.
#'
#' @param batch_drc_results Output from \code{\link{batch_drc_analysis}}. Can be either
#'   the complete wrapper object or the extracted \code{drc_results} component.
#' @param output_dir Character string specifying the main output directory where
#'   plots will be saved. Default is "DRC_Plots".
#' @param organize_by Character string specifying how to organize the output files.
#'   Options are:
#'   \itemize{
#'     \item{"plate" (default):}{ Create subfolders for each plate.}
#'     \item{"compound":}{ Create subfolders for each compound.}
#'     \item{"flat":}{ Save all plots in the same directory.}
#'   }
#' @param compounds_to_plot Optional character vector of compound names to plot.
#'   If NULL (default), all compounds are plotted.
#' @param plates_to_plot Optional character vector of plate names to plot.
#'   If NULL (default), all plates are processed.
#' @param format Character string specifying the output format for plots.
#'   Options include "png", "pdf", "svg", "jpeg". Default is "png".
#' @param width Numeric value specifying the plot width in inches. Default is 10.
#' @param height Numeric value specifying the plot height in inches. Default is 10.
#' @param dpi Numeric value specifying the resolution for raster formats (png, jpeg).
#'   Default is 600.
#' @param point_color Character string specifying the color of data points.
#'   Default is "black".
#' @param point_size Numeric value specifying the size of data points.
#'   Default is 2.
#' @param show_ic50_line Logical value indicating whether to show a vertical line
#'   at the IC50 value. Default is FALSE.
#' @param plot_title Controls the plot title. Can be:
#'   \itemize{
#'     \item{\code{FALSE} (default):}{ No title displayed.}
#'     \item{\code{TRUE}:}{ Automatic title with compound name.
#'       If the model failed to converge, adds "(Model failed)" to the title.}
#'     \item{character:}{ Custom title text applied to all plots.}
#'   }
#' @param verbose Logical value indicating whether to display progress messages.
#'   Default is TRUE.
#' @param ... Additional arguments passed to \code{\link{plot_dose_response}}.
#'   See that function's documentation for all available parameters
#'   (e.g., \code{line_color}, \code{show_legend}, \code{y_limits},
#'   \code{axis_label_size}, etc.).
#'
#' @return An invisible list containing:
#'   \itemize{
#'     \item{\code{total}:}{ Total number of compounds processed.}
#'     \item{\code{successes}:}{ Number of successfully generated plots.}
#'     \item{\code{failures}:}{ Number of failed plots.}
#'     \item{\code{failed_compounds}:}{ Character vector of compounds that failed.}
#'     \item{\code{error_messages}:}{ List of error messages for failed plots.}
#'     \item{\code{output_dir}:}{ Path to the output directory.}
#'     \item{\code{organization}:}{ Organization mode used.}
#'     \item{\code{point_color}:}{ Color used for data points.}
#'     \item{\code{timestamp}:}{ Timestamp of when the function was run.}
#'   }
#'
#' @details
#' This function scans all plates in the batch analysis result, identifies
#' compounds with successful dose-response fits, and generates individual plots
#' for each compound using \code{\link{plot_dose_response}}. The plots are saved
#' with filenames based on compound names and organized according to the
#' \code{organize_by} parameter.
#'
#' The function automatically handles different naming formats:
#' \itemize{
#'   \item{"Construct | Compound" format}
#'   \item{"Construct:Compound" format}
#'   \item{Simple compound names}
#' }
#'
#' File naming convention:
#' \itemize{
#'   \item{When \code{organize_by = "plate"}:}{ \code{compound_name.format} saved in plate subfolder}
#'   \item{When \code{organize_by = "compound"}:}{ \code{plate_name_compound_name.format} saved in compound subfolder}
#'   \item{When \code{organize_by = "flat"}:}{ \code{plate_name_compound_name.format} in main directory}
#' }
#'
#' @note
#' This function requires the \pkg{ggplot2} package. If a plot fails to generate,
#' the error is captured and the function continues processing remaining compounds.
#'
#' @seealso
#' \code{\link{plot_dose_response}} for the underlying plotting function and
#' all available aesthetic parameters;
#' \code{\link{batch_drc_analysis}} for generating the input data.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#' # Basic usage - saves all plots organized by plate
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_results,
#'   output_dir = "All_Plots"
#' )
#'
#' # Save only specific compounds, organized by compound name
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_results,
#'   compounds_to_plot = c("MDCV001", "GZD824", "MLI-2"),
#'   organize_by = "compound",
#'   output_dir = "Selected_Compounds"
#' )
#'
#' # Customize plot appearance
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_results,
#'   output_dir = "Custom_Plots",
#'   point_color = "blue",
#'   point_size = 3,
#'   show_ic50_line = TRUE,
#'   plot_title = TRUE,  # Show automatic titles
#'   y_limits = c(0, 100),  # Passed to plot_dose_response via ...
#'   axis_label_size = 16,
#'   axis_text_size = 14
#' )
#'
#' # Save specific plates with custom title
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_results,
#'   plates_to_plot = c("plate_01", "plate_03"),
#'   output_dir = "Selected_Plates",
#'   plot_title = "Dose-Response Curve"
#' )
#'
#' # Process only compounds from a specific construct
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_results,
#'   compounds_to_plot = c("MDCV001", "MDCV002", "MDCV003"),
#'   output_dir = "LRRK2_Compounds"
#' )
#' }
#'
#' @export

batch_save_all_drc_plots <- function(batch_drc_results,
                                     output_dir = "DRC_Plots",
                                     organize_by = "plate",
                                     compounds_to_plot = NULL,
                                     plates_to_plot = NULL,
                                     format = "png",
                                     width = 10,
                                     height = 10,
                                     dpi = 600,
                                     point_color = "black",
                                     verbose = TRUE,
                                     show_ic50_line = FALSE,
                                     plot_title = FALSE,
                                     point_size = 2,
                                     ...) {

  # ============================================================================
  # 1. VALIDATION AND SETUP
  # ============================================================================

  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  # Helper function for safe filename generation
  safe_filename <- function(string) {
    if (is.null(string) || is.na(string)) return("unknown")
    # Remove any characters that could cause filesystem issues
    s <- gsub("[^[:alnum:]._-]", "_", string)
    s <- gsub("_+", "_", s)
    s <- gsub("^_|_$", "", s)
    if (nchar(s) == 0) return("unknown")
    return(s)
  }

  # Helper function to extract compound name properly
  extract_compound_name <- function(compound_string) {
    if (is.null(compound_string)) return("Unknown")

    # Remove replicate suffix if present (.1, .2, etc.)
    name <- gsub("\\.\\d+$", "", compound_string)

    # Handle "Construct | Compound" format
    if (grepl(" \\| ", name)) {
      parts <- strsplit(name, " \\| ")[[1]]
      return(trimws(parts[2]))
    }

    # Handle "Construct:Compound" format
    if (grepl(":", name)) {
      parts <- strsplit(name, ":")[[1]]
      return(trimws(parts[2]))
    }

    return(name)
  }

  # Helper function to extract construct name
  extract_construct_name <- function(compound_string) {
    if (is.null(compound_string)) return("Unknown")

    # Remove replicate suffix if present (.1, .2, etc.)
    name <- gsub("\\.\\d+$", "", compound_string)

    # Handle "Construct | Compound" format
    if (grepl(" \\| ", name)) {
      parts <- strsplit(name, " \\| ")[[1]]
      return(trimws(parts[1]))
    }

    # Handle "Construct:Compound" format
    if (grepl(":", name)) {
      parts <- strsplit(name, ":")[[1]]
      return(trimws(parts[1]))
    }

    return("Unknown")
  }

  # Extract drc_results if batch_drc_results is the wrapper object
  if (is.list(batch_drc_results)) {
    if ("drc_results" %in% names(batch_drc_results)) {
      if (verbose) message("Detected batch_drc_analysis wrapper. Extracting drc_results...")
      drc_results <- batch_drc_results$drc_results
    } else {
      drc_results <- batch_drc_results
    }
  } else {
    stop("batch_drc_results must be a list")
  }

  # Get plate names
  plate_names <- names(drc_results)
  if (is.null(plate_names) || length(plate_names) == 0) {
    stop("No plates found in drc_results")
  }

  # Filter plates if specified
  if (!is.null(plates_to_plot)) {
    plate_names <- intersect(plate_names, plates_to_plot)
    if (length(plate_names) == 0) {
      stop("No valid plates specified")
    }
  }

  if (verbose) {
    message("Found ", length(plate_names), " plates to process")
    message("Output directory: ", output_dir)
  }

  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ============================================================================
  # 2. SCAN FOR VALID COMPOUNDS ACROSS ALL PLATES
  # ============================================================================

  if (verbose) message("\nScanning for valid compounds...")

  compounds_list <- list()

  for (plate_name in plate_names) {
    plate <- drc_results[[plate_name]]

    # Check if plate has drc_result
    if (is.null(plate$drc_result)) {
      if (verbose > 1) message("  Skipping ", plate_name, ": no drc_result")
      next
    }

    # Get detailed results
    detailed <- plate$drc_result$detailed_results
    if (is.null(detailed) || !is.list(detailed)) {
      if (verbose > 1) message("  Skipping ", plate_name, ": no detailed_results")
      next
    }

    # For each compound in the plate
    for (i in seq_along(detailed)) {
      result <- detailed[[i]]

      # Check if fit was successful
      if (!isTRUE(result$success)) next

      # Get compound and construct names
      compound_name <- extract_compound_name(result$compound)
      construct_name <- extract_construct_name(result$compound)

      # Filter by compound if specified
      if (!is.null(compounds_to_plot) && !compound_name %in% compounds_to_plot) next

      # Store compound info
      compounds_list <- append(compounds_list, list(list(
        plate = plate_name,
        construct = construct_name,
        compound = compound_name,
        compound_full = result$compound,
        index = i,
        results_obj = plate$drc_result  # Pass the full results object for plot_dose_response
      )))
    }
  }

  if (length(compounds_list) == 0) {
    stop("No valid compounds found to plot")
  }

  if (verbose) {
    message("Found ", length(compounds_list), " valid compounds")
    message("  - Plates: ", paste(unique(sapply(compounds_list, function(x) x$plate)), collapse = ", "))
    message("  - Compounds: ", length(unique(sapply(compounds_list, function(x) x$compound))))
  }

  # ============================================================================
  # 3. CREATE DIRECTORY STRUCTURE
  # ============================================================================

  if (organize_by == "plate") {
    # Create subfolders for each plate
    for (plate_name in unique(sapply(compounds_list, function(x) x$plate))) {
      plate_dir <- file.path(output_dir, safe_filename(plate_name))
      if (!dir.exists(plate_dir)) {
        dir.create(plate_dir, recursive = TRUE)
      }
    }
  } else if (organize_by == "compound") {
    # Create subfolders for each compound
    for (compound_name in unique(sapply(compounds_list, function(x) x$compound))) {
      compound_dir <- file.path(output_dir, safe_filename(compound_name))
      if (!dir.exists(compound_dir)) {
        dir.create(compound_dir, recursive = TRUE)
      }
    }
  }

  # ============================================================================
  # 4. GENERATE ALL PLOTS
  # ============================================================================

  if (verbose) message("\nGenerating plots...")

  total <- length(compounds_list)
  successes <- 0
  failures <- 0
  failed_list <- character()
  error_messages <- list()

  # Progress bar
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }

  for (i in seq_along(compounds_list)) {
    info <- compounds_list[[i]]

    # Determine output path based on organization
    if (organize_by == "plate") {
      # plate/compound.format
      filename <- paste0(safe_filename(info$compound), ".", format)
      output_path <- file.path(output_dir, safe_filename(info$plate), filename)
    } else if (organize_by == "compound") {
      # compound/plate.format
      filename <- paste0(safe_filename(info$plate), "_", safe_filename(info$compound), ".", format)
      output_path <- file.path(output_dir, safe_filename(info$compound), filename)
    } else {
      # flat: plate_compound.format
      filename <- paste0(safe_filename(info$plate), "_", safe_filename(info$compound), ".", format)
      output_path <- file.path(output_dir, filename)
    }

    # Create directory if it doesn't exist
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

    # Generate plot
    tryCatch({
      # Call plot_dose_response with the correct parameters
      # Set point_color to black (default) and let other parameters pass through
      plot_dose_response(
        results = info$results_obj,
        compound_index = info$index,
        save_plot = output_path,
        plot_width = width,
        plot_height = height,
        plot_dpi = dpi,
        point_color = point_color,  # Passando a cor dos pontos (padrão = "black")
        show_ic50_line = show_ic50_line,
        verbose = FALSE,
        plot_title = plot_title,
        point_size = point_size,
        ...  # Pass any additional arguments to plot_dose_response
      )
      successes <- successes + 1
    }, error = function(e) {
      failures <- failures + 1
      failed_list <<- c(failed_list, paste(info$plate, info$compound, sep = "/"))
      error_messages[[length(error_messages) + 1]] <<- paste(info$plate, info$compound, ":", e$message)
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # ============================================================================
  # 5. SUMMARY AND RETURN
  # ============================================================================

  if (verbose) {
    message("\n")
    message("========================================")
    message("PLOT GENERATION COMPLETE")
    message("========================================")
    message("Total compounds: ", total)
    message("Successful: ", successes)
    message("Failed: ", failures)
    message("Point color: ", point_color)  # Mostrar a cor usada
    message("Output directory: ", normalizePath(output_dir))

    if (failures > 0) {
      message("\nFailed compounds:")
      for (f in failed_list) {
        message("  - ", f)
      }

      if (verbose > 1) {
        message("\nError details:")
        for (err in error_messages) {
          message("  - ", err)
        }
      }
    }

    # Show directory structure
    message("\nDirectory structure:")
    if (organize_by == "plate") {
      for (plate in unique(sapply(compounds_list, function(x) x$plate))) {
        plate_files <- list.files(file.path(output_dir, safe_filename(plate)), pattern = paste0("\\.", format, "$"))
        message("  ", plate, "/ : ", length(plate_files), " files")
      }
    } else if (organize_by == "compound") {
      for (compound in unique(sapply(compounds_list, function(x) x$compound))) {
        compound_files <- list.files(file.path(output_dir, safe_filename(compound)), pattern = paste0("\\.", format, "$"))
        message("  ", compound, "/ : ", length(compound_files), " files")
      }
    } else {
      all_files <- list.files(output_dir, pattern = paste0("\\.", format, "$"))
      message("  Flat structure: ", length(all_files), " files in root")
    }
  }

  # ============================================================================
  # 6. RETURN INVISIBLE SUMMARY
  # ============================================================================

  invisible(list(
    total = total,
    successes = successes,
    failures = failures,
    failed_compounds = failed_list,
    error_messages = error_messages,
    output_dir = output_dir,
    organization = organize_by,
    point_color = point_color,
    timestamp = Sys.time()
  ))
}

