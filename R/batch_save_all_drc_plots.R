#' Save All Dose-Response Curves from Batch Analysis Results
#'
#' @description
#' This function automatically generates and saves Dose-Response Curve (DRC) plots for all
#' successful fits found in batch analysis results. It creates a systematic file structure
#' and handles naming, organization, and error management.
#'
#' @param batch_drc_results A named list containing the results of the DRC batch analysis.
#'   Structure expected: \code{list(plate_name = list(drc_result = list(detailed_results = ...)))}.
#'   Can also accept a wrapper object from \code{batch_drc_analysis()} which contains a \code{drc_results} element.
#' @param output_dir Character. Directory where plots will be saved (default: "DRC_Plots").
#'   Created if it doesn't exist.
#' @param create_subfolders Logical. If \code{TRUE} (default), creates separate subfolders for each plate.
#'   If \code{FALSE}, all plots are saved directly in \code{output_dir}.
#' @param file_prefix Character. Prefix for all saved files (default: "DRC").
#' @param file_extension Character. File format extension (default: "png").
#'   Other options: "pdf", "jpg", "tiff", "svg".
#' @param overwrite Logical. If \code{TRUE}, overwrites existing files. If \code{FALSE} (default),
#'   skips files that already exist.
#' @param plates_to_process Character vector. Specific plate names to process.
#'   If \code{NULL} (default), processes all plates.
#' @param compounds_to_exclude Character vector. Compound names to exclude from plotting.
#'   Useful for removing controls or failed compounds.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages to the console.
#' @param show_legend Logical. If \code{TRUE}, displays legend in plots. If \code{FALSE} (default),
#'   omits legend for cleaner individual plots.
#' @param ... Additional arguments passed to \code{\link{plot_drc_batch}}.
#'   Common options: \code{y_limits}, \code{colors}, \code{plot_width}, \code{plot_height}, etc.
#'
#' @return An invisible list containing:
#'   \itemize{
#'     \item \code{total}: Total number of plots attempted
#'     \item \code{success}: Number of successfully saved plots
#'     \item \code{failed}: Number of failed plots
#'     \item \code{output_dir}: Path where plots were saved
#'   }
#'
#' @details
#' This function scans through all batch DRC results, identifies successful fits,
#' and generates individual plots for each construct-compound combination.
#'
#' \strong{Filename Generation:} Files are named using the pattern:
#' \code{[prefix]_[construct]_[compound].[extension]} (with plate subfolders if enabled).
#' Special characters are replaced with underscores for compatibility.
#'
#' \strong{File Organization:}
#' \itemize{
#'   \item With subfolders: \code{output_dir/plate_name/DRC_construct_compound.png}
#'   \item Without subfolders: \code{output_dir/DRC_construct_compound_plate.png}
#' }
#'
#' \strong{Error Handling:} The function continues processing even if individual plots fail,
#' recording failures in the return object and optionally printing warnings.
#'
#' @examples
#' \dontrun{
#' # 1. Save all plots with default settings
#' batch_save_all_drc_plots(batch_drc_results)
#'
#' # 2. Save only specific plates as PDFs
#' batch_save_all_drc_plots(
#'   batch_drc_results,
#'   plates_to_process = c("Plate1", "Plate2"),
#'   file_extension = "pdf",
#'   overwrite = TRUE
#' )
#'
#' # 3. Customize plot appearance and organization
#' batch_save_all_drc_plots(
#'   batch_drc_results,
#'   output_dir = "Figures/DRC_Curves",
#'   create_subfolders = FALSE,
#'   file_prefix = "Curve",
#'   y_limits = c(0, 150),
#'   colors = "viridis",
#'   plot_width = 10,
#'   plot_height = 6,
#'   plot_dpi = 300
#' )
#'
#' # 4. Exclude specific compounds and show legend
#' batch_save_all_drc_plots(
#'   batch_drc_results,
#'   compounds_to_exclude = c("DMSO", "Control"),
#'   show_legend = TRUE,
#'   verbose = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link{plot_drc_batch}} for individual plot generation
#' \code{\link{batch_drc_analysis}} for generating batch DRC results
#'
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr bind_rows
#' @export
batch_save_all_drc_plots <- function(batch_drc_results,
                                     output_dir = "DRC_Plots",
                                     create_subfolders = TRUE,
                                     file_prefix = "DRC",
                                     file_extension = "png",
                                     overwrite = FALSE,
                                     plates_to_process = NULL,
                                     compounds_to_exclude = NULL,
                                     verbose = TRUE,
                                     show_legend = FALSE,
                                     ...) {

  # ============================================================================
  # 1. SETUP AND DATA EXTRACTION
  # ============================================================================
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required")
  }

  # Extract drc_results if a wrapper object is provided
  if (is.list(batch_drc_results) && "drc_results" %in% names(batch_drc_results)) {
    if (verbose) message("Detected analysis wrapper object. Extracting 'drc_results'...")
    batch_drc_results <- batch_drc_results$drc_results
  }

  # Helper function for safe filename generation
  make_safe_filename <- function(string) {
    if (is.null(string) || is.na(string)) return("unknown")
    s <- gsub("[^[:alnum:]]+", "_", string)
    s <- gsub("^_|_$", "", s)
    return(s)
  }

  # Helper to extract detailed results from plate objects
  get_detailed_results <- function(plate_obj) {
    if (is.null(plate_obj$drc_result)) return(NULL)
    # Try different possible locations for results
    res <- plate_obj$drc_result$detailed_results
    if (is.null(res)) res <- plate_obj$drc_result$curve_results
    if (is.null(res)) res <- plate_obj$drc_result$fits
    if (is.list(res) && length(res) > 0) return(res)
    return(NULL)
  }

  # ============================================================================
  # 2. EXTRACT VALID CONSTRUCT-COMPOUND COMBINATIONS
  # ============================================================================
  extract_combinations <- function(batch_results) {
    combos_list <- list()
    all_plates <- names(batch_results)

    # Filter plates if specified
    target_plates <- if (!is.null(plates_to_process)) {
      intersect(all_plates, plates_to_process)
    } else {
      all_plates
    }

    for (plate_name in target_plates) {
      plate_res_list <- get_detailed_results(batch_results[[plate_name]])
      if (is.null(plate_res_list)) next

      for (i in seq_along(plate_res_list)) {
        res <- plate_res_list[[i]]

        # Check if fit was successful
        has_success <- !is.null(res$success) && isTRUE(res$success)
        if (!has_success) next

        # Parse construct and compound names
        info <- tryCatch({
          raw_name <- res$compound
          if (is.null(raw_name)) raw_name <- paste0("Unknown_", i)

          clean <- trimws(gsub("\\.\\d+$", "", raw_name))

          # Handle "Construct | Compound" format
          if (grepl(" \\| ", clean)) {
            parts <- strsplit(clean, " \\| ")[[1]]
            clean <- parts[1]
          }

          parts <- strsplit(clean, ":")[[1]]
          parts <- trimws(parts)

          if (length(parts) >= 2) {
            list(construct = parts[1], compound = parts[2])
          } else {
            list(construct = parts[1], compound = parts[1])
          }
        }, error = function(e) {
          list(construct = "Unknown", compound = "Unknown")
        })

        # Skip excluded compounds
        if (!is.null(compounds_to_exclude) && info$compound %in% compounds_to_exclude) {
          next
        }

        combos_list[[length(combos_list) + 1]] <- data.frame(
          plate = plate_name,
          construct = info$construct,
          compound = info$compound,
          construct_compound = paste(info$construct, info$compound, sep = ":"),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(combos_list) == 0) return(NULL)
    return(dplyr::bind_rows(combos_list))
  }

  if (verbose) message("Scanning results for successful fits...")
  combos_df <- extract_combinations(batch_drc_results)

  if (is.null(combos_df) || nrow(combos_df) == 0) {
    warning("No valid/successful DRC results found to plot.")
    return(invisible(NULL))
  }

  # ============================================================================
  # 3. SETUP OUTPUT DIRECTORY STRUCTURE
  # ============================================================================
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (create_subfolders) {
    for (plate in unique(combos_df$plate)) {
      plate_path <- file.path(output_dir, make_safe_filename(plate))
      if (!dir.exists(plate_path)) {
        dir.create(plate_path, recursive = TRUE)
      }
    }
  }

  # ============================================================================
  # 4. GENERATE AND SAVE PLOTS
  # ============================================================================
  total_plots <- nrow(combos_df)
  successes <- 0
  failures <- 0

  if (verbose) {
    message("Starting generation of ", total_plots, " plots...")
    pb <- txtProgressBar(min = 0, max = total_plots, style = 3)
  }

  for (i in seq_len(total_plots)) {
    combo <- combos_df[i, ]

    safe_construct <- make_safe_filename(combo$construct)
    safe_compound <- make_safe_filename(combo$compound)
    safe_plate <- make_safe_filename(combo$plate)

    # Generate filename
    if (create_subfolders) {
      filename <- sprintf("%s_%s_%s.%s", file_prefix, safe_construct, safe_compound, file_extension)
      output_path <- file.path(output_dir, safe_plate, filename)
    } else {
      filename <- sprintf("%s_%s_%s_%s.%s", file_prefix, safe_construct, safe_compound, safe_plate, file_extension)
      output_path <- file.path(output_dir, filename)
    }

    # Skip if file exists and overwrite is FALSE
    if (file.exists(output_path) && !overwrite) {
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    # Generate and save plot
    tryCatch({
      suppressMessages({
        plot_drc_batch(
          batch_drc_results = batch_drc_results,
          construct_compound = combo$construct_compound,
          save_plot = output_path,
          verbose = FALSE,
          show_legend = show_legend,
          plot_title = combo$compound,
          ...
        )
      })
      successes <- successes + 1
    }, error = function(e) {
      failures <- failures + 1
      if (verbose) {
        warning(sprintf("Failed to plot %s: %s", combo$construct_compound, e$message))
      }
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("\nComplete! Successful: ", successes, " | Failed: ", failures)
    message("Plots saved to: ", normalizePath(output_dir))
  }

  # ============================================================================
  # 5. RETURN SUMMARY
  # ============================================================================
  return(invisible(list(
    total = total_plots,
    success = successes,
    failed = failures,
    output_dir = output_dir
  )))
}
