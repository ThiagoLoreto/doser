#' Batch Save All DRC Plots
#'
#' Automatically iterates through a batch of Dose-Response Curve (DRC) analysis results,
#' generates plots for all valid construct-compound combinations, and saves them to a specified directory.
#'
#' @description
#' This function extracts all successful models from the provided `batch_drc_results` list,
#' parses the construct and compound names, sanitizes them for file creation, and saves
#' the corresponding plots. It supports creating subfolders per plate, filtering specific plates,
#' and excluding specific compounds.
#'
#' @param batch_drc_results A list containing the results of a batch DRC analysis.
#'   The structure is expected to be a list of plates, where each plate contains
#'   `drc_result$detailed_results`.
#' @param output_dir Character string. The main directory where plots will be saved.
#'   Defaults to "DRC_Plots". The directory is created if it does not exist.
#' @param create_subfolders Logical. If \code{TRUE} (default), creates a subfolder
#'   for each plate within \code{output_dir}. If \code{FALSE}, all plots are saved
#'   directly in \code{output_dir} with the plate name appended to the filename to prevent duplicates.
#' @param file_prefix Character string. A prefix added to the start of every filename.
#'   Defaults to "DRC".
#' @param file_extension Character string. The file format/extension for the saved plots
#'   (e.g., "png", "pdf", "svg", "tiff"). Defaults to "png".
#' @param overwrite Logical. If \code{TRUE}, existing files with the same name will be overwritten.
#'   If \code{FALSE} (default), existing files are skipped to save time.
#' @param plates_to_process Character vector (optional). A list of specific plate names
#'   (keys in `batch_drc_results`) to process. If \code{NULL} (default), all plates are processed.
#' @param compounds_to_exclude Character vector (optional). A list of compound names
#'   to exclude from plotting. Matches against the parsed compound name.
#' @param verbose Logical. If \code{TRUE} (default), displays a progress bar and
#'   summary statistics (success/failure counts) in the console.
#' @param show_legend Logical. Whether to display the legend in the generated plots.
#'   Defaults to \code{FALSE} for batch processing to save space or keep plots clean.
#'   Can be overridden to \code{TRUE} if needed.
#' @param ... Additional arguments passed directly to the underlying \code{plot_drc_batch} function
#'   and/or the plot saving mechanism. Common useful arguments include:
#'   \itemize{
#'     \item \code{plot_width}: Width of the output image (in inches).
#'     \item \code{plot_height}: Height of the output image (in inches).
#'     \item \code{plot_dpi}: Resolution (e.g., 300).
#'     \item \code{y_limits}: A numeric vector of length 2 (e.g., \code{c(0, 1.5)}).
#'     \item \code{colors}: A vector of colors for the curves.
#'     \item \code{show_IC50}: Logical, to display the IC50 value on the plot.
#'   }
#'
#' @return An invisible list containing processing statistics:
#'   \itemize{
#'     \item \code{total}: Total number of plots attempted.
#'     \item \code{success}: Number of successfully saved plots.
#'     \item \code{failed}: Number of failed plots.
#'   }
#'
#' @details
#' The function uses a robust sanitization method to ensure filenames do not contain
#' illegal characters (converting special characters to underscores).
#' It assumes the compound names in the results are formatted as "Construct | Compound"
#' or "Construct:Compound".
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' batch_save_all_drc_plots(my_batch_results)
#'
#' # Customize output with higher resolution, no subfolders, and specific colors
#' batch_save_all_drc_plots(
#'   batch_drc_results = my_batch_results,
#'   output_dir = "Final_Plots",
#'   create_subfolders = FALSE,
#'   plot_dpi = 600,
#'   colors = c("black", "red"),
#'   show_IC50 = TRUE
#' )
#' }
#'
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
  # 1. SETUP & EXTRACTION
  # ============================================================================
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")

  # --- AJUSTE 1: Extração robusta do objeto de análise (Wrapper) ---
  if (is.list(batch_drc_results) && "drc_results" %in% names(batch_drc_results)) {
    if (verbose) message("Detected analysis wrapper object. Extracting 'drc_results'...")
    # Trabalhamos apenas com a lista de resultados interna
    batch_drc_results <- batch_drc_results$drc_results
  }

  # Helper: Nome de arquivo seguro
  make_safe_filename <- function(string) {
    if (is.null(string) || is.na(string)) return("unknown")
    s <- gsub("[^[:alnum:]]+", "_", string)
    s <- gsub("^_|_$", "", s)
    return(s)
  }

  # --- AJUSTE 2: Helper para buscar resultados onde quer que estejam ---
  get_detailed_results <- function(plate_obj) {
    if (is.null(plate_obj$drc_result)) return(NULL)
    res <- plate_obj$drc_result$detailed_results
    if (is.null(res)) res <- plate_obj$drc_result$curve_results
    if (is.null(res)) res <- plate_obj$drc_result$fits
    if (is.list(res) && length(res) > 0) return(res)
    return(NULL)
  }

  # ============================================================================
  # 2. EXTRACT COMBINATIONS TO PLOT
  # ============================================================================
  extract_combinations <- function(batch_results) {
    combos_list <- list()
    all_plates <- names(batch_results)

    target_plates <- if (!is.null(plates_to_process)) intersect(all_plates, plates_to_process) else all_plates

    for (plate_name in target_plates) {
      # Usa o helper robusto aqui
      plate_res_list <- get_detailed_results(batch_results[[plate_name]])

      if (is.null(plate_res_list)) next

      for (i in seq_along(plate_res_list)) {
        res <- plate_res_list[[i]]

        # Verifica sucesso (compatível com nova estrutura)
        has_success <- !is.null(res$success) && isTRUE(res$success)
        if (!has_success) next

        # Parsing seguro do nome
        info <- tryCatch({
          raw_name <- res$compound
          if (is.null(raw_name)) raw_name <- paste0("Unknown_", i)

          clean <- trimws(gsub("\\.\\d+$", "", raw_name))

          # Lógica para Target:Compound ou Target | Compound
          if (grepl(" \\| ", clean)) {
            parts <- strsplit(clean, " \\| ")[[1]]
            clean <- parts[1]
          }

          parts <- strsplit(clean, ":")[[1]]
          parts <- trimws(parts)

          if(length(parts) >= 2) {
            list(const = parts[1], comp = parts[2])
          } else {
            list(const = parts[1], comp = parts[1])
          }
        }, error = function(e) list(const = "Unknown", comp = "Unknown"))

        if (!is.null(compounds_to_exclude) && info$comp %in% compounds_to_exclude) next

        combos_list[[length(combos_list) + 1]] <- data.frame(
          plate = plate_name,
          construct = info$const,
          compound = info$comp,
          # Chave única para o plot_drc_batch encontrar o dado específico
          construct_compound = paste(info$const, info$comp, sep = ":"),
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
  # 3. SETUP DIRECTORIES
  # ============================================================================
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (create_subfolders) {
    for (p in unique(combos_df$plate)) {
      p_path <- file.path(output_dir, make_safe_filename(p))
      if (!dir.exists(p_path)) dir.create(p_path, recursive = TRUE)
    }
  }

  # ============================================================================
  # 4. GENERATE PLOTS
  # ============================================================================
  total_plots <- nrow(combos_df)
  successes <- 0
  failures <- 0

  if (verbose) {
    message("Starting generation of ", total_plots, " plots...")
    pb <- txtProgressBar(min = 0, max = total_plots, style = 3)
  }

  for (i in 1:total_plots) {
    combo <- combos_df[i, ]

    safe_construct <- make_safe_filename(combo$construct)
    safe_compound  <- make_safe_filename(combo$compound)
    safe_plate     <- make_safe_filename(combo$plate)

    # Nome do arquivo
    fname <- sprintf("%s_%s_%s.%s", file_prefix, safe_construct, safe_compound, file_extension)

    # Caminho completo
    if (create_subfolders) {
      output_path <- file.path(output_dir, safe_plate, fname)
    } else {
      fname_flat <- sprintf("%s_%s_%s_%s.%s", file_prefix, safe_construct, safe_compound, safe_plate, file_extension)
      output_path <- file.path(output_dir, fname_flat)
    }

    # Verifica overwrite
    if (file.exists(output_path) && !overwrite) {
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    # Tenta plotar
    tryCatch({
      suppressMessages({
        # Passamos a lista 'unwrapped' batch_drc_results aqui
        # Como plot_drc_batch já é robusta, ela vai funcionar bem
        plot_drc_batch(
          batch_drc_results = batch_drc_results,
          target_compound = combo$construct_compound, # Passa a chave "Target:Compound"
          save_plot = output_path,
          verbose = FALSE,
          show_legend = show_legend,
          plot_title = combo$compound, # Força o título simples (Nome do composto)
          ...
        )
      })
      successes <- successes + 1
    }, error = function(e) {
      failures <- failures + 1
      # Opcional: imprimir erro se falhar muito
      # warning(paste("Failed:", output_path, "-", e$message))
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("\nDone! Success: ", successes, " | Failed: ", failures)
    message("Plots saved in: ", output_dir)
  }

  return(invisible(list(total = total_plots, success = successes, failed = failures)))
}
