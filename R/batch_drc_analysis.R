#' Batch Dose–Response Curve (DRC) Analysis for Multiple Plates
#'
#' @description
#' `batch_drc_analysis()` performs automated dose–response curve (DRC) fitting
#' across multiple plates previously processed by `batch_ratio_analysis()`.
#'
#' It extracts the `modified_ratio_table` from each plate, applies a
#' 3-parameter logistic regression via [`fit_drc_3pl()`], and generates:
#'
#' * Per-plate DRC result files (optional)
#' * A consolidated batch DRC report (`batch_drc_analysis_report.xlsx`)
#'
#' The function computes fitting success rates, R-squared metrics, curve
#' qualities, and compiles them into structured summary sheets.
#'
#'
#' @param batch_results A named list containing results from
#'   [`batch_ratio_analysis()`]. Each element must contain the field
#'   `result$modified_ratio_table`, which will be used for the DRC fit.
#' @param normalize Logical. Whether to normalize responses inside `fit_drc_3pl()`.
#' @param enforce_bottom_threshold Logical. If `TRUE`, forces the lower plateau of
#'   the fitted curve to stay above `bottom_threshold`.
#' @param bottom_threshold Numeric. Minimum acceptable bottom asymptote value.
#' @param r_sqr_threshold Minimum acceptable R-squared for accepting a curve fit.
#' @param output_dir Directory where individual plate results and consolidated
#'   batch reports will be saved. Defaults to the working directory.
#' @param generate_reports Logical. If `TRUE` (default), generates a consolidated
#'   Excel report summarizing all DRC fits.
#' @param verbose Logical. If `TRUE`, prints progress details.
#'
#'
#' @details
#' For each plate in `batch_results`, the function:
#'
#' 1. Extracts the `modified_ratio_table` produced during ratio normalization.
#' 2. Ensures the table contains valid data (non-empty and with column names).
#' 3. Performs a 3-parameter logistic fit via [`fit_drc_3pl()`], where:
#'    * The first column is assumed to be log(inhibitor concentration)
#'    * Remaining columns are compound responses
#' 4. Stores:
#'    * Full DRC results for each compound
#'    * Summary tables, final tables, and curve quality metrics
#'    * Optional Excel result files per plate
#'
#' After processing all plates, if `generate_reports = TRUE`, a consolidated
#' Excel file (`batch_drc_analysis_report.xlsx`) is created, containing:
#'
#' * **Summary** sheet — plate-level statistics
#' * **All_Results** — merged results for all compounds across plates
#' * **Curve_Quality** — QC-oriented summary (R², slope, quality flag)
#' * A **_summary** and **_final_summary** sheet per plate
#'
#'
#' @return
#' A named list of DRC analysis results for each plate.
#' Each entry contains:
#'
#' \describe{
#'   \item{plate_info}{A list with metadata: `data_file`, `info_sheet`, `sheet_number`.}
#'   \item{drc_result}{Full result object from `fit_drc_3pl()`, including fitted parameters,
#'   summary tables, final summary tables, QC metrics, and counts.}
#'   \item{timestamp}{Time when the plate was processed.}
#' }
#'
#' If `generate_reports = TRUE`, also saves:
#'
#' * `batch_drc_analysis_report.xlsx`
#' * Individual Excel files: `drc_results_<plate>.xlsx`
#'
#'
#' @examples
#' \dontrun{
#'
#' # Assuming you already ran batch_ratio_analysis()
#' ratio_results <- batch_ratio_analysis("experiment_folder")
#'
#' # Run DRC for all plates
#' drc_results <- batch_drc_analysis(ratio_results)
#'
#' # Specify thresholds and output directory
#' drc_results <- batch_drc_analysis(
#'   batch_results = ratio_results,
#'   r_sqr_threshold = 0.9,
#'   bottom_threshold = 50,
#'   output_dir = "drc_output"
#' )
#'
#' # Skip generating the Excel report
#' batch_drc_analysis(
#'   batch_results = ratio_results,
#'   generate_reports = FALSE
#' )
#' }
#'
#'
#' @seealso
#' * [`fit_drc_3pl()`] — Fits the dose–response curve for a single plate.
#' * [`batch_ratio_analysis()`] — Preprocessing step generating modified ratio tables.
#' * `openxlsx` — Excel manipulation used for reporting.
#'
#'
#' @export


batch_drc_analysis <- function(batch_results,
                               normalize = FALSE,
                               enforce_bottom_threshold = FALSE,
                               bottom_threshold = 60,
                               r_sqr_threshold = 0.8,
                               output_dir = NULL,
                               generate_reports = TRUE,
                               verbose = TRUE) {

  # ============================================================================
  # 0. SETUP & DEPENDENCIES
  # ============================================================================
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (generate_reports && !requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required.")

  # Helper: Operador %||% para fallback seguro
  `%||%` <- function(a, b) {
    if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
  }

  # Helper: Sanitização de nome de arquivo
  sanitize_filename <- function(name, max_len = 50) {
    clean <- gsub("[^A-Za-z0-9_.-]", "_", name)
    clean <- gsub("_+", "_", clean)
    if (nchar(clean) > max_len) clean <- substr(clean, 1, max_len)
    return(clean)
  }

  # Validação de Input
  if (!is.list(batch_results) || length(batch_results) == 0) {
    stop("batch_results must be a non-empty list.")
  }

  # Definição de Diretórios
  if (is.null(output_dir)) output_dir <- getwd()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Cria subpasta APENAS se gerar relatórios
  detailed_dir <- file.path(output_dir, "Detailed_Reports")
  if (generate_reports && !dir.exists(detailed_dir)) {
    dir.create(detailed_dir, recursive = TRUE)
  }

  # ============================================================================
  # 1. FUNÇÕES INTERNAS
  # ============================================================================

  generate_drc_batch_report <- function(drc_results, main_dir, sub_dir, verbose = TRUE) {

    # Helper para nome de abas seguro
    get_safe_sheet_name <- function(base_name, suffix = "", existing_names = c()) {
      max_len <- 31 - nchar(suffix) - 3
      clean_base <- gsub("[^A-Za-z0-9_]", "_", base_name)
      candidate_base <- substr(clean_base, 1, max_len)

      final_name <- paste0(candidate_base, suffix)
      counter <- 1
      while (final_name %in% existing_names) {
        final_name <- paste0(candidate_base, "_", counter, suffix)
        counter <- counter + 1
      }
      return(final_name)
    }

    path_pharma <- file.path(main_dir, "Pharmacology_Summary.xlsx")
    path_details <- file.path(sub_dir, "Batch_Analysis_Details.xlsx")

    wb_pharma <- openxlsx::createWorkbook()
    wb_details <- openxlsx::createWorkbook()

    # Listas acumuladoras
    summary_list <- list()
    all_results_list <- list()
    quality_list <- list()
    pharm_list <- list()

    used_sheet_names <- c("Summary", "All_Results", "Curve_Quality")

    # --- LOOP DE GERAÇÃO DO REPORT ---
    for (plate_name in names(drc_results)) {
      plate_res_obj <- drc_results[[plate_name]]

      # Validação de Estrutura
      if (is.null(plate_res_obj$drc_result) || is.null(plate_res_obj$drc_result$summary_table)) next
      drc_summary <- plate_res_obj$drc_result$summary_table
      if (nrow(drc_summary) == 0) next

      # 1. Standard Summary
      n_compounds <- nrow(drc_summary)
      successful_fits <- sum(!is.na(drc_summary$IC50))

      summary_list[[length(summary_list) + 1]] <- data.frame(
        Plate_Name = plate_name,
        Data_File = plate_res_obj$plate_info$data_file,
        N_Compounds = n_compounds,
        Successful_Fits = successful_fits,
        Success_Rate = round(successful_fits / n_compounds * 100, 1),
        Avg_R_squared = round(mean(drc_summary$R_squared, na.rm = TRUE), 3),
        Good_Curves = sum(grepl("Good curve", drc_summary$Curve_Quality, ignore.case = TRUE)),
        stringsAsFactors = FALSE
      )

      # 2. All Results (Sem pipe)
      drc_summary_copy <- drc_summary
      drc_summary_copy$Plate <- plate_name

      # Reordena colunas usando subsetting básico
      cols_order <- c("Plate", setdiff(names(drc_summary_copy), "Plate"))
      drc_summary_copy <- drc_summary_copy[, cols_order, drop = FALSE]

      all_results_list[[length(all_results_list) + 1]] <- drc_summary_copy

      # 3. Quality List (Evita colunas duplicadas explicitamente)
      desired_cols <- c("Plate", "Compound", "Curve_Quality", "R_squared", "Max_Slope",
                        "Ideal_Hill_Slope", "Bottom", "Top", "LogIC50", "IC50")
      existing_cols <- intersect(desired_cols, names(drc_summary_copy))
      quality_list[[length(quality_list) + 1]] <- drc_summary_copy[, existing_cols, drop = FALSE]

      # 4. PHARMACOLOGY SUMMARY (Robusto)
      res_root <- plate_res_obj$drc_result
      detailed_res <- res_root$detailed_results %||% res_root$curve_results %||% res_root$fits
      if (!is.list(detailed_res)) detailed_res <- list()

      if (length(detailed_res) > 0) {

        for (i in seq_along(detailed_res)) {
          res <- detailed_res[[i]]
          # Verifica se fit existiu
          if (is.null(res$success) || !isTRUE(res$success)) next

          # Parsing Name
          clean_name <- gsub("\\.\\d+$", "", res$compound %||% "Unknown")
          parts <- strsplit(clean_name, " \\| |:")[[1]]
          parts <- trimws(parts)
          construct_name <- parts[1]
          compound_name  <- if(length(parts) > 1) parts[2] else parts[1]

          # --- pIC50 ---
          log_ic50 <- NA_real_
          if (!is.null(res$parameters) && length(res$parameters$Value) >= 3) {
            log_ic50 <- res$parameters$Value[3]
          }
          pic50 <- if (!is.na(log_ic50)) -log_ic50 else NA_real_

          # --- CI e DELTA POSITIVO ---
          ci_log_lower_bound <- NA_real_
          ci_log_upper_bound <- NA_real_

          if (!is.null(res$confidence_intervals) && !is.null(res$confidence_intervals$LogIC50)) {
            ci_log_lower_bound <- res$confidence_intervals$LogIC50[1]
            ci_log_upper_bound <- res$confidence_intervals$LogIC50[2]
          }

          pic50_diff_upper <- NA_real_
          pic50_diff_lower <- NA_real_

          if (!is.na(pic50) && !is.na(ci_log_lower_bound) && !is.na(ci_log_upper_bound)) {
            abs_pic50_upper <- -ci_log_lower_bound
            abs_pic50_lower <- -ci_log_upper_bound

            # Deltas Positivos (Magnitude)
            pic50_diff_upper <- abs_pic50_upper - pic50
            pic50_diff_lower <- pic50 - abs_pic50_lower
          }

          # --- NORMALIZED SPAN (Fórmula Solicitada) ---
          norm_span <- NA_real_
          if (!is.null(res$data) && nrow(res$data) >= 2) {
            # Ordenação segura
            d_ord <- res$data[order(res$data$log_inhibitor), ]

            min_c <- d_ord$log_inhibitor[1]
            max_c <- d_ord$log_inhibitor[nrow(d_ord)]

            # Identifica índices para média
            idx_min <- which(d_ord$log_inhibitor == min_c)
            idx_max <- which(d_ord$log_inhibitor == max_c)

            if (length(idx_min) > 0 && length(idx_max) > 0) {
              c0   <- mean(d_ord$response[idx_min], na.rm = TRUE)
              c100 <- mean(d_ord$response[idx_max], na.rm = TRUE)

              denom <- c100 - c0

              # span_val geralmente é parâmetro 5 (Span) vindo do fit_drc_3pl
              span_val <- NA_real_
              if (length(res$parameters$Value) >= 5) span_val <- res$parameters$Value[5]

              # FÓRMULA: 1 - (span / (c100 - c0))
              # OBSERVAÇÃO: Com esta fórmula:
              # - Se span = denom (cobre toda janela) → norm_span = 0
              # - Se span = 0.5*denom (cobre metade) → norm_span = 0.5
              # - Se span = 0 (não cobre) → norm_span = 1
              # Portanto, norm_span < 0.5 significa que a curva cobre MAIS da metade da janela
              if (!is.na(denom) && abs(denom) > 1e-6 && !is.na(span_val)) {
                norm_span <- 1 - (span_val / denom)
              }
            }
          }

          ideal_hill <- res$ideal_hill_slope %||% NA_real_

          # --- FLAGS ---
          flag_collector <- character()

          # 1. CI indefinido
          if (is.na(pic50_diff_lower) || is.na(pic50_diff_upper)) {
            flag_collector <- c(flag_collector, "Undefined CI")
          } else {
            # 2. CI muito amplo (> 0.4642 log units ≈ 1.5-fold em escala linear)
            if (pic50_diff_lower > 0.4642) flag_collector <- c(flag_collector, "Unstable pIC50 (Lower CI > 0.4642)")
            if (pic50_diff_upper > 0.4642) flag_collector <- c(flag_collector, "Unstable pIC50 (Upper CI > 0.4642)")
          }

          # 3. Hill slope fora do ideal
          if (!is.na(ideal_hill) && (ideal_hill < 0.5 || ideal_hill > 1.5)) {
            flag_collector <- c(flag_collector, "Hill Slope < 0.5 or > 1.5")
          }

          # 4. Normalized Span < 0.5 (CRITÉRIO SOLICITADO)
          # Com a fórmula acima, isso significa que a curva cobre MAIS que 50% da janela dinâmica
          if (!is.na(norm_span) && norm_span < 0.5) {
            flag_collector <- c(flag_collector, "Normalized Span < 0.5")
          }

          final_flags <- if(length(flag_collector) > 0) {
            paste(flag_collector, collapse = "; ")
          } else {
            "OK"
          }

          pharm_list[[length(pharm_list) + 1]] <- data.frame(
            Plate = plate_name,
            Construct = construct_name,
            Compound = compound_name,
            pIC50 = round(pic50, 3),
            CI_95_Upper = round(pic50_diff_upper, 3),
            CI_95_Lower = round(pic50_diff_lower, 3),
            Ideal_Hill_Slope = round(ideal_hill, 3),
            Normalized_Span = round(norm_span, 3),
            Flags = final_flags,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # Consolidação (com verificação de listas vazias)
    if (length(summary_list) == 0) {
      summary_data <- data.frame()
    } else {
      summary_data <- dplyr::bind_rows(summary_list)
    }

    if (length(all_results_list) == 0) {
      all_results_combined <- data.frame()
    } else {
      all_results_combined <- dplyr::bind_rows(all_results_list)
    }

    if (length(quality_list) == 0) {
      quality_combined <- data.frame()
    } else {
      quality_combined <- dplyr::bind_rows(quality_list)
    }

    if (length(pharm_list) == 0) {
      pharm_combined <- data.frame()
    } else {
      pharm_combined <- dplyr::bind_rows(pharm_list)
    }

    # --- ESCRITA EXCEL ---

    # Arquivo 1: Pharmacology Summary
    openxlsx::addWorksheet(wb_pharma, "Pharmacology_Summary")
    if (nrow(pharm_combined) > 0) {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary", pharm_combined)
    } else {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary",
                          data.frame(Note = "No pharmacology data available"))
    }
    openxlsx::saveWorkbook(wb_pharma, path_pharma, overwrite = TRUE)

    # Arquivo 2: Detailed Report
    openxlsx::addWorksheet(wb_details, "Summary")
    openxlsx::writeData(wb_details, "Summary", summary_data)

    openxlsx::addWorksheet(wb_details, "All_Results")
    if (nrow(all_results_combined) > 0) {
      openxlsx::writeData(wb_details, "All_Results", all_results_combined)
    }

    openxlsx::addWorksheet(wb_details, "Curve_Quality")
    if (nrow(quality_combined) > 0) {
      openxlsx::writeData(wb_details, "Curve_Quality", quality_combined)
    }

    # Worksheets individuais para cada placa
    for (plate_name in names(drc_results)) {
      plate_res_obj <- drc_results[[plate_name]]
      if (is.null(plate_res_obj$drc_result)) next
      if (is.null(plate_res_obj$drc_result$summary_table)) next

      sheet_sum <- get_safe_sheet_name(plate_name, "_sum", used_sheet_names)
      used_sheet_names <- c(used_sheet_names, sheet_sum)

      openxlsx::addWorksheet(wb_details, sheet_sum)
      openxlsx::writeData(wb_details, sheet_sum, plate_res_obj$drc_result$summary_table)
    }

    openxlsx::saveWorkbook(wb_details, path_details, overwrite = TRUE)

    if (verbose) {
      message("Reports generated successfully:")
      message("  1. ", path_pharma)
      message("  2. ", path_details)
    }
  }

  # Função robusta para extração de dados
  extract_data_for_drc <- function(plate_result) {
    # Busca hierárquica por dados válidos
    search_locations <- list()

    # Nível 1: Via $result (estrutura mais comum)
    if (!is.null(plate_result$result)) {
      search_locations <- c(search_locations, list(
        plate_result$result$modified_ratio_table,
        plate_result$result$processed_data,
        plate_result$result$ratio_table,
        plate_result$result$normalized_data,
        plate_result$result$data
      ))
    }

    # Nível 2: Diretamente nos atributos
    search_locations <- c(search_locations, list(
      plate_result$modified_ratio_table,
      plate_result$processed_data,
      plate_result$ratio_table,
      plate_result$normalized_data,
      plate_result$data
    ))

    # Nível 3: A própria plate_result se for data.frame
    if (is.data.frame(plate_result) && nrow(plate_result) > 0) {
      search_locations <- c(search_locations, list(plate_result))
    }

    # Busca sequencial
    for (dt in search_locations) {
      if (!is.null(dt) && is.data.frame(dt) && nrow(dt) > 0 && ncol(dt) >= 3) {
        # Verifica se primeira coluna parece ser numérica (concentração)
        first_col <- dt[[1]]
        if (is.numeric(first_col) ||
            (is.character(first_col) &&
             all(!is.na(suppressWarnings(as.numeric(first_col[!is.na(first_col)])))))) {
          return(dt)
        }
      }
    }
    return(NULL)
  }

  # Função simples de preparação (sem criação de réplicas artificiais)
  prepare_drc_data <- function(data_table) {
    if (is.null(data_table)) return(NULL)
    # Retorna os dados sem modificação
    return(data_table)
  }

  # ============================================================================
  # 2. MAIN EXECUTION
  # ============================================================================

  if (verbose) {
    message("==========================================================")
    message("STARTING BATCH DOSE-RESPONSE ANALYSIS")
    message("==========================================================")
    message("Main Output: ", output_dir)
  }

  drc_results <- list()
  failed_plates <- character()
  total_plates <- length(batch_results)

  # Loop principal com contador
  for (i in seq_along(batch_results)) {
    plate_name <- names(batch_results)[i]
    if (verbose) {
      message(sprintf("\nProcessing %d/%d: %s", i, total_plates, plate_name))
    }

    tryCatch({
      proc_start <- Sys.time()

      # 1. Extração de dados
      data_table <- extract_data_for_drc(batch_results[[plate_name]])
      if (is.null(data_table)) {
        stop("No valid data table found in plate result.")
      }

      # Validação básica do formato
      if (verbose) {
        message(sprintf("  Data shape: %d rows × %d columns",
                        nrow(data_table), ncol(data_table)))

        # Verifica se primeira coluna é numérica
        if (!is.numeric(data_table[[1]])) {
          message("  Note: First column may not be numeric (check concentration values)")
        }
      }

      # 2. Preparação (Sem duplicatas artificiais)
      drc_data <- prepare_drc_data(data_table)

      # 3. Caminho do arquivo individual
      output_file <- NULL
      if (generate_reports) {
        clean_name <- sanitize_filename(plate_name)
        output_file <- file.path(detailed_dir, paste0("drc_", clean_name, ".xlsx"))
      }

      # 4. Ajuste da Curva com tratamento de erro interno
      plate_drc_result <- tryCatch({
        fit_drc_3pl(
          data = drc_data,
          output_file = output_file,
          normalize = normalize,
          verbose = FALSE,
          enforce_bottom_threshold = enforce_bottom_threshold,
          bottom_threshold = bottom_threshold,
          r_sqr_threshold = r_sqr_threshold
        )
      }, error = function(e) {
        # Retorna estrutura padrão em caso de erro interno
        return(list(
          successful_fits = 0,
          n_compounds = ncol(drc_data) - 1,
          summary_table = data.frame(),
          detailed_results = list(),
          error = e$message
        ))
      })

      # 5. Metadados
      d_file <- batch_results[[plate_name]]$data_file %||% "unknown"
      i_sheet <- batch_results[[plate_name]]$info_sheet %||% "unknown"
      s_num <- batch_results[[plate_name]]$sheet_number %||% "unknown"

      drc_results[[plate_name]] <- list(
        plate_info = list(
          original_name = plate_name,
          data_file = d_file,
          info_sheet = i_sheet,
          sheet_number = s_num
        ),
        drc_result = plate_drc_result,
        processing_timestamp = Sys.time(),
        processing_time = as.numeric(difftime(Sys.time(), proc_start, units = "secs"))
      )

      if (verbose) {
        succ <- plate_drc_result$successful_fits %||% 0
        tot <- plate_drc_result$n_compounds %||% 0
        proc_time <- difftime(Sys.time(), proc_start, units = "secs")
        message(sprintf("  -> Success: %d/%d compounds (%.1f sec)", succ, tot, proc_time))

        if (!is.null(plate_drc_result$error)) {
          message(sprintf("  -> DRC warning: %s", plate_drc_result$error))
        }
      }

    }, error = function(e) {
      warning(sprintf("Failed to process plate '%s': %s", plate_name, e$message))
      failed_plates <- c(failed_plates, plate_name)
      if (verbose) message(sprintf("  -> Error: %s", e$message))
    })
  }

  # ============================================================================
  # 3. REPORTING & RETURN
  # ============================================================================

  report_info <- NULL
  if (length(drc_results) > 0 && generate_reports) {
    if (verbose) {
      message("\n", paste(rep("=", 50), collapse = ""))
      message("Generating consolidated reports...")
    }

    tryCatch({
      report_info <- generate_drc_batch_report(drc_results, output_dir, detailed_dir, verbose)
    }, error = function(e) {
      warning("Failed to generate master reports: ", e$message)
      if (verbose) message("  -> Report generation failed: ", e$message)
    })
  }

  # ============================================================================
  # 4. FINAL SUMMARY
  # ============================================================================

  if (verbose) {
    message("\n", paste(rep("=", 50), collapse = ""))
    message("BATCH ANALYSIS COMPLETE")
    message(paste(rep("-", 50), collapse = ""))
    message(sprintf("Total Processed:    %d", total_plates))
    message(sprintf("Successful Plates:  %d", length(drc_results)))
    message(sprintf("Failed Plates:      %d", length(failed_plates)))

    if (length(failed_plates) > 0) {
      message("\nFailed Plates List:")
      if (length(failed_plates) <= 10) {
        for (fp in failed_plates) message(sprintf("  - %s", fp))
      } else {
        message(sprintf("  - First 10: %s", paste(failed_plates[1:10], collapse = ", ")))
      }
    }

    if (generate_reports) {
      message("\nOutput Files:")
      message(sprintf("  • Pharmacology Summary: %s", file.path(output_dir, "Pharmacology_Summary.xlsx")))
      message(sprintf("  • Detailed Reports:      %s", detailed_dir))
    }
    message(paste(rep("=", 50), collapse = ""))
  }

  # ============================================================================
  # 5. RETURN RESULTS
  # ============================================================================

  return(invisible(list(
    drc_results = drc_results,
    metadata = list(
      total_plates = total_plates,
      success_count = length(drc_results),
      failed_plates = failed_plates,
      output_directory = output_dir,
      detailed_reports_dir = if (generate_reports) detailed_dir else NULL,
      timestamp = Sys.time(),
      processing_summary = if (length(drc_results) > 0) {
        # Calcula tempo médio de processamento
        times <- sapply(drc_results, function(x) x$processing_time %||% 0)
        list(
          avg_time = mean(times),
          total_time = sum(times)
        )
      } else NULL
    ),
    settings = list(
      normalize = normalize,
      bottom_threshold = bottom_threshold,
      r_sqr_threshold = r_sqr_threshold,
      enforce_bottom = enforce_bottom_threshold,
      generate_reports = generate_reports
    ),
    report_info = report_info
  )))
}
