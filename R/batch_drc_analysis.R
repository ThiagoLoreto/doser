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
  
  # Internal function to generate DRC batch report
  generate_drc_batch_report <- function(drc_results, output_dir, verbose = TRUE) {
    
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for generating reports.")
    }
    
    report_path <- file.path(output_dir, "batch_drc_analysis_report.xlsx")
    wb <- openxlsx::createWorkbook()
    
    # Summary sheet
    openxlsx::addWorksheet(wb, "Summary")
    
    summary_data <- data.frame(
      Plate_Name = character(),
      Data_File = character(),
      N_Compounds = integer(),
      Successful_Fits = integer(),
      Success_Rate = numeric(),
      Avg_R_squared = numeric(),
      stringsAsFactors = FALSE
    )
    
    # All results consolidated sheet
    openxlsx::addWorksheet(wb, "All_Results")
    all_results_list <- list()
    
    # Process each plate for summary statistics
    for (plate_name in names(drc_results)) {
      plate_result <- drc_results[[plate_name]]
      drc_summary <- plate_result$drc_result$summary_table
      
      # Calculate summary statistics
      n_compounds <- nrow(drc_summary)
      successful_fits <- sum(!is.na(drc_summary$IC50))
      success_rate <- round(successful_fits / n_compounds * 100, 1)
      avg_r_squared <- round(mean(drc_summary$R_squared, na.rm = TRUE), 3)
      
      summary_data <- rbind(summary_data, data.frame(
        Plate_Name = plate_name,
        Data_File = plate_result$plate_info$data_file,
        N_Compounds = n_compounds,
        Successful_Fits = successful_fits,
        Success_Rate = success_rate,
        Avg_R_squared = avg_r_squared,
        stringsAsFactors = FALSE
      ))
      
      # Add plate results with identification
      plate_results_with_id <- cbind(Plate = plate_name, drc_summary)
      all_results_list[[plate_name]] <- plate_results_with_id
    }
    
    # Combine all results
    all_results_combined <- do.call(rbind, all_results_list)
    rownames(all_results_combined) <- NULL
    
    # Write data to sheets
    openxlsx::writeData(wb, "Summary", summary_data)
    openxlsx::writeData(wb, "All_Results", all_results_combined)
    
    # Add individual sheets for each plate
    for (plate_name in names(drc_results)) {
      plate_result <- drc_results[[plate_name]]
      
      # Plate summary table
      openxlsx::addWorksheet(wb, paste0(plate_name, "_summary"))
      openxlsx::writeData(wb, paste0(plate_name, "_summary"), plate_result$drc_result$summary_table)
      
      # Final summary table (transposed)
      openxlsx::addWorksheet(wb, paste0(plate_name, "_final_summary"))
      openxlsx::writeData(wb, paste0(plate_name, "_final_summary"), 
                          plate_result$drc_result$final_summary_table, rowNames = TRUE)
    }
    
    # Curve quality sheet
    openxlsx::addWorksheet(wb, "Curve_Quality")
    quality_data <- all_results_combined[, c("Plate", "Compound", "Curve_Quality", "R_squared", "Max_Slope")]
    openxlsx::writeData(wb, "Curve_Quality", quality_data)
    
    openxlsx::saveWorkbook(wb, report_path, overwrite = TRUE)
    
    if (verbose) {
      message("Batch DRC report saved: ", report_path)
    }
  }
  
  # ========== MAIN FUNCTION BODY ==========
  
  # Validate input
  if (length(batch_results) == 0) {
    stop("Batch results list is empty")
  }
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (verbose) {
    message("Starting batch DRC analysis for ", length(batch_results), " plates")
    message("Output directory: ", output_dir)
    message("Generate reports: ", generate_reports)
  }
  
  # Store DRC results for each plate
  drc_results <- list()
  
  # Process each plate
  for (plate_name in names(batch_results)) {
    
    if (verbose) message("\nProcessing DRC for: ", plate_name)
    
    tryCatch({
      # Extract modified_ratio_table from plate results
      modified_table <- batch_results[[plate_name]]$result$modified_ratio_table
      
      # Validate table exists and has data
      if (is.null(modified_table) || nrow(modified_table) == 0 || ncol(modified_table) == 0) {
        warning("No modified_ratio_table found for ", plate_name, " - skipping")
        next
      }
      
      if (verbose) {
        message("  Data dimensions: ", nrow(modified_table), " rows x ", ncol(modified_table), " columns")
        message("  Column names: ", paste(colnames(modified_table), collapse = ", "))
      }
      
      # Prepare data for DRC analysis
      # First column should be log(inhibitor), remaining columns are responses
      drc_data <- modified_table
      
      # Run dose-response analysis
      plate_drc_result <- fit_drc_3pl(
        data = drc_data,
        output_file = if(generate_reports) file.path(output_dir, paste0("drc_results_", plate_name, ".xlsx")) else NULL,
        normalize = normalize,
        verbose = verbose,
        enforce_bottom_threshold = enforce_bottom_threshold,
        bottom_threshold = bottom_threshold,
        r_sqr_threshold = r_sqr_threshold
      )
      
      # Store results
      drc_results[[plate_name]] <- list(
        plate_info = batch_results[[plate_name]][c("data_file", "info_sheet", "sheet_number")],
        drc_result = plate_drc_result,
        timestamp = Sys.time()
      )
      
      if (verbose) message("+ Successfully processed DRC for ", plate_name)
      
    }, error = function(e) {
      warning("Failed to process DRC for ", plate_name, ": ", e$message)
      if (verbose) message("X Failed to process DRC for ", plate_name)
    })
  }
  
  # Generate consolidated report if requested
  if (length(drc_results) > 0 && generate_reports) {
    generate_drc_batch_report(drc_results, output_dir, verbose)
  } else if (length(drc_results) > 0) {
    if (verbose) message("Skipping DRC report generation as requested")
  } else {
    warning("No plates were successfully processed for DRC analysis")
  }
  
  # Final statistics
  if (verbose) {
    successful_plates <- length(drc_results)
    total_compounds <- sum(sapply(drc_results, function(x) x$drc_result$n_compounds))
    successful_fits <- sum(sapply(drc_results, function(x) x$drc_result$successful_fits))
    
    message("\n" , strrep("=", 60))
    message("BATCH DRC ANALYSIS COMPLETED")
    message(strrep("=", 60))
    message("Plates processed: ", successful_plates, "/", length(batch_results))
    message("Total compounds analyzed: ", total_compounds)
    message("Successful fits: ", successful_fits, " (", 
            round(successful_fits/total_compounds * 100, 1), "%)")
    message("Reports generated: ", generate_reports)
    message(strrep("=", 60))
  }
  
  return(drc_results)
}