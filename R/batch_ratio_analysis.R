#' Batch Processing of Ratio-Based Dose-Response Analyses
#'
#' @description
#' `batch_ratio_analysis()` automates the batch processing of multiple
#' ratio-based dose–response experiments.  
#' It scans an info table file (typically `info_tables.xlsx`) containing one
#' sheet per plate/experiment, identifies sheets whose names end with a numeric
#' suffix (e.g., `"Plate_01"`, `"Exp_03"`, `"05"`), and automatically pairs each
#' sheet with the corresponding raw data file.  
#'
#' **Important:**  
#' The numeric suffix of the sheet name **must match** the numeric identifier in
#' the raw data filename. Example:
#'
#' • Sheet name: `"Plate_01"` → numeric ID: `"01"`  
#' • Raw data file: `"data_01.xlsx"`  
#'
#' This numeric matching is required for the batch processing to correctly link
#' each info sheet with its corresponding data file.
#'
#' The function runs [`ratio_dose_response()`] for each matched pair, saves
#' individual results, and optionally generates a consolidated Excel report.
#'
#'
#' @param directory Path to the folder containing all raw data files and the
#'   info table file. Defaults to the current working directory.
#' @param control_0perc Value or identifier of the **0% control** (minimum).
#'   Passed directly to `ratio_dose_response()`.
#' @param control_100perc Value or identifier of the **100% control** (maximum).
#'   Passed directly to `ratio_dose_response()`.
#' @param split_replicates Logical. If `TRUE` (default), replicates are processed
#'   individually inside `ratio_dose_response()`.
#' @param low_value_threshold Minimum acceptable signal value for filtering
#'   low-quality data.
#' @param info_file Name of the Excel file containing info sheets for each plate.
#'   Defaults to `"info_tables.xlsx"`.
#' @param data_pattern Regular expression used to find raw data files. The default
#'   (`"_\\d+\\.xlsx$"`) matches filenames such as `"data_01.xlsx"` and `"plate_12.xlsx"`.
#' @param output_dir Directory where processed result files and consolidated
#'   reports will be saved. If `NULL`, defaults to `directory`.
#' @param generate_reports If `TRUE` (default), creates a consolidated
#'   `batch_analysis_report.xlsx` file containing a summary and quality metrics
#'   for all processed plates.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Loads the info table file and identifies all sheet names ending with a
#'    numeric suffix.  
#' 2. Extracts the numeric portion from each sheet name (e.g., `"Plate_03"` → `"03"`).  
#' 3. Searches for a raw data file in `directory` whose filename contains the
#'    *same* numeric identifier.  
#'
#'    **This numeric matching is required** for the function to determine which
#'    raw data file corresponds to which info sheet.  
#'    If the sheet is `"Exp_07"`, the data file must contain `"_07.xlsx"`.  
#'
#' 4. Reads the raw data and the corresponding info table.
#' 5. Executes [`ratio_dose_response()`] for each matched pair.
#' 6. Saves per-plate results into Excel files (`results_<number>.xlsx`) unless
#'    `generate_reports = FALSE`.
#' 7. Optionally generates a consolidated report containing:
#'    - A "Summary" sheet with plate status, number of targets, data file used,
#'      and overall quality classification.  
#'    - Individual sheets containing interval mean tables (`interval_means`) for
#'      each processed plate.
#'
#'
#' @return
#' Returns a named list where each element corresponds to a successfully processed
#' plate/experiment.  
#'
#' Each entry contains:
#' \describe{
#'   \item{data_file}{Raw data filename associated with the sheet.}
#'   \item{info_sheet}{Name of the info sheet used.}
#'   \item{sheet_number}{Numeric identifier extracted from the sheet name.}
#'   \item{result}{Full output of `ratio_dose_response()`.}
#' }
#'
#' If `generate_reports = TRUE`, an Excel file  
#' `"batch_analysis_report.xlsx"`  
#' is also created in `output_dir`.
#'
#'
#' @examples
#' \dontrun{
#'
#' # Basic usage (using current working directory):
#' batch_ratio_analysis()
#'
#' # Running with specified controls and output directory:
#' batch_ratio_analysis(
#'   directory = "experiment_data",
#'   control_0perc = "DMSO",
#'   control_100perc = "Positive_Control",
#'   output_dir = "processed_results"
#' )
#'
#' # Running without generating the consolidated report:
#' batch_ratio_analysis(
#'   directory = "raw_data",
#'   generate_reports = FALSE
#' )
#'
#' }
#'
#'
#' @seealso
#' * [`ratio_dose_response()`] – The function executed for each matched plate.  
#' * `openxlsx` – Package used for reading and writing Excel files.
#'
#'
#' @export




batch_ratio_analysis <- function(directory = getwd(),
                                 control_0perc = NULL,
                                 control_100perc = NULL,
                                 split_replicates = TRUE,
                                 low_value_threshold = 1000,
                                 info_file = "info_tables.xlsx",
                                 data_pattern = "_\\d+\\.xlsx$",
                                 output_dir = NULL,
                                 generate_reports = TRUE,
                                 verbose = TRUE) {
  
  # Internal function to generate batch report
  generate_batch_report <- function(results, directory, verbose = TRUE) {
    
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required.")
    }
    
    report_path <- file.path(directory, "batch_analysis_report.xlsx")
    wb <- openxlsx::createWorkbook()
    
    # Summary sheet
    openxlsx::addWorksheet(wb, "Summary")
    
    summary_data <- data.frame(
      Sheet_Name = character(),
      Sheet_Number = character(),
      Data_File = character(),
      Status = character(),
      N_Targets = integer(),
      Overall_Quality = character(),
      stringsAsFactors = FALSE
    )
    
    # Collect summary data from all plates
    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result
      
      n_targets <- 0
      overall_quality <- "Not assessed"
      
      if (!is.null(plate_result$interval_means)) {
        n_targets <- ncol(plate_result$interval_means)
        if ("Overall_Quality" %in% rownames(plate_result$interval_means)) {
          qualities <- as.character(plate_result$interval_means["Overall_Quality", ])
          overall_quality <- paste(unique(qualities), collapse = ", ")
        }
      }
      
      summary_data <- rbind(summary_data, data.frame(
        Sheet_Name = plate_sheet,
        Sheet_Number = results[[plate_sheet]]$sheet_number,
        Data_File = results[[plate_sheet]]$data_file,
        Status = "Completed",
        N_Targets = n_targets,
        Overall_Quality = overall_quality,
        stringsAsFactors = FALSE
      ))
    }
    
    openxlsx::writeData(wb, "Summary", summary_data)
    
    # Add individual plate quality sheets
    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result
      
      if (!is.null(plate_result$interval_means)) {
        openxlsx::addWorksheet(wb, plate_sheet)
        
        quality_data <- as.data.frame(t(plate_result$interval_means))
        quality_data <- cbind(Target = rownames(quality_data), quality_data)
        rownames(quality_data) <- NULL
        
        openxlsx::writeData(wb, plate_sheet, quality_data)
      }
    }
    
    openxlsx::saveWorkbook(wb, report_path, overwrite = TRUE)
    
    if (verbose) {
      message("Batch report saved: ", report_path)
      message("Successfully processed ", length(results), " plates")
    }
  }
  
  # ========== MAIN FUNCTION BODY ==========
  
  # Use current directory if not specified
  if (directory == getwd()) {
    if (verbose) message("Using current working directory: ", directory)
  }
  
  # Create output directory if specified
  if (is.null(output_dir)) {
    output_dir <- directory
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) message("Created output directory: ", output_dir)
    }
  }
  
  # Validate directory exists
  if (!dir.exists(directory)) {
    stop("Directory not found: ", directory)
  }
  
  # Check if info file exists
  info_path <- file.path(directory, info_file)
  if (!file.exists(info_path)) {
    stop("Info file not found: ", info_path)
  }
  
  # Load available sheets from info file
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Please install it.")
  }
  
  info_sheets <- openxlsx::getSheetNames(info_path)
  
  # Find sheets that end with numbers (e.g., "plate_01", "exp_01", "01")
  number_sheets <- info_sheets[grepl("\\d+$", info_sheets)]
  
  if (length(number_sheets) == 0) {
    stop("No valid number sheets found in info file. Sheets should end with numbers (ex: 'plate_01', 'exp_01', '01')")
  }
  
  if (verbose) {
    message("Found ", length(number_sheets), " number sheets: ", 
            paste(number_sheets, collapse = ", "))
    message("Generate reports: ", generate_reports)
  }
  
  # Find data files matching the pattern
  data_files <- list.files(directory, pattern = data_pattern, full.names = FALSE)
  data_files <- data_files[!grepl("^~\\$", data_files)] # Exclude Excel temp files
  
  if (verbose) {
    message("Found ", length(data_files), " data files matching pattern")
    message("Output directory: ", output_dir)
  }
  
  results <- list()
  
  # Process each info sheet
  for (info_sheet in number_sheets) {
    
    # Extract number from sheet name (ignoring any prefix)
    sheet_number <- gsub("^.*?(\\d+)$", "\\1", info_sheet)
    
    # Find matching data files for this number
    pattern <- paste0("_", sheet_number, "\\.xlsx$")
    matching_files <- data_files[grepl(pattern, data_files)]
    
    if (length(matching_files) == 0) {
      if (verbose) message("No data file found for sheet ", info_sheet, " (number ", sheet_number, "), skipping")
      next
    }
    
    # Use first matching file (should be only one)
    data_filename <- matching_files[1]
    
    if (length(matching_files) > 1) {
      warning("Multiple files match sheet ", info_sheet, ". Using: ", data_filename)
    }
    
    data_path <- file.path(directory, data_filename)
    
    if (verbose) message("Processing sheet '", info_sheet, "' (number ", sheet_number, ") with data: ", data_filename)
    
    tryCatch({
      # Load data and info table
      data <- openxlsx::read.xlsx(data_path)
      info_table <- openxlsx::read.xlsx(info_path, sheet = info_sheet)
      
      # Run ratio dose-response analysis
      result <- ratio_dose_response(
        data = data,
        control_0perc = control_0perc,
        control_100perc = control_100perc,
        split_replicates = split_replicates,
        info_table = info_table,
        low_value_threshold = low_value_threshold,
        verbose = verbose,
        save_to_excel = if(generate_reports) file.path(output_dir, paste0("results_", sheet_number, ".xlsx")) else NULL
      )
      
      # Store results
      results[[info_sheet]] <- list(
        data_file = data_filename,
        info_sheet = info_sheet,
        sheet_number = sheet_number,
        result = result
      )
      
      if (verbose) message("+ Successfully processed sheet ", info_sheet)
      
    }, error = function(e) {
      warning("Failed to process sheet ", info_sheet, " (", data_filename, "): ", e$message)
      if (verbose) message("X Failed to process sheet ", info_sheet)
    })
  }
  
  # Generate consolidated report if requested
  if (length(results) > 0 && generate_reports) {
    generate_batch_report(results, output_dir, verbose)
  } else if (length(results) > 0) {
    if (verbose) message("Skipping report generation as requested")
  } else {
    warning("No plates were successfully processed")
  }
  
  return(results)
}
