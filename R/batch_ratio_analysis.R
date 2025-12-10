#' Batch Processing of Ratio-Based Dose-Response Analyses with Version Control
#'
#' @description
#' `batch_ratio_analysis()` automates the batch processing of multiple
#' ratio-based dose–response experiments with support for two different
#' processing algorithms. It scans an info table file (typically `info_tables.xlsx`)
#' containing one sheet per plate/experiment, identifies sheets whose names end
#' with a numeric suffix, and automatically pairs each sheet with the corresponding
#' raw data file.
#'
#' **Key Features:**
#' - **Two processing algorithms**: Choose between original (`v1`) and enhanced (`v2`)
#' - **Automatic file matching**: Pairs info sheets with data files using numeric IDs
#' - **Comprehensive reporting**: Generates detailed Excel reports for each plate
#' - **Batch summary**: Creates consolidated report with quality metrics
#'
#' **Important:**
#' The numeric suffix of the sheet name **must match** the numeric identifier in
#' the raw data filename. Example:
#'
#' • Sheet name: `"Plate_01"` → numeric ID: `"01"`
#' • Raw data file: `"data_01.xlsx"`, `"experiment_01.xlsx"`
#'
#' This numeric matching is required for the batch processing to correctly link
#' each info sheet with its corresponding data file.
#'
#' @param directory Path to the folder containing all raw data files and the
#'   info table file. Defaults to the current working directory.
#' @param control_0perc Value or identifier of the **0% control** (minimum).
#'   Interpretation depends on `function_version`:
#'   - **For `function_version = "v1"`**: Must be a single column name (e.g., `"DMSO"`)
#'   - **For `function_version = "v2"`**: Can be either:
#'     * A single numeric value (e.g., `16`) for fixed-value 0% control
#'     * A single column name (e.g., `"DMSO_Column"`) for column-based control
#' @param control_100perc Value or identifier of the **100% control** (maximum).
#'   Interpretation depends on `function_version`:
#'   - **For `function_version = "v1"`**: Must be a single column name (e.g., `"Staurosporine"`)
#'   - **For `function_version = "v2"`**: Can be either:
#'     * A numeric vector of column positions (e.g., `c(12, 24)`) for multiple controls
#'     * A character vector of column names (e.g., `c("Ctrl_A", "Ctrl_B")`)
#' @param split_replicates Logical. If `TRUE` (default), technical replicates are
#'   separated into individual columns with `.2` suffix.
#' @param low_value_threshold Minimum acceptable signal value for filtering
#'   low-quality data. Values below this threshold are replaced with `NA`.
#'   Defaults to `1000`.
#' @param info_file Name of the Excel file containing info sheets for each plate.
#'   Defaults to `"info_tables.xlsx"`.
#' @param data_pattern Regular expression used to find raw data files. The default
#'   (`"_\\d+\\.xlsx$"`) matches filenames such as `"data_01.xlsx"` and `"plate_12.xlsx"`.
#' @param output_dir Directory where processed result files and consolidated
#'   reports will be saved. If `NULL`, defaults to `directory`.
#' @param generate_reports Logical. If `TRUE` (default):
#'   - Creates individual Excel files for each plate: `results_<number>.xlsx`
#'   - Creates consolidated report: `batch_analysis_report.xlsx`
#' @param function_version Character specifying which processing algorithm to use:
#'   - `"v1"` (default): Original `ratio_dose_response()` function
#'   - `"v2"`: Enhanced `ratio_dose_response_v2()` function with flexible controls
#' @param verbose Logical. If `TRUE` (default), prints detailed progress messages.
#'
#' @details
#' ## Function Versions
#'
#' ### Version 1 (`"v1"` - Default)
#' Uses the original `ratio_dose_response()` function with traditional control specification:
#' - `control_0perc`: Single column name
#' - `control_100perc`: Single column name
#' - **Best for**: Legacy data and compatibility with previous analyses
#'
#' ### Version 2 (`"v2"`)
#' Uses the enhanced `ratio_dose_response_v2()` function with new features:
#' - **Fixed-value controls**: Specify 0% control as a numeric value (e.g., `16`)
#' - **Multiple 100% controls**: Average multiple control columns
#' - **Column positions**: Identify controls by column number
#' - **New control columns**: Creates `Fixed_0perc` and `Mean_100perc` columns
#' - **Better quality reporting**: Enhanced control information
#' - **Best for**: New experiments and advanced control strategies
#'
#' ## Processing Pipeline
#'
#' The function performs the following steps:
#'
#' 1. **Sheet Identification**: Loads info file and finds sheets ending with numbers
#' 2. **File Matching**: Extracts numeric IDs and matches with data files
#' 3. **Data Processing**: Executes selected version of ratio analysis
#' 4. **Result Saving**: Creates individual Excel files for each plate
#' 5. **Report Generation**: Optionally creates consolidated batch report
#'
#' ## File Matching Logic
#'
#' **This numeric matching is required** for the function to determine which
#' raw data file corresponds to which info sheet.
#' The matching follows this pattern:
#' ```
#' Sheet Name → Extracted Number → Data File Pattern
#' "Plate_01" → "01" → "*_01.xlsx"
#' "Exp_07"   → "07" → "*_07.xlsx"
#' "15"       → "15" → "*_15.xlsx"
#' ```
#'
#' ## Output Files
#'
#' When `generate_reports = TRUE`, the following files are created:
#'
#' ### Individual Plate Results (`results_<number>.xlsx`)
#' Each file contains these sheets:
#' - `Modified_Ratio_Table`: Processed data ready for curve fitting
#' - `Original_Ratio_Table`: Raw calculated ratios
#' - `General_Means`: Overall control means (when available)
#' - `Interval_Means`: Construct-specific quality metrics with row names preserved
#' - `Control_Info`: Control processing details (v2 only)
#' - `Processing_Info`: Metadata about the analysis
#'
#' ### Consolidated Report (`batch_analysis_report.xlsx`)
#' Contains:
#' - **Summary sheet**: Overview of all processed plates
#' - **Individual quality sheets**: `interval_means` for each plate
#'
#' ## Quality Assessment
#'
#' When info tables are provided, the function calculates:
#' - **Luciferase Signal**: Signal intensity categories (high/medium/low/insufficient)
#' - **Z-Score**: Assay robustness metric
#' - **Assay Window**: Dynamic range calculation
#' - **Overall Quality**: Composite assessment
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
#'   \item{function_version}{Which processing algorithm was used (`"v1"` or `"v2"`).}
#'   \item{control_0perc}{Control 0% specification as provided.}
#'   \item{control_100perc}{Control 100% specification as provided.}
#'   \item{result}{Full output of the selected ratio analysis function.}
#' }
#'
#' If `generate_reports = TRUE`, two types of Excel files are created:
#' 1. Individual plate results: `results_<number>.xlsx`
#' 2. Consolidated batch report: `batch_analysis_report.xlsx`
#'
#' @examples
#' \dontrun{
#' # ===================================================================
#' # EXAMPLE 1: Default settings (v1 with column names)
#' # ===================================================================
#' batch_ratio_analysis(
#'   directory = "experiment_data",
#'   control_0perc = "DMSO_Control",
#'   control_100perc = "Staurosporine_100pct",
#'   # function_version = "v1" is default
#' )
#'
#' # ===================================================================
#' # EXAMPLE 2: Enhanced v2 with fixed-value controls
#' # ===================================================================
#' batch_ratio_analysis(
#'   directory = "new_experiment_data",
#'   control_0perc = 16,                # Fixed value for 0% control
#'   control_100perc = c(12, 24),       # Column positions for duplicate 100% controls
#'   function_version = "v2",
#'   output_dir = "processed_v2_results",
#'   verbose = TRUE
#' )
#'
#' # ===================================================================
#' # EXAMPLE 3: v2 with column names (backward compatible)
#' # ===================================================================
#' batch_ratio_analysis(
#'   directory = "mixed_data",
#'   control_0perc = "DMSO_Column",
#'   control_100perc = "Stauro_Column",
#'   function_version = "v2",          # v2 works with column names too
#'   split_replicates = FALSE,
#'   generate_reports = TRUE
#' )
#'
#' # ===================================================================
#' # EXAMPLE 4: Minimal usage (current directory, default everything)
#' # ===================================================================
#' # Assumes: info_tables.xlsx and data_01.xlsx, data_02.xlsx, etc. exist
#' results <- batch_ratio_analysis()
#'
#' # Access results for first plate
#' plate1 <- results[["Plate_01"]]$result
#' analysis_ready_data <- plate1$modified_ratio_table
#' quality_metrics <- plate1$interval_means
#' }
#'
#' @section Version Selection Guide:
#'
#' **When to use `v1` (default):**
#' - Analyzing legacy data from previous experiments
#' - Need compatibility with existing analysis pipelines
#' - Using traditional column-based controls
#' - Preferring the original, well-tested algorithm
#'
#' **When to use `v2`:**
#' - Starting new experiments with fixed 0% control values
#' - Have duplicate 100% control columns that need averaging
#' - Want more detailed control information in reports
#' - Need to specify controls by column position instead of name
#' - Require enhanced quality reporting
#'
#' @section Troubleshooting:
#'
#' **Common issues and solutions:**
#'
#' 1. **"No data file found for sheet..."**
#'    - Check that sheet names end with numbers (e.g., `"Plate_01"`, not `"Plate_One"`)
#'    - Verify data files contain matching numbers (e.g., `"_01.xlsx"`)
#'    - Ensure `data_pattern` matches your filename convention
#'
#' 2. **Control specification errors**
#'    - For `v1`: Both controls must be single column names
#'    - For `v2`: Check if using correct type (numeric vs character)
#'    - Verify control columns exist in your data
#'
#' 3. **Excel file issues**
#'    - Ensure `openxlsx` package is installed
#'    - Close Excel files before running analysis
#'    - Check file permissions in output directory
#'
#' @seealso
#' * [`ratio_dose_response()`] – Original processing function (v1)
#' * [`ratio_dose_response_v2()`] – Enhanced processing function (v2)
#' * `openxlsx` – Package used for reading and writing Excel files
#' * [`fit_drc_3pl()`] – For subsequent dose-response curve fitting
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
                                 function_version = "v1",
                                 verbose = TRUE) {

  # Validate function_version argument
  valid_versions <- c("v1", "v2")
  if (!function_version %in% valid_versions) {
    stop("function_version must be either 'v1' or 'v2'")
  }

  # Show informative message about version differences
  if (verbose) {
    message("=== BATCH ANALYSIS SETTINGS ===")
    message("Using function version: ", function_version)
    if (function_version == "v1") {
      message("  • Original ratio_dose_response function")
      message("  • control_0perc: must be column name (e.g., 'DMSO_Control')")
      message("  • control_100perc: must be single column name (e.g., 'Stauro_Control')")
    } else {
      message("  • Enhanced ratio_dose_response_v2 function")
      message("  • control_0perc: can be fixed value (e.g., 16) OR column name")
      message("  • control_100perc: can be column positions (e.g., c(12, 24)) OR column names")
    }
    message("================================")
  }

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
      Function_Version = character(),
      Status = character(),
      N_Constructs = integer(),
      Overall_Quality = character(),
      Control_Structure = character(),
      stringsAsFactors = FALSE
    )

    # Collect summary data from all plates
    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result

      n_constructs <- 0
      overall_quality <- "Not assessed"
      control_structure <- "Unknown"

      # Get control structure information
      if (!is.null(plate_result$control_info)) {
        if ("new_columns_created" %in% names(plate_result$control_info) &&
            !is.null(plate_result$control_info$new_columns_created)) {
          control_structure <- paste(plate_result$control_info$new_columns_created, collapse = ", ")
        } else {
          control_structure <- "Original columns"
        }
      }

      if (!is.null(plate_result$interval_means)) {
        n_constructs <- ncol(plate_result$interval_means)
        if ("Overall_Quality" %in% rownames(plate_result$interval_means)) {
          qualities <- as.character(plate_result$interval_means["Overall_Quality", ])
          overall_quality <- paste(unique(qualities), collapse = ", ")
        }
      }

      summary_data <- rbind(summary_data, data.frame(
        Sheet_Name = plate_sheet,
        Sheet_Number = results[[plate_sheet]]$sheet_number,
        Data_File = results[[plate_sheet]]$data_file,
        Function_Version = results[[plate_sheet]]$function_version,
        Status = "Completed",
        N_Constructs = n_constructs,
        Overall_Quality = overall_quality,
        Control_Structure = control_structure,
        stringsAsFactors = FALSE
      ))
    }

    openxlsx::writeData(wb, "Summary", summary_data)

    # Style header row
    header_style <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD",
                                          halign = "center", textDecoration = "Bold",
                                          border = "TopBottom", borderColour = "#4F81BD")
    openxlsx::addStyle(wb, "Summary", header_style, rows = 1, cols = 1:ncol(summary_data))

    # Add individual plate quality sheets
    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result

      if (!is.null(plate_result$interval_means)) {
        openxlsx::addWorksheet(wb, plate_sheet)

        # Add plate info as first rows
        info_text <- c(
          paste("Data file:", results[[plate_sheet]]$data_file),
          paste("Info sheet:", results[[plate_sheet]]$info_sheet),
          paste("Function version:", results[[plate_sheet]]$function_version),
          paste("Control 0%:", ifelse(is.numeric(results[[plate_sheet]]$control_0perc),
                                      paste("Fixed value", results[[plate_sheet]]$control_0perc),
                                      results[[plate_sheet]]$control_0perc)),
          paste("Control 100%:", ifelse(is.numeric(results[[plate_sheet]]$control_100perc),
                                        paste("Positions", paste(results[[plate_sheet]]$control_100perc, collapse = ", ")),
                                        paste(results[[plate_sheet]]$control_100perc, collapse = ", "))),
          ""
        )

        openxlsx::writeData(wb, plate_sheet, info_text, startRow = 1)

        quality_data <- as.data.frame(t(plate_result$interval_means))
        quality_data <- cbind(Construct = rownames(quality_data), quality_data)
        rownames(quality_data) <- NULL

        openxlsx::writeData(wb, plate_sheet, quality_data, startRow = length(info_text) + 1)

        # Style header for quality data
        openxlsx::addStyle(wb, plate_sheet, header_style,
                           rows = length(info_text) + 1,
                           cols = 1:ncol(quality_data))
      }
    }

    openxlsx::saveWorkbook(wb, report_path, overwrite = TRUE)

    if (verbose) {
      message("\n=== BATCH REPORT SUMMARY ===")
      message("Report saved: ", report_path)
      message("Successfully processed ", length(results), " plates")
      message("\nFunction versions used:")
      versions_used <- unique(sapply(results, function(x) x$function_version))
      for (ver in versions_used) {
        n_ver <- sum(sapply(results, function(x) x$function_version == ver))
        message("  • ", ver, ": ", n_ver, " plate(s)")
      }
      message("=============================\n")
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

      # Select function based on version
      if (function_version == "v1") {
        # Use original function
        result <- ratio_dose_response(
          data = data,
          control_0perc = control_0perc,
          control_100perc = control_100perc,
          split_replicates = split_replicates,
          info_table = info_table,
          low_value_threshold = low_value_threshold,
          verbose = verbose,
          save_to_excel = NULL
        )
      } else {
        # Use enhanced v2 function
        result <- ratio_dose_response_v2(
          data = data,
          control_0perc = control_0perc,
          control_100perc = control_100perc,
          split_replicates = split_replicates,
          info_table = info_table,
          low_value_threshold = low_value_threshold,
          verbose = verbose,
          save_to_excel = NULL
        )
      }

      # --- SAVE INDIVIDUAL PLATE RESULTS TO EXCEL ---
      if (generate_reports) {
        excel_path <- file.path(output_dir, paste0("results_", sheet_number, ".xlsx"))

        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          stop("Package 'openxlsx' is required to save Excel files.")
        }

        wb <- openxlsx::createWorkbook()

        # 1. Modified Ratio Table sheet
        openxlsx::addWorksheet(wb, "Modified_Ratio_Table")
        modified_data <- cbind(RowNames = rownames(result$modified_ratio_table),
                               result$modified_ratio_table)
        openxlsx::writeData(wb, "Modified_Ratio_Table", modified_data, rowNames = FALSE)

        # 2. Original Ratio Table sheet
        openxlsx::addWorksheet(wb, "Original_Ratio_Table")
        original_data <- cbind(RowNames = rownames(result$original_ratio_table),
                               result$original_ratio_table)
        openxlsx::writeData(wb, "Original_Ratio_Table", original_data, rowNames = FALSE)

        # 3. General Means sheet (if exists)
        if (!is.null(result$general_means)) {
          openxlsx::addWorksheet(wb, "General_Means")
          openxlsx::writeData(wb, "General_Means", result$general_means)
        }

        # 4. Interval Means sheet
        if (!is.null(result$interval_means)) {
          openxlsx::addWorksheet(wb, "Interval_Means")

          interval_df <- as.data.frame(result$interval_means)

          interval_with_rownames <- cbind(
            Metric = rownames(interval_df),
            interval_df
          )

          rownames(interval_with_rownames) <- NULL

          # Writing Excel
          openxlsx::writeData(wb, "Interval_Means", interval_with_rownames, rowNames = FALSE)

          header_style <- openxlsx::createStyle(
            fontColour = "#FFFFFF",
            fgFill = "#4F81BD",
            halign = "center",
            textDecoration = "Bold",
            border = "TopBottom",
            borderColour = "#4F81BD"
          )

          openxlsx::addStyle(wb, "Interval_Means", header_style,
                             rows = 1, cols = 1:ncol(interval_with_rownames))

          alt_row_style <- openxlsx::createStyle(
            fgFill = "#F2F2F2"
          )

          n_rows <- nrow(interval_with_rownames)
          if (n_rows > 1) {
            for (row in seq(2, n_rows, by = 2)) {
              openxlsx::addStyle(wb, "Interval_Means", alt_row_style,
                                 rows = row, cols = 1:ncol(interval_with_rownames))
            }
          }
        }

        # 5. Control Info sheet (if exists in v2)
        if (!is.null(result$control_info)) {
          openxlsx::addWorksheet(wb, "Control_Info")

          control_df <- data.frame(
            Parameter = names(result$control_info),
            Value = sapply(result$control_info, function(x) {
              if (is.null(x)) return("NULL")
              if (length(x) == 0) return("Empty")
              if (is.character(x)) return(paste(x, collapse = ", "))
              if (is.numeric(x)) return(paste(x, collapse = ", "))
              if (is.logical(x)) return(as.character(x))
              return(as.character(x))
            }),
            stringsAsFactors = FALSE
          )

          openxlsx::writeData(wb, "Control_Info", control_df)

          openxlsx::addStyle(wb, "Control_Info", header_style, rows = 1, cols = 1:2)
        }

        # 6. Processing Info sheet
        openxlsx::addWorksheet(wb, "Processing_Info")

        processing_info <- data.frame(
          Parameter = c(
            "Data File",
            "Info Sheet",
            "Function Version",
            "Control 0%",
            "Control 100%",
            "Split Replicates",
            "Low Value Threshold",
            "Processing Date"
          ),
          Value = c(
            data_filename,
            info_sheet,
            function_version,
            ifelse(is.numeric(control_0perc),
                   paste("Fixed value:", control_0perc),
                   ifelse(is.null(control_0perc), "Not specified", control_0perc)),
            ifelse(is.numeric(control_100perc),
                   paste("Positions:", paste(control_100perc, collapse = ", ")),
                   ifelse(is.null(control_100perc), "Not specified",
                          paste(control_100perc, collapse = ", "))),
            split_replicates,
            low_value_threshold,
            as.character(Sys.time())
          ),
          stringsAsFactors = FALSE
        )

        openxlsx::writeData(wb, "Processing_Info", processing_info)
        openxlsx::addStyle(wb, "Processing_Info", header_style, rows = 1, cols = 1:2)

        # Salvar arquivo Excel
        openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)

        if (verbose) {
          message("  • Excel file saved: ", basename(excel_path))
        }
      }

      # Store results with additional metadata
      results[[info_sheet]] <- list(
        data_file = data_filename,
        info_sheet = info_sheet,
        sheet_number = sheet_number,
        function_version = function_version,
        control_0perc = control_0perc,
        control_100perc = control_100perc,
        result = result
      )

      if (verbose) message("✓ Successfully processed sheet ", info_sheet, " with ", function_version)

    }, error = function(e) {
      warning("Failed to process sheet ", info_sheet, " (", data_filename, ") with ", function_version, ": ", e$message)
      if (verbose) message("✗ Failed to process sheet ", info_sheet)
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
