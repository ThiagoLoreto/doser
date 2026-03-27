#' Batch NanoBRET Ratio Dose-Response Analysis
#'
#' Processes an entire directory of raw BMG PHERAstar Excel files in one
#' call: reads each plate, computes BRET ratios via
#' \code{\link{ratio_dose_response_v2}} (or the legacy
#' \code{ratio_dose_response}), and optionally writes per-plate and
#' summary Excel reports.
#'
#' @description
#' This is the \strong{top-level entry point} for the NanoBRET processing
#' pipeline.  It pairs each numbered sheet in the info file with the
#' matching raw data file (matched by the trailing plate number in the
#' filename), calls the chosen ratio function, and collects results into a
#' named list suitable for passing directly to
#' \code{\link{rout_outliers_batch}}.
#'
#' @section File matching:
#' The info file (default \code{info_tables.xlsx}) must contain one sheet
#' per plate, named with a trailing integer (e.g. \code{Sheet1},
#' \code{Plate_01}, \code{2}).  Each sheet is matched to a raw data file
#' whose name ends with \code{_<number>.xlsx} (controlled by
#' \code{data_pattern}).  If multiple files match, the first is used and a
#' warning is issued.
#'
#' @section Raw file parsing:
#' Raw files are parsed by the internal helper \code{.read_nanobret_raw()},
#' which detects the two emission tables by content (searching for rows
#' matching \code{"\\d+\\.\\s*Raw Data"} in column 2) rather than by
#' fixed row position.  This makes the function robust to files that
#' contain an extra \code{ID1:} metadata row (produced when more than one
#' plate is measured in the same run), as well as any other variation in
#' the number of header rows above the tables.  After detection, the
#' metadata block is normalised so that the column-number header always
#' lands at row 9 -- the position expected by both
#' \code{ratio_dose_response} and \code{\link{ratio_dose_response_v2}}.
#' Plate row labels A--P are validated before the data is passed downstream.
#'
#' @param directory Character string.  Directory containing the raw data
#'   files and the info file (default: current working directory).
#'
#' @param control_0perc 0\% control specification passed to the ratio
#'   function.  For \code{function_version = "v1"}: a column name string.
#'   For \code{function_version = "v2"}: either a single numeric value
#'   (fixed-value mode) or a column name string.
#'
#' @param control_100perc 100\% control specification passed to the ratio
#'   function.  For \code{function_version = "v1"}: a column name string.
#'   For \code{function_version = "v2"}: a numeric vector of column
#'   positions or a column name string.
#'
#' @param split_replicates Logical.  Split technical replicates into
#'   separate \code{.2} columns (default \code{TRUE}).  Passed to the
#'   ratio function.
#'
#' @param low_value_threshold Numeric.  Donor-channel values below this
#'   threshold are replaced with \code{NA} before ratio calculation
#'   (default \code{1000}).  Passed to the ratio function.
#'
#' @param info_file Character string.  Name of the info Excel file
#'   (default \code{"info_tables.xlsx"}).  Must reside in \code{directory}.
#'
#' @param data_pattern Regular expression used to identify raw data files
#'   in \code{directory} (default \code{"_\\\\d+\\.xlsx$"}).  Files
#'   whose names begin with \code{~$} (temporary Excel lock files) are
#'   always excluded.
#'
#' @param output_dir Character string.  Directory for output files.
#'   Defaults to \code{directory}.  Created automatically if it does not
#'   exist.
#'
#' @param generate_reports Logical.  Save per-plate Excel files
#'   (\code{results_<N>.xlsx}) and a consolidated
#'   \code{batch_analysis_report.xlsx} (default \code{TRUE}).
#'
#' @param function_version Character.  \code{"v1"} to call the legacy
#'   \code{ratio_dose_response()} or \code{"v2"} (default) to call
#'   \code{\link{ratio_dose_response_v2}}.
#'
#' @param verbose Logical.  Print progress messages (default \code{TRUE}).
#'
#' @param selected_columns Integer vector of 1-based data column indices to
#'   pass to the ratio function (\code{NULL} uses all 24 columns).
#'
#' @return A named list with one element per successfully processed plate.
#'   The element name is the info-sheet name (e.g. \code{"Sheet1"}).
#'   Each element is a list with:
#'   \describe{
#'     \item{\code{data_file}}{Filename of the raw data file used.}
#'     \item{\code{info_sheet}}{Name of the info sheet used.}
#'     \item{\code{sheet_number}}{Plate number extracted from the sheet
#'       name.}
#'     \item{\code{function_version}}{\code{"v1"} or \code{"v2"}.}
#'     \item{\code{control_0perc}}{Value passed as \code{control_0perc}.}
#'     \item{\code{control_100perc}}{Value passed as
#'       \code{control_100perc}.}
#'     \item{\code{selected_columns}}{Value passed as
#'       \code{selected_columns}.}
#'     \item{\code{result}}{The full return value of the ratio function
#'       (see \code{\link{ratio_dose_response_v2}} for the structure).}
#'   }
#'   Plates that fail processing are omitted from the list and a warning is
#'   issued.  If no plates succeed, a warning is issued and an empty list
#'   is returned.
#'
#' @section Output files (when generate_reports = TRUE):
#' \describe{
#'   \item{\code{results_<N>.xlsx}}{Per-plate workbook with sheets:
#'     \code{Modified_Ratio_Table}, \code{Original_Ratio_Table},
#'     \code{General_Means}, \code{Interval_Means},
#'     \code{Control_Info}, \code{Selected_Columns_Info},
#'     \code{Processing_Info}.}
#'   \item{\code{batch_analysis_report.xlsx}}{Consolidated workbook with a
#'     \code{Summary} sheet (one row per plate) and one quality-metrics
#'     sheet per plate.}
#' }
#'
#' @examples
#' \dontrun{
#' # --- Fixed-value mode (v2, recommended) ---
#' results <- batch_ratio_analysis(
#'   directory        = "data/experiment_01/",
#'   control_0perc    = 16,
#'   control_100perc  = c(12, 24),
#'   function_version = "v2",
#'   split_replicates = TRUE
#' )
#'
#' # Access the processed table for plate 1
#' results$Sheet1$result$modified_ratio_table
#'
#' # --- Legacy column-name mode (v1) ---
#' results_v1 <- batch_ratio_analysis(
#'   directory        = "data/legacy/",
#'   control_0perc    = "DMSO_Control",
#'   control_100perc  = "Staurosporine",
#'   function_version = "v1"
#' )
#'
#' # --- Pipe into outlier detection ---
#' results_clean <- rout_outliers_batch(results, Q = 0.01)
#' }
#'
#' @seealso
#' \code{\link{ratio_dose_response_v2}} for the single-plate ratio
#' function called internally.
#'
#' \code{\link{rout_outliers_batch}} for the next step in the
#' pipeline.
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
                                 verbose = TRUE,
                                 selected_columns = NULL) {

  # -- Validate function_version ----------------------------------------------
  valid_versions <- c("v1", "v2")
  if (!function_version %in% valid_versions)
    stop("function_version must be either 'v1' or 'v2'")

  # -- Internal helper: parse a raw BMG PHERAstar Excel file ------------------
  #
  # Both ratio_dose_response() (v1) and ratio_dose_response_v2() use hard-coded
  # row indices: data[9,] for colnames, data[10:25,] for table 1 (450-80 B),
  # data[28:43,] for table 2 (610-LP A).  These indices are only correct when
  # the column-number header row lands at position 9 of the data frame.
  #
  # When more than one plate is measured in the same run the reader inserts an
  # ID1: row in the metadata block, shifting everything down by one and breaking
  # all three indices.  This helper detects the tables by content, then
  # normalises the metadata block so the column-number header always lands at
  # row 9 - regardless of how many metadata rows the file contains.
  # It also validates that plate row labels A-P are in the expected positions.
  .read_nanobret_raw <- function(path, function_version = "v2",
                                 validate = TRUE, verbose = FALSE) {

    if (!file.exists(path))
      stop("File not found: ", path)

    # 1. Read entire sheet with no layout assumptions
    raw <- openxlsx::read.xlsx(path, sheet = 1L,
                                colNames      = FALSE,
                                skipEmptyRows = FALSE,
                                skipEmptyCols = FALSE)

    # 2. Locate table-title rows - "1. Raw Data (450-80 B)" etc. in column 2
    col2 <- as.character(raw[, 2L])
    title_rows <- grep("^\\d+\\.\\s*Raw Data", col2)

    if (length(title_rows) == 0L)
      stop("No 'Raw Data' table titles found in column 2 of: ", basename(path),
           "\nExpected rows like '1. Raw Data (450-80 B)'.")

    first_title_row <- title_rows[1L]
    col_hdr_row_raw <- first_title_row + 1L   # col-number header immediately follows

    if (verbose)
      message(sprintf(
        "  [read] found %d table(s) in '%s' at raw row(s) %s",
        length(title_rows), basename(path),
        paste(title_rows, collapse = ", ")))

    # 3. Collect non-empty metadata rows above the first title
    meta_raw    <- raw[seq_len(first_title_row - 1L), , drop = FALSE]
    is_nonempty <- apply(meta_raw, 1L, function(r) any(!is.na(r)))
    meta_ne     <- meta_raw[is_nonempty, , drop = FALSE]

    # Target: exactly 9 non-empty metadata rows so that after
    # colNames=TRUE / skipEmptyRows=TRUE the col-header lands at data row 9.
    target <- 9L
    n_meta  <- nrow(meta_ne)

    if (n_meta < target) {
      # Pad with non-NA placeholder rows (empty string in col 1) so they
      # survive the all-NA row filter applied below.
      n_pad <- target - n_meta
      pad   <- as.data.frame(matrix(NA_character_, nrow = n_pad, ncol = ncol(raw)),
                              stringsAsFactors = FALSE)
      colnames(pad) <- colnames(meta_ne)
      pad[, 1L] <- ""
      meta_ne <- rbind(meta_ne, pad)
    } else if (n_meta > target) {
      meta_ne <- meta_ne[seq_len(target), , drop = FALSE]
    }

    # 4. Assemble normalised sheet: trimmed metadata + col-header + data
    tail_raw      <- raw[col_hdr_row_raw:nrow(raw), , drop = FALSE]
    normalised    <- rbind(meta_ne, tail_raw)
    rownames(normalised) <- NULL

    # 5. Mimic read.xlsx(colNames=TRUE, skipEmptyRows=TRUE)
    colnames(normalised) <- as.character(normalised[1L, ])
    normalised <- normalised[-1L, , drop = FALSE]
    all_na     <- apply(normalised, 1L, function(r) all(is.na(r)))
    result     <- normalised[!all_na, , drop = FALSE]
    rownames(result) <- NULL

    # 6. Validate plate row labels A-P at the positions v1/v2 will read
    if (validate) {
      expected <- LETTERS[1:16]
      tmp      <- result
      colnames(tmp) <- as.character(tmp[9L, ])
      t1 <- as.character(tmp[10:25, 1L])
      t2 <- as.character(tmp[28:43, 1L])
      bad1 <- which(t1 != expected)
      bad2 <- which(t2 != expected)

      if (length(bad1) > 0L || length(bad2) > 0L) {
        msgs <- character(0L)
        if (length(bad1) > 0L)
          msgs <- c(msgs, sprintf(
            "Table 1 (450-80 B): position(s) %s have label(s) '%s' (expected '%s')",
            paste(bad1, collapse = ", "),
            paste(t1[bad1], collapse = ", "),
            paste(expected[bad1], collapse = ", ")))
        if (length(bad2) > 0L)
          msgs <- c(msgs, sprintf(
            "Table 2 (610-LP A): position(s) %s have label(s) '%s' (expected '%s')",
            paste(bad2, collapse = ", "),
            paste(t2[bad2], collapse = ", "),
            paste(expected[bad2], collapse = ", ")))
        fn_name <- if (function_version == "v1") "ratio_dose_response()"
                   else                          "ratio_dose_response_v2()"
        stop(
          "Row label validation failed in '", basename(path), "'.\n",
          "The normalised layout does not have A-P in the expected positions\n",
          "required by ", fn_name, " (data[10:25,1] and data[28:43,1]).\n",
          "This usually means the file has an unexpected number of metadata rows\n",
          "or the table structure differs from the standard BMG PHERAstar export.\n",
          paste(msgs, collapse = "\n")
        )
      }

      if (verbose)
        message("  [read] row label validation passed (A-P in both subtables)")
    }

    result
  }   # end .read_nanobret_raw

  # -- Internal helper: generate consolidated Excel report --------------------
  .generate_batch_report <- function(results, directory, verbose = TRUE) {

    report_path <- file.path(directory, "batch_analysis_report.xlsx")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary")

    summary_data <- data.frame(
      Sheet_Name = character(), Sheet_Number = character(),
      Data_File = character(), Function_Version = character(),
      Selected_Columns = character(), Status = character(),
      N_Constructs = integer(), Overall_Quality = character(),
      Control_Structure = character(), stringsAsFactors = FALSE
    )

    for (plate_sheet in names(results)) {
      plate_result      <- results[[plate_sheet]]$result
      n_constructs      <- 0L
      overall_quality   <- "Not assessed"
      control_structure <- "Unknown"
      selected_cols_info <- "All"

      if (!is.null(plate_result$control_info)) {
        if ("new_columns_created" %in% names(plate_result$control_info) &&
            !is.null(plate_result$control_info$new_columns_created))
          control_structure <- paste(plate_result$control_info$new_columns_created, collapse = ", ")
        else
          control_structure <- "Original columns"
      }

      if (!is.null(plate_result$selected_columns_info))
        selected_cols_info <- plate_result$selected_columns_info$description

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
        Selected_Columns = selected_cols_info,
        Status = "Completed",
        N_Constructs = n_constructs,
        Overall_Quality = overall_quality,
        Control_Structure = control_structure,
        stringsAsFactors = FALSE
      ))
    }

    openxlsx::writeData(wb, "Summary", summary_data)
    header_style <- openxlsx::createStyle(
      fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center",
      textDecoration = "Bold", border = "TopBottom", borderColour = "#4F81BD"
    )
    openxlsx::addStyle(wb, "Summary", header_style, rows = 1, cols = 1:ncol(summary_data))

    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result
      if (!is.null(plate_result$interval_means)) {
        openxlsx::addWorksheet(wb, plate_sheet)
        info_text <- c(
          paste("Data file:",        results[[plate_sheet]]$data_file),
          paste("Info sheet:",       results[[plate_sheet]]$info_sheet),
          paste("Function version:", results[[plate_sheet]]$function_version),
          paste("Control 0%:",  ifelse(is.numeric(results[[plate_sheet]]$control_0perc),
                                       paste("Fixed value", results[[plate_sheet]]$control_0perc),
                                       results[[plate_sheet]]$control_0perc)),
          paste("Control 100%:", ifelse(is.numeric(results[[plate_sheet]]$control_100perc),
                                        paste("Positions", paste(results[[plate_sheet]]$control_100perc, collapse = ", ")),
                                        paste(results[[plate_sheet]]$control_100perc, collapse = ", "))),
          paste("Selected columns:", ifelse(!is.null(results[[plate_sheet]]$selected_columns),
                                            paste(results[[plate_sheet]]$selected_columns, collapse = ", "),
                                            "All")),
          ""
        )
        openxlsx::writeData(wb, plate_sheet, info_text, startRow = 1)
        quality_data <- as.data.frame(t(plate_result$interval_means))
        quality_data <- cbind(Construct = rownames(quality_data), quality_data)
        rownames(quality_data) <- NULL
        openxlsx::writeData(wb, plate_sheet, quality_data, startRow = length(info_text) + 1)
        openxlsx::addStyle(wb, plate_sheet, header_style,
                           rows = length(info_text) + 1, cols = 1:ncol(quality_data))
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
        message(" . ", ver, ": ", n_ver, " plate(s)")
      }
      message("=============================\n")
    }
  }   # end .generate_batch_report

  # ========== MAIN FUNCTION BODY ==========

  if (verbose) {
    message("=== BATCH ANALYSIS SETTINGS ===")
    message("Using function version: ", function_version)
    if (function_version == "v1") {
      message(" . Original ratio_dose_response function")
      message(" . control_0perc: must be column name (e.g., 'DMSO_Control')")
      message(" . control_100perc: must be single column name (e.g., 'Stauro_Control')")
    } else {
      message(" . Enhanced ratio_dose_response_v2 function")
      message(" . control_0perc: can be fixed value (e.g., 16) OR column name")
      message(" . control_100perc: can be column positions (e.g., c(12, 24)) OR column names")
    }
    if (!is.null(selected_columns)) {
      message(" . Selected columns: ", paste(selected_columns, collapse = ", "))
      message("   (User view: 1 = first data column, excludes row names column)")
    } else {
      message(" . Selected columns: All data columns (1:24)")
    }
    message("================================")
  }

  if (directory == getwd() && verbose)
    message("Using current working directory: ", directory)

  if (is.null(output_dir)) {
    output_dir <- directory
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }

  if (!dir.exists(directory))
    stop("Directory not found: ", directory)

  info_path <- file.path(directory, info_file)
  if (!file.exists(info_path))
    stop("Info file not found: ", info_path)

  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required. Please install it.")

  info_sheets   <- openxlsx::getSheetNames(info_path)
  number_sheets <- info_sheets[grepl("\\d+$", info_sheets)]

  if (length(number_sheets) == 0L)
    stop("No valid number sheets found in info file.")

  if (verbose) {
    message("Found ", length(number_sheets), " number sheets: ",
            paste(number_sheets, collapse = ", "))
    message("Generate reports: ", generate_reports)
  }

  data_files <- list.files(directory, pattern = data_pattern, full.names = FALSE)
  data_files <- data_files[!grepl("^~\\$", data_files)]

  if (verbose) {
    message("Found ", length(data_files), " data files matching pattern")
    message("Output directory: ", output_dir)
  }

  results <- list()

  for (info_sheet in number_sheets) {

    sheet_number    <- gsub("^.*?(\\d+)$", "\\1", info_sheet)
    pattern         <- paste0("_", sheet_number, "\\.xlsx$")
    matching_files  <- data_files[grepl(pattern, data_files)]

    if (length(matching_files) == 0L) {
      if (verbose) message("No data file found for sheet ", info_sheet,
                           " (number ", sheet_number, "), skipping")
      next
    }

    data_filename <- matching_files[1L]
    if (length(matching_files) > 1L)
      warning("Multiple files match sheet ", info_sheet, ". Using: ", data_filename)

    data_path <- file.path(directory, data_filename)
    if (verbose) message("Processing sheet '", info_sheet, "' (number ",
                         sheet_number, ") with data: ", data_filename)

    tryCatch({

      info_table <- openxlsx::read.xlsx(info_path, sheet = info_sheet)

      # Parse raw data - robust to ID1 row and any other extra metadata rows;
      # validates A-P row labels before passing to ratio_dose_response.
      data <- .read_nanobret_raw(data_path,
                                 function_version = function_version,
                                 validate         = TRUE,
                                 verbose          = verbose)

      if (function_version == "v1") {
        result <- ratio_dose_response(
          data                = data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          save_to_excel       = NULL,
          selected_columns    = selected_columns
        )
      } else {
        result <- ratio_dose_response_v2(
          data                = data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          save_to_excel       = NULL,
          selected_columns    = selected_columns
        )
      }

      # Save individual plate results to Excel
      if (generate_reports) {
        excel_path <- file.path(output_dir, paste0("results_", sheet_number, ".xlsx"))
        wb <- openxlsx::createWorkbook()

        # BUG_7 fix: define header_style once here so it is always in scope,
        # regardless of whether result$interval_means is NULL.
        header_style <- openxlsx::createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center",
          textDecoration = "Bold", border = "TopBottom", borderColour = "#4F81BD"
        )

        openxlsx::addWorksheet(wb, "Modified_Ratio_Table")
        openxlsx::writeData(wb, "Modified_Ratio_Table",
                            cbind(RowNames = rownames(result$modified_ratio_table),
                                  result$modified_ratio_table),
                            rowNames = FALSE)

        openxlsx::addWorksheet(wb, "Original_Ratio_Table")
        openxlsx::writeData(wb, "Original_Ratio_Table",
                            cbind(RowNames = rownames(result$original_ratio_table),
                                  result$original_ratio_table),
                            rowNames = FALSE)

        if (!is.null(result$general_means)) {
          openxlsx::addWorksheet(wb, "General_Means")
          openxlsx::writeData(wb, "General_Means", result$general_means)
        }

        if (!is.null(result$interval_means)) {
          openxlsx::addWorksheet(wb, "Interval_Means")
          interval_df <- as.data.frame(result$interval_means)
          interval_with_rownames <- cbind(Metric = rownames(interval_df), interval_df)
          rownames(interval_with_rownames) <- NULL
          openxlsx::writeData(wb, "Interval_Means", interval_with_rownames, rowNames = FALSE)

          openxlsx::addStyle(wb, "Interval_Means", header_style,
                             rows = 1, cols = 1:ncol(interval_with_rownames))

          alt_row_style <- openxlsx::createStyle(fgFill = "#F2F2F2")
          n_rows <- nrow(interval_with_rownames)
          if (n_rows > 1L)
            for (row in seq(2L, n_rows, by = 2L))
              openxlsx::addStyle(wb, "Interval_Means", alt_row_style,
                                 rows = row, cols = 1:ncol(interval_with_rownames))
        }

        if (!is.null(result$control_info)) {
          openxlsx::addWorksheet(wb, "Control_Info")
          control_df <- data.frame(
            Parameter = names(result$control_info),
            Value = sapply(result$control_info, function(x) {
              if (is.null(x))        return("NULL")
              if (length(x) == 0L)   return("Empty")
              if (is.character(x))   return(paste(x, collapse = ", "))
              if (is.numeric(x))     return(paste(x, collapse = ", "))
              if (is.logical(x))     return(as.character(x))
              return(as.character(x))
            }),
            stringsAsFactors = FALSE
          )
          openxlsx::writeData(wb, "Control_Info", control_df)
          openxlsx::addStyle(wb, "Control_Info", header_style, rows = 1, cols = 1:2)
        }

        if (!is.null(result$selected_columns_info)) {
          openxlsx::addWorksheet(wb, "Selected_Columns_Info")
          selected_cols_df <- data.frame(
            Parameter = names(result$selected_columns_info),
            Value = sapply(result$selected_columns_info, function(x) {
              if (is.null(x))      return("NULL")
              if (length(x) == 0L) return("Empty")
              if (is.character(x)) return(x)
              if (is.numeric(x))   return(paste(x, collapse = ", "))
              return(as.character(x))
            }),
            stringsAsFactors = FALSE
          )
          openxlsx::writeData(wb, "Selected_Columns_Info", selected_cols_df)
          openxlsx::addStyle(wb, "Selected_Columns_Info", header_style, rows = 1, cols = 1:2)
        }

        openxlsx::addWorksheet(wb, "Processing_Info")
        processing_info <- data.frame(
          Parameter = c("Data File", "Info Sheet", "Function Version",
                        "Control 0%", "Control 100%", "Selected Columns",
                        "Split Replicates", "Low Value Threshold", "Processing Date"),
          Value = c(
            data_filename, info_sheet, function_version,
            ifelse(is.numeric(control_0perc),
                   paste("Fixed value:", control_0perc),
                   ifelse(is.null(control_0perc), "Not specified", control_0perc)),
            ifelse(is.numeric(control_100perc),
                   paste("Positions:", paste(control_100perc, collapse = ", ")),
                   ifelse(is.null(control_100perc), "Not specified",
                          paste(control_100perc, collapse = ", "))),
            ifelse(is.null(selected_columns), "All columns",
                   paste("Data columns:", paste(selected_columns, collapse = ", "))),
            split_replicates, low_value_threshold, as.character(Sys.time())
          ),
          stringsAsFactors = FALSE
        )
        openxlsx::writeData(wb, "Processing_Info", processing_info)
        openxlsx::addStyle(wb, "Processing_Info", header_style, rows = 1, cols = 1:2)

        openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
        if (verbose) message(" . Excel file saved: ", basename(excel_path))
      }

      results[[info_sheet]] <- list(
        data_file        = data_filename,
        info_sheet       = info_sheet,
        sheet_number     = sheet_number,
        function_version = function_version,
        control_0perc    = control_0perc,
        control_100perc  = control_100perc,
        selected_columns = selected_columns,
        result           = result
      )

      if (verbose) message("OK - Successfully processed sheet ", info_sheet,
                           " with ", function_version)

    }, error = function(e) {
      warning("Failed to process sheet ", info_sheet, " (", data_filename,
              ") with ", function_version, ": ", e$message)
      if (verbose) message("X Failed to process sheet ", info_sheet)
    })
  }

  if (length(results) > 0L && generate_reports) {
    .generate_batch_report(results, output_dir, verbose)
  } else if (length(results) > 0L) {
    if (verbose) message("Skipping report generation as requested")
  } else {
    warning("No plates were successfully processed")
  }

  return(results)
}
