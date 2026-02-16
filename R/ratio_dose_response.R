#' Process and Analyze Raw Dose-Response Ratio Data
#'
#' Processes raw experimental data from dose-response assays to calculate BRET ratios,
#' perform quality assessment, and prepare data for downstream curve fitting analysis.
#' This function handles the initial data processing pipeline from raw plate reader
#' data to analysis-ready formatted data.
#'
#' @param data Data frame containing **raw dose-response experimental data** with specific
#'   structure typically exported from plate readers. Must have at least 43 rows with
#'   column names in row 9.
#' @param control_0perc Character specifying the column name for 0% control (background control,
#'   typically vehicle-treated samples like DMSO).
#' @param control_100perc Character specifying the column name for 100% control (positive control,
#'   typically maximum inhibition samples).
#' @param split_replicates Logical indicating whether to split experimental replicates
#'   into separate columns (default: TRUE).
#' @param info_table Data frame with experimental metadata containing at least 4 columns:
#'   log(inhibitor), Plate_Row, Construct, and Compound information.
#' @param save_to_excel Character string specifying Excel file path for saving processed results
#'   (default: NULL, no saving).
#' @param verbose Logical indicating whether to display progress messages
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{original_ratio_table}: Original calculated ratio table from raw data
#'   \item \code{modified_ratio_table}: Processed and formatted table ready for analysis
#'   \item \code{general_means}: Data frame with general control means across all rows
#'   \item \code{interval_means}: Data frame with construct-specific means and quality metrics
#'   \item \code{construct_intervals}: List mapping constructs to row intervals
#' }
#'
#' @details
#' This function performs the **initial processing of raw dose-response data** from BRET assays
#' or similar experimental formats. It transforms raw experimental measurements into analysis-ready
#' data through the following pipeline:
#'
#' \strong{Raw Data Processing Pipeline:}
#' \itemize{
#'   \item \strong{Raw Data Input}: Accepts direct output from plate readers or experimental data systems
#'   \item \strong{Data Extraction}: Separates measurement subtables (typically donor and acceptor channels)
#'   \item \strong{Quality Filtering}: Removes low-intensity signals (<1000) that may represent failed wells
#'   \item \strong{Ratio Calculation}: Computes BRET ratios as (subtable2 / subtable1) * 1000
#'   \item \strong{Control Processing}: Handles control wells separately for normalization
#'   \item **Quality Assessment**: Calculates assay performance metrics (Z-score, Assay Window)
#'   \item **Data Formatting**: Transposes and structures data for downstream analysis
#' }
#'
#' \strong{Typical Raw Data Structure:}
#' The function expects **raw experimental data** in the following format:
#' \preformatted{
#'   Rows 1-8:     Instrument headers, plate layout, experimental metadata (ignored)
#'   Row 9:        Column names (well identifiers or sample names)
#'   Rows 10-25:   First measurement subtable (e.g., donor channel or timepoint 1)
#'   Rows 26-27:   Separator or additional headers (ignored)
#'   Rows 28-43:   Second measurement subtable (e.g., acceptor channel or timepoint 2)
#' }
#'
#' \strong{Typical info_table structure:}
#' \preformatted{
#' Column 1: log(inhibitor) - Numeric values of inhibitor concentrations in log scale
#' Column 2: Plate_Row - Row identifiers matching the ratio table (e.g., "A", "B", "C", etc.)
#' Column 3: Construct - Construct protein identifiers (e.g., "BRD4", "EGFR", "KRAS")
#' Column 4: Compound - Compound identifiers (e.g., "JQ1", "Gefitinib", "ARS-1620")
#' }
#'
#' \strong{Quality Metrics for Raw Data Assessment:}
#' \itemize{
#'   \item \strong{Luciferase Signal}: Assesses raw signal intensity from experimental measurements
#'   \item \strong{Z-Score}: Evaluates assay robustness from raw control data variability
#'   \item \strong{Assay Window}: Calculates dynamic range from raw control measurements
#'   \item \strong{Overall Quality}: Determines if raw data quality supports further analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Typical workflow starting with raw experimental data
#' # Load raw data directly from plate reader export
#' raw_plate_data <- read.csv("plate_reader_export.csv")
#'
#' # Process raw data with control definitions
#' processed_data <- ratio_dose_response(
#'   data = raw_plate_data,  # Raw experimental data
#'   control_0perc = "DMSO_Ctrl",      # Background control from raw data
#'   control_100perc = "Stauro_Ctrl",  # Positive control from raw data
#'   info_table = sample_design,       # Experimental design metadata
#'   save_to_excel = "processed_data.xlsx"
#' )
#'
#' # The output is now ready for curve fitting analysis
#' analysis_ready_data <- processed_data$modified_ratio_table
#'
#' # Proceed to dose-response curve fitting
#' dr_results <- fit_drc_3pl(analysis_ready_data, normalize = TRUE)
#' }
#'
#' @section Raw Data Requirements:
#' The input data should be **raw experimental measurements** with:
#' \itemize{
#'   \item \strong{Direct instrument output}: Minimal preprocessing required
#'   \item \strong{Consistent structure}: Fixed row positions for data extraction
#'   \item \strong{Proper controls}: Clearly identified control wells in column names
#'   \item \strong{Metadata alignment}: Info table matching experimental design
#' }
#'
#' @section Expected Raw Data Sources:
#' \itemize{
#'   \item Plate reader exports (Tecan, BMG Labtech, PerkinElmer)
#'   \item BRET assay raw measurements
#'   \item Luminescence or fluorescence intensity data
#'   \item High-throughput screening raw data
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for the next step in the analysis pipeline
#' \code{\link{plot_dose_response}} for visualization of processed data
#'
#' @export
#'
#' @references
#' For raw data processing in dose-response assays:
#' \itemize{
#'   \item "BRET Assay Development Guide" (PerkinElmer)
#'   \item "High-Throughput Screening Data Analysis" (Inglese et al.)
#'   \item Journal of Biomolecular Screening raw data standards
#' }




ratio_dose_response <- function(data,
                                control_0perc = NULL, control_100perc = NULL,
                                split_replicates = TRUE, info_table = NULL,
                                save_to_excel = NULL, verbose = TRUE,
                                low_value_threshold = 1000,
                                selected_columns = NULL) {

  # --- DATA VALIDATION ---
  if (nrow(data) < 43) {
    stop("Data must have at least 43 rows")
  }

  colnames(data) <- as.character(data[9, ])

  # Extract complete subtables first
  subtable1_full <- data[10:25, 1:25]  # First measurement (luciferase)
  subtable2_full <- data[28:43, 1:25]  # Second measurement (normalizer)

  final_rownames <- subtable1_full[, 1]

  # --- COLUMN SELECTION LOGIC ---
  if (!is.null(selected_columns)) {
    # Validate selected_columns are numeric indices
    if (!is.numeric(selected_columns)) {
      stop("selected_columns must be numeric indices (e.g., c(2:23))")
    }

    # Adjust indices
    data_columns_indices <- selected_columns + 1

    # Validate indices are within bounds
    if (max(data_columns_indices) > ncol(subtable1_full)) {
      stop("Selected column index ", max(selected_columns),
           " is out of bounds. Maximum allowed: ", ncol(subtable1_full) - 1)
    }
    if (min(data_columns_indices) < 2) {
      stop("Selected column indices must be >= 1")
    }

    # Validate number of data columns is even
    if (length(selected_columns) %% 2 != 0) {
      warning("Number of selected data columns is not even (", length(selected_columns),
              " columns selected). This may cause issues with split_replicates.")
    }

    # Always keep column 1 (row names) + selected data columns
    columns_to_keep <- c(1, data_columns_indices)

    # Apply column selection to both subtables
    subtable1 <- subtable1_full[, columns_to_keep, drop = FALSE]
    subtable2 <- subtable2_full[, columns_to_keep, drop = FALSE]

  } else {
    # If no selection, use all columns
    subtable1 <- subtable1_full
    subtable2 <- subtable2_full
  }

  # Convert to numeric data frames
  convert_to_numeric_df <- function(df, rownames_vec) {
    num_df <- as.data.frame(apply(df[, -1], 2, as.numeric))
    rownames(num_df) <- rownames_vec
    return(num_df)
  }

  subtable1_num <- convert_to_numeric_df(subtable1, final_rownames)
  subtable2_num <- convert_to_numeric_df(subtable2, final_rownames)

  # --- CONTROL COLUMN MAPPING ---
  map_control_column <- function(control_spec) {
    if (is.null(control_spec)) return(NULL)

    if (is.numeric(control_spec)) {
      # User provided a data column index
      # Need to map to actual column index
      actual_col_index <- control_spec + 1

      # Check if exists
      if (actual_col_index > ncol(data)) {
        stop("Control column index ", control_spec, " is out of bounds. ",
             "Maximum allowed: ", ncol(data) - 1)
      }

      # Get column name
      col_name <- as.character(data[9, actual_col_index])
      return(list(index = actual_col_index, name = col_name, user_index = control_spec))
    } else {
      # User provided a column name directly
      return(list(index = which(colnames(data) == control_spec)[1],
                  name = control_spec,
                  user_index = NA))
    }
  }

  # Map controls
  control_0_info <- map_control_column(control_0perc)
  control_100_info <- map_control_column(control_100perc)

  # Check if controls are in selected columns
  if (!is.null(control_0_info)) {
    if (!control_0_info$name %in% colnames(subtable1_num)) {
      stop("Control column '", control_0_info$name, "' (data column ",
           ifelse(is.na(control_0_info$user_index), "named", control_0_info$user_index),
           ") not found in selected columns.")
    }
  }

  if (!is.null(control_100_info)) {
    if (!control_100_info$name %in% colnames(subtable1_num)) {
      stop("Control column '", control_100_info$name, "' (data column ",
           ifelse(is.na(control_100_info$user_index), "named", control_100_info$user_index),
           ") not found in selected columns.")
    }
  }

  # --- LOW VALUE FILTERING ---
  replace_low_values <- function(df, threshold = low_value_threshold) {
    for (col in colnames(df)) {
      low_vals <- df[[col]] < threshold & !is.na(df[[col]])
      if (any(low_vals)) {
        n_replaced <- sum(low_vals, na.rm = TRUE)
        warning("Replaced ", n_replaced, " value(s) < ", threshold, " with NA in column '", col, "'")
        df[[col]][low_vals] <- NA
      }
    }
    return(df)
  }

  subtable1_num <- replace_low_values(subtable1_num)

  if (any(subtable1_num == 0, na.rm = TRUE)) {
    warning("Division by zero detected in ratio calculation - replacing with NA")
    subtable1_num[subtable1_num == 0] <- NA
  }

  # --- CORE RATIO CALCULATION ---
  ratio <- (subtable2_num / subtable1_num) * 1000
  ratio_modified <- ratio  # Copy for modifications

  result <- list()

  # --- INFO TABLE PROCESSING ---
  # Process metadata for biological replicates
  if (!is.null(info_table)) {
    if (ncol(info_table) < 4) {
      stop("Info table must have at least 4 columns: log(inhibitor), Plate_Row, Construct, and Compound")
    }

    base_id_values <- paste(info_table[[3]], info_table[[4]], sep = ":")
    info_table$Base_ID <- base_id_values

    id_counts <- table(base_id_values)
    duplicate_ids <- names(id_counts)[id_counts > 1]

    # Distinguish biological replicates by adding suffixes
    if (length(duplicate_ids) > 0) {
      suffix_counter <- setNames(rep(1, length(duplicate_ids)), duplicate_ids)
      new_construct_values <- info_table[[3]]

      for (i in seq_along(base_id_values)) {
        current_base_id <- base_id_values[i]
        if (current_base_id %in% duplicate_ids) {
          if (suffix_counter[current_base_id] > 1) {
            new_construct_values[i] <- paste0(info_table[[3]][i], "_", suffix_counter[current_base_id])
          }
          suffix_counter[current_base_id] <- suffix_counter[current_base_id] + 1
        }
      }

      info_table$Construct_Modified <- new_construct_values

      if (verbose) {
        message("Found and automatically distinguished ", length(duplicate_ids),
                " biological replicate(s): ", paste(duplicate_ids, collapse = ", "))
      }
    } else {
      info_table$Construct_Modified <- info_table[[3]]
    }

    info_table$ID <- paste(info_table$Construct_Modified, info_table[[4]], sep = ":")
  }

  # --- QUALITY CONTROL CALCULATIONS ---
  # Calculate assay quality metrics if metadata and controls are provided
  if (!is.null(info_table) && !is.null(control_0_info) && !is.null(control_100_info)) {

    plate_row_values <- info_table[[2]]
    construct_values <- info_table$Construct_Modified

    plate_row_to_index <- setNames(seq_along(plate_row_values), plate_row_values)

    # Group constructs by plate rows
    construct_groups <- list()
    unique_constructs <- unique(construct_values)

    for (construct in unique_constructs) {
      construct_indices <- which(construct_values == construct)
      construct_groups[[construct]] <- plate_row_values[construct_indices]
    }

    row_intervals <- lapply(construct_groups, function(plate_rows) {
      indices <- plate_row_to_index[plate_rows]
      as.numeric(indices[!is.na(indices)])
    })

    row_intervals <- row_intervals[sapply(row_intervals, length) > 0]

    mean_columns <- c(control_0_info$name, control_100_info$name)
    existing_columns <- mean_columns[mean_columns %in% colnames(ratio)]

    if (length(existing_columns) > 0 && length(row_intervals) > 0) {
      # Calculate general means for control columns
      general_means <- colMeans(ratio[, existing_columns, drop = FALSE], na.rm = TRUE)

      result$general_means <- data.frame(
        Type = "General",
        Column = names(general_means),
        Mean = as.numeric(general_means),
        Rows = "All (1-16)",
        stringsAsFactors = FALSE
      )

      # Helper function to determine overall quality from multiple metrics
      get_lowest_comment <- function(luciferase_comment, assay_window_comment, z_score_comment) {
        priority_order <- c("insufficient", "low", "medium", "high")

        luciferase_level <- strsplit(luciferase_comment, " ")[[1]][1]
        assay_level <- strsplit(assay_window_comment, " ")[[1]][1]
        z_level <- strsplit(z_score_comment, " ")[[1]][1]

        luciferase_score <- match(luciferase_level, priority_order)
        assay_score <- match(assay_level, priority_order)
        z_score <- match(z_level, priority_order)

        scores_to_consider <- c(luciferase_score, assay_score)
        if (!is.na(z_score) && z_level != "insufficient") {
          scores_to_consider <- c(scores_to_consider, z_score)
        }

        if (all(is.na(scores_to_consider))) {
          return("insufficient")
        }

        min_score <- min(scores_to_consider, na.rm = TRUE)
        return(priority_order[min_score])
      }

      # Calculate quality metrics for each construct
      calculate_construct_means <- function() {
        interval_means_list <- list()

        for (construct_name in unique_constructs) {
          if (construct_name %in% names(row_intervals)) {
            interval_rows <- row_intervals[[construct_name]]
            valid_rows <- interval_rows[interval_rows >= 1 & interval_rows <= nrow(ratio)]

            if (length(valid_rows) > 0) {
              # Luciferase signal analysis
              subtable1_interval_data <- as.matrix(subtable1_num[valid_rows, ])
              mean_subtable1_all <- mean(subtable1_interval_data, na.rm = TRUE)

              luciferase_comment <- if (is.na(mean_subtable1_all)) {
                "insufficient luciferase signal"
              } else if (mean_subtable1_all > 100000) {
                "high (>100000)"
              } else if (mean_subtable1_all > 10000) {
                "medium (10000<x<100000)"
              } else if (mean_subtable1_all > 1000) {
                "low (1000<x<10000)"
              } else {
                "insufficient luciferase signal"
              }

              # Initialize quality metrics
              z_score <- NA
              assay_window <- NA
              assay_window_comment <- NA
              z_score_comment <- NA

              # Calculate Z-score and assay window if we have both controls
              if (length(existing_columns) == 2) {
                control_0_data <- ratio[valid_rows, control_0_info$name]
                control_100_data <- ratio[valid_rows, control_100_info$name]

                mean_0 <- mean(control_0_data, na.rm = TRUE)
                mean_100 <- mean(control_100_data, na.rm = TRUE)
                sd_0 <- sd(control_0_data, na.rm = TRUE)
                sd_100 <- sd(control_100_data, na.rm = TRUE)

                z_score <- if (!is.na(mean_100) && !is.na(mean_0) && (mean_100 - mean_0) != 0) {
                  1 - (3 * (sd_100 + sd_0) / (mean_100 - mean_0))
                } else {
                  NA
                }

                assay_window <- if (!is.na(mean_100) && !is.na(mean_0) && mean_0 != 0) {
                  mean_100 / mean_0
                } else {
                  NA
                }

                # Quality comments based on calculated values
                assay_window_comment <- if (is.na(assay_window)) {
                  "insufficient"
                } else if (assay_window > 3) {
                  "high (>3)"
                } else if (assay_window > 2) {
                  "medium (2<x<3)"
                } else if (assay_window > 1.5) {
                  "low (<2)"
                } else {
                  "insufficient"
                }

                z_score_comment <- if (is.na(z_score)) {
                  "insufficient"
                } else if (z_score > 0.7) {
                  "high (>0.7)"
                } else if (z_score > 0.5) {
                  "medium (0.5<x<0.7)"
                } else if (z_score > 0.25) {
                  "low (<0.5)"
                } else {
                  "insufficient"
                }
              }

              lowest_comment <- get_lowest_comment(luciferase_comment, assay_window_comment, z_score_comment)

              # Calculate control means for this construct
              mean_background <- if (length(existing_columns) >= 1) mean(ratio[valid_rows, control_0_info$name], na.rm = TRUE) else NA
              mean_positive <- if (length(existing_columns) >= 2) mean(ratio[valid_rows, control_100_info$name], na.rm = TRUE) else NA

              # Replace individual control values with construct means
              if (length(existing_columns) >= 1) {
                ratio_modified[valid_rows, control_0_info$name] <<- mean_background
              }
              if (length(existing_columns) >= 2) {
                ratio_modified[valid_rows, control_100_info$name] <<- mean_positive
              }

              interval_means_list[[construct_name]] <- data.frame(
                Type = "Construct_Interval",
                Construct = construct_name,
                Average_Background = mean_background,
                SD_Background = if (length(existing_columns) >= 1) sd(ratio[valid_rows, control_0_info$name], na.rm = TRUE) else NA,
                Average_Positive_Ctrl = mean_positive,
                SD_Positive_Ctrl = if (length(existing_columns) >= 2) sd(ratio[valid_rows, control_100_info$name], na.rm = TRUE) else NA,
                Average_luciferase_signal = mean_subtable1_all,
                Luciferase_signal_comment = luciferase_comment,
                Z_Score = z_score,
                Assay_z_Comment = z_score_comment,
                Assay_Window = assay_window,
                Assay_window_Comment = assay_window_comment,
                Overall_Quality = lowest_comment,
                Rows = paste0(construct_name, " (rows ", paste(range(valid_rows), collapse = "-"), ")"),
                Rows_Count = length(valid_rows),
                stringsAsFactors = FALSE
              )
            }
          }
        }
        return(interval_means_list)
      }

      interval_means_list <- calculate_construct_means()

      if (length(interval_means_list) > 0) {
        result$interval_means <- do.call(rbind, interval_means_list)
        rownames(result$interval_means) <- NULL

        # Reorder columns for better presentation
        result$interval_means <- result$interval_means[, c(
          "Type", "Construct",
          "Average_Positive_Ctrl", "SD_Positive_Ctrl",
          "Average_Background", "SD_Background",
          "Average_luciferase_signal", "Z_Score", "Assay_Window",
          "Luciferase_signal_comment", "Assay_window_Comment",
          "Assay_z_Comment", "Overall_Quality",
          "Rows", "Rows_Count"
        )]

        # Transpose for better Excel presentation
        if (nrow(result$interval_means) > 0) {
          interval_means_clean <- result$interval_means[, -c(1, 2)]
          rownames(interval_means_clean) <- result$interval_means$Construct
          interval_means_transposed <- as.data.frame(t(interval_means_clean))
          result$interval_means <- interval_means_transposed
        }

        result$construct_intervals <- row_intervals
      }
    }
  }

  # --- COLUMN REORGANIZATION ---
  # Place control columns at beginning and end
  if (!is.null(control_0_info) && !is.null(control_100_info)) {
    control_cols <- c(control_0_info$name, control_100_info$name)

    missing_controls <- control_cols[!control_cols %in% colnames(ratio_modified)]
    if (length(missing_controls) > 0) {
      stop("Control columns not found: ", paste(missing_controls, collapse = ", "))
    }

    other_columns <- setdiff(colnames(ratio_modified), control_cols)
    ratio_modified <- ratio_modified[, c(control_0_info$name, other_columns, control_100_info$name)]
  }

  # Transpose the table
  ratio_modified_transposed <- as.data.frame(t(ratio_modified))
  colnames(ratio_modified_transposed) <- rownames(ratio_modified)

  # --- COLUMN RENAMING WITH INFO TABLE ---
  # Replace column names with IDs from info_table
  if (!is.null(info_table)) {
    mapping <- setNames(info_table$ID, info_table[[2]])
    new_colnames <- mapping[colnames(ratio_modified_transposed)]

    colnames(ratio_modified_transposed) <- ifelse(
      is.na(new_colnames),
      colnames(ratio_modified_transposed),
      new_colnames
    )
  }

  # --- TECHNICAL REPLICATE SPLITTING ---
  # Split experimental concentrations into two technical replicates
  final_table <- if (split_replicates) {
    split_replicates_func <- function(df) {
      n_rows <- nrow(df)
      if (n_rows < 3) {
        warning("Not enough rows to split replicates. Returning original table.")
        return(df)
      }

      # First and last rows are controls
      control_rows <- c(1, n_rows)
      exp_rows <- 2:(n_rows - 1)  # Middle rows are experimental concentrations
      split_pt <- floor(length(exp_rows) / 2)

      rep1_rows <- exp_rows[1:split_pt]
      rep2_rows <- exp_rows[(split_pt + 1):length(exp_rows)]

      new_table <- data.frame()

      for (col in colnames(df)) {
        rep1 <- df[c(control_rows[1], rep1_rows, control_rows[2]), col]
        rep2 <- df[c(control_rows[1], rep2_rows, control_rows[2]), col]

        if (ncol(new_table) == 0) {
          new_table <- data.frame(rep1, rep2)
          colnames(new_table) <- c(col, paste0(col, ".2"))
        } else {
          new_table[[col]] <- rep1
          new_table[[paste0(col, ".2")]] <- rep2
        }
      }

      new_rownames <- c(
        rownames(df)[control_rows[1]],
        rownames(df)[rep1_rows],
        rownames(df)[control_rows[2]]
      )

      rownames(new_table) <- new_rownames
      return(new_table)
    }
    split_replicates_func(ratio_modified_transposed)
  } else {
    ratio_modified_transposed
  }

  # --- ADD LOG(INHIBITOR) COLUMN ---
  # Add log(inhibitor) concentration as first column
  if (!is.null(info_table)) {
    log_col <- info_table[[1]]
    n_needed <- nrow(final_table)

    if (length(log_col) > 0 && !is.na(log_col[1])) {
      warning("First row of log(inhibitor) is not NA. Shifting values and setting first row to NA.")

      adjusted_log_col <- c(NA, log_col[1:min(length(log_col), n_needed - 1)])

      if (length(adjusted_log_col) < n_needed) {
        adjusted_log_col <- c(adjusted_log_col, rep(NA, n_needed - length(adjusted_log_col)))
      }
    } else {
      adjusted_log_col <- if (length(log_col) >= n_needed) {
        log_col[1:n_needed]
      } else {
        c(log_col, rep(NA, n_needed - length(log_col)))
      }
    }

    final_table <- cbind(adjusted_log_col, final_table)
    colnames(final_table)[1] <- colnames(info_table)[1]
  }

  # --- EXCEL EXPORT ---
  if (!is.null(save_to_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required to save Excel files. Please install it.")
    }

    tryCatch({
      wb <- openxlsx::createWorkbook()

      openxlsx::addWorksheet(wb, "Modified_Ratio_Table")
      openxlsx::writeData(wb, "Modified_Ratio_Table",
                          cbind(RowNames = rownames(final_table), final_table),
                          rowNames = FALSE)

      if (!is.null(result$original_ratio_table)) {
        openxlsx::addWorksheet(wb, "Original_Ratio_Table")
        original_with_rownames <- cbind(RowNames = rownames(result$original_ratio_table),
                                        result$original_ratio_table)
        openxlsx::writeData(wb, "Original_Ratio_Table", original_with_rownames, rowNames = FALSE)
      }

      if (!is.null(result$general_means)) {
        openxlsx::addWorksheet(wb, "General_Means")
        openxlsx::writeData(wb, "General_Means", result$general_means)
      }

      if (!is.null(result$interval_means)) {
        openxlsx::addWorksheet(wb, "Interval_Means")
        openxlsx::writeData(wb, "Interval_Means", result$interval_means)
      }

      openxlsx::saveWorkbook(wb, save_to_excel, overwrite = TRUE)
      if (verbose) message("Excel file saved successfully: ", save_to_excel)

    }, error = function(e) {
      warning("Failed to save Excel file: ", e$message)
    })
  }

  # --- FINAL RESULT ASSEMBLY ---
  result$original_ratio_table <- ratio
  result$modified_ratio_table <- final_table

  # Add selected columns information
  result$selected_columns_info <- if (!is.null(selected_columns)) {
    list(
      user_indices = selected_columns,
      actual_indices = selected_columns + 1,
      description = paste("Data columns", paste(selected_columns, collapse = ", "))
    )
  } else {
    list(
      user_indices = 1:24,
      actual_indices = 2:25,
      description = "All data columns (1:24)"
    )
  }

  return(result)
}
