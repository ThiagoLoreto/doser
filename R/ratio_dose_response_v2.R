#' Process and Analyze Dose-Response Ratio Data with Flexible Control Options
#'
#' Processes experimental data from dose-response assays to calculate ratios,
#' perform quality assessment, and prepare data for downstream curve fitting analysis.
#' This enhanced version supports both traditional column-based controls and new fixed-value
#' controls with multiple 100% control columns.
#'
#' @param data Data frame containing **raw dose-response experimental data** with specific
#'   structure typically exported from plate readers. Must have at least 43 rows with
#'   column names in row 9.
#' @param control_0perc **Either** a single numeric value (e.g., 16) for fixed 0% control 
#'   **or** a character string specifying the column name for 0% control (background control).
#'   When using fixed value, a new column 'Fixed_0perc' is created with this value.
#' @param control_100perc **Either** a numeric vector of column positions (e.g., c(12, 24)) 
#'   **or** a character vector of column names for 100% control(s). When using multiple
#'   columns, their row-wise means are calculated to create a new 'Mean_100perc' column.
#' @param split_replicates Logical indicating whether to split experimental replicates
#'   into separate columns (default: TRUE). Creates technical replicates with ".2" suffix.
#' @param info_table Data frame with experimental metadata containing at least 4 columns:
#'   log(inhibitor), Plate_Row, Construct, and Compound information.
#' @param save_to_excel Character string specifying Excel file path for saving processed results
#'   (default: NULL, no saving).
#' @param verbose Logical indicating whether to display progress messages (default: TRUE).
#' @param low_value_threshold Numeric threshold for filtering low luciferase signals 
#'   (default: 3000). Values below this are replaced with NA.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{original_ratio_table}: Original calculated ratio table from data
#'   \item \code{modified_ratio_table}: Processed and formatted table ready for analysis
#'   \item \code{general_means}: Data frame with general control means (when controls provided)
#'   \item \code{interval_means}: Data frame with construct-specific means and quality metrics
#'   \item \code{construct_intervals}: List mapping constructs to row intervals
#'   \item \code{control_info}: Information about how controls were processed
#' }
#'
#' @details
#' This function processes dose-response data through a comprehensive pipeline. The enhanced
#' version supports two control specification modes:
#'
#' \strong{New Fixed-Value Mode (Recommended for New Experiments):}
#' \itemize{
#'   \item \code{control_0perc = 16}: Uses fixed value 16 for all 0% control measurements
#'   \item \code{control_100perc = c(12, 24)}: Uses columns 12 and 24 as duplicate 100% controls
#'   \item Creates new columns: 'Fixed_0perc' (all values = 16) and 'Mean_100perc' (row-wise means)
#'   \item Removes original 100% control columns after calculating means
#'   \item Reorganizes: Fixed_0perc (first), experimental columns (middle), Mean_100perc (last)
#' }
#'
#' \strong{Traditional Column Name Mode (Backward Compatible):}
#' \itemize{
#'   \item \code{control_0perc = "DMSO"}: Uses column named "DMSO" for 0% control
#'   \item \code{control_100perc = "Staurosporine"}: Uses single column for 100% control
#'   \item Maintains original column structure and positions
#' }
#'
#' \strong{Key Features of Enhanced Version:}
#' \itemize{
#'   \item **Flexible Control Specification**: Accepts both fixed values and column names/positions
#'   \item **Multiple 100% Controls**: Automatically averages duplicate 100% control columns
#'   \item **Automatic Biological Replicate Handling**: Distinguishes replicates with suffixes
#'   \item **Quality Control Metrics**: Calculates Z-score, assay window, and signal quality
#'   \item **Technical Replicate Splitting**: Separates replicates into distinct columns
#'   \item **Metadata Integration**: Maps experimental metadata to data columns
#' }
#'
#' \strong{Data Processing Pipeline:}
#' 1. **Data Validation**: Checks minimum row requirements and structure
#' 2. **Subtable Extraction**: Separates measurement subtables (rows 10-25 and 28-43)
#' 3. **Low Value Filtering**: Replaces values below threshold with NA
#' 4. **Ratio Calculation**: Computes (subtable2 / subtable1) * 1000
#' 5. **Control Processing**: Implements fixed-value or column-based control strategy
#' 6. **Quality Assessment**: Calculates metrics for each construct
#' 7. **Data Transformation**: Transposes and renames columns based on metadata
#' 8. **Replicate Splitting**: Separates technical replicates if requested
#' 9. **Log(Inhibitor) Addition**: Adds concentration information as first column
#'
#' \strong{Data Structure Requirements:}
#' \preformatted{
#'   Rows 1-8:     Headers and metadata (ignored)
#'   Row 9:        Column names for data extraction
#'   Rows 10-25:   First measurement (e.g., donor channel)
#'   Rows 26-27:   Separator (ignored)
#'   Rows 28-43:   Second measurement (e.g., acceptor channel)
#' }
#'
#' \strong{Quality Metrics (When info_table Provided):}
#' \itemize{
#'   \item **Luciferase Signal**: >100000 (high), 10000-100000 (medium), 1000-10000 (low)
#'   \item **Assay Window**: >3 (high), 2-3 (medium), 1.5-2 (low)
#'   \item **Z-Score**: >0.7 (high), 0.5-0.7 (medium), 0.25-0.5 (low)
#'   \item **Overall Quality**: Determined by lowest-performing metric
#' }
#'
#' @examples
#' \dontrun{
#' # --- NEW FIXED-VALUE MODE (Recommended) ---
#' # Using fixed 0% control and column positions for 100% control
#' processed_data <- ratio_dose_response_v2(
#'   data = raw_experimental_data,
#'   control_0perc = 16,                # Fixed value for 0% control
#'   control_100perc = c(12, 24),       # Column positions for duplicate 100% controls
#'   info_table = experiment_design,
#'   split_replicates = TRUE,
#'   save_to_excel = "processed_results.xlsx"
#' )
#'
#' # Access the processed table ready for curve fitting
#' analysis_ready <- processed_data$modified_ratio_table
#'
#' # View quality metrics for each construct
#' quality_metrics <- processed_data$interval_means
#'
#' # --- TRADITIONAL MODE (Backward Compatible) ---
#' # Using column names for both controls
#' processed_old <- ratio_dose_response_v2(
#'   data = raw_experimental_data,
#'   control_0perc = "DMSO_Control",
#'   control_100perc = "Staurosporine_100pct",
#'   info_table = experiment_design
#' )
#'
#' # --- WITHOUT METADATA (Basic Processing) ---
#' # Just calculate ratios with new control structure
#' basic_processed <- ratio_dose_response_v2(
#'   data = raw_data,
#'   control_0perc = 16,
#'   control_100perc = c(12, 24),
#'   info_table = NULL,
#'   split_replicates = FALSE
#' )
#' }
#'
#' @section Control Processing Details:
#' \strong{When using fixed-value mode (control_0perc is numeric):}
#' \enumerate{
#'   \item Creates new column 'Fixed_0perc' with all values equal to control_0perc
#'   \item Calculates row-wise means of specified 100% control columns
#'   \item Creates new column 'Mean_100perc' with these averages
#'   \item Removes original 100% control columns from the data
#'   \item Reorganizes columns: Fixed_0perc → experimental columns → Mean_100perc
#'   \item After transposition: Fixed_0perc becomes first row, Mean_100perc becomes last row
#' }
#'
#' \strong{When using column name mode (control_0perc is character):}
#' \enumerate{
#'   \item Uses specified columns directly as controls
#'   \item Reorganizes columns: control_0perc column → experimental columns → control_100perc column
#'   \item No new columns created, original structure maintained
#' }
#'
#' @section Biological Replicate Handling:
#' The function automatically detects and distinguishes biological replicates
#' (same Construct:Compound combination) by adding numeric suffixes:
#' \itemize{
#'   \item First occurrence: "PIP4K2C:AZ-3458"
#'   \item Second occurrence: "PIP4K2C:AZ-3458_2"
#'   \item Third occurrence: "PIP4K2C:AZ-3458_3"
#' }
#' This allows proper analysis of replicate experiments within the same dataset.
#'
#' @section Output Structure:
#' The modified_ratio_table has the following structure:
#' \preformatted{
#'   Column 1:    log(inhibitor) concentration values (first row = NA)
#'   Columns 2-n: Experimental constructs (with .2 suffix for second technical replicate)
#'   
#'   Row structure:
#'   - Row 1:      0% control values (Fixed_0perc if using fixed-value mode)
#'   - Rows 2-(n-1): Experimental concentration measurements
#'   - Row n:      100% control values (Mean_100perc if using fixed-value mode)
#' }
#'
#' @section Note on Control Specifications:
#' \itemize{
#'   \item For **new experiments**, use fixed-value mode: \code{control_0perc = 16, control_100perc = c(12, 24)}
#'   \item For **legacy data compatibility**, use column name mode: \code{control_0perc = "ColumnName", control_100perc = "ColumnName"}
#'   \item You can mix modes: fixed 0% with named 100% control or vice versa
#'   \item Column positions are 1-indexed (first data column after row names = position 1)
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for dose-response curve fitting
#' \code{\link{plot_dose_response}} for visualization
#' \code{\link{calculate_assay_metrics}} for additional quality assessments
#'
#' @export
#'
#' @references
#' For dose-response data processing and quality control:
#' \itemize{
#'   \item "Dose-Response Data Analysis in Drug Discovery" (Motulsky & Christopoulos, 2004)
#'   \item "Assay Guidance Manual: Quantitative Biology and Pharmacology in Preclinical Drug Discovery" (NIH)
#'   \item "Best Practices in Dose-Response Assay Development" (Journal of Biomolecular Screening)
#' }



ratio_dose_response_v2 <- function(data,
                                   control_0perc = NULL,
                                   control_100perc = NULL,
                                   split_replicates = TRUE,
                                   info_table = NULL,
                                   save_to_excel = NULL,
                                   verbose = TRUE,
                                   low_value_threshold = 3000) {
  
  # --- DATA VALIDATION AND PREPARATION ---
  if (nrow(data) < 43) {
    stop("Data must have at least 43 rows")
  }
  
  # Set column names from row 9 of data
  colnames(data) <- as.character(data[9, ])
  
  # Extract the two main subtables (luciferase and normalizer readings)
  subtable1 <- data[10:25, 1:25]  # First measurement
  subtable2 <- data[28:43, 1:25]  # Second measurement
  
  final_rownames <- subtable1[, 1]
  
  convert_to_numeric_df <- function(df, rownames_vec) {
    num_df <- as.data.frame(apply(df[, -1], 2, as.numeric))
    rownames(num_df) <- rownames_vec
    return(num_df)
  }
  
  subtable1_num <- convert_to_numeric_df(subtable1, final_rownames)
  subtable2_num <- convert_to_numeric_df(subtable2, final_rownames)
  
  # Filter out low luciferase values (potential background noise)
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
  
  # --- CORE CALCULATION: RATIO ---
  # Calculate ratio: (subtable2 / subtable1) * 1000
  ratio <- (subtable2_num / subtable1_num) * 1000
  
  ratio_modified <- ratio  # Copy for modifications
  
  result <- list()
  
  # --- CONTROL ARGUMENT PROCESSING ---
  # Determine if control_0perc is a fixed value or column name
  control_0_is_value <- FALSE
  control_0_value <- NULL
  
  if (!is.null(control_0perc)) {
    if (is.numeric(control_0perc) && length(control_0perc) == 1) {
      control_0_is_value <- TRUE
      control_0_value <- control_0perc
      if (verbose) message("Using fixed value ", control_0_value, " for 0% control")
    } else if (is.character(control_0perc) && length(control_0perc) == 1) {
      control_0_is_value <- FALSE
      if (verbose) message("Using column '", control_0perc, "' for 0% control")
    } else {
      stop("control_0perc must be either a single numeric value or a single column name")
    }
  }
  
  # Determine if control_100perc are column positions or names
  control_100_positions <- NULL
  control_100_names <- NULL
  using_positions <- FALSE
  
  if (!is.null(control_100perc)) {
    if (is.numeric(control_100perc)) {
      control_100_positions <- control_100perc
      using_positions <- TRUE
      if (verbose) {
        message("Using column positions ", paste(control_100_positions, collapse = ", "), 
                " for 100% control")
      }
    } else if (is.character(control_100perc)) {
      control_100_names <- control_100perc
      if (verbose && length(control_100_names) > 1) {
        message("Using column names for 100% control: ", 
                paste(control_100_names, collapse = ", "))
      }
    } else {
      stop("control_100perc must be either numeric positions or character column names")
    }
  }
  
  # --- INFO TABLE PROCESSING ---
  # Process metadata for biological replicates and construct identification
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
  if (!is.null(info_table) && !is.null(control_0perc) && !is.null(control_100perc)) {
    
    plate_row_values <- info_table[[2]]
    construct_values <- info_table$Construct_Modified
    
    plate_row_to_index <- setNames(seq_along(plate_row_values), plate_row_values)
    
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
    
    # Convert positions to column names if using positions
    if (using_positions) {
      all_colnames <- colnames(ratio_modified)
      valid_positions <- control_100_positions[control_100_positions >= 1 & 
                                                 control_100_positions <= length(all_colnames)]
      control_100_colnames <- all_colnames[valid_positions]
    } else {
      control_100_colnames <- control_100_names
    }
    
    existing_100_cols <- control_100_colnames[control_100_colnames %in% colnames(ratio)]
    
    if (length(existing_100_cols) > 0 && length(row_intervals) > 0) {
      
      # General means calculation for 100% controls
      if (length(existing_100_cols) > 0) {
        general_mean_100 <- mean(as.matrix(ratio[, existing_100_cols, drop = FALSE]), na.rm = TRUE)
        
        result$general_means <- data.frame(
          Type = "General",
          Control_Type = "100%_Control",
          Mean = general_mean_100,
          Columns_Used = paste(existing_100_cols, collapse = ", "),
          stringsAsFactors = FALSE
        )
      }
      
      # Helper function to determine overall quality from multiple metrics
      get_lowest_comment <- function(luciferase_comment, assay_window_comment, z_score_comment) {
        priority_order <- c("insufficient", "low", "medium", "high")
        
        extract_level <- function(comment) {
          if (is.null(comment) || is.na(comment) || comment == "") {
            return("insufficient")
          }
          parts <- strsplit(as.character(comment), " ")[[1]]
          if (length(parts) > 0) {
            return(parts[1])
          } else {
            return("insufficient")
          }
        }
        
        luciferase_level <- extract_level(luciferase_comment)
        assay_level <- extract_level(assay_window_comment)
        z_level <- extract_level(z_score_comment)
        
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
              
              # Initialize variables
              z_score <- NA
              assay_window <- NA
              assay_window_comment <- NA
              z_score_comment <- NA
              mean_background <- NA
              mean_positive <- NA
              sd_background <- NA
              sd_positive <- NA
              
              # Calculate metrics if we have fixed 0% value and 100% controls
              if (control_0_is_value && length(existing_100_cols) > 0) {
                mean_background <- control_0_value
                
                # Calculate mean and SD of 100% controls for this construct
                control_100_data <- ratio[valid_rows, existing_100_cols, drop = FALSE]
                mean_positive <- mean(as.matrix(control_100_data), na.rm = TRUE)
                sd_positive <- sd(as.matrix(control_100_data), na.rm = TRUE)
                
                # Z-score and assay window calculations
                if (!is.na(mean_positive) && !is.na(mean_background)) {
                  sd_background <- 0  # Fixed value has no variability
                  
                  z_score <- if ((mean_positive - mean_background) != 0) {
                    1 - (3 * (sd_positive + sd_background) / (mean_positive - mean_background))
                  } else {
                    NA
                  }
                  
                  assay_window <- if (mean_background != 0) {
                    mean_positive / mean_background
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
                } else {
                  assay_window_comment <- "insufficient"
                  z_score_comment <- "insufficient"
                }
              }
              
              lowest_comment <- get_lowest_comment(luciferase_comment, assay_window_comment, z_score_comment)
              
              interval_means_list[[construct_name]] <- data.frame(
                Type = "Construct_Interval",
                Construct = construct_name,
                Average_Background = mean_background,
                SD_Background = sd_background,
                Average_Positive_Ctrl = mean_positive,
                SD_Positive_Ctrl = sd_positive,
                Average_luciferase_signal = mean_subtable1_all,
                Luciferase_signal_comment = luciferase_comment,
                Z_Score = z_score,
                Assay_z_Comment = z_score_comment,
                Assay_Window = assay_window,
                Assay_window_Comment = assay_window_comment,
                Overall_Quality = lowest_comment,
                Rows = paste0(construct_name, " (rows ", paste(range(valid_rows), collapse = "-"), ")"),
                Rows_Count = length(valid_rows),
                Background_Type = ifelse(control_0_is_value, "Fixed_Value", "Column"),
                Background_Value = ifelse(control_0_is_value, control_0_value, NA),
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
          "Rows", "Rows_Count", "Background_Type", "Background_Value"
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
  
  # --- NEW CONTROL COLUMN IMPLEMENTATION ---
  # Convert positions to column names if needed
  if (using_positions) {
    all_colnames <- colnames(ratio_modified)
    valid_positions <- control_100_positions[control_100_positions >= 1 & 
                                               control_100_positions <= length(all_colnames)]
    control_100_colnames <- all_colnames[valid_positions]
  } else {
    control_100_colnames <- control_100_names
  }
  
  existing_100_cols <- control_100_colnames[control_100_colnames %in% colnames(ratio_modified)]
  
  # When using fixed 0% value and multiple 100% control columns
  if (control_0_is_value && length(existing_100_cols) > 0) {
    
    # Calculate row-wise mean of 100% control columns
    if (length(existing_100_cols) == 1) {
      mean_100perc <- ratio_modified[, existing_100_cols]
    } else {
      mean_100perc <- rowMeans(ratio_modified[, existing_100_cols, drop = FALSE], na.rm = TRUE)
    }
    
    # 1. Create new column for fixed 0% control value
    ratio_modified$Fixed_0perc <- control_0_value
    
    # 2. Create new column for 100% control (mean of specified columns)
    ratio_modified$Mean_100perc <- mean_100perc
    
    # 3. Remove original 100% control columns
    ratio_modified <- ratio_modified[, !colnames(ratio_modified) %in% existing_100_cols]
    
    # 4. Reorder columns: Fixed_0perc first, experimental columns middle, Mean_100perc last
    all_cols <- colnames(ratio_modified)
    control_cols <- c("Fixed_0perc", "Mean_100perc")
    other_cols <- setdiff(all_cols, control_cols)
    
    ratio_modified <- ratio_modified[, c("Fixed_0perc", other_cols, "Mean_100perc")]
    
    if (verbose) {
      message("Created new control structure with fixed 0% = ", control_0_value, 
              " and mean 100% from ", length(existing_100_cols), " column(s)")
    }
    
  } else if (!control_0_is_value && length(existing_100_cols) == 1) {
    control_cols <- c(control_0perc, existing_100_cols)
    missing_controls <- control_cols[!control_cols %in% colnames(ratio_modified)]
    if (length(missing_controls) > 0) {
      stop("Control columns not found: ", paste(missing_controls, collapse = ", "))
    }
    
    other_columns <- setdiff(colnames(ratio_modified), control_cols)
    ratio_modified <- ratio_modified[, c(control_0perc, other_columns, existing_100_cols)]
  }
  
  # --- TABLE TRANSPOSITION ---
  # Transpose so constructs become columns and concentrations become rows
  ratio_modified_transposed <- as.data.frame(t(ratio_modified))
  colnames(ratio_modified_transposed) <- rownames(ratio_modified)
  
  # --- COLUMN RENAMING WITH INFO_TABLE ---
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
      
      # Ensure both replicates have same number of rows
      if (length(rep1_rows) != length(rep2_rows)) {
        min_len <- min(length(rep1_rows), length(rep2_rows))
        rep1_rows <- rep1_rows[1:min_len]
        rep2_rows <- rep2_rows[1:min_len]
      }
      
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
        rownames(df)[control_rows[1]],  # 0% control
        rownames(df)[rep1_rows],        # First technical replicate
        rownames(df)[control_rows[2]]   # 100% control
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
    
    if (length(log_col) > 0) {
      if (!is.na(log_col[1])) {
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
  
  # Control information for user reference
  result$control_info <- list(
    control_0_type = ifelse(control_0_is_value, "Fixed_Value", "Column"),
    control_0_value = if (control_0_is_value) control_0_value else control_0perc,
    control_100_input = if (using_positions) control_100_positions else control_100_names,
    control_100_actual_columns = existing_100_cols,
    new_columns_created = if (control_0_is_value && length(existing_100_cols) > 0) 
      c("Fixed_0perc", "Mean_100perc") else NULL
  )
  
  return(result)
}
