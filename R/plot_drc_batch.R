#' Plot Dose-Response Curves from Batch DRC Analysis
#'
#' @description
#' Creates publication-quality dose-response curve plots from batch DRC analysis results.
#' The function supports plotting single or multiple compounds across different plates,
#' with intelligent title generation, flexible formatting options, and automatic
#' handling of experimental duplicates.
#'
#' @details
#' This function extracts compound information from batch DRC results and creates
#' customizable dose-response curve plots. Key features include:
#' \itemize{
#'   \item Intelligent compound selection using target:compound syntax or position
#'   \item Automatic generation of descriptive plot titles
#'   \item Support for experimental duplicates with automatic labeling
#'   \item Flexible color and shape customization
#'   \item Publication-quality formatting with configurable themes
#'   \item Option to differentiate duplicates by both colors and shapes
#'   \item Error bars for data points with standard deviation
#'   \item Multiple output formats and saving options
#' }
#'
#' The function automatically detects the plotting strategy based on the number
#' of selected compounds and plates:
#' \itemize{
#'   \item **Single compound, multiple plates**: Colors (and optionally shapes) represent different plates/duplicates
#'   \item **Multiple compounds, single plate**: Colors represent different compounds
#'   \item **Multiple compounds, multiple plates**: Colors represent compounds, grouping by plates
#' }
#'
#' @param batch_drc_results List containing batch DRC analysis results from
#'   \code{\link{batch_drc_analysis}}. Each element should contain a \code{drc_result}
#'   with \code{detailed_results}.
#' @param target_compound Character string specifying which compound(s) to plot.
#'   Can be in several formats:
#'   \itemize{
#'     \item \code{"target:compound"} - Exact target:compound pair (e.g., "MEK:PD0325901")
#'     \item \code{"target"} - All compounds for a specific target (e.g., "MEK")
#'     \item \code{"compound"} - All compounds with specific compound name (e.g., "PD0325901")
#'   }
#'   Case-insensitive partial matching is supported.
#' @param position Integer specifying the position/index of the compound to plot.
#'   Useful for quickly plotting specific compounds without knowing their names.
#' @param y_limits Numeric vector of length 2 specifying Y-axis limits (default: c(0, 150)).
#' @param colors Character vector of colors for plot elements. If NULL, default
#'   ggplot2 colors are used. Length should match the number of groups/duplicates.
#' @param point_shapes Numeric vector of point shapes for duplicates (when
#'   \code{shape_by_duplicate = TRUE}). If NULL, default shapes are used.
#' @param show_error_bars Logical indicating whether to show error bars for
#'   data points (default: TRUE).
#' @param legend_position Character specifying legend position. One of: "none",
#'   "left", "right", "bottom", "top" (default: "right").
#' @param show_legend Logical indicating whether to show the legend (default: TRUE).
#' @param shape_by_duplicate Logical indicating whether to use different point
#'   shapes for different duplicates (in addition to colors). Only applies when
#'   plotting a single compound across multiple plates (default: FALSE).
#' @param show_grid Logical indicating whether to show plot grid lines (default: FALSE).
#' @param save_plot Character string specifying file path to save the plot.
#'   If NULL, plot is displayed but not saved. Supported formats include
#'   .png, .pdf, .jpg, .tiff.
#' @param plot_width Numeric specifying plot width in inches (default: 12).
#' @param plot_height Numeric specifying plot height in inches (default: 8).
#' @param plot_dpi Numeric specifying plot resolution in dots per inch (default: 600).
#' @param plot_title Character string for custom plot title. If NULL, an
#'   intelligent title is generated based on the selected compounds.
#' @param legend_title Character string for custom legend title. If NULL,
#'   an appropriate title is generated automatically.
#' @param y_axis_title Character string for custom Y-axis title. If NULL,
#'   title is set based on normalization status.
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#' @param axis_text_color Character specifying axis text color (default: "black").
#' @param axis_text_size Numeric specifying axis text size (default: 12).
#' @param axis_title_color Character specifying axis title color (default: "black").
#' @param axis_title_size Numeric specifying axis title size (default: 14).
#'
#' @return
#' Returns a ggplot2 object containing the dose-response curve plot. The plot
#' object has an additional "metadata" attribute containing:
#' \itemize{
#'   \item \code{selected_groups}: Names of selected compound groups
#'   \item \code{n_groups}: Number of compound groups
#'   \item \code{n_plates}: Number of plates
#'   \item \code{unique_targets}: Unique target names
#'   \item \code{unique_compounds}: Unique compound names
#'   \item \code{match_type}: Type of match used for compound selection
#'   \item \code{input_used}: Input that produced the match
#'   \item \code{is_normalized}: Whether data was normalized
#'   \item \code{plot_title}: Final plot title used
#'   \item \code{y_axis_title}: Final Y-axis title used
#'   \item \code{duplicate_labels}: Ordered duplicate labels
#'   \item \code{shape_by_duplicate}: Whether shapes were used for duplicates
#'   \item \code{show_grid}: Whether grid was shown
#' }
#'
#' @section Intelligent Title Generation:
#' The function automatically generates descriptive plot titles based on the
#' selected compounds:
#' \itemize{
#'   \item Single compound: Uses compound name
#'   \item Single target: Uses target name
#'   \item Multiple compounds with same target: Uses target name
#'   \item Multiple different compounds: "Multiple Compounds (n)"
#'   \item Custom title: Used if provided via \code{plot_title}
#' }
#'
#' @section Duplicate Handling:
#' Experimental plates are automatically labeled as "Duplicate 01", "Duplicate 02", etc.
#' The function extracts plate numbers from plate names and orders them numerically.
#' When \code{shape_by_duplicate = TRUE}, different point shapes are used in addition
#' to colors to visually distinguish between duplicates.
#'
#' @section Error Bars:
#' Error bars show mean +- standard deviation for each concentration point.
#' They are only displayed when replicate measurements are available and
#' \code{show_error_bars = TRUE}.
#'
#' @examples
#' \dontrun{
#' # Load example batch DRC results
#' data("example_batch_drc_results")
#'
#' # Plot specific compound by target:compound syntax
#' p1 <- plot_drc_batch(
#'   batch_drc_results = example_batch_drc_results,
#'   target_compound = "MEK:PD0325901"
#' )
#'
#' # Plot with shapes to differentiate duplicates
#' p2 <- plot_drc_batch(
#'   batch_drc_results = example_batch_drc_results,
#'   target_compound = "MEK:PD0325901",
#'   shape_by_duplicate = TRUE
#' )
#'
#' # Plot all compounds for a specific target
#' p3 <- plot_drc_batch(
#'   batch_drc_results = example_batch_drc_results,
#'   target_compound = "MEK"
#' )
#'
#' # Plot with custom colors and save to file
#' p4 <- plot_drc_batch(
#'   batch_drc_results = example_batch_drc_results,
#'   target_compound = "MEK:PD0325901",
#'   colors = c("#E41A1C", "#377EB8", "#4DAF4A"),
#'   save_plot = "drc_plot.png",
#'   plot_width = 10,
#'   plot_height = 6
#' )
#'
#' # Plot with grid and custom theme
#' p5 <- plot_drc_batch(
#'   batch_drc_results = example_batch_drc_results,
#'   target_compound = "MEK:PD0325901",
#'   show_grid = TRUE,
#'   axis_text_color = "darkblue",
#'   axis_title_color = "darkred"
#' )
#'
#' # Access plot metadata
#' metadata <- attr(p1, "metadata")
#' print(metadata$n_plates)
#' print(metadata$unique_compounds)
#' }
#'
#' @seealso
#' \code{\link{batch_drc_analysis}} for generating batch DRC results
#' \code{\link{plot_biological_replicates}} for comparing biological replicates
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @import ggplot2
#' @import scales
#' @export



plot_drc_batch <- function(batch_drc_results,
                                     target_compound = NULL,
                                     position = NULL,
                                     y_limits = c(0, 150),
                                     colors = NULL,
                                     point_shapes = NULL,
                                     show_error_bars = TRUE,
                                     legend_position = "right",
                                     show_legend = TRUE,
                                     shape_by_duplicate = FALSE,
                                     show_grid = FALSE,
                                     save_plot = NULL,
                                     plot_width = 12,
                                     plot_height = 8,
                                     plot_dpi = 600,
                                     plot_title = NULL,
                                     legend_title = NULL,
                                     y_axis_title = NULL,
                                     verbose = TRUE,
                                     axis_text_color = "black",
                                     axis_text_size = 12,
                                     axis_title_color = "black",
                                     axis_title_size = 14) {
  
  # ============================================================================
  # 1. DEPENDENCY CHECK AND INITIAL VALIDATION
  # ============================================================================
  
  # Check required packages
  required_packages <- c("ggplot2", "scales")
  missing_packages <- sapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  })
  
  if (any(missing_packages)) {
    stop("The following packages are required: ", 
         paste(required_packages[missing_packages], collapse = ", "))
  }
  
  # Load packages
  library(ggplot2)
  library(scales)
  
  # Validate input
  if (length(batch_drc_results) == 0) {
    stop("Batch DRC results list is empty")
  }
  
  if (verbose) {
    message("Processing DRC batch results for plotting...")
    message("Number of plates: ", length(batch_drc_results))
  }
  
  # ============================================================================
  # 2. HELPER FUNCTIONS
  # ============================================================================
  
  # Detect input type (target:compound format or single entry)
  detect_input_type <- function(input) {
    if (is.null(input)) return("none")
    if (grepl(":", input)) {
      parts <- strsplit(input, ":")[[1]]
      if (length(parts) == 2 && nchar(parts[1]) > 0 && nchar(parts[2]) > 0) {
        return("both")
      }
    }
    return("unknown")
  }
  
  # Search for matches in index lists (exact or partial)
  find_matches <- function(pattern, index_list, exact_first = TRUE) {
    # Exact match
    if (exact_first && pattern %in% names(index_list)) {
      return(list(matches = index_list[[pattern]], type = "exact"))
    }
    
    # Partial match (case-insensitive)
    matches_idx <- grepl(pattern, names(index_list), ignore.case = TRUE)
    if (any(matches_idx)) {
      matching_items <- unlist(index_list[matches_idx], recursive = FALSE)
      return(list(matches = matching_items, type = "partial"))
    }
    
    return(NULL)
  }
  
  # ============================================================================
  # 3. EXTRACT COMPOUND INFORMATION FROM ALL PLATES
  # ============================================================================
  
  extract_all_compounds_info <- function(batch_drc_results) {
    all_compounds <- list()
    
    for (plate_name in names(batch_drc_results)) {
      plate_result <- batch_drc_results[[plate_name]]$drc_result$detailed_results
      
      if (verbose && length(plate_result) > 0) {
        message("  Extracting compounds from plate: ", plate_name)
      }
      
      for (i in seq_along(plate_result)) {
        result <- plate_result[[i]]
        
        # Skip failed fits
        if (is.null(result$model) || !isTRUE(result$success)) next
        
        # Extract and parse compound name
        compound_name <- result$compound
        
        # Parse name format (target:compound | target:compound.2)
        if (grepl(" \\| ", compound_name)) {
          main_part <- strsplit(compound_name, " \\| ")[[1]][1]
          name_parts <- strsplit(main_part, ":")[[1]]
          target <- ifelse(length(name_parts) >= 1, name_parts[1], main_part)
          compound <- ifelse(length(name_parts) >= 2, name_parts[2], target)
        } else if (grepl(":", compound_name)) {
          name_parts <- strsplit(compound_name, ":")[[1]]
          target <- name_parts[1]
          compound <- ifelse(length(name_parts) > 1, name_parts[2], target)
        } else {
          target <- compound_name
          compound <- compound_name
        }
        
        # Clean names (remove .2 suffix)
        target_clean <- trimws(gsub("\\.2$", "", target))
        compound_clean <- trimws(gsub("\\.2$", "", compound))
        
        # Get experimental data
        data <- result$data
        if (is.null(data)) next
        
        valid_data <- data[is.finite(data$log_inhibitor) & is.finite(data$response), ]
        if (nrow(valid_data) < 2) next
        
        tryCatch({
          # Generate fitted curve points
          x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
          x_seq <- seq(x_range[1], x_range[2], length.out = 100)
          predicted <- predict(result$model, newdata = data.frame(log_inhibitor = x_seq))
          
          curve_data <- data.frame(
            log_inhibitor = x_seq,
            response = predicted,
            plate = plate_name,
            target = target_clean,
            compound = compound_clean,
            target_compound = paste(target_clean, compound_clean, sep = ":"),
            original_index = i
          )
          
          # Calculate point statistics (mean, SD, n)
          conc_levels <- unique(valid_data$log_inhibitor)
          point_stats <- do.call(rbind, lapply(conc_levels, function(conc) {
            conc_data <- valid_data[valid_data$log_inhibitor == conc, ]
            if (nrow(conc_data) > 0) {
              data.frame(
                log_inhibitor = conc,
                mean_response = mean(conc_data$response, na.rm = TRUE),
                sd_response = sd(conc_data$response, na.rm = TRUE),
                n_points = nrow(conc_data),
                plate = plate_name,
                target = target_clean,
                compound = compound_clean,
                target_compound = paste(target_clean, compound_clean, sep = ":"),
                original_index = i
              )
            }
          }))
          
          # Store compound information
          compound_info <- list(
            curve_data = curve_data,
            point_data = if (!is.null(point_stats) && nrow(point_stats) > 0) point_stats else data.frame(),
            target = target_clean,
            compound = compound_clean,
            target_compound = paste(target_clean, compound_clean, sep = ":"),
            plate = plate_name
          )
          
          key <- paste(target_clean, compound_clean, plate_name, i, sep = "_")
          all_compounds[[key]] <- compound_info
          
        }, error = function(e) {
          if (verbose) message("    Error processing compound: ", e$message)
        })
      }
    }
    
    return(all_compounds)
  }
  
  # Extract all compounds from batch results
  all_compounds <- extract_all_compounds_info(batch_drc_results)
  
  if (length(all_compounds) == 0) {
    stop("No compounds with valid data found in batch results")
  }
  
  if (verbose) {
    message("Total valid compounds extracted: ", length(all_compounds))
  }
  
  # ============================================================================
  # 4. CREATE SEARCH INDEX FOR COMPOUNDS
  # ============================================================================
  
  create_search_index <- function(all_compounds) {
    index <- list(by_target = list(), by_compound = list(), by_target_compound = list())
    
    for (key in names(all_compounds)) {
      info <- all_compounds[[key]]
      
      # Index by target
      if (!info$target %in% names(index$by_target)) {
        index$by_target[[info$target]] <- list()
      }
      index$by_target[[info$target]][[key]] <- info
      
      # Index by compound
      if (!info$compound %in% names(index$by_compound)) {
        index$by_compound[[info$compound]] <- list()
      }
      index$by_compound[[info$compound]][[key]] <- info
      
      # Index by target:compound
      tc <- info$target_compound
      if (!tc %in% names(index$by_target_compound)) {
        index$by_target_compound[[tc]] <- list()
      }
      index$by_target_compound[[tc]][[key]] <- info
    }
    
    return(index)
  }
  
  search_index <- create_search_index(all_compounds)
  
  if (verbose) {
    message("Search index created:")
    message("  Unique targets: ", length(search_index$by_target))
    message("  Unique compounds: ", length(search_index$by_compound))
    message("  Unique target:compound pairs: ", length(search_index$by_target_compound))
  }
  
  # ============================================================================
  # 5. SELECT COMPOUNDS BASED ON USER CRITERIA
  # ============================================================================
  
  selected_compounds <- list()
  match_type <- "unknown"
  input_used <- NULL
  
  # Selection by position
  if (!is.null(position)) {
    if (position > length(all_compounds)) {
      stop("Position ", position, " is out of range. Only ", length(all_compounds), " compounds available.")
    }
    selected_key <- names(all_compounds)[position]
    selected_compounds[[selected_key]] <- all_compounds[[selected_key]]
    match_type <- "position"
    input_used <- paste("Position", position)
    
    if (verbose) {
      info <- all_compounds[[selected_key]]
      message("Selected by position ", position, ": ", info$target_compound, " (plate: ", info$plate, ")")
    }
  } 
  # Selection by target_compound
  else if (!is.null(target_compound)) {
    input_type <- detect_input_type(target_compound)
    input_used <- target_compound
    
    if (verbose) message("Input '", target_compound, "' detected as type: ", input_type)
    
    if (input_type == "both") {
      # Search for target:compound
      result <- find_matches(target_compound, search_index$by_target_compound)
      if (!is.null(result)) {
        selected_compounds <- result$matches
        match_type <- ifelse(result$type == "exact", "target_compound", "target_compound_partial")
      } else {
        # Separate and search individually
        parts <- strsplit(target_compound, ":")[[1]]
        target_part <- parts[1]
        compound_part <- ifelse(length(parts) > 1, parts[2], "")
        
        # Try target
        result <- find_matches(target_part, search_index$by_target)
        if (!is.null(result)) {
          selected_compounds <- result$matches
          match_type <- ifelse(result$type == "exact", "target", "target_partial")
          input_used <- target_part
        } 
        # Try compound
        else if (compound_part != "") {
          result <- find_matches(compound_part, search_index$by_compound)
          if (!is.null(result)) {
            selected_compounds <- result$matches
            match_type <- ifelse(result$type == "exact", "compound", "compound_partial")
            input_used <- compound_part
          }
        }
      }
    } 
    else {
      # Search first as compound, then as target
      result <- find_matches(target_compound, search_index$by_compound)
      if (!is.null(result)) {
        selected_compounds <- result$matches
        match_type <- ifelse(result$type == "exact", "compound", "compound_partial")
      } else {
        result <- find_matches(target_compound, search_index$by_target)
        if (!is.null(result)) {
          selected_compounds <- result$matches
          match_type <- ifelse(result$type == "exact", "target", "target_partial")
        }
      }
    }
    
    if (length(selected_compounds) == 0) {
      stop("No match found for '", target_compound, "'.")
    }
  } 
  # No criteria provided - use first compound
  else {
    if (verbose) message("No selection criteria provided. Using first compound.")
    first_key <- names(all_compounds)[1]
    selected_compounds[[first_key]] <- all_compounds[[first_key]]
    match_type <- "first"
  }
  
  if (length(selected_compounds) == 0) {
    stop("No compounds selected for plotting.")
  }
  
  if (verbose) {
    message("Selected ", length(selected_compounds), " compounds (match type: ", match_type, ")")
    message("Show grid: ", show_grid)
  }
  
  # ============================================================================
  # 6. GROUP COMPOUNDS FOR PLOTTING
  # ============================================================================
  
  selected_groups <- list()
  for (key in names(selected_compounds)) {
    info <- selected_compounds[[key]]
    group_key <- info$target_compound
    
    if (!group_key %in% names(selected_groups)) {
      selected_groups[[group_key]] <- list()
    }
    selected_groups[[group_key]][[info$plate]] <- info
  }
  
  n_groups <- length(selected_groups)
  
  # ============================================================================
  # 7. PREPARE PLOT DATA
  # ============================================================================
  
  prepare_plot_data <- function(selected_groups) {
    process_compound_data <- function(group_key, plate_name, compound_info) {
      # Create duplicate label (e.g., "Duplicate 01", "Duplicate 02")
      plate_num <- gsub(".*?([0-9]+).*", "\\1", plate_name)
      duplicate_label <- if (grepl("^[0-9]+$", plate_num)) {
        paste("Duplicate", sprintf("%02d", as.numeric(plate_num)))
      } else {
        paste("Duplicate", plate_name)
      }
      
      plot_group <- paste(group_key, plate_name, sep = "_")
      
      # Process curve data
      curve_df <- compound_info$curve_data
      curve_df$group_key <- group_key
      curve_df$plate <- plate_name
      curve_df$duplicate_label <- duplicate_label
      curve_df$plot_group <- plot_group
      
      # Process point data
      point_df <- compound_info$point_data
      if (nrow(point_df) > 0) {
        point_df$group_key <- group_key
        point_df$plate <- plate_name
        point_df$duplicate_label <- duplicate_label
        point_df$plot_group <- plot_group
      }
      
      return(list(curves = curve_df, points = point_df))
    }
    
    all_curves <- list()
    all_points <- list()
    
    for (group_key in names(selected_groups)) {
      for (plate_name in names(selected_groups[[group_key]])) {
        compound_info <- selected_groups[[group_key]][[plate_name]]
        processed <- process_compound_data(group_key, plate_name, compound_info)
        
        all_curves[[length(all_curves) + 1]] <- processed$curves
        if (nrow(processed$points) > 0) {
          all_points[[length(all_points) + 1]] <- processed$points
        }
      }
    }
    
    if (length(all_curves) == 0) {
      stop("No curve data available for plotting")
    }
    
    list(
      curves = do.call(rbind, all_curves),
      points = if (length(all_points) > 0) do.call(rbind, all_points) else data.frame()
    )
  }
  
  plot_data <- prepare_plot_data(selected_groups)
  n_plates <- length(unique(plot_data$curves$plate))
  
  # ============================================================================
  # 8. ANALYZE DATA FOR TITLE AND LEGEND GENERATION
  # ============================================================================
  
  unique_targets <- unique(sapply(selected_compounds, function(x) x$target))
  unique_compounds <- unique(sapply(selected_compounds, function(x) x$compound))
  
  if (verbose) {
    message("Analysis for title generation:")
    message("  Unique targets: ", paste(unique_targets, collapse = ", "))
    message("  Unique compounds: ", paste(unique_compounds, collapse = ", "))
    message("  Match type: ", match_type)
    message("  Input used: ", ifelse(is.null(input_used), "None", input_used))
    message("  Shape by duplicate: ", shape_by_duplicate)
    message("Plot data prepared: ", nrow(plot_data$curves), " curve points")
    if (nrow(plot_data$points) > 0) {
      message("  and ", nrow(plot_data$points), " data points")
    }
  }
  
  # Check if data was normalized
  is_normalized <- FALSE
  if (length(batch_drc_results) > 0) {
    first_plate <- names(batch_drc_results)[1]
    if (!is.null(batch_drc_results[[first_plate]]$drc_result$normalized)) {
      is_normalized <- batch_drc_results[[first_plate]]$drc_result$normalized
    }
  }
  
  if (verbose) {
    message("Data normalization status: ", ifelse(is_normalized, "Normalized", "Not normalized"))
  }
  
  # Set Y axis title
  y_axis_title_final <- if (is.null(y_axis_title)) {
    ifelse(is_normalized, "Normalized BRET ratio [%]", "BRET ratio")
  } else {
    y_axis_title
  }
  
  # ============================================================================
  # 9. GENERATE INTELLIGENT PLOT TITLE
  # ============================================================================
  
  generate_intelligent_title <- function(match_type, input_used, unique_targets, unique_compounds, n_groups) {
    
    # Custom title provided by user
    if (!is.null(plot_title)) {
      if (verbose) message("Using custom plot title: ", plot_title)
      return(plot_title)
    }
    
    # Single compound (all have same compound name)
    if (length(unique_compounds) == 1 && unique_compounds[1] != "") {
      title <- unique_compounds[1]
      if (verbose) message("Title strategy: Single compound -> ", title)
    }
    # Single target (all have same target name)
    else if (length(unique_targets) == 1 && unique_targets[1] != "") {
      title <- unique_targets[1]
      if (verbose) message("Title strategy: Single target -> ", title)
    }
    # Based on match type
    else if (match_type %in% c("compound", "compound_partial")) {
      title <- input_used
      if (verbose) message("Title strategy: Matched by compound -> ", title)
    }
    else if (match_type %in% c("target", "target_partial")) {
      title <- input_used
      if (verbose) message("Title strategy: Matched by target -> ", title)
    }
    else if (match_type == "target_compound") {
      # Extract only compound from target:compound format
      if (grepl(":", input_used)) {
        parts <- strsplit(input_used, ":")[[1]]
        if (length(parts) >= 2) {
          title <- parts[2]  # Use only compound name
        } else {
          title <- input_used
        }
      } else {
        title <- input_used
      }
      if (verbose) message("Title strategy: Matched by target:compound -> ", title)
    }
    # Default descriptive title
    else {
      if (n_groups > 1) {
        title <- paste("Multiple Compounds (", n_groups, ")", sep = "")
      } else {
        title <- "Dose-Response Curve"
      }
      if (verbose) message("Title strategy: Default -> ", title)
    }
    
    return(title)
  }
  
  plot_title_final <- generate_intelligent_title(match_type, input_used, unique_targets, 
                                                 unique_compounds, n_groups)
  
  # ============================================================================
  # 10. SET LEGEND TITLE AND ORDER DUPLICATES
  # ============================================================================
  
  # Set legend title
  if (is.null(legend_title)) {
    legend_title_final <- if (n_groups == 1 && n_plates > 1) {
      "Duplicate"
    } else if (n_groups > 1 && n_plates == 1) {
      "Compound"
    } else if (n_groups > 1 && n_plates > 1) {
      "Group"
    } else {
      ifelse(length(unique_compounds) == 1, "Compound", "Target")
    }
  } else {
    legend_title_final <- legend_title
  }
  
  # Order duplicates numerically
  get_duplicate_number <- function(label) {
    num <- gsub(".*?([0-9]+).*", "\\1", label)
    ifelse(grepl("^[0-9]+$", num), as.numeric(num), 999)
  }
  
  unique_duplicates <- unique(plot_data$curves$duplicate_label)
  sorted_duplicates <- unique_duplicates[order(sapply(unique_duplicates, get_duplicate_number))]
  
  # Convert to factor to maintain order
  plot_data$curves$duplicate_label <- factor(plot_data$curves$duplicate_label, levels = sorted_duplicates)
  if (nrow(plot_data$points) > 0) {
    plot_data$points$duplicate_label <- factor(plot_data$points$duplicate_label, levels = sorted_duplicates)
  }
  
  # ============================================================================
  # 11. CONFIGURE SHAPES FOR DUPLICATES (IF REQUESTED)
  # ============================================================================
  
  if (shape_by_duplicate && n_plates > 1 && n_groups == 1) {
    default_shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10, 11, 12, 13, 14)
    shape_values <- if (!is.null(point_shapes) && length(point_shapes) >= n_plates) {
      point_shapes[1:n_plates]
    } else {
      default_shapes[1:n_plates]
    }
    shape_mapping <- setNames(shape_values, sorted_duplicates)
  }
  
  # ============================================================================
  # 12. CREATE THE PLOT
  # ============================================================================
  
  p <- ggplot2::ggplot()
  
  # Configure colors
  if (is.null(colors)) {
    colors <- if (n_groups == 1) scales::hue_pal()(n_plates) else scales::hue_pal()(n_groups)
  }
  
  # Function to define plot elements based on strategy
  add_plot_elements <- function() {
    if (n_groups == 1) {
      # Single compound group, multiple plates
      if (shape_by_duplicate && n_plates > 1) {
        # Use both colors and shapes for duplicates
        return(list(
          geom_line = ggplot2::aes(x = log_inhibitor, y = response, 
                                   color = duplicate_label, group = plot_group),
          geom_point = ggplot2::aes(x = log_inhibitor, y = mean_response, 
                                    color = duplicate_label, shape = duplicate_label),
          scale_color = ggplot2::scale_color_manual(
            values = setNames(colors[1:n_plates], sorted_duplicates), 
            name = legend_title_final, breaks = sorted_duplicates),
          scale_shape = ggplot2::scale_shape_manual(
            values = shape_mapping, name = legend_title_final, breaks = sorted_duplicates),
          error_aes = ggplot2::aes(x = log_inhibitor, 
                                   ymin = mean_response - sd_response, 
                                   ymax = mean_response + sd_response, 
                                   color = duplicate_label)
        ))
      } else {
        # Use only colors for duplicates
        return(list(
          geom_line = ggplot2::aes(x = log_inhibitor, y = response, 
                                   color = duplicate_label, group = plot_group),
          geom_point = ggplot2::aes(x = log_inhibitor, y = mean_response, 
                                    color = duplicate_label),
          scale_color = ggplot2::scale_color_manual(
            values = setNames(colors[1:n_plates], sorted_duplicates), 
            name = legend_title_final, breaks = sorted_duplicates),
          scale_shape = NULL,
          error_aes = ggplot2::aes(x = log_inhibitor, 
                                   ymin = mean_response - sd_response, 
                                   ymax = mean_response + sd_response, 
                                   color = duplicate_label)
        ))
      }
    } else {
      # Multiple compound groups
      return(list(
        geom_line = ggplot2::aes(x = log_inhibitor, y = response, 
                                 color = group_key, group = plot_group),
        geom_point = ggplot2::aes(x = log_inhibitor, y = mean_response, 
                                  color = group_key),
        scale_color = ggplot2::scale_color_manual(
          values = setNames(colors[1:n_groups], names(selected_groups)), 
          name = legend_title_final),
        scale_shape = NULL,
        error_aes = ggplot2::aes(x = log_inhibitor, 
                                 ymin = mean_response - sd_response, 
                                 ymax = mean_response + sd_response, 
                                 color = group_key)
      ))
    }
  }
  
  # Get plot elements and add to plot
  plot_elements <- add_plot_elements()
  
  p <- p + 
    ggplot2::geom_line(data = plot_data$curves, plot_elements$geom_line, 
                       linewidth = 1, alpha = 0.8)
  
  if (nrow(plot_data$points) > 0) {
    p <- p + ggplot2::geom_point(data = plot_data$points, plot_elements$geom_point, 
                                 size = 3)
  }
  
  p <- p + plot_elements$scale_color
  
  if (!is.null(plot_elements$scale_shape)) {
    p <- p + plot_elements$scale_shape
  }
  
  # Add error bars if requested
  if (show_error_bars && nrow(plot_data$points) > 0 && "sd_response" %in% colnames(plot_data$points)) {
    p <- p + ggplot2::geom_errorbar(data = plot_data$points, plot_elements$error_aes, 
                                    width = 0.05, linewidth = 0.5, alpha = 0.6)
  }
  
  # ============================================================================
  # 13. CONFIGURE PLOT THEME
  # ============================================================================
  
  base_theme <- ggplot2::theme_minimal() + 
    ggplot2::theme(
      legend.position = ifelse(show_legend, legend_position, "none"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = ggplot2::element_text(color = axis_text_color, size = axis_text_size),
      axis.title = ggplot2::element_text(color = axis_title_color, size = axis_title_size, face = "bold"),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.x.top = ggplot2::element_blank(),
      axis.line.y.right = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      panel.border = ggplot2::element_blank()
    )
  
  # Remove grid if requested
  if (!show_grid) {
    base_theme <- base_theme + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  }
  
  # Apply theme and labels
  p <- p + 
    ggplot2::labs(x = expression(paste("Log"[10], " Concentration [M]")),
                  y = y_axis_title_final,
                  title = plot_title_final) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    base_theme
  
  # Adjust legend when using shapes
  if (shape_by_duplicate && n_groups == 1 && n_plates > 1) {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(shape = shape_mapping[sorted_duplicates]))
    )
  }
  
  # ============================================================================
  # 14. DISPLAY AND SAVE PLOT
  # ============================================================================
  
  print(p)
  
  if (!is.null(save_plot)) {
    ggplot2::ggsave(save_plot, plot = p, width = plot_width, 
                    height = plot_height, dpi = plot_dpi)
    if (verbose) message("Plot saved as: ", save_plot)
  }
  
  # ============================================================================
  # 15. RETURN METADATA
  # ============================================================================
  
  metadata <- list(
    selected_groups = names(selected_groups),
    n_groups = n_groups,
    n_plates = n_plates,
    unique_targets = unique_targets,
    unique_compounds = unique_compounds,
    match_type = match_type,
    input_used = input_used,
    is_normalized = is_normalized,
    plot_title = plot_title_final,
    y_axis_title = y_axis_title_final,
    duplicate_labels = sorted_duplicates,
    shape_by_duplicate = shape_by_duplicate,
    show_grid = show_grid
  )
  
  attr(p, "metadata") <- metadata
  
  return(p)
}
