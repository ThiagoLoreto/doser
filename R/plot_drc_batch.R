#' Plot Dose-Response Curves (DRC) from Batch Analysis Results
#'
#' @description
#' This function generates publication-quality Dose-Response Curve (DRC) plots from a list of
#' batch processing results. It handles data extraction, statistical aggregation, and visualization
#' in a single step. It supports advanced features like scientific color palettes, faceting,
#' automatic IC50/EC50 annotation, and custom styling.
#'
#' @param batch_drc_results A named list containing the results of the DRC batch analysis.
#'   Structure expected: \code{list(plate_name = list(drc_result = list(detailed_results = ...)))}.
#' @param target_compound Character string. The specific target/compound to plot.
#'   Supports two formats:
#'   \itemize{
#'     \item \strong{"Target:Compound"}: Performs an exact match (e.g., "PAK3:IPA-3").
#'     \item \strong{"SearchString"}: Performs a fuzzy search across targets and compounds (e.g., "IPA-3").
#'   }
#'   Ignored if \code{position} is provided.
#' @param position Integer. Selects a compound by its numerical index in the list of unique results.
#'   Useful for iterating through compounds blindly.
#' @param y_limits Numeric vector of length 2 (e.g., \code{c(0, 150)}). Sets the Y-axis limits.
#'   If \code{NULL} (default), limits are automatic.
#'   \emph{Note:} Ignored if \code{facet_scales} contains "free".
#' @param colors Character vector or string. Controls the line/point colors.
#'   \itemize{
#'     \item \strong{Manual}: A vector of color codes/names (e.g., \code{c("red", "blue")}).
#'     \item \strong{Scientific Palettes (ggsci)}: "npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb", "d3".
#'     \item \strong{Perceptual Palettes (viridis)}: "viridis", "magma", "plasma", "inferno", "cividis", "turbo".
#'     \item \strong{RColorBrewer}: Names like "Set1", "Dark2", "Paired".
#'   }
#' @param point_shapes Integer vector. Custom shapes for points (0-25). Used only if \code{shape_by_duplicate = TRUE}.
#' @param show_error_bars Logical. If \code{TRUE} (default), adds error bars (mean +/- SD) to the points.
#' @param legend_position Character. Position of the legend ("right", "bottom", "top", "none").
#' @param show_legend Logical. Whether to display the legend.
#' @param shape_by_duplicate Logical. If \code{TRUE} and plotting a single compound group,
#'   different replicates (plates) will have different point shapes.
#' @param show_grid Logical. If \code{TRUE}, displays grid lines. If \code{FALSE} (default),
#'   produces a clean "classic" scientific style (axis lines only).
#' @param font_family Character. Font family for text elements.
#'   Common options: "sans" (Arial/Helvetica), "serif" (Times New Roman), "mono" (Courier).
#'   For custom fonts, use the \code{extrafont} package.
#' @param facet_by Character. Variable to split the plot into panels.
#'   Options: \code{"compound"}, \code{"target"}, \code{"plate"}. Default is \code{NULL} (single plot).
#' @param facet_ncol Integer. Number of columns for the faceted grid.
#' @param facet_scales Character. Scaling of axes in facets.
#'   Options: \code{"fixed"} (default), \code{"free"}, \code{"free_y"}, \code{"free_x"}.
#' @param show_IC50 Logical. If \code{TRUE}, extracts the IC50/EC50 from the model and displays it
#'   in the plot legend. Also prints a summary table to the console.
#' @param save_plot Character string. File path to save the plot (e.g., "plot.png", "plot.pdf").
#'   If \code{NULL}, the plot is not saved automatically.
#' @param plot_width Numeric. Width of the saved plot in inches.
#' @param plot_height Numeric. Height of the saved plot in inches.
#' @param plot_dpi Integer. Resolution of the saved plot (default: 600).
#' @param plot_title Character. Custom main title. If \code{NULL}, generates a smart title based on selection.
#' @param legend_title Character. Custom legend title.
#' @param y_axis_title Character. Custom Y-axis title. Default is "BRET ratio".
#' @param verbose Logical. If \code{TRUE}, prints progress messages and IC50 tables to the console.
#' @param axis_text_color Character. Color of axis tick labels.
#' @param axis_text_size Numeric. Size of axis tick labels.
#' @param axis_title_color Character. Color of axis titles.
#' @param axis_title_size Numeric. Size of axis titles.
#'
#' @return A \code{ggplot2} object. This allows further modification using standard ggplot syntax
#'   (e.g., \code{plot + theme_dark()}).
#'
#' @details
#' \strong{Data Extraction:} The function iterates through the nested list structure of \code{batch_drc_results},
#' extracting raw data points and model predictions for valid fits.
#'
#' \strong{IC50 Extraction:} The function attempts to extract the IC50/EC50 parameter (usually 'e')
#' from the underlying \code{drc} model object. This value is displayed in the console and optionally in the legend.
#'
#' \strong{Faceting vs. Limits:} If \code{facet_scales} is set to "free" or "free_y", the global \code{y_limits}
#' argument will be ignored to allow each panel to scale independently.
#'
#' @examples
#' \dontrun{
#' # 1. Basic usage: Plot specific compound by name
#' plot_drc_batch(results, target_compound = "PAK3:IPA-3")
#'
#' # 2. Search fuzzy match and use JAMA colors
#' plot_drc_batch(results, target_compound = "IPA-3", colors = "jama")
#'
#' # 3. Plot by position (e.g., the 5th compound found)
#' plot_drc_batch(results, position = 5)
#'
#' # 4. Facet by compound (compare multiple compounds side-by-side)
#' plot_drc_batch(results,
#'                target_compound = "PAK3",
#'                facet_by = "compound",
#'                facet_ncol = 3)
#'
#' # 5. Show IC50 in legend and use Times New Roman
#' plot_drc_batch(results,
#'                target_compound = "TargetA:CompB",
#'                show_IC50 = TRUE,
#'                font_family = "serif")
#' }
#'
#' @seealso
#' \code{\link{batch_drc_analysis}} for generating batch DRC results
#' \code{\link{plot_biological_replicates}} for comparing biological replicates
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import scales
#' @importFrom RColorBrewer brewer.pal
#' @export





plot_drc_batch <- function(batch_drc_results,
                           target_compound = NULL,
                           position = NULL,
                           y_limits = NULL,
                           colors = NULL,
                           point_shapes = NULL,
                           show_error_bars = TRUE,
                           legend_position = "right",
                           show_legend = TRUE,
                           shape_by_duplicate = FALSE,
                           show_grid = FALSE,
                           font_family = "sans",
                           facet_by = NULL,
                           facet_ncol = NULL,
                           facet_scales = "fixed",
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
  # 1. SETUP
  # ============================================================================
  required_packages <- c("ggplot2", "scales", "dplyr", "tidyr", "tibble")
  missing_packages <- sapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  })

  if (any(missing_packages)) {
    stop("The following packages are required: ", paste(required_packages[missing_packages], collapse = ", "))
  }

  suppressPackageStartupMessages({
    library(ggplot2); library(scales); library(dplyr); library(tidyr); library(tibble)
  })

  # Extração robusta do objeto de análise
  if (is.list(batch_drc_results) && "drc_results" %in% names(batch_drc_results)) {
    if (verbose) message("Detected analysis wrapper object. Extracting 'drc_results'...")
    batch_drc_results <- batch_drc_results$drc_results
  }

  if (length(batch_drc_results) == 0) stop("The DRC results list is empty.")
  if (verbose) message("Processing DRC batch results for plotting...")

  # ============================================================================
  # 2. HELPER FUNCTIONS
  # ============================================================================
  parse_compound_name <- function(name) {
    if (is.null(name)) return(list(target = "Unknown", compound = "Unknown"))
    clean_name <- trimws(gsub("\\.\\d+$", "", name))

    if (grepl(" \\| ", clean_name)) {
      parts <- strsplit(clean_name, " \\| ")[[1]]
      clean_name <- parts[1]
    }

    if (grepl(":", clean_name)) {
      parts <- strsplit(clean_name, ":")[[1]]
      return(list(target = trimws(parts[1]), compound = if(length(parts)>1) trimws(parts[2]) else trimws(parts[1])))
    } else {
      return(list(target = clean_name, compound = clean_name))
    }
  }

  generate_palette <- function(user_input, n_needed) {
    if (is.null(user_input)) return(scales::hue_pal()(n_needed))
    if (length(user_input) > 1) return(rep(user_input, length.out = n_needed))

    pal_name <- tolower(user_input)
    if (pal_name %in% c("viridis", "magma", "plasma", "inferno", "cividis", "turbo")) {
      if (requireNamespace("viridis", quietly = TRUE)) return(viridis::viridis_pal(option = pal_name)(n_needed))
    }
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      if (pal_name %in% tolower(rownames(RColorBrewer::brewer.pal.info))) {
        return(colorRampPalette(RColorBrewer::brewer.pal(8, pal_name))(n_needed))
      }
    }
    return(rep(user_input, length.out = n_needed))
  }

  get_detailed_results <- function(plate_obj) {
    if (is.null(plate_obj$drc_result)) return(NULL)
    res <- plate_obj$drc_result$detailed_results
    if (is.null(res)) res <- plate_obj$drc_result$curve_results
    if (is.null(res)) res <- plate_obj$drc_result$fits
    if (is.list(res) && length(res) > 0) return(res)
    return(NULL)
  }

  # ============================================================================
  # 3. DATA EXTRACTION (SEM IC50)
  # ============================================================================
  raw_data_list <- list(); curve_data_list <- list(); counter <- 0

  for (plate_name in names(batch_drc_results)) {
    plate_res_list <- get_detailed_results(batch_drc_results[[plate_name]])
    if (is.null(plate_res_list)) next

    for (i in seq_along(plate_res_list)) {
      res <- plate_res_list[[i]]
      has_success <- !is.null(res$success) && isTRUE(res$success)
      has_data <- !is.null(res$data) && nrow(res$data) > 0

      if (!has_success || !has_data) next

      comp_name <- res$compound
      if(is.null(comp_name)) comp_name <- paste0("Cmpd_", i)

      info <- parse_compound_name(comp_name)
      tc_key <- paste(info$target, info$compound, sep = ":")
      unique_id <- paste(plate_name, i, sep = "_")
      counter <- counter + 1

      # (REMOVIDO: Extração de IC50)

      valid_data <- res$data %>%
        filter(is.finite(log_inhibitor) & is.finite(response)) %>%
        mutate(plate = plate_name, target = info$target, compound = info$compound,
               target_compound = tc_key, unique_id = unique_id)

      if (nrow(valid_data) >= 2) {
        raw_data_list[[counter]] <- valid_data
        model_obj <- res$model
        if (!is.null(model_obj)) {
          try({
            x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
            x_span <- diff(x_range)
            x_seq <- seq(x_range[1] - 0.1*x_span, x_range[2] + 0.1*x_span, length.out = 100)
            pred_y <- predict(model_obj, newdata = data.frame(log_inhibitor = x_seq))
            curve_data_list[[counter]] <- data.frame(
              log_inhibitor = x_seq, response = pred_y, plate = plate_name,
              target = info$target, compound = info$compound, target_compound = tc_key, unique_id = unique_id
            )
          }, silent = TRUE)
        }
      }
    }
  }

  if (length(raw_data_list) == 0) stop("No valid data found.")

  df_raw_master <- bind_rows(raw_data_list)
  df_curve_master <- bind_rows(curve_data_list)
  # (REMOVIDO: df_ec50_master)

  # ============================================================================
  # 4. FILTERING
  # ============================================================================
  unique_combinations <- df_raw_master %>% select(target_compound, target, compound, unique_id) %>% distinct()
  selected_ids <- NULL; match_desc <- ""

  if (!is.null(position)) {
    unique_tc <- unique(unique_combinations$target_compound)
    if (position > length(unique_tc)) stop("Invalid position index.")
    target_tc <- unique_tc[position]
    selected_ids <- unique_combinations %>% filter(target_compound == target_tc) %>% pull(unique_id)
    match_desc <- paste("Pos:", position)

  } else if (!is.null(target_compound)) {
    search_term <- target_compound
    if (grepl(":", search_term)) {
      parts <- strsplit(search_term, ":")[[1]]
      target_query <- trimws(parts[1]); compound_query <- if(length(parts) > 1) trimws(parts[2]) else ""
      matches <- unique_combinations %>% filter(target == target_query, compound == compound_query)
      if (nrow(matches) == 0) matches <- unique_combinations %>% filter(tolower(target) == tolower(target_query), tolower(compound) == tolower(compound_query))
    } else {
      matches <- unique_combinations %>% filter(grepl(search_term, compound, ignore.case=TRUE) | grepl(search_term, target, ignore.case=TRUE))
    }
    if (nrow(matches) == 0) stop("No match found for: ", search_term)
    selected_ids <- matches$unique_id; match_desc <- search_term

  } else {
    first_tc <- unique_combinations$target_compound[1]
    selected_ids <- unique_combinations %>% filter(target_compound == first_tc) %>% pull(unique_id)
    match_desc <- "First compound"
  }

  plot_raw <- df_raw_master %>% filter(unique_id %in% selected_ids)
  plot_curve <- df_curve_master %>% filter(unique_id %in% selected_ids)
  # (REMOVIDO: selected_ec50s)

  # ============================================================================
  # 5. DATA AGGREGATION & LEGEND
  # ============================================================================
  plot_stats <- plot_raw %>%
    group_by(plate, target, compound, target_compound, unique_id, log_inhibitor) %>%
    summarise(mean_response = mean(response, na.rm=T),
              sd_response = sd(response, na.rm=T), .groups="drop")

  n_groups <- length(unique(plot_raw$target_compound))
  n_distinct_targets <- length(unique(plot_raw$target))

  plate_levels <- unique(plot_raw$plate)
  nums <- suppressWarnings(as.numeric(gsub("\\D", "", plate_levels)))
  if (!any(is.na(nums)) && length(nums) == length(plate_levels)) plate_levels <- plate_levels[order(nums)] else plate_levels <- sort(plate_levels)

  create_cols <- function(df) {
    # (REMOVIDO: left_join com ec50)
    df %>% mutate(
      plate = factor(plate, levels = plate_levels),
      dup_lbl = ifelse(grepl("\\d", plate), paste("Plate", gsub("\\D", "", plate)), as.character(plate)),

      base_label = if(n_groups > 1) {
        if (n_distinct_targets == 1) compound else target_compound
      } else {
        dup_lbl
      },

      legend_group = factor(base_label)
    )
  }

  plot_curve <- create_cols(plot_curve)
  plot_stats <- create_cols(plot_stats)

  # ============================================================================
  # 6. GGPLOT
  # ============================================================================
  if (is.null(plot_title)) {
    if (n_groups == 1) {
      comp_name <- unique(plot_raw$compound)[1]; targ_name <- unique(plot_raw$target)[1]
      if (targ_name == "Unknown" || targ_name == "" || targ_name == comp_name) final_title <- comp_name
      else final_title <- paste(targ_name, comp_name, sep = " : ")
    } else {
      if (!is.null(facet_by)) final_title <- "Comparative DRC Analysis"
      else final_title <- match_desc
    }
  } else { final_title <- plot_title }

  n_legend_items <- length(unique(plot_curve$legend_group))
  final_colors <- generate_palette(colors, n_legend_items)

  use_shapes <- shape_by_duplicate && n_groups == 1
  if (use_shapes) {
    if (is.null(point_shapes)) final_shapes <- 15:(15 + n_legend_items - 1) else final_shapes <- rep(point_shapes, length.out = n_legend_items)
    names(final_shapes) <- levels(plot_curve$legend_group)
  }

  yt <- if(is.null(y_axis_title)) "Response" else y_axis_title
  lt <- if(!is.null(legend_title)) legend_title else (if(n_groups > 1) "Compound" else "Replicate")

  p <- ggplot() +
    geom_line(data = plot_curve,
              aes(x = log_inhibitor, y = response, color = legend_group, group = unique_id),
              linewidth = 1, alpha = 0.8) +
    geom_point(data = plot_stats,
               aes(x = log_inhibitor, y = mean_response, color = legend_group,
                   shape = if(use_shapes) legend_group else NULL),
               size = 3)

  if (show_error_bars) {
    p <- p + geom_errorbar(data = plot_stats,
                           aes(x = log_inhibitor, ymin = mean_response - sd_response, ymax = mean_response + sd_response, color = legend_group),
                           width = 0.05, linewidth = 0.5, alpha = 0.6)
  }

  p <- p + scale_color_manual(values = final_colors, name = lt)

  if (use_shapes) p <- p + scale_shape_manual(values = final_shapes, name = lt)
  else p <- p + guides(shape = "none")

  if (!is.null(facet_by)) {
    if (!facet_by %in% names(plot_curve)) warning("Facet column '", facet_by, "' not found.")
    else p <- p + facet_wrap(as.formula(paste("~", facet_by)), ncol = facet_ncol, scales = facet_scales)
  }

  base_theme <- if(show_grid) theme_minimal(base_family = font_family) else theme_classic(base_family = font_family)

  p <- p +
    labs(title = final_title, x = expression(paste("Log"[10], " Concentration [M]")), y = yt) +
    base_theme +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = element_text(color = axis_text_color, size = axis_text_size),
      axis.title = element_text(color = axis_title_color, size = axis_title_size, face = "bold"),
      legend.position = if(show_legend) legend_position else "none",
      legend.title = element_text(face="bold"),
      legend.text = element_text(size = 10),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text = element_text(face = "bold", size = 11)
    )

  if (!show_grid) p <- p + theme(panel.grid = element_blank())

  if (!is.null(y_limits)) {
    if (!is.null(facet_by) && grepl("free", facet_scales)) {
      if(verbose) message("Warning: 'y_limits' ignored due to free facet scales.")
    } else {
      p <- p + coord_cartesian(ylim = y_limits)
    }
  }

  # ============================================================================
  # 7. OUTPUT
  # ============================================================================
  print(p)

  if (!is.null(save_plot)) {
    ggsave(save_plot, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi)
    if(verbose) message("Plot saved to: ", save_plot)
  }

  return(invisible(p))
}

