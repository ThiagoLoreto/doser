#' Fit a 3-Parameter Logistic Dose-Response Model
#'
#' The `fit_drc_3pl()` function fits a three-parameter logistic (3PL) dose-response
#' model to experimental data. It supports duplicate measurements, normalization,
#' and optional filtering of low bottom values. The function computes fitted parameters,
#' model statistics, and diagnostic metrics, providing a concise summary of model performance.
#'
#' @param data A data frame containing dose-response data. The first column should
#'   contain log10(inhibitor concentration) values, followed by pairs of response
#'   columns (duplicate measurements).
#' @param output_file Optional character string specifying the output file path.
#'   Supports both Excel (.xlsx) and CSV formats.
#' @param normalize Logical indicating whether to normalize response data.
#'   If TRUE, responses are normalized from 0-100% based on the first and last values.
#'   Default is FALSE.
#' @param verbose Logical indicating whether to display progress messages and
#'   summary statistics. Default is TRUE
#' @param enforce_bottom_threshold Logical indicating whether to exclude IC50 values
#'   for curves where the Bottom parameter exceeds a threshold. Default is FALSE.
#' @param bottom_threshold Numeric value specifying the Bottom threshold when
#'   `enforce_bottom_threshold = TRUE`. Curves with Bottom <= this value will have
#'   IC50 values set to NA. Default is 60.
#' @param r_squared_threshold Numeric value between 0 and 1 specifying the R2 cutoff
#'   for curve quality assessment. Curves with R2 below this value will be flagged
#'   as "Low R2" in the quality assessment (default: 0.8).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{summary_table}: Data frame with fitted parameters for all compounds
#'   \item \code{detailed_results}: List containing detailed fitting results for each compound
#'   \item \code{n_compounds}: Number of compounds analyzed
#'   \item \code{successful_fits}: Number of successful curve fits
#'   \item \code{normalized}: Logical indicating if data was normalized
#'   \item \code{original_data}: Original input data
#'   \item \code{processed_data}: Processed data (normalized if requested)
#'   \item \code{threshold_settings}: List with threshold settings and affected compounds (if applicable)
#' }
#'
#' @details
#' This function performs comprehensive dose-response analysis using the following steps:
#'
#' \strong{Data Structure:}
#' - First column: log10(inhibitor concentration) values
#' - Subsequent columns: Pairs of response measurements (technical duplicates)
#' - Example: Column1 = log_conc, Column2 = Response1, Column3 = Response2, Column4 = Response3, Column5 = Response4, etc.
#'
#' \strong{Model Fitting:}
#' - Uses 3-parameter logistic model with fixed Hill slope of -1
#' - Implements multiple starting value strategies for robust convergence
#' - Includes biological plausibility checks and automatic parameter corrections
#' - Calculates ideal Hill slope using 4-parameter model for comparison
#'
#' \strong{Output Parameters:}
#' - Bottom: Minimum response asymptote
#' - Top: Maximum response asymptote
#' - LogIC50: log10(IC50) value
#' - IC50: Half-maximal inhibitory concentration
#' - Span: Top - Bottom (response range)
#' - Confidence intervals for all parameters
#' - Goodness-of-fit metrics (R2, Syx, etc.)
#' - Curve quality assessment
#' - Maximum slope and ideal Hill slope
#'
#' @examples
#' \dontrun{
#' # Basic usage with example data
#' results <- fit_dose_response(dose_response_data)
#'
#' # With normalization and output file
#' results <- fit_dose_response(
#'   data = my_data,
#'   output_file = "results.xlsx",
#'   normalize = TRUE
#' )
#'
#' # With bottom threshold enforcement
#' results <- fit_dose_response(
#'   data = my_data,
#'   enforce_bottom_threshold = TRUE,
#'   bottom_threshold = 60
#' )
#'
#' # Access results
#' summary_table <- results$summary_table
#' successful_fits <- results$successful_fits
#'
#' # View curve quality distribution
#' table(results$summary_table$Curve_Quality)
#' }
#'
#' @section Data Format:
#' The input data should be structured as follows:
#' \preformatted{
#'   log_conc  Targ/Comp1_rep1  Targ/Comp1_rep2  Targ/Comp2_rep1  Targ/Comp2_rep2 ...
#'     NA      25.1              25.1            25.1            25.1
#'   -8.0      30.5              31.2            22.1            23.4
#'   -7.0      45.2              46.8            35.6            36.9
#'   ...       ...            ...             ...             ...
#' }
#'
#' @section Curve Quality Assessment:
#' Curves are automatically assessed for quality based on:
#' \itemize{
#'   \item Maximum slope (shallow curves flagged)
#'   \item Response span (<20% flagged)
#'   \item R-squared values (<0.5 flagged)
#'   \item Biological plausibility of parameters
#'   \item logIC50 range (>0.666 flagged)
#' }
#'
#' @export
#' @seealso
#' Useful related functions:
#' \itemize{
#'   \item [plot_dose_response()] for plotting results
#'   \item plot_multiple_compounds() for exporting results
#'   \item save_multiple_sheets() for saving results to Excel
#'   \item \code{\link[stats]{nls}} for nonlinear regression details
#' }
#'
#' @export
#'
#' @references
#' For methodological details on dose-response analysis:
#' \itemize{
#'   \item GraphPad Prism Curve Fitting Guide
#'   \item NIH Assay Guidance Manual
#' }
#'
#' @importFrom stats median predict sd setNames
#' @importFrom dplyr %>% mutate across
#' @export

fit_drc_3pl <- function(data, output_file = NULL, normalize = FALSE, verbose = TRUE,
                        enforce_bottom_threshold = FALSE, bottom_threshold = 60,
                        r_sqr_threshold = 0.8) {

  # Check and load required packages
  required_packages <- c("dplyr", "stats", "graphics", "grDevices")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Required packages are not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install using: install.packages(c(",
         paste0("\"", missing_packages, "\"", collapse = ", "), "))")
  }

  library(dplyr, warn.conflicts = FALSE)

  # Constants
  PARAM_NAMES <- c("Bottom", "Top", "LogIC50", "IC50", "Span")

  # NULL coalescing operator
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

  # 3PL model for inhibition (Hill Slope = -1)
  three_param_model_inhibition <- function(log_inhibitor, Bottom, Top, LogIC50) {
    Bottom + (Top - Bottom) / (1 + 10^(log_inhibitor - LogIC50))
  }

  # 3PL model for activation (Hill Slope = +1)
  three_param_model_activation <- function(log_inhibitor, Bottom, Top, LogIC50) {
    Bottom + (Top - Bottom) / (1 + 10^(LogIC50 - log_inhibitor))
  }

  # 4-parameter logistic model for ideal Hill Slope calculation
  four_param_model <- function(log_inhibitor, Bottom, Top, LogIC50, HillSlope) {
    Bottom + (Top - Bottom) / (1 + 10^((log_inhibitor - LogIC50) * HillSlope))
  }

  # Detect curve type based on data pattern
  detect_curve_type <- function(data) {
    responses <- data$response[!is.na(data$response)]
    log_inhibitor <- data$log_inhibitor[!is.na(data$response)]

    if (length(responses) < 4) return("unknown")

    sorted_data <- data[order(data$log_inhibitor), ]
    sorted_responses <- sorted_data$response[!is.na(sorted_data$response)]

    if (length(sorted_responses) < 4) return("unknown")

    initial_avg <- mean(head(sorted_responses, 3), na.rm = TRUE)
    final_avg <- mean(tail(sorted_responses, 3), na.rm = TRUE)

    if (initial_avg > final_avg + 15) {
      return("inhibition")
    } else if (final_avg > initial_avg + 15) {
      return("activation")
    } else {
      return("flat")
    }
  }

  # Get appropriate model based on curve type
  get_model_for_fitting <- function(data) {
    curve_type <- detect_curve_type(data)
    if (curve_type == "activation") {
      return(three_param_model_activation)
    } else {
      return(three_param_model_inhibition)
    }
  }

  # Correct Hill Slope sign based on curve type
  correct_hill_slope <- function(hill_slope, data) {
    if (is.na(hill_slope)) return(NA)

    curve_type <- detect_curve_type(data)

    if (curve_type == "inhibition" && hill_slope > 0) {
      return(-abs(hill_slope))
    } else if (curve_type == "activation" && hill_slope < 0) {
      return(abs(hill_slope))
    } else {
      return(hill_slope)
    }
  }

  # Create empty result structure for failed fits
  create_empty_result <- function(comp_name, reason = "Model failed") {
    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = rep(NA, 5), stringsAsFactors = FALSE),
      confidence_intervals = list(
        Bottom = c(NA, NA), Top = c(NA, NA), LogIC50 = c(NA, NA), IC50 = c(NA, NA),
        Bottom_Lower = NA, Bottom_Upper = NA, Top_Lower = NA, Top_Upper = NA
      ),
      goodness_of_fit = list(
        R_squared = NA, Syx = NA, Sum_of_Squares = NA,
        Total_Sum_of_Squares = NA, Regression_Sum_of_Squares = NA, Degrees_of_Freedom = NA
      ),
      curve_quality = paste("Fit", tolower(reason)),
      max_slope = NA,
      ideal_hill_slope = NA,
      model = NULL,
      success = FALSE,
      compound = comp_name
    )
  }

  # Recalculate dependent parameters after corrections
  recalculate_dependent_params <- function(params, corrections) {
    corrected_params <- params

    if ("Bottom" %in% names(corrections)) corrected_params[1] <- corrections$Bottom
    if ("Top" %in% names(corrections)) corrected_params[2] <- corrections$Top
    if ("LogIC50" %in% names(corrections)) corrected_params[3] <- corrections$LogIC50

    c(corrected_params[1:3], 10^corrected_params[3], corrected_params[2] - corrected_params[1])
  }

  # Recalculate confidence intervals after parameter corrections
  recalculate_ci <- function(original_ci, corrections) {
    ci <- original_ci

    if ("Bottom" %in% names(corrections)) {
      ci$Bottom <- ci$Bottom_Lower <- ci$Bottom_Upper <- NA
    }
    if ("Top" %in% names(corrections)) {
      ci$Top <- ci$Top_Lower <- ci$Top_Upper <- NA
    }

    if ("Bottom" %in% names(corrections) || "Top" %in% names(corrections) || "LogIC50" %in% names(corrections)) {
      ci$LogIC50 <- ci$IC50 <- c(NA, NA)
    } else if (!any(is.na(ci$LogIC50))) {
      ci$IC50 <- 10^ci$LogIC50
    }

    ci
  }

  # Check biological plausibility and apply corrections if needed
  check_biological_plausibility <- function(params, data) {
    responses <- data$response[!is.na(data$response)]
    exp_min <- min(responses, na.rm = TRUE)
    exp_max <- max(responses, na.rm = TRUE)

    curve_type <- detect_curve_type(data)

    bottom_limits <- if (curve_type == "activation") c(-100, Inf) else c(-100, 600)
    top_limits <- c(0, 700)
    logIC50_limits <- c(-20, 5)

    corrections <- list()
    reasons <- list()

    if (params[1] < bottom_limits[1] || params[1] > bottom_limits[2] || !is.finite(params[1])) {
      corrections$Bottom <- max(bottom_limits[1], min(bottom_limits[2], exp_min))
      reasons$Bottom <- sprintf("Biologically implausible (%.2f). Using experimental minimum: %.2f",
                                params[1], exp_min)
    }

    if (params[2] < top_limits[1] || params[2] > top_limits[2] || !is.finite(params[2])) {
      corrections$Top <- max(top_limits[1], min(top_limits[2], exp_max))
      reasons$Top <- sprintf("Biologically implausible (%.2f). Using experimental maximum: %.2f",
                             params[2], exp_max)
    }

    if (params[3] < logIC50_limits[1] || params[3] > logIC50_limits[2] || !is.finite(params[3])) {
      fallback <- median(data$log_inhibitor, na.rm = TRUE)
      corrections$LogIC50 <- max(logIC50_limits[1], min(logIC50_limits[2], fallback))
      reasons$LogIC50 <- sprintf("Biologically implausible (%.2f). Using median concentration: %.2f",
                                 params[3], fallback)
    }

    if (length(corrections) > 0) {
      return(list(
        corrected_params = recalculate_dependent_params(params, corrections),
        corrections_applied = corrections,
        correction_reasons = reasons,
        needs_correction = TRUE
      ))
    }

    list(needs_correction = FALSE)
  }

  # Prepare and validate data for analysis
  prepare_and_validate_data <- function(pair_data) {
    log_inhibitor <- pair_data[, 1]
    response <- as.numeric(unlist(pair_data[, -1]))

    df <- data.frame(log_inhibitor = rep(log_inhibitor, 2), response = response)
    df_clean <- stats::na.omit(df)

    if (nrow(df_clean) < 4 || stats::sd(df_clean$response, na.rm = TRUE) < 1e-6) {
      return(list(valid = FALSE, df_clean = df_clean))
    }

    min_resp <- min(df_clean$response, na.rm = TRUE)
    max_resp <- max(df_clean$response, na.rm = TRUE)

    approx_ic50 <- tryCatch({
      suppressWarnings(
        stats::approx(df_clean$response, df_clean$log_inhibitor,
                      xout = (min_resp + max_resp) / 2)$y
      )
    }, error = function(e) NULL) %||% stats::median(df_clean$log_inhibitor, na.rm = TRUE)

    list(valid = TRUE, df_clean = df_clean,
         min_response = min_resp, max_response = max_resp, approx_ic50 = approx_ic50)
  }

  # Robust nonlinear fitting with model selection
  try_robust_fit <- function(df_clean, min_resp, max_resp, approx_ic50) {
    if (nrow(df_clean) < 4) return(NULL)

    control_configs <- list(
      default = stats::nls.control(maxiter = 500, tol = 1e-04, minFactor = 1/4096, warnOnly = TRUE),
      relaxed = stats::nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024, warnOnly = TRUE)
    )

    # Determine which model to use based on curve type
    model_func <- get_model_for_fitting(df_clean)
    curve_type <- detect_curve_type(df_clean)

    if (verbose) {
      hill_direction <- if (curve_type == "activation") "+1" else "-1"
      cat("  - Detected curve type:", curve_type, "-> Using Hill Slope =", hill_direction, "\n")
    }

    start_strategies <- list(
      list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50),
      list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50),
      list(Bottom = min_resp * 1.2, Top = max_resp * 0.8, LogIC50 = approx_ic50),
      list(Bottom = max(0, min_resp - 10), Top = min(150, max_resp + 10), LogIC50 = approx_ic50)
    )

    for (start_vals in start_strategies) {
      fit <- tryCatch({
        stats::nls(
          response ~ model_func(log_inhibitor, Bottom, Top, LogIC50),
          data = df_clean, start = start_vals, control = control_configs$default, algorithm = "port"
        )
      }, error = function(e) NULL)

      if (!is.null(fit)) return(fit)
    }

    tryCatch({
      stats::nls(
        response ~ model_func(log_inhibitor, Bottom, Top, LogIC50),
        data = df_clean, start = start_strategies[[1]], control = control_configs$relaxed, algorithm = "port"
      )
    }, error = function(e) NULL)
  }

  # Calculate ideal Hill Slope using 4-parameter model
  calculate_ideal_hill_slope <- function(df_clean, initial_params) {
    if (nrow(df_clean) < 4) return(NA)

    tryCatch({
      start_vals <- list(
        Bottom = initial_params[1],
        Top = initial_params[2],
        LogIC50 = initial_params[3],
        HillSlope = if (detect_curve_type(df_clean) == "activation") 1 else -1
      )

      fit_4param <- stats::nls(
        response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
        data = df_clean, start = start_vals,
        control = stats::nls.control(maxiter = 200, tol = 1e-04, warnOnly = TRUE),
        algorithm = "port"
      )

      hill_slope <- stats::coef(fit_4param)["HillSlope"]
      hill_slope_value <- round(as.numeric(hill_slope), 3)

      corrected_hill_slope <- correct_hill_slope(hill_slope_value, df_clean)

      return(corrected_hill_slope)

    }, error = function(e) {
      return(NA)
    })
  }

  # Calculate goodness of fit metrics
  calculate_goodness_of_fit <- function(fit, df_clean) {
    tryCatch({
      residuals <- stats::resid(fit)
      mean_resp <- mean(df_clean$response, na.rm = TRUE)
      total_ss <- sum((df_clean$response - mean_resp)^2)
      residual_ss <- sum(residuals^2)
      regression_ss <- total_ss - residual_ss
      n_obs <- nrow(df_clean)
      dof <- max(1, n_obs - length(stats::coef(fit)))

      list(
        R_squared = if (total_ss > 0) regression_ss / total_ss else 0,
        Syx = sqrt(residual_ss / dof),
        Sum_of_Squares = residual_ss,
        Degrees_of_Freedom = dof
      )
    }, error = function(e) {
      list(R_squared = NA, Syx = NA, Sum_of_Squares = NA, Degrees_of_Freedom = NA)
    })
  }

  # Calculate confidence intervals safely
  calculate_ci <- function(fit) {
    safe_ci <- function(param_name) {
      tryCatch({
        prof <- stats::profile(fit, which = param_name, alpha = 0.05)
        unname(stats::confint(prof, level = 0.95))
      }, error = function(e) {
        tryCatch({
          if (requireNamespace("lmtest", quietly = TRUE)) {
            unname(lmtest::confint2(fit, level = 0.95)[param_name, ])
          } else {
            se <- summary(fit)$coefficients[param_name, "Std. Error"]
            est <- stats::coef(fit)[param_name]
            z <- stats::qnorm(0.975)
            c(est - z * se, est + z * se)
          }
        }, error = function(e2) c(NA, NA))
      })
    }

    bottom_ci <- safe_ci("Bottom")
    top_ci <- safe_ci("Top")
    logIC50_ci <- safe_ci("LogIC50")

    list(
      Bottom = bottom_ci, Top = top_ci, LogIC50 = logIC50_ci, IC50 = 10^logIC50_ci,
      Bottom_Lower = bottom_ci[1], Bottom_Upper = bottom_ci[2],
      Top_Lower = top_ci[1], Top_Upper = top_ci[2]
    )
  }

  # Calculate curve quality metrics
  calculate_curve_quality <- function(params, gof_results, plausibility_check = NULL, logIC50_ci = NULL) {
    tryCatch({
      span <- params[5]
      max_slope <- -span * log(10) / 4

      quality_flags <- character()
      if (abs(max_slope) < 5) quality_flags <- c(quality_flags, "Very shallow slope")
      else if (abs(max_slope) < 15) quality_flags <- c(quality_flags, "Shallow slope")
      if (abs(span) < 20) quality_flags <- c(quality_flags, "Small span")
      if (gof_results$R_squared < r_sqr_threshold) quality_flags <- c(quality_flags, "Low R²")
      if (!is.null(logIC50_ci) && !any(is.na(logIC50_ci))) {
        ci_range <- abs(logIC50_ci[2] - logIC50_ci[1])
        if (ci_range > 0.666) {
          quality_flags <- c(quality_flags, "Wide logIC50 CI range")
        }
      }

      if (!is.null(plausibility_check) && plausibility_check$needs_correction) {
        quality_flags <- c(quality_flags, "Parameters corrected (biologically implausible)")
      }

      list(
        quality = if (length(quality_flags) == 0) "Good curve" else paste(quality_flags, collapse = "; "),
        max_slope = max_slope
      )
    }, error = function(e) {
      list(quality = "Error in quality assessment", max_slope = NA)
    })
  }

  # Main analysis function for each compound pair
  analyze_single_pair <- function(pair_data, comp_name) {
    prepared <- prepare_and_validate_data(pair_data)
    if (!prepared$valid) {
      if (verbose) warning("Data validation failed for ", comp_name)
      return(create_empty_result(comp_name, "Validation failed"))
    }

    # Detect curve type before fitting
    curve_type <- detect_curve_type(prepared$df_clean)

    fit <- try_robust_fit(prepared$df_clean, prepared$min_response,
                          prepared$max_response, prepared$approx_ic50)

    if (is.null(fit)) {
      if (verbose) warning("Could not fit model for ", comp_name)
      result <- create_empty_result(comp_name, "Approximate (model failed)")
      result$parameters$Value <- c(prepared$min_response, prepared$max_response,
                                   prepared$approx_ic50, 10^prepared$approx_ic50,
                                   prepared$max_response - prepared$min_response)
      result$curve_quality <- "Model failed - approximate parameters"
      result$success = TRUE
      result$curve_type <- curve_type
      return(result)
    }

    params <- tryCatch(unname(stats::coef(fit)), error = function(e) NULL)
    if (is.null(params) || any(is.na(params))) {
      if (verbose) warning("Error extracting coefficients for ", comp_name)
      return(create_empty_result(comp_name, "Coefficient extraction failed"))
    }

    # Process results
    initial_params <- c(params, 10^params[3], params[2] - params[1])
    ci_results <- calculate_ci(fit)
    gof_results <- calculate_goodness_of_fit(fit, prepared$df_clean)
    plausibility_check <- check_biological_plausibility(params, prepared$df_clean)

    if (plausibility_check$needs_correction) {
      if (verbose) {
        cat("Biological plausibility correction for ", comp_name, ":\n", sep = "")
        for (param in names(plausibility_check$correction_reasons)) {
          cat("  - ", param, ": ", plausibility_check$correction_reasons[[param]], "\n", sep = "")
        }
      }
      final_params <- plausibility_check$corrected_params
      ci_results <- recalculate_ci(ci_results, plausibility_check$corrections_applied)
    } else {
      final_params <- initial_params
    }

    curve_quality_info <- calculate_curve_quality(final_params, gof_results, plausibility_check, ci_results$LogIC50)

    ideal_hill_slope <- calculate_ideal_hill_slope(prepared$df_clean, final_params[1:3])

    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = final_params, stringsAsFactors = FALSE),
      confidence_intervals = ci_results,
      goodness_of_fit = gof_results,
      curve_quality = curve_quality_info$quality,
      max_slope = curve_quality_info$max_slope,
      ideal_hill_slope = ideal_hill_slope,
      model = fit,
      success = TRUE,
      compound = comp_name,
      biological_plausibility_check = plausibility_check,
      data = prepared$df_clean,
      curve_type = curve_type
    )
  }

  # Data validation and normalization
  if (ncol(data) < 3) stop("Data must have at least 3 columns")
  if ((ncol(data) - 1) %% 2 != 0) stop("Number of response columns must be even")

  original_data <- data

  if (normalize) {
    if (verbose) cat("Applying normalization to data...\n")
    data <- data %>% mutate(across(-1, ~ {
      values <- as.numeric(as.character(.x))
      values_clean <- values[!is.na(values)]

      if (length(values_clean) < 2) return(rep(NA, length(values)))

      first_val <- values_clean[1]
      last_val <- values_clean[length(values_clean)]

      if (is.na(first_val) || is.na(last_val) || (last_val - first_val) == 0) {
        return(rep(NA, length(values)))
      }

      (values - first_val) / (last_val - first_val) * 100
    }))
  }

  # Main analysis loop
  if (verbose) cat("Starting dose-response analysis...\n")
  n_pairs <- (ncol(data) - 1) %/% 2

  all_results <- lapply(seq_len(n_pairs), function(pair) {
    if (verbose) cat("Processing duplicate pair", pair, "/", n_pairs, "...\n")

    col1 <- 1 + (pair * 2 - 1)
    col2 <- 1 + (pair * 2)
    col_names <- colnames(data)

    comp_name <- if (!is.null(col_names) && length(col_names) >= col2) {
      paste(col_names[col1], col_names[col2], sep = " | ")
    } else {
      paste("Compound", pair)
    }

    analyze_single_pair(data[, c(1, col1, col2)], comp_name)
  })

  # Create summary table
  summary_table <- do.call(rbind, lapply(all_results, function(result) {
    if (isTRUE(result$success) && !all(is.na(result$parameters$Value))) {
      params <- unname(result$parameters$Value)
      gof <- result$goodness_of_fit
      ci <- result$confidence_intervals

      # Apply threshold ONLY for inhibition curves
      apply_threshold <- FALSE
      if (enforce_bottom_threshold && !is.na(params[1]) && params[1] >= bottom_threshold) {
        curve_type <- result$curve_type %||% detect_curve_type(result$data)
        if (curve_type %in% c("inhibition", "flat")) {
          apply_threshold <- TRUE
        }
      }

      data.frame(
        Compound = strsplit(result$compound, " \\| ")[[1]][1],
        Bottom = round(params[1], 3),
        Top = round(params[2], 3),
        LogIC50 = if (!apply_threshold) round(params[3], 3) else NA,
        IC50 = if (!apply_threshold) format(params[4], scientific = TRUE) else NA,
        Bottom_Lower_95CI = if (!is.na(ci$Bottom_Lower)) sprintf("%.3f", ci$Bottom_Lower) else NA,
        Bottom_Upper_95CI = if (!is.na(ci$Bottom_Upper)) sprintf("%.3f", ci$Bottom_Upper) else NA,
        Top_Lower_95CI = if (!is.na(ci$Top_Lower)) sprintf("%.3f", ci$Top_Lower) else NA,
        Top_Upper_95CI = if (!is.na(ci$Top_Upper)) sprintf("%.3f", ci$Top_Upper) else NA,
        LogIC50_Lower_95CI = if (!apply_threshold && !is.na(ci$LogIC50[1])) round(ci$LogIC50[1], 3) else NA,
        LogIC50_Upper_95CI = if (!apply_threshold && !is.na(ci$LogIC50[2])) round(ci$LogIC50[2], 3) else NA,
        IC50_Lower_95CI = if (!apply_threshold && !is.na(ci$IC50[1])) format(ci$IC50[1], scientific = TRUE) else NA,
        IC50_Upper_95CI = if (!apply_threshold && !is.na(ci$IC50[2])) format(ci$IC50[2], scientific = TRUE) else NA,
        Span = round(params[5], 3),
        R_squared = round(gof$R_squared, 3),
        Syx = round(gof$Syx, 3),
        Sum_of_Squares = round(gof$Sum_of_Squares, 3),
        Degrees_of_Freedom = gof$Degrees_of_Freedom,
        Max_Slope = round(result$max_slope %||% NA, 3),
        Ideal_Hill_Slope = result$ideal_hill_slope %||% NA,
        Curve_Quality = result$curve_quality %||% "Not assessed",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Compound = strsplit(result$compound, " \\| ")[[1]][1],
        Bottom = NA, Top = NA, LogIC50 = NA, IC50 = NA,
        Bottom_Lower_95CI = NA, Bottom_Upper_95CI = NA,
        Top_Lower_95CI = NA, Top_Upper_95CI = NA,
        LogIC50_Lower_95CI = NA, LogIC50_Upper_95CI = NA,
        IC50_Lower_95CI = NA, IC50_Upper_95CI = NA, Span = NA,
        R_squared = NA, Syx = NA, Sum_of_Squares = NA, Degrees_of_Freedom = NA,
        Max_Slope = NA, Ideal_Hill_Slope = NA,
        Curve_Quality = "Fit failed",
        stringsAsFactors = FALSE
      )
    }
  }))

  # Create final_summary_table (transposed version)
  if (nrow(summary_table) > 0) {
    compound_names <- summary_table$Compound
    transposed_data <- as.data.frame(t(summary_table[, -1]))
    colnames(transposed_data) <- compound_names
    final_summary_table <- transposed_data

  } else {
    final_summary_table <- data.frame()
  }

  # Identify compounds affected by threshold ONLY for inhibition
  threshold_affected <- character()
  if (enforce_bottom_threshold) {
    for (result in all_results) {
      if (isTRUE(result$success) &&
          !is.na(result$parameters$Value[1]) &&
          result$parameters$Value[1] >= bottom_threshold &&
          result$curve_type == "inhibition") {
        comp_name <- strsplit(result$compound, " \\| ")[[1]][1]
        threshold_affected <- c(threshold_affected, comp_name)
      }
    }
  }

  # Print summary statistics
  if (verbose) {
    successful <- sum(!is.na(summary_table$IC50))
    total <- nrow(summary_table)
    success_rate <- round(successful / total * 100, 1)
    threshold_count <- length(threshold_affected)

    cat("\n", strrep("=", 50), "\n", sep = "")
    cat("DOSE-RESPONSE ANALYSIS COMPLETED SUCCESSFULLY!\n")
    cat(strrep("=", 50), "\n")
    cat("SUMMARY STATISTICS:\n")
    cat("  • Compounds analyzed: ", total, "\n")
    cat("  • Successful fits: ", successful, " (", success_rate, "%)\n", sep = "")
    cat("  • Failed fits: ", total - successful, "\n")

    if (enforce_bottom_threshold && threshold_count > 0) {
      cat("  • IC50 values excluded (Bottom ≥", bottom_threshold, "): ", threshold_count, "\n", sep = "")

      cat("\nCOMPOUNDS WITH EXCLUDED IC50 VALUES:\n")
      for (i in seq_along(threshold_affected)) {
        bottom_val <- summary_table$Bottom[summary_table$Compound == threshold_affected[i]]
        cat("  • ", threshold_affected[i], " (Bottom = ",
            round(bottom_val, 1), ")\n", sep = "")
      }
    }

    cat("  • Problematic curves: ", sum(grepl("shallow|Low R²|Small span", summary_table$Curve_Quality)), "\n")

    if ("Curve_Quality" %in% names(summary_table)) {
      cat("\nCURVE QUALITY DISTRIBUTION:\n")
      quality_counts <- table(summary_table$Curve_Quality)
      for (quality in names(sort(quality_counts, decreasing = TRUE))) {
        count <- quality_counts[quality]
        cat("  • ", sprintf("%-35s: %d (%.1f%%)", quality, count, count/total*100), "\n")
      }
    }
    cat(strrep("=", 50), "\n\n")
  }

  # Save results to file
  if (!is.null(output_file)) {
    file_ext <- tolower(tools::file_ext(output_file))

    if (file_ext == "xlsx" && requireNamespace("openxlsx", quietly = TRUE)) {
      # Criar workbook manualmente para controlar melhor a formatação
      wb <- openxlsx::createWorkbook()

      # Sheet 1: Summary (não tem row names, então podemos usar writeData normalmente)
      openxlsx::addWorksheet(wb, "Summary")
      openxlsx::writeData(wb, "Summary", summary_table, rowNames = FALSE)

      # Sheet 2: Final_Summary (PRECISA preservar row names)
      openxlsx::addWorksheet(wb, "Final_Summary")

      if (!is.null(final_summary_table) && nrow(final_summary_table) > 0) {
        # Converter row names para uma coluna chamada "Parameter"
        final_summary_with_rownames <- data.frame(
          Parameter = rownames(final_summary_table),
          final_summary_table,
          row.names = NULL,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        openxlsx::writeData(wb, "Final_Summary", final_summary_with_rownames, rowNames = FALSE)
      } else {
        # Se a tabela estiver vazia
        openxlsx::writeData(wb, "Final_Summary",
                            data.frame(Note = "No final summary data available"),
                            rowNames = FALSE)
      }

      # Salvar o workbook
      openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

      if (verbose) {
        message("Results saved to Excel with two sheets: ", output_file)
      }
    } else {
      # Fallback para CSV se não tiver openxlsx
      if (file_ext == "xlsx") {
        output_file <- sub("\\.xlsx$", ".csv", output_file, ignore.case = TRUE)
        if (verbose) warning("openxlsx not available. Falling back to CSV format...")
      }
      utils::write.csv(summary_table, output_file, row.names = FALSE)
      if (verbose) message("Results saved to: ", output_file)
    }
  }

  # Return results object
  list(
    summary_table = summary_table,
    final_summary_table = final_summary_table,
    detailed_results = all_results,
    n_compounds = n_pairs,
    successful_fits = sum(!is.na(summary_table$IC50)),
    normalized = normalize,
    original_data = original_data,
    processed_data = data,
    threshold_settings = if (enforce_bottom_threshold) {
      list(
        enforce_bottom_threshold = enforce_bottom_threshold,
        bottom_threshold = bottom_threshold,
        affected_compounds = threshold_affected
      )
    } else {
      NULL
    }
  )
}
