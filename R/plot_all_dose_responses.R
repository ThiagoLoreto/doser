plot_all_dose_responses <- function(results, 
                                    compounds = "all",
                                    output_dir = "dose_response_plots",
                                    file_prefix = "dose_response",
                                    file_extension = "png",
                                    ...) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine which compounds to plot
  n_compounds <- length(results$detailed_results)
  
  if (identical(compounds, "all")) {
    compounds_to_plot <- 1:n_compounds
  } else if (is.numeric(compounds)) {
    # Validate compound indices
    invalid_compounds <- compounds[compounds < 1 | compounds > n_compounds]
    if (length(invalid_compounds) > 0) {
      stop("Invalid compound indices: ", paste(invalid_compounds, collapse = ", "))
    }
    compounds_to_plot <- compounds
  } else {
    stop("compounds must be 'all' or a numeric vector of indices")
  }
  
  cat("Generating", length(compounds_to_plot), "dose-response plots...\n")
  
  # Generate plots for each compound
  for (i in compounds_to_plot) {
    compound_name <- results$detailed_results[[i]]$compound
    compound_name_parts <- strsplit(compound_name, " \\| ")[[1]]
    compound_name_display <- compound_name_parts[1]
    
    # Create safe filename
    safe_filename <- gsub("[^a-zA-Z0-9_-]", "_", compound_name_display)
    filename <- file.path(output_dir, 
                          paste0(file_prefix, "_", safe_filename, ".", file_extension))
    
    # Generate and save plot
    plot_dose_response(results = results, 
                       compound_index = i, 
                       save_plot = filename,
                       ...)
    
    cat("  - Plot", i, "of", length(compounds_to_plot), ":", compound_name_display, "\n")
  }
  
  cat("All plots saved in:", output_dir, "\n")
  return(invisible(compounds_to_plot))
}