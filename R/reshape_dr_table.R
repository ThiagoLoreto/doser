reshape_dr_table <- function(results_table, output_file = NULL, decimal_comma = FALSE) {
  # Check if table has expected structure
  required_cols <- c("Compound", "Bottom", "Top", "LogIC50", "IC50", "Span", 
                     "R_squared", "Syx", "Max_Slope", "Curve_Quality", "Degrees_of_Freedom")
  
  if (!all(required_cols %in% colnames(results_table))) {
    stop("Table does not have the expected dose-response analysis structure")
  }
  
  # Transpose table - compounds become columns
  transposed <- as.data.frame(t(results_table[, -1]))
  colnames(transposed) <- results_table$Compound
  
  # Define all sections and parameters in order
  sections <- list(
    list(header = "log(inhibitor) vs. response (three parameters)", params = NULL),
    list(header = "Best-fit values", params = c("Bottom", "Top", "LogIC50", "IC50", "Span")),
    list(header = "Lower 95% conf. limit (profile likelihood)", 
         params = c("Bottom_Lower_95CI", "Top_Lower_95CI", "LogIC50_Lower_95CI", "IC50_Lower_95CI")),
    list(header = "Upper 95% conf. limit (profile likelihood)", 
         params = c("Bottom_Upper_95CI", "Top_Upper_95CI", "LogIC50_Upper_95CI", "IC50_Upper_95CI")),
    list(header = "Goodness of Fit", 
         params = c("Degrees_of_Freedom", "R_squared", "Sum_of_Squares", "Syx")),
    list(header = NULL, params = c("Max_Slope", "Curve_Quality"))
  )
  
  # Build final table structure
  final_table <- data.frame(matrix(NA, nrow = 0, ncol = ncol(transposed)))
  colnames(final_table) <- colnames(transposed)
  
  # Populate table with sections and parameters
  for (section in sections) {
    # Add section header if it exists
    if (!is.null(section$header)) {
      header_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(transposed)))
      colnames(header_row) <- colnames(transposed)
      rownames(header_row) <- section$header
      final_table <- rbind(final_table, header_row)
    }
    
    # Add section parameters
    for (param in section$params) {
      if (param %in% rownames(transposed)) {
        row_data <- transposed[param, , drop = FALSE]
        final_table <- rbind(final_table, row_data)
      }
    }
  }
  
  # Save to Excel if output_file specified
  if (!is.null(output_file)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("The 'openxlsx' package is required to save Excel files. Install using: install.packages('openxlsx')")
    }
    
    # Prepare table for Excel export
    excel_table <- final_table
    excel_table$Parameter <- rownames(excel_table)
    excel_table <- excel_table[, c("Parameter", colnames(final_table))]
    
    # Convert decimal points to commas if requested
    if (decimal_comma) {
      for (col in 2:ncol(excel_table)) {
        excel_table[[col]] <- as.character(excel_table[[col]])
        excel_table[[col]] <- gsub("\\.", ",", excel_table[[col]])
        excel_table[[col]][is.na(excel_table[[col]])] <- NA
      }
    }
    
    # Save to Excel
    openxlsx::write.xlsx(excel_table, output_file, rowNames = FALSE)
    cat("Table saved to:", output_file, "\n")
  }
  
  return(final_table)
}
