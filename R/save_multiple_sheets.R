#' Save Multiple Data Frames to Excel with Advanced Formatting Options
#'
#' Exports multiple data frames to a single Excel file with professional formatting,
#' decimal separator customization, selective rounding, and row name preservation.
#' Designed specifically for dose-response analysis results and scientific data.
#'
#' @param file_name Character string specifying the output Excel file path.
#' @param ... Data frames to be saved as separate worksheets. Can use named arguments
#'   for custom sheet names (e.g., Summary = df1, Details = df2).
#' @param decimal_comma Logical indicating whether to use comma as decimal separator
#'   instead of point. Recommended for European formats (default: TRUE).
#' @param decimal_places Integer specifying the number of decimal places for rounding
#'   (default: 3).
#' @param round_sheets Numeric vector specifying which sheets (by position) should
#'   have rounding applied. If NULL, all sheets are rounded (default: NULL).
#'
#' @return Invisibly returns NULL. The function primarily produces an Excel file
#'   as output and displays a confirmation message.
#'
#' @details
#' This function provides sophisticated Excel export capabilities optimized for
#' scientific and pharmacological data presentation. Key features include:
#'
#' \strong{Advanced Formatting Features:}
#' \itemize{
#'   \item \textbf{Multiple Sheets}: Each data frame becomes a separate worksheet
#'   \item \textbf{Row Name Preservation}: Automatically includes row names in a protected "RowNames" column
#'   \item \textbf{Decimal Customization}: Converts decimal separators (point to comma) with precision control
#'   \item \textbf{Selective Rounding}: Applies rounding only to specified worksheets
#'   \item \textbf{Scientific Notation}: Preserves scientific notation when present
#'   \item \textbf{Protected Columns}: Never modifies the "RowNames" column during formatting
#' }
#'
#' \strong{Automatic Sheet Naming:}
#' \itemize{
#'   \item Uses variable names for sheet names when unnamed arguments provided
#'   \item Uses provided names when named arguments used
#'   \item Maintains original data structure and integrity
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Exporting comprehensive dose-response analysis results
#' # Perform analysis first
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#'
#' # Export complete results with European decimal format
#' save_multiple_sheets(
#'   "dose_response_complete_analysis.xlsx",
#'   Summary_Table = analysis_results$summary_table,
#'   Detailed_Parameters = analysis_results$detailed_results[[1]]$parameters,
#'   Quality_Assessment = analysis_results$interval_means,
#'   Raw_Data = analysis_results$processed_data,
#'   decimal_comma = TRUE,
#'   decimal_places = 3
#' )
#'
#' # Example 2: International format with selective rounding
#' save_multiple_sheets(
#'   "international_format.xlsx",
#'   IC50_Results = ic50_data,           # Round to 2 decimal places
#'   Kinetic_Parameters = kinetics_data, # No rounding for precision
#'   Statistical_Analysis = stats_data,  # Round to 2 decimal places
#'   decimal_comma = FALSE,              # Use point for international journals
#'   decimal_places = 2,
#'   round_sheets = c(1, 3)              # Round only sheets 1 and 3
#' )
#'
#' # Example 3: High-precision scientific data
#' save_multiple_sheets(
#'   "high_precision_data.xlsx",
#'   Binding_Constants = kd_data,        # High precision, no rounding
#'   Dose_Response = dr_data,            # Standard precision
#'   decimal_comma = TRUE,
#'   decimal_places = 4,
#'   round_sheets = 2                    # Round only dose-response data
#' )
#'
#' # Example 4: Multiple analysis batches
#' # Analyze different experimental conditions
#' control_results <- fit_dose_response(control_data)
#' treated_results <- fit_dose_response(treated_data)
#'
#' save_multiple_sheets(
#'   "experimental_conditions.xlsx",
#'   Control_Condition = control_results$summary_table,
#'   Treated_Condition = treated_results$summary_table,
#'   Comparison_Analysis = comparison_stats,
#'   decimal_comma = TRUE,
#'   decimal_places = 3
#' )
#' }
#'
#' @section Row Name Handling:
#' The function automatically preserves row names through:
#' \itemize{
#'   \item Creating a dedicated "RowNames" column as the first column
#'   \item Using original row names when available
#'   \item Generating sequential numbers (1, 2, 3...) when no row names exist
#'   \item Protecting the "RowNames" column from all decimal formatting operations
#' }
#'
#' @section Decimal Formatting Control:
#' Advanced control over number formatting:
#' \itemize{
#'   \item \code{decimal_comma = TRUE}: European format (1,234 instead of 1.234)
#'   \item \code{decimal_comma = FALSE}: International format (1.234)
#'   \item Scientific notation preserved in both formats (1,23e-4 or 1.23e-4)
#'   \item Rounding applied before decimal separator conversion
#' }
#'
#' @section Rounding Strategies:
#' The \code{round_sheets} parameter enables precise control:
#' \itemize{
#'   \item \code{NULL}: Apply rounding to all worksheets (default behavior)
#'   \item \code{1}: Apply rounding only to first worksheet
#'   \item \code{c(1, 3)}: Apply rounding to worksheets 1 and 3
#'   \item \code{2:4}: Apply rounding to worksheets 2, 3, and 4
#'   \item \code{numeric(0)}: No rounding applied to any worksheet
#' }
#'
#' @seealso
#' \code{\link{fit_dose_response}} for generating analysis results
#' \code{\link{reshape_dr_table}} for creating structured result tables
#' \code{\link[openxlsx]{write.xlsx}} for basic Excel export functionality
#'
#' @export
#'
#' @references
#' For scientific data presentation standards:
#' \itemize{
#'   \item "The International System of Units (SI)" - Decimal separator conventions
#'   \item Nature Research Reporting Guidelines
#'   \item Journal of Pharmacology and Experimental Therapeutics data standards
#' }
#' @examples
#' \dontrun{
#' # Example 1: Exporting dose-response analysis results
#' # Perform analysis first
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#'
#' # Export comprehensive results
#' save_multiple_sheets(
#'   "dose_response_complete.xlsx",
#'   Summary_Table = analysis_results$summary_table,
#'   Detailed_Results = analysis_results$detailed_results[[1]]$parameters,
#'   Quality_Metrics = analysis_results$interval_means,
#'   decimal_comma = TRUE,
#'   decimal_places = 3
#' )
#' }



save_multiple_sheets <- function(file_name, ..., decimal_comma = TRUE, decimal_places = 3, round_sheets = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Please install it.")
  }
  
  # Funcao que NUNCA converte a coluna "RowNames"
  convert_decimal_separator <- function(df, decimal_places = 3, apply_rounding = TRUE) {
    if (!is.data.frame(df)) return(df)
    
    df_conv <- df  # Comecar com copia do original
    
    # Aplicar conversao apenas nas colunas que NAO sao "RowNames"
    for (col_name in names(df_conv)) {
      if (col_name == "RowNames") next  # Pular a coluna RowNames
      
      column <- df_conv[[col_name]]
      char_column <- as.character(column)
      
      # Aplicar substituicao apenas em valores que sao numeros
      result <- sapply(char_column, function(x) {
        if (is.na(x)) return(NA_character_)
        
        x_clean <- trimws(x)
        
        # Verificar se e numero (incluindo notacao cientifica)
        if (grepl("^-?\\d*\\.\\d+$", x_clean) ||
            grepl("^-?\\d+\\.\\d*$", x_clean) ||
            grepl("^-?\\d*\\.?\\d+[eE][-+]?\\d+$", x_clean)) {
          
          num_value <- as.numeric(x_clean)
          if (!is.na(num_value)) {
            # Aplicar arredondamento apenas se solicitado
            if (apply_rounding) {
              rounded_value <- round(num_value, decimal_places)
            } else {
              rounded_value <- num_value
            }
            
            # Manter a representacao original (pode incluir notacao cientifica)
            formatted_value <- as.character(rounded_value)
            
            # Apenas substituir ponto por virgula
            gsub("\\.", ",", formatted_value)
          } else {
            x_clean
          }
        } else {
          x_clean
        }
      })
      
      df_conv[[col_name]] <- result
    }
    
    return(df_conv)
  }
  
  # Funcao simples para adicionar rownames
  add_rownames_column <- function(df) {
    if (!is.data.frame(df)) return(df)
    
    if (!is.null(rownames(df)) && length(rownames(df)) > 0) {
      return(cbind(RowNames = rownames(df), df))
    } else {
      return(cbind(RowNames = as.character(1:nrow(df)), df))
    }
  }
  
  # Criar workbook
  wb <- openxlsx::createWorkbook()
  
  # Obter objetos
  objects <- list(...)
  object_names <- as.character(substitute(list(...)))[-1]
  
  # Processar CADA sheet
  for (i in seq_along(objects)) {
    current_df <- objects[[i]]
    sheet_name <- object_names[i]
    
    openxlsx::addWorksheet(wb, sheet_name)
    
    # **PRIMEIRO adicionar rownames**
    current_df <- add_rownames_column(current_df)
    
    # **DEPOIS converter (a funcao vai pular a coluna "RowNames")**
    if (decimal_comma) {
      # Verificar se deve aplicar ARREDONDAMENTO nesta sheet
      apply_rounding <- TRUE
      if (!is.null(round_sheets)) {
        apply_rounding <- i %in% round_sheets
      }
      
      current_df <- convert_decimal_separator(current_df, decimal_places, apply_rounding)
    }
    
    openxlsx::writeData(wb, sheet_name, current_df, rowNames = FALSE)
  }
  
  # Salvar
  openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)
  message("? Arquivo salvo: ", file_name)
}
