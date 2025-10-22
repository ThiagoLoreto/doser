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
