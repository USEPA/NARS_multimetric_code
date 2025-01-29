#' @title Metric redundancy test
#' @description
#' This function evaluates redundancy of metrics by comparing Pearson 
#' correlations based on a user supplied
#' cutoff. 
#' @param df Wide data frame with one row per sample. It is assumed to only 
#' contain one sample per sites, so if there are revisits, keep only the one 
#' visit. It is also assumed that all metrics are numeric, so remove any that 
#' are not. 
#' @param id_vars This is a character vector with any non-metric variables in 
#' the data frame.
#' @param cutoff Minimum absolute correlation value at which metrics are 
#' flagged as redundant
#' @returns Output dataset has `metric` with metric name, and `redund_met` is a 
#' string containing all metrics with an absolute correlation of at least the 
#' cutoff value, separated by commas. 

redund_test <- function(df, id_vars, cutoff){

  # Melt input dataset by all variables in data frame that are not numeric metrics (or remove those variables)
  met_names <- names(df)[!names(df) %in% id_vars]
  df_in <- df |> 
    mutate(across(all_of(met_names), as.numeric))

  corr_in <- tidyr::pivot_longer(df_in,
                                 cols = met_names, 
                                 names_to='variable',
                                 values_drop_na=TRUE) |>
    filter(!is.infinite(value)) |>
    mutate(value = as.numeric(value))
  # Create empty data frame 
  corr_out <- data.frame(metric = character(), 
                         redund_met = character())
  
  for(i in 1:length(met_names)){
    print(i)
    # Make sure to alter to exclude the id variables used in melting the input data
    cor_met <- cor(subset(corr_in, variable==met_names[i], select='value'),
                  subset(df_in, select=!(names(df_in) %in% c(id_vars, met_names[i]))),
                  method="pearson")
    red_met <- data.frame(metric = attr(cor_met,"dimnames")[[2]],
                         r = cor_met[1:length(cor_met)]) |> 
      subset(abs(r)>=cutoff, select='metric')
    red_list <- data.frame(metric = met_names[i],
                          redund_met = paste(red_met$metric, collapse=","))
    corr_out <- rbind(corr_out, red_list)
  }
  return(corr_out)
}

