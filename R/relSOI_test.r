#' @title Relative scope of impairment
#' @description
#' This function calculate the relative scope of impairment of a metric, based 
#' on the reference designations provided by the user. Relative
#' scope of impairment (SOI) is calculated as the ratio of the range
#' of possible impairment (values beyond the poorer quality 25th
#' percentile) and the interquartile range (IQR) of reference sites, a measure of
#' variability. Relative SOI values smaller than 1 are more desirable,
#' with values greater than 1 indicating too much variability among
#' reference sites compared to the range of impairment.
#' @param df Wide data frame with one row per sample. It is assumed to only 
#' contain one sample per sites, so if there are revisits, keep only the one 
#' visit. It is also assumed that all metrics are numeric, so remove any that 
#' are not. This data frame should include variables to identify samples and a 
#' single variable for level of disturbance. The data frame should include only 
#' the most and least disturbed classes. Some consolidation of classes may be 
#' desirable, depending on the situation.
#' @param id_vars This is a character vector with variables to identify samples
#' in the data frame.
#' @param ref_var String identifying the variable with disturbance condition 
#' for each site
#' @param least String representing value of ref_var indicating least 
#' disturbed condition
#' @param most String representing value of ref_var indicating most 
#' disturbed condition
#' @returns The output data frame includes metric, least disturbed 25th 
#' percentile (L_q1), least disturbed 75th percentile (L_q3), least disturbed
#' median (L_median), most disturbed median (M_median), maximum value (max_val),
#' metric direction (POS/NEG, metdir), and  relative scope of impairment (soi). 
#' Desirable values of SOI are 1 or less. Higher values indicate excessive 
#' reference site variability.
#' @references 
#' Blocksom, KA, and BR Johnson. 2009. Development of a regional
#' macroinvertebrate index for large river bioassessment. Ecological 
#' Indicators 9:313-328.
#' 
#' Klemm, DJ, KA Blocksom, WT Thoeny, FA Fulk, AT Herlihy, PR Kaufmann, and
#' SM Cormier. 2002. Methods development and use of macroinvertebrates as 
#' indicators of ecological condition for streams in the Mid-Atlantic 
#' Highlands region. Environmental Monitoring and Assessment 78:169-212.
#' 
#' USEPA, 1998. Lake and Reservoir Bioassessment and Biocriteria Technical 
#' Guidance Document, EPA/841/B-98/007, U.S. Environmental Protection Agency, 
#' Office of Wetlands, Oceans, and Watersheds, Office of Science and 
#' Technology, Office of Water, Washington, D.C.
relSOI_test <- function(df, id_vars, ref_var, least, most){
  met_names <- names(df)[!names(df) %in% c(ref_var, id_vars)]
  
  names(df)[names(df) == ref_var] <- 'ref_vals'
  
  df_in <- df |> 
    mutate(across(all_of(met_names), as.numeric)) |> 
    pivot_longer(cols = all_of(met_names), 
                 names_to='metric', 
                 values_drop_na=TRUE) |> 
    filter(ref_vals %in% c(least, most)) |> 
    mutate(ref_vals = case_when(
      ref_vals == least ~ 'L',
      ref_vals == most ~ 'M'
    ))
  
  soi_out <- data.frame(data.frame(metric = character(),
                                   L_q1 = numeric(),
                                   L_q3 = numeric(), 
                                   L_median = numeric(),
                                   M_median = numeric(),
                                   max_val = numeric(),
                                   metdir = character(),
                                   roi = numeric()))

  for(i in 1:length(met_names)){
  ## Use interquartile ranges to simulate comparison of boxplots and use scoring from Tetra Tech
    quants <- df_in |> 
      filter(metric == met_names[i]) |> 
      group_by(metric, ref_vals) |> 
      summarise(q1 = quantile(value, probs = 0.25),
                q3 = quantile(value, probs = 0.75),
                median = quantile(value, probs = 0.50),
                .groups = 'drop') |> 
      pivot_longer(cols = q1:median) |> 
      pivot_wider(id_cols = metric, 
                  names_from = c(ref_vals, name),
                  names_sep = '_') |> 
      select(-M_q1, -M_q3)
    
    maxes <- df_in |> 
      filter(metric == met_names[i]) |> 
      summarise(max_val = max(value))
    # Determine direction of each metric: positive means value increases
    # with better condition, negative means value increases as 
    # condition decreases - i.e. most disturbed median > least disturbed median
    metdir <- quants |> 
      mutate(metdir = if_else(L_median > M_median, 'POS', 'NEG')) |> 
      merge(maxes)
    
    # Now combine into one data frame to calculate ROI
    soi <- metdir |> 
      mutate(soi = case_when(
        metdir == 'POS' ~ round((L_q3 - L_q1)/L_q1, 3),
        metdir == 'NEG' ~ round((L_q3 - L_q1)/(max_val - L_q3), 3)
      ))
    
    soi_out <- bind_rows(soi_out, soi)
  }
  return(soi_out)
}
