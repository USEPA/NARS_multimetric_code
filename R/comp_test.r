#' @title Metric comparison test
#' @description
#' This function compares the most and least disturbed using a Kruskal-Wallis 
#' test because metrics are typically somewhat skewed. It also performs a box 
#' plot comparison using scoring developed by Tetra Tech for Florida SCI.

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

comp_test <- function(df, id_vars, ref_var, least, most){
  met_names <- names(df)[!names(df) %in% c(ref_var, id_vars)]
  
  df_in <- df |> 
    mutate(across(all_of(met_names), as.numeric)) |> 
    pivot_longer(cols = met_names, 
                        names_to='variable', 
                        values_drop_na=TRUE) |> 
    dplyr::rename(ref_vals = ref_var) |> 
    filter(ref_vals %in% c(least, most)) |> 
    mutate(ref_vals = case_when(
      ref_vals == least ~ 'L',
      ref_vals == most ~ 'M'
    ))

  # Create empty data frame to accept test output
  comp_out <- data.frame(metric = character(),
                        KW_stat = numeric(),
                        KW_pval = numeric(), 
                        box_score = integer())
  
  # For each metric, first run a Kruskal-Wallis test comparing best and worst categories
  # Then determine quantiles to simulate boxplot comparisons as used by Tetra Tech
  for(i in 1:length(met_names)){
    # print(metList[i])
    ## Obtain K-W p-value
    kw <- kruskal.test(value~factor(eval(as.name(ref_var))),
                       data = subset(df_in, variable==met_names[i]))
    
    ## Use interquartile ranges to simulate comparison of boxplots and use scoring from Tetra Tech
    quants <- df_in |> 
      filter(variable == met_names[i]) |> 
      group_by(ref_var) |> 
      summarise(p25 = quantile(value, probs = 0.25),
                p75 = quantile(value, probs = 0.75),
                median = quantile(value, probs = 0.50)) |> 
      ungroup()
    
    quants_1 <- data.frame(metric = met_names[i],
                           pivot_longer(quants, 
                                        cols=names(quants)[!names(quants) %in% ref_var], 
                                        names_to='variable')) |> 
      mutate(variable = paste(eval(as.name(ref_var)), variable, sep='_'))  |> 
      pivot_wider(id_cols='metric', names_from='variable')

    quants_2 <- mutate(quants_1,
                       munder=ifelse(M_median < L_p25, 1, 0),
                       mover=ifelse(M_median > L_p75, 1, 0),
                       lunder=ifelse(L_median < M_p25, 1, 0),
                       lover=ifelse(L_median > M_p75, 1, 0),
                       overlap=case_when(
                         munder==1 & lover==1 & (M_p75 < L_p25) ~ 1,
                         mover==1 & lunder==1 & (L_p75 < M_p25) ~ 1,
                         .default = 0),
                       box_score = sum(munder, mover, lunder, lover, overlap))
    
    temp_comp <- data.frame(metric = met_names[i],
                           KW_stat = round(kw$statistic, 2),
                           KW_pval = round(kw$p.value, 4),
                           box_score=quants.2$box_score)
    comp_out <- rbind(comp_out,temp_comp)
  }
  return(comp_out)
}
