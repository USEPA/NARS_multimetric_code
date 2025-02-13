#' @title Signal-to-noise test
#' @description
#' This function calculates a version of the signal-to-noise approach used
#' in NARS and in EMAP. It is intended to use variation due to within-cycle 
#' revisits to estimate noise and variance across sites to estimate signal. 
#' The model used is a random effects model and assumes that sites are selected 
#' using a probability-based design. 
#' @param df This data frame is in wide format with one row per sample and 
#' assumed to contain only numeric metrics. It should also contain all visits 
#' to each site, and at least a subset of sites must have multiple visits. 
#' If there are calibration and validation subsets, this df should only contain 
#' calibration samples. Only identifying variables and metrics should be 
#' included in this data frame, and only numeric metrics can be included.
#' @param id_vars_samp A character vector containing variables that identify
#' individual samples (e.g., UID)
#' @param id_vars_site A string containing a variable name that 
#' identifies sites (e.g., SITE_ID). This cannot be the same or a subset of
#' variables in `id_vars_samp`. This should not be a variable that occurs across
#' survey cycles (e.g., UNIQUE_ID), because this would assume that variation 
#' across cycles is part of the noise, rather than a true change over time.
#' @param year A string containing the name of the year variable if sites are
#' revisited across years (as well as within year). The default value is NULL.
#' @returns Output is a data frame with the metric, signal, noise, and 
#' sn_ration, and comment. If the model is undefined for a given metric, the 
#' values for signal, noise, and sn_ratio will be NA with a comment. Otherwise,
#' each of these values will be provided.  
#' 
sn_test <- function(df,
                   id_vars_samp,
                   id_vars_site,
                   year = NULL){
  
  options(warn=2)
    # Do some error checking first
  if(id_vars_site %in% id_vars_samp){
    return(print("id_vars_samp CANNOT be a part of id_vars_samp"))
  }
  # Make sure there are some repeated sites
  if(nrow(df[duplicated(df[, id_vars_site]),])<1){
    return(print("You need multiple visits for at least 1 site"))
  }  
  print(c('Number of revisits:', nrow(df[duplicated(df[, id_vars_site]),])))

  met_names <- names(df)[!names(df) %in% c(id_vars_samp, id_vars_site, year)]
  if(!is.null(year)){
    df_in <- pivot_longer(df, 
                           cols=all_of(met_names),
                           names_to='variable', 
                           values_drop_na=TRUE) %>%
      filter(!is.infinite(value)) %>%
      mutate(variable = as.character(variable), 
             value = as.numeric(value)) |> 
      dplyr::rename(c(site = all_of(id_vars_site),
                      year = year))
  }else{
    df_in <- pivot_longer(df, 
                          cols=all_of(met_names),
                          names_to='variable', 
                          values_drop_na=TRUE) %>%
      filter(!is.infinite(value)) %>%
      mutate(variable = as.character(variable),
             value = as.numeric(value)) |> 
      dplyr::rename(c(site = all_of(id_vars_site)))
  }
  
  ## Signal-to-noise ratio
  # Create empty data frame to accept output for each metric
  sn_out <- data.frame(metric = character(),
                      signal = numeric(),
                      noise = numeric(),
                      sn_ratio = numeric(),
                      com = character())
  
  # For each metric in parList, run a linear mixed-effects model with SITE_ID as a random effect
  for(i in 1:length(met_names)){
    in_met <- subset(df_in, variable == met_names[i])
    # print(parList[i])					
    # Run model
    if(!is.null(year)){
      sn <- try(lmer(value~year + (1|year:site), 
                     in_met),
                silent = TRUE)
      
      if(class(sn)=='try-error'){
        sn <- try(lmer(value~year + (1|year:site),
                       in_met, 
                       control= lmerControl(optimizer = "bobyqa", 
                                            optCtrl = list(maxfun = 100000))),
                  silent=TRUE)
      }
    }else{
      sn <- try(lmer(value~(1|site), in_met), 
                silent=TRUE)
    }
    
    # If model output is error, send to output data frame
    if(class(sn)=='try-error'){
      s_to_n <- data.frame(metric = met_names[i],
                           signal = NA,
                           noise = NA,
                           sn_ratio = NA,
                           com = sn[1])
    }else{
      # If model output value, determine signal and noise by extracting variance components due to site and error
      var_comp <- VarCorr(sn)
      if(!is.null(year)){
        s_to_n <- data.frame(metric = met_names[i],
                           signal = round(var_comp$'year:site'[1], 2), 
                           noise = round(attr(var_comp, "sc")^2, 2),
                           sn_ratio = round(var_comp$'year:site'[1]/(attr(var_comp,"sc")^2), 2),
                           com = NA)
      }else{
        s_to_n <- data.frame(metric = met_names[i], 
                             signal = round(var_comp$site[1], 2), 
                             noise = round(attr(var_comp,"sc")^2, 2),
                             sn_ratio = round(var_comp$site[1]/(attr(var_comp,"sc")^2), 2),
                             com = NA)
      }
    }			
    sn_out <- rbind(sn_out, s_to_n)
  }
  return(sn_out)
  
}