#' @title Metric range test
#' @description
#' This function evaluates the range for a set of numeric metrics as 
#' a screening step for building a multimetric indicator. Specifically, 
#' this function is intended to identify metrics with a limited range or 
#' those that are highly skewed, such that most of the values are a single
#' or very few values. It allows user-specified parameters to perform 
#' these evaluations, with default values provided. 
#' @param df Input data frame in wide format with one row per sample and 
#' assumed to contain only numeric metrics. It should also contain only the 
#' one visit to each site. If there are calibration and validation subsets,
#' this data frame should only contain calibration samples.
#' @param perc_vars A character vector containing names or 
#' partial names (e.g., PTAX, PIND) that clearly identify the metrics that 
#' are percentages.
#' @param id_vars A character vector containing any variables that 
#' identify samples and are not metrics, making it wise to drop any variables 
#' from the input data frame that are not necessary.
#' @param quant_zero This value is the maximum allowable proportion of 
#' samples equal to either the minimum or the maximum metric value. If the
#' proportion exceeds this value, the metric fails the range test. 
#' Proportion should be a number between 0 and 1 and default is 0.75.
#' @param pass_ If the possibility of a partial pass is desired, 
#' provide an alternate value to quant_zero which is lower than quant_zero. 
#' If the proportion of samples is equal to the minimum or maximum value is 
#' between pass_ and quant_zero, a value of PASS- is assigned. If left blank, 
#' this part is not performed.
#' @param quant_range This value determines what upper proportion of sites 
#' is used to evaluate metric range. For example, a value of 0.8 considers the 
#' range of values in the max - 20th percentile. To examine the whole range, 
#' use a value of 1. Value should be a number between 0 and 1, but the default 
#' value is 0.8. This value particularly affects non-percentage metrics 
#' because they do not have a set maximum. The function will fail metrics 
#' for which most of the values are compressed to a small range. An example 
#' would be a metric with a range of 0-30, but with only 25\% of values 
#' representing range 0-27 and 75\% of values between 28-30. If this input is
#' set to 1, non-percentage metrics will always pass, but percentage metrics
#' will fail if the full range of values is < \strong{quant_range_perc}.
#' @param quant_range_perc This numeric value represents a discrete value 
#' to which ranges of percentage metrics are compared. For example, if 
#' \strong{quant_range} = 0.75 and \strong{quant_range_perc} = 15, and 
#' the range of the upper 75\% of values is less than 15, the metric fails.
#' @param quant_range_oth Used for non-percentage metrics (e.g., richness), 
#' this numeric proportion of the full range is used as a point of comparison 
#' for the range determined using \strong{quant_range}. For example, the default value 
#' of 1/3 multiplies the full range by 1/3. For a quant_range of 0.75 and 
#' quant_range_oth of 1/3, the difference between the 25th percentile and the 
#' max is compared to the max/3. If the range is less than max/3, the metric 
#' fails. Value should be a number between 0 and 1.
#' @returns For each metric, provides metric name, overall RANGE_TEST result (PASS/ 
#' PASS-/FAIL), skewness test result (zero_test, PASS/PASS-/FAIL), and
#' range test result (rg_lim, PASS/FAIL). Along with these, the quantiles used
#' to assess range and skewness are also provided: p0 = minimum, p100 = maximum,
#' plower = the (1 - \strong{quant_zero}) quantile, pupper = the \strong{quant_zero} quantile,
#' prob_lower_rg = the (1 - \strong{quant_range}) quantile, pmid = median value.
range_test <- function(df, 
                      perc_vars,
                      id_vars,
                      quant_zero = 0.75,
                      pass_ = NULL,
                      quant_range = 0.8,
                      quant_range_perc = 15, 
                      quant_range_oth = 1/3){

  met_names <- names(df)[!names(df) %in% id_vars]
  perc_vars_1 <- paste(perc_vars, collapse='|')
  perc_names <- grep(perc_vars_1, names(df), value=TRUE)
  
  # Check ranges of values of input arguments
  if(quant_zero < 0|quant_zero > 1){
    return(print("quant_zero value needs to be between 0 and 1"))    
  }
  if(quant_range < 0|quant_range > 1){
    return(print("quant_range value needs to be between 0 and 1"))
  }
  if(quant_range_oth < 0|quant_range_oth > 1){
    return(print("quant_range_oth value needs to be between 0 and 1"))
  }
  if(!is.null(pass_) & pass_ < 0|pass_ > quant_zero){
    return(print("pass_ value needs to be > 0 and < quant_zero"))
  }
  
  # Now melt input data frame
  in_long <- mutate(df, across(all_of(met_names), as.numeric)) |> 
    tidyr::pivot_longer(cols = all_of(met_names), names_to='variable',
                                values_drop_na=TRUE) |>
    dplyr::filter(!is.infinite(value)) 
  
  # First get quantiles based on inputs
  rg_q_all <- in_long |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(
      p0 = min(value, na.rm = T),
      p100 = max(value, na.rm = T),
      plower = quantile(value, probs = (1-quant_zero), na.rm=T),
      pupper = quantile(value, probs = (quant_zero), na.rm=T),
      prob_lower_rg = quantile(value, probs=(1 - quant_range), na.rm=T)) |> 
    dplyr::ungroup()
    
  if(!is.null(pass_)){
    rg_q_pass_ <- in_long |> 
      dplyr::group_by(variable) |> 
      dplyr::summarise(pmid = quantile(value, probs = pass_, na.rm=T)) |> 
      dplyr::ungroup()
    
    rg_q_all <- merge(rg_q_all, rg_q_pass_, by='variable')
  }else{
    rg_q_all <- dplyr::mutate(pmid=NA)
  }
  
  # Range test - apply percentiles calculated above to test range of metrics
  rg_test_all <- mutate(rg_q_all,
                        zero_test = dplyr::case_when(
                          pupper == p0 ~ 'FAIL', # large proportion equal to minimum
                          plower == p100 ~ 'FAIL', # large proportion equal to maximum
                          pupper == 0 ~ 'FAIL', # all zeros
                          !is.null(pass_) & pmid == p0 ~ 'PASS-', # pass- proportion equal to minimum
                          .default = 'PASS'),
                      rg_lim = dplyr::case_when(
                        (variable %in% perc_names) & ((p100-prob_lower_rg) < quant_range_perc) ~ 'FAIL', # for % variable, examined range is less than specified min range
                        (!(variable %in% perc_names) & ((p100-prob_lower_rg) < ((p100-p0)*quant_range_oth))) ~ 'FAIL', # if not % variable, examined range smaller than proportion of full range
                        .default = 'PASS'))
  
  rg_out_all <- mutate(rg_test_all,
                       RANGE_TEST = dplyr::case_when(
                         zero_test=='FAIL'|rg_lim=='FAIL' ~ 'FAIL', # if either part of test fails, FAIL
                         zero_test == 'PASS-' ~ 'PASS-', # if first test gets PASS- but rg_lim = PASS
                         .default = 'PASS'), # otherwise PASS
                       METRIC = as.character(variable)) |>
    select(-variable)
  
  return(rg_out_all)
  
}

