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
#' are percentages
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
#' use a value of 1. Value should be a number between 0 and 1.
#' @param quant_range_perc This numeric value represents a discrete value 
#' to which ranges of percentage metrics are compared. For example, if 
#' quant_range = 0.75 and quant_range_perc = 15, and the range of the upper 75% 
#' of values is less than 15, the metric fails.
#' @param quant_range_oth Used for non-percentage metrics (e.g., richness), 
#' this numeric proportion of the full range is used as a point of comparison 
#' for the range determined using quant_range. For example, the default value 
#' of 1/3 multiplies the full range by 1/3. For a quant_range of 0.75 and 
#' quant_range_oth of 1/3, the difference between the 25th percentile and the 
#' max is compared to the max/3. If the range is less than max/3, the metric 
#' fails. Value should be a number between 0 and 1.

range_test <- function(df, 
                      perc_vars,
                      id_vars,
                      quant_zero = 0.75,
                      pass_ = NULL,
                      quant_range = 1,
                      quant_range_perc = 15, 
                      quant_range_oth = 1/3){

  metNames <- names(dfIn)[names(dfIn) %nin% idVars]
  percVars.1 <- paste(percVars,collapse='|')
  percNames <- grep(percVars.1,names(dfIn),value=TRUE)
  
  # Check ranges of values of input arguments
  if(quantZero<0|quantZero>1){
    return(print("quantZero value needs to be between 0 and 1"))    
  }
  if(quantRange<0|quantRange>1){
    return(print("quantRange value needs to be between 0 and 1"))
  }
  if(quantRange.oth<0|quantRange.oth>1){
    return(print("quantRange.oth value needs to be between 0 and 1"))
  }
  if(!is.null(pass_) & pass_<0|pass_>quantZero){
    return(print("pass_ value needs to be > 0 and < quantZero"))
  }
  
  # Now melt input data frame
  inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% idVars], names_to='variable',
                                values_drop_na=TRUE) %>%
    dplyr::filter(!is.infinite(value)) %>% 
    mutate(value=as.numeric(value))
  
  # First get quantiles based on inputs
  rgQAll <- ddply(inLong,c('variable'),summarise,p0=min(value,na.rm=T),p100=max(value,na.rm=T)
                  ,plower=quantile(value,probs=(1-quantZero),na.rm=T),pupper=quantile(value,probs=(quantZero),na.rm=T),
                  prob.lower.rg=quantile(value, probs=(1-quantRange), na.rm=T))
  
  if(!is.null(pass_)){
    rgQpass_ <- ddply(inLong,c('variable'),summarise,pmid=quantile(value,probs=pass_,na.rm=T))
    rgQAll <- merge(rgQAll,rgQpass_,by='variable')
  }else{
    rgQAll <- mutate(pmid=NA)
  }
  
  # Range test - apply percentiles calculated above to test range of metrics
  rgTestAll <- mutate(rgQAll,zeroTest=ifelse(pupper==p0|plower==p100|pupper==0,'FAIL',
                                             ifelse(!is.null(pass_) & pmid==p0,'PASS-','PASS')),
                      rgLim=ifelse((variable %in% percNames & (p100-prob.lower.rg)<quantRange.perc)|
                                     (variable %nin% percNames & (p100-prob.lower.rg)<((p100-p0)*quantRange.oth)),'FAIL','PASS'))
  
  rgOutAll <- mutate(rgTestAll,RANGE_TEST=ifelse(zeroTest=='FAIL'|rgLim=='FAIL','FAIL'
                                                 ,ifelse(zeroTest=='PASS-','PASS-','PASS')),METRIC=as.character(variable)) %>%
    select(-variable)
  
  return(rgOutAll)
  
}

