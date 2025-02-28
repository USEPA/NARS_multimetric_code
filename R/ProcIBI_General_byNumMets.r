#' @title Construct and evaluate MMI based on various numbers of metrics, 
#' regardless of metric type.
#' @description
#' This function allows the user to automate the process of constructing and 
#' comparing multimetric indices built from various numbers of metrics. 
#' This particular function allows the user to 
#' decide on the number of iterations to compare, with a default of 1000. 
#' This version of the function uses an input of metric names and
#' uses this information to construct random combinations of metrics built 
#' using specified numbers of metric. 
#' Several tests are carried out and results are provided for 
#' each constructed MMI as the output of this function. The tests include the 
#' mean and standard deviation of reference sites, proportion of intermediate 
#' and most disturbed sites statistically below reference (interval test), 
#' and mean and max correlations among metrics within a given MMI. A signal-to-
#' noise statistic is also calculated for each constructed MMI. 
#' This code can be customized for the type of data available to the user, so 
#' some tests can be dropped and others added. The primary benefit of this
#' code is the set up of randomly selected combinations of metrics. 
#' @param df Input data frame with each sample in one row, with sample-identifying 
#' variables and metrics as columns. Assumed to have metrics already scored or 
#' rescaled (in whatever manner is desired). Ideally, the MMI is 
#' built on only calibration sites, and it is assumed any sites in the input 
#' dataset will be included in MMI development, so exclude any sites that 
#' shouldn't be used ahead of time. It is assumed that there are multiple visits 
#' for some sites and that the `indexvis` argument identifies the visit to use 
#' for most tests (except signal-to-noise). The default value for the Site 
#' identifier is `SITE_ID`.
#' @param metList Vector of metric names to be used in selecting metrics to
#' combine into candidate multimetric indices (MMIs). If scored metrics have 
#' different names than raw metric values, be sure these 
#' names are the ones in this data frame.
#' @param idVars Variable or combination of variables in `df` used to 
#' identify individual samples. 
#' @param siteVar String containing name of variable that identifies unique 
#' site ID (across years if site was sampled in multiple years). 
#' Default value is SITE_ID.
#' @param refVar String containing name of reference variable in the dataset
#' @param year String containing the name of the year variable in the dataset. 
#' The default value is NULL. If NULL, signal-to-noise calculation is based
#' only on within year revisits.
#' @param indexvis String containing name of variable to identify whether a 
#' visit is the index visit to a site, with values of Yes/No. If a record has
#' an indexvis value of No, it will only be included in signal-to-noise 
#' calculations.
#' @param least Value of `refVar` for least disturbed sites in dataset.
#' @param most Value of `refVar` for most disturbed sites in dataset.
#' @param nummets Numeric vector with numbers of metrics to test 
#' (i.e., c(4,6) tests combinations of 4 metrics and 6 metrics).
#' @param nsamp Number of random combinations to test for each number of 
#' metrics, with a default of 1000.
#' @param seed Random seed to supply to random number generator used to 
#' randomly select metric combinations. Default is 20160310.
#' @returns A data frame containing the iteration number, 
#' and for each iteration, the combination of metrics in the MMI,
#' percent of most disturbed sites falling statistically below reference 
#' (interval test, pct_MMI_M), mean MMI score among reference sites (mn_mmi_ref),
#' standard deviation among reference sites (sd_mmi_ref), 
#' signal-to-noise ratio, maximum correlation among
#' metrics (max_corr), and mean correlation among metrics (mean_corr).
#' 
prIBI_byNumMets <- function(df, 
                  metList,
                  idVars, 
                  siteVar = 'SITE_ID', 
                  refVar, 
                  year = NULL, 
                  indexvis, 
                  least,
                  most,
                  nummets,
                  nsamp = 1000, 
                  seed = 20160310){

  # Step 1 - Prepare datasets from inputs
  df_1 <- base::subset(df, select=names(df) %in% 
                         c(idVars, refVar, siteVar, indexvis, 
                           year, metList)) 
  
  if(!is.null(year)){
    names(df_1)[names(df_1)==year] <- 'year'
  }
  
  df_1 <- dplyr::rename(df_1, site_use = all_of(siteVar))
  
  sum(complete.cases(df_1[, metList])); #check for missing metric data;
  # delete any sites having incomplete metric data for all candidates;
  df_2 <- df_1[complete.cases(df_1[, metList]),]
  # Table of sites by REF_NWCA
  table(df_2[, refVar])
  
  # Datasets are ready for random MMI construction
  # Step 2 - Calculate critical values of interval-test statistics
  # First, critical values for MMI tests
  type_cnt <- with(df_2, table(eval(as.name(refVar)), eval(as.name(indexvis)))) # count sites, each disturbance class (L-I-M)
  print(type_cnt)
  
  n_ref <- type_cnt[least, 'Yes'] # number of reference sites that are index visits
  print(n_ref)
  
  # Univariate F test at 0.05 level. See Kilgour;
  # Use Kilgour notation for 0.05, even though it seems backwards to obtain non-central F
  F_uni_05 <- qf(p = 0.95, 
                 df1 = 1, 
                 df2 = n_ref - 1, 
                 ncp = n_ref*3.841)
  
  # Step 3 - Create largest metric base combinations and create data frames corresponding
  #   to the nummets argument
  # For a single base combo, all subsets are nested
  
  # Select a random set of the max value in nummets metrics, then randomly select next highest
  # metrics from that max set. The set of the next highest value is then selected from that, and so on.
  nummets <- nummets[order(desc(nummets))]
  numnums <- length(nummets)
  maxnum <- max(nummets)
  
  # define empty data frames for base combo and added metrics;
  matOut <- vector("list", numnums)
  for(k in 1:numnums){
    df_name <- paste("met_samp.", nummets[k], sep='')
    assign(df_name, matrix(rep(' ', nummets[k]*nsamp), ncol=nummets[k]))
    matOut[[k]] <- eval(as.name(df_name))
  }
  
  # loop over base samples,, and populate the added-metrics in data frames;
  set.seed(seed)
  for(k in 1:nsamp) {
    for(i in 1:numnums){
      if(i==1){
        matOut[[i]][k,1:nummets[i]] <- sample(metList, size=nummets[i])
      }else{
        matOut[[i]][k,1:nummets[i]] <- sample(matOut[[i-1]][k, 1:nummets[i-1]], 
                                              size=nummets[i], replace=FALSE)
      }
    }
  }
  
  
  # Step 5 - Calculate MMIs (Loop over all metric combos and cases, calculate number of sites outside for each)**
  # Set up results data frame
  ntrial <- length(nummets)*nsamp # total number trials = base*numcases
  print(ntrial)
  
  nmet_res <- expand.grid(nummets=nummets, baseind = seq(1:nsamp))
  
  #performance metrics are percent outside reference, for MMI for Most (M) versus intermediate (I) disturbed;
  nmet_res <- data.frame(grp = rep(seq(1:numnums), nsamp),
                         nmet_res,
                         pct_MMI_M = rep(NA, ntrial),
                         mets = rep(NA, ntrial),
                         mn_mmi_ref = rep(NA, ntrial),
                         sd_mmi_ref = rep(NA, ntrial),
                         sn_mmi = rep(NA, ntrial),
                         max_corr = rep(NA, ntrial), 
                         mean_corr = rep(NA, ntrial))
  
  
  #loop over base samples and their add-metric cases
  start_time=proc.time()
  for(j in 1:ntrial) {
    if(j%%1000==0) print(j)
    
    # first create and report the current metric set
    for(indx in 1:numnums){
      if(nmet_res$nummets[j]==nummets[indx]){
        curmets <- matOut[[indx]][nmet_res$baseind[j],]
      }
    }
    nmet_res$mets[j] <- paste(as.character(curmets),collapse=", ")
    
    # Next are interval tests for all test sites
    # Calculate vector of MMI at all sites, equal to sum of scored metrics,
    # Scaled by 10/nummet, to put on 100-point scale
    IBI_cur <- rowSums(mets[, curmets])*10/nmet_res$nummet[j]
    # mean of IBI and and its squared SE, for ref sites and first visit only
    # mean of IBI and and its squared SE, for ref sites and first visit only
    ibi_ref_mn <- mean(IBI_cur[df_2[, refVar]==least & df_2[, indexvis]=='Yes'])
    nmet_res$mn_mmi_ref[j] <- ibi_ref_mn
    
    ibi_ref_se2 <- var(IBI_cur[df_2[, refVar]==least & df_2[, indexvis]=='Yes'])/n_ref
    nmet_res$sd_mmi_ref[j] <- sqrt(ibi_ref_se2*n_ref)
    
    # next extract the subset of test-site MMI scores, for either most disturbed (M) or Intermediate 
    # distubance (I)
    IBI_sub_M <- IBI_cur[df_2[, refVar]==most & df_2[, indexvis]=='Yes'] 
    
    # Vector of 1-sample F-scores, all test sites;
    F_ibi_M <- ((ibi_ref_mn-IBI_sub_M)^2)/ibi_ref_se2
    
    # Percent impaired sites. Under IBI (=MMI), impaired if F.IBI > F.uni.05. Kilgour p.546   
    nmet_res$pct_MMI_M[j] <- 100*sum(F_ibi_M > F_uni_05)/length(F_ibi_M)
    
    # S/N ratio for MMI 
    if(!is.null(year)){
      cur_temp <- data.frame(sitevar = df_2$site_use, 
                             year = df_2$year,
                             IBI_cur)
    }else{
      cur_temp <- data.frame(sitevar = df_2$site_use,
                             IBI_cur)
    }
    
    if(!is.null(year)){
      sn <- try(lmer(IBI_cur~year + (1|year:sitevar), cur_temp), silent=TRUE)
      
      if(class(sn)=='try-error'){
        sn <- try(lmer(value~year + (1|year:sitevar), inMet, 
                       control= lmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 100000))), 
                  silent=TRUE)
        
        if(class(sn)=='try-error'){
          sn <- try(lmer(value~year + (1|sitevar), inMet, 
                         control = lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 100000))), 
                    silent=TRUE)
        }
      }
    }else{
      sn <- try(lmer(IBI_cur~(1|sitevar), cur_temp),silent=TRUE)
    }
    
    if(class(sn)=='try-error'){
      nmet_res$sn_mmi[j] <- NA
    }else{
      # If model output value, determine signal and noise by extracting variance components 
      # due to sitevar and error
      varcomp <- VarCorr(sn)
      
      if(!is.null(year)){
        nmet_res$sn_mmi[j] <- round(varcomp$'year:sitevar'[1]/(attr(varcomp,"sc")^2),digits=2)
      }else{
        nmet_res$sn_mmi[j] <- round(varcomp$sitevar[1]/(attr(varcomp,"sc")^2),digits=2)
      }
      
    }  		      
    
    # Max and mean correlations among metrics in MMI
    corr_cur <- as.dist(cor(mets[, curmets], method='pearson'))
    nmet_res$max_corr[j] <- round(max(corr_cur), digits=3)
    nmet_res$mean_corr[j] <- round(mean(corr_cur), digits=3)
  } # end of trial loop
  
  elaps <- proc.time()-start_time
  print(c("elapsed time = ",elaps))
  
  return(nmet_res)
}

