# ProcIBI_General_VariousMetNum_forDistribution.r
# Purpose: Construction and Evaluation of MMIs Based on Randomly Selected Metrics to Identify Optimum Metric Number - 
#       No metric types
# 3/9/2016 Karen Blocksom altered from existing code for Vegetation MMI for NWCA 2011 to work as a function
# 10/27/2020 Updated by Karen Blocksom to add multi-year samples for S/N and also to add arguments.


# Load Packages for Analysis
library(Hmisc)
library(plyr)
library(lme4)


prIBI <- function(metList, dfIn, idVars, siteVar='SITE_ID', refVar, year=NULL, indexvis, least,
                  most, nummets, nsamp=1000, seed=20160310, grpName='National',dirName){
#  function(metList,dfIn,idVars,refVar,least,most,nummets,nsamp=1000,seed=20160309){
  # dfIn - Input data frame with each sample in one row, with sample identifying variables and metrics as columsn.
  #       Assumed to have metrics already scored or rescaled (in whatever manner is desired). Ideally, the MMI is 
  #       built on only calibration sites, and it is assumed any sites in the input dataset will be included in 
  #       MMI development, so exclude any sites that shouldn't be used ahead of time.It is assumed that there
  #       are multiple visits for some sites and that the visit variable is 'VISIT_NO'. Site identifier is assumed 
  #       to have the name SITE_ID.
  #
  # metList - if scored metrics have different names than raw metric values, be sure these names are the ones 
  #       in this character vector
  #
  # idVars - Variables in dataset used to identify individual samples. This should include 
  #   variables for individual samples
  #         as well as sites, and visit number should be identified using variable VISIT_NO
  #
  # siteVar - string containing name of variable that identifies unique site ID (across years if 
  #   site was sampled in multiple years). Default is SITE_ID.
  #
  # refVar - string containing name of reference variable
  #
  # year - string containing name of Year variable if sites are revisited across years (as well as within year), 
  #   default is 'YEAR'. Set to NULL (default) if no samples across years.
  #
  # indexvis = string containing name of variable to identify whether a visit is the index visit to a
  #     site, with values of Y/N (Yes/No)
  #
  # least - value of refVar for least disturbed sites
  #
  # most - value of refVar for most disturbed sites
  #
  # nummets - numeric vector with numbers of metrics to test (i.e., c(4,6) tests combinations of 4 
  #   metrics and 6 metrics)
  #
  # nsamp - Number of random combinations to test for each number of metrics, with default value of 1000
  #
  # seed - Seed for random number generator. The default value is 20160309. 
  #
  # grpName - String for subset of data being used to generate VMMI. Default is 'National' and will
  #   be added to the file name containing the scored metrics.
  ###########################################################################################################

  # Step 1 - Prepare datasets from inputs
  dfIn.1 <- subset(dfIn,select=names(dfIn) %in% c(idVars,refVar,siteVar,indexvis,year,metList))  

  if(!is.null(year)){
    names(dfIn.1)[names(dfIn.1)==year] <- 'year'
  }
  
  dfIn.1$sitevar <- dfIn.1[,siteVar]
  
  sum(complete.cases(dfIn.1[,metList])); #check for missing metric data;
  # delete any sites having incomplete metric data for all candidates;
  dfIn.2 <- dfIn.1[complete.cases(dfIn.1[,metList]),]
  # Table of sites by refVar
  table(dfIn.2[,refVar])

  # Datasets are ready for random MMI construction
  # Step 2 - Setup Tasks
  # scored.mets<-subset(dfIn.2,select=names(dfIn.2) %in% c(metList))
  
  scoremet <- function(met){
    quants <- quantile(dfIn.2[,met],probs=c(0.05,0.95),na.rm = TRUE) # quantiles for floor/ceiling based on all samples;

    # check metric response direction by comparing metric means at ref and trashed
    mndiff <- mean(dfIn.2[dfIn.2[, refVar]==least & dfIn.2[,indexvis]=='Y', met])- mean(dfIn.2[dfIn.2[,refVar]==most & dfIn.2[,indexvis]=='Y', met]) # negative value indicates negative metric
    flush.console()
    
    print(c(met,mndiff))    #print the mean difference just FYI;
    # score the metric. Does a linear interpolation of metric to 0-10 range, direction depends on meandiff
    if(mndiff>0){
      zz <- approx(x=quants,y=c(0,10),xout=dfIn.2[,met],method='linear',yleft=0,yright=10)$y
    }else{
      zz <- approx(x=quants,y=c(10,0),xout=dfIn.2[,met],method='linear',yleft=10,yright=0)$y
    }
    zz
  }  #End of function. return a vector of "scored" metric values for all sites
  
  # Now apply the above function to get scored values for a set of metrics, at all sites, and put results into new data frame
  scored.mets <- as.data.frame(sapply(metList,scoremet)) %>% subset(select=metList)
  
  scored.mets.out <- data.frame(dfIn.2[,c(idVars,refVar,siteVar,indexvis,'year')], scored.mets, stringsAsFactors=F)
  write.csv(scored.mets.out, paste(paste(dirName, "NWCA2016_VarNum_Calib_Scoredmetrics_",sep='/'), grpName, ".csv", sep=""), row.names=F)

  # Step 3 - Calculate critical values of interval-test statistics
  # First, critical values for MMI tests
  type.cnt <- table(dfIn.2[,refVar],dfIn.2[,indexvis]) # count sites, each disturbance class (L-I-M)
  print(type.cnt)
  
  n.ref <- type.cnt[least,'Y'] # number of reference sites; remove visit 2 data for reference sites
  print(n.ref)
  
  # Univariate F test at 0.05 level. See Kilgour;
  # Use Kilgour notation for 0.05, even though it seems backwards to obtain non-central F
  F.uni.05 <- qf(p=.95,df1=1,df2=n.ref-1,ncp=n.ref*3.841)
  
  # Step 4 - Create largest metric base combinations and create data frames corresponding
  #   to the nummets argument
  # For a single base combo, all subsets are nested

  # Select a random set of the max value in nummets metrics, then randomly select next highest
  # metrics from that max set. The set of the next highest value is then selected from that, and so on.
  nummets <- nummets[order(desc(nummets))]
  numnums <- length(nummets)
  maxnum <- max(nummets)
  
  # define empty data frames for base combo and added metrics;
  matOut <- vector("list",numnums)
  for(k in 1:numnums){
    df.name <- paste("met.samp.",nummets[k],sep='')
    assign(df.name, matrix(rep(' ',nummets[k]*nsamp),ncol=nummets[k]))
    matOut[[k]] <- eval(as.name(df.name))
  }
  
  # loop over base samples,, and populate the added-metrics in data frames;
  set.seed(seed)
  for(k in 1:nsamp) {
    for(i in 1:numnums){
      if(i==1){
        matOut[[i]][k,1:nummets[i]] <- sample(metList,size=nummets[i])
      }else{
        matOut[[i]][k,1:nummets[i]] <- sample(matOut[[i-1]][k,1:nummets[i-1]],size=nummets[i],replace=F)
      }
    }
  }
  
  
  # Step 5 - Calculate MMIs (Loop over all metric combos and cases, calculate number of sites outside for each)**
  # Set up results data frame
  ntrial <- length(nummets)*nsamp # total number trials = base*numcases
  print(ntrial)
  
  nmet.res <- expand.grid(nummets=nummets,baseind=seq(1:nsamp))
  
  #performance metrics are percent outside reference, for MMI for Most (M) versus intermediate (I) disturbed;
  nmet.res <- data.frame(grp=rep(seq(1:numnums),nsamp),nmet.res,pct.MMI.M=rep(NA,ntrial),
                      mets=rep(NA,ntrial),
                      mn.mmi.ref=rep(NA,ntrial),sd.mmi.ref=rep(NA,ntrial),sn.mmi=rep(NA,ntrial),
                      max.corr=rep(NA,ntrial),mean.corr=rep(NA,ntrial),stringsAsFactors=F)
  
       
   #loop over base samples and their add-metric cases
   start.time=proc.time()
      for(j in 1:ntrial) {
        if(j%%1000==0) print(j)
  
      # first create and report the current metric set
      for(indx in 1:numnums){
        if(nmet.res$nummets[j]==nummets[indx]){
          curmets<- matOut[[indx]][nmet.res$baseind[j],]
        }
      }
      nmet.res$mets[j] <- paste(as.character(curmets),collapse=", ")
       
    # Next are interval tests for all test sites
    # Calculate vector of MMI at all sites, equal to sum of scored metrics,
    # Scaled by 10/nummet, to put on 100-point scale
      IBI.cur <- rowSums(scored.mets[,curmets])*10/nmet.res$nummet[j]
    # mean of IBI and and its squared SE, for ref sites and first visit only
      ibi.ref.mn <- mean(IBI.cur[dfIn.2[,refVar]==least & dfIn.2[,indexvis]=='Y'])
      nmet.res$mn.mmi.ref[j] <- ibi.ref.mn
  
      ibi.ref.se2 <- var(IBI.cur[dfIn.2[,refVar]==least & dfIn.2[,indexvis]=='Y'])/n.ref
      nmet.res$sd.mmi.ref[j] <- sqrt(ibi.ref.se2*n.ref)
  
    # next extract the subset of test-site MMI scores, for either most disturbed (M) or Intermediate 
    # distubance (I)
      IBI.sub.M <- IBI.cur[as.character(dfIn.2[,refVar])==most & dfIn.2[,indexvis]=='Y'] 
        
    # Vector of 1-sample F-scores, all test sites;
      F.ibi.M <- ((ibi.ref.mn-IBI.sub.M)^2)/ibi.ref.se2

    # Percent impaired sites. Under IBI (=MMI), impaired if F.IBI > F.uni.05. Kilgour p.546   
     nmet.res$pct.MMI.M[j] <- 100*sum(F.ibi.M > F.uni.05)/length(F.ibi.M)
  
    # S/N ratio for MMI 
     if(!is.null(year)){
       cur.temp <- data.frame(sitevar = as.character(dfIn.2$sitevar), year = dfIn.2$year,
                              IBI.cur, stringsAsFactors=FALSE)
     }else{
       cur.temp <- data.frame(sitevar=as.character(dfIn.2$sitevar),IBI.cur,
                              stringsAsFactors=FALSE)
     }
     # cur.temp <- data.frame(sitevar=as.character(dfIn.2[,sitevar]),IBI.cur,stringsAsFactors=FALSE)
      
      if(!is.null(year)){
        sn <- try(lmer(IBI.cur~year + (1|year:sitevar),cur.temp),silent=TRUE)
        
        if(class(sn)=='try-error'){
          sn <- try(lmer(value~year + (1|year:site),inMet, control = lmerControl(optimizer = "bobyqa", 
                                                        optCtrl = list(maxfun = 100000))), silent=TRUE)
          if(class(sn)=='try-error'){
            sn <- try(lmer(value~year + (1|site), inMet, control = lmerControl(optimizer = "bobyqa", 
                                                        optCtrl = list(maxfun = 100000))), silent=TRUE)
          }
          
        }
      }else{
        sn <- try(lmer(IBI.cur~(1|sitevar),cur.temp),silent=TRUE)
      }
#      sn <- try(lmer(IBI.cur~(1|SITE_ID),cur.temp),silent=TRUE) 
      if(class(sn)=='try-error'){
        nmet.res$sn.mmi[j] <- NA
      }else{
        # If model output value, determine signal and noise by extracting variance components 
        #  due to sitevar and error
        varcomp <- VarCorr(sn)
        if(!is.null(year)){
          nmet.res$sn.mmi[j] <- round(varcomp$'year:site'[1]/(attr(varcomp,"sc")^2),digits=2)
        }else{
          nmet.res$sn.mmi[j] <- round(varcomp$sitevar[1]/(attr(varcomp,"sc")^2),digits=2)
        }
        #nmet.res$sn.mmi[j] <- round(varcomp$SITE_ID[1]/(attr(varcomp,"sc")^2),digits=2)
      }  		      
      
    # Max and mean correlations among metrics in MMI
      corr.cur <- as.dist(cor(scored.mets[,curmets],method='pearson'))
      nmet.res$max.corr[j] <- round(max(corr.cur),digits=3)
      nmet.res$mean.corr[j] <- round(mean(corr.cur),digits=3)
    } # end of trial loop

    elaps <- proc.time()-start.time
    print(c("elapsed time = ",elaps))

  return(nmet.res)
}
  

