library(tidyverse)

source(paste0(here::here(), "/r/range_test.r"))

veg_met <- read.delim("c:/users/kblockso/onedrive - environmental protection agency (EPA)/documents/nars/rscripts/metric eval vegetation testing/VegMetrics.tab", sep = '\t')

sites <- read.delim("c:/users/kblockso/onedrive - environmental protection agency (EPA)/documents/nars/rscripts/metric eval vegetation testing/SiteInfo.tab", sep = '\t') |> 
  select(UID,REF_NWCA)

# Make sure all metrics in this data frame are numeric and remove any in this step that are not
v_in_first <- filter(veg_met, SITE_USE != 'NWCA_REVISITS') |> 
  select(-c(DOM_SANDT, LITTER_TYPE, 
                           PUBLICATION_DATE, VISIT_NO, SITE_USE, STATE, SITE_ID))
  

rg_test_data <- read_csv(paste0(here::here(), "/Range_test_testdata.csv"))

test_rg <- range_test(rg_test_data, 
                      perc_vars=c('PCTN','XCOV','XABCOV','RFREQ','FREQ','RIMP','IMP','XRCOV'),
                      id_vars='UID', quant_zero=0.75,
                      pass_=0.5, quant_range=1,
                      quant_range_perc=15, quant_range_oth=1/3)

# Save test_rg as expected test results
# write_csv(test_rg, paste0(here::here(), "/Range_test_output.csv"))
# Compare output to expected results 
exp_rg <- read_csv(paste0(here::here(), "/Range_test_output.csv"))

test_rg_long <- pivot_longer(test_rg, cols = p0:pmid)
exp_rg_long <- pivot_longer(exp_rg, cols = p0:pmid)

# Compare numeric values
inner_join(test_rg_long, exp_rg_long, by = c('METRIC', 'name')) |> 
  filter(value.x!=value.y) # Should be 0

inner_join(test_rg, exp_rg, by = 'METRIC') |> 
  filter(zero_test.x != zero_test.y|
           rg_lim.x!=rg_lim.y|
           RANGE_TEST.x!=RANGE_TEST.y) # Should be 0
