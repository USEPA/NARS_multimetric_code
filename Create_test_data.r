library(tidyverse)
library(tibble)

# Create data for the range test that will trigger each type of fail or pass
veg_met <- read.delim("c:/users/kblockso/onedrive - environmental protection agency (EPA)/documents/nars/rscripts/metric eval vegetation testing/VegMetrics.tab", sep = '\t')

sites <- read.delim("c:/users/kblockso/onedrive - environmental protection agency (EPA)/documents/nars/rscripts/metric eval vegetation testing/SiteInfo.tab", sep = '\t') |> 
  select(UID,REF_NWCA)

# Make sure all metrics in this data frame are numeric and remove any in this step that are not
v_in_first <- filter(veg_met, SITE_USE != 'NWCA_REVISITS') |> 
  select(-c(DOM_SANDT, LITTER_TYPE, 
            PUBLICATION_DATE, VISIT_NO, SITE_USE, STATE, SITE_ID))

# Data for redundancy
v_in_redund <- select(v_in_first, UID, N_MONOCOT, N_MONOCOTS_NAT, N_MONOCOTS_ALIEN,
                      FQAI_ALL, FQAI_FREQ_ALL, FQAI_COV_ALL, FQAI_COV_NAT)

write_csv(v_in_redund, paste0(here::here(), "/Redund_test_testdata.csv"))

# Data for comparison
v_in_comp <- v_in_redund |> 
  inner_join(sites, by = 'UID')

write_csv(v_in_comp, paste0(here::here(), "/Boxplot_comp_test_testdata.csv"))

# Data for signal-to-noise
v_in_sn <- veg_met |> 
  select(UID, SITE_ID, TOTN_SPP, TOTN_NATSPP, FQAI_ALL, PCTN_HSEN, TOTN_ADVSPP)

write_csv(v_in_sn, paste0(here::here(), "/SN_test_testdata.csv"))
