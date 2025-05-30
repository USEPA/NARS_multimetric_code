library(tidyverse)
library(lme4)

source(paste0(here::here(), "/r/range_test.r"))
 # Range test
rg_test_data <- read_csv(paste0(here::here(), "/examples/Range_test_testdata.csv"))

test_rg <- range_test(rg_test_data, 
                      perc_vars=c('PCTN','XCOV','XABCOV','RFREQ','FREQ','RIMP','IMP','XRCOV'),
                      id_vars='UID', quant_zero=0.75,
                      pass_=0.5, quant_range=1,
                      quant_range_perc=15, quant_range_oth=1/3)


# Redundancy test
source(paste0(here::here(), "/r/redund_test.r"))

redund_test_data <- read_csv(paste0(here::here(), "/examples/Redund_test_testdata.csv"))
test_redund <- redund_test(redund_test_data, id_vars = 'UID', cutoff = 0.70)


# Comparison test
source(paste0(here::here(), "/r/comp_test.r"))

comp_test_data <- read_csv(paste0(here::here(), "/examples/Boxplot_comp_test_testdata.csv")) |> 
  mutate(REF_ALT = case_when(
    REF_NWCA == 'L' ~ 'R',
    REF_NWCA == 'M' ~ 'T',
    .default = REF_NWCA
  )) |> 
  select(-REF_NWCA)

test_comp <- comp_test(comp_test_data, 
                       id_vars = 'UID',
                       ref_var = 'REF_ALT', 
                       least = 'R', 
                       most = 'T')

comp_test_long <- pivot_longer(comp_test_data, cols = N_MONOCOT:FQAI_COV_NAT) |> 
  filter(REF_ALT != 'I')
ggplot(comp_test_long, aes(x = REF_ALT, y = value)) +
  geom_boxplot() +
  facet_wrap(~ name, scales = 'free_y')

# Signal-to-noise
source(paste0(here::here(), "/r/sn_test.r"))

sn_test_data <- read_csv(paste0(here::here(), "/examples/SN_test_testdata.csv"))

test_sn <- sn_test(df = sn_test_data,
                   id_vars_samp = 'UID',
                   id_vars_site = 'SITE_ID') 

# ProcIBI by type code
source(paste0(here::here(), "/r/ProcIBI_General_byMetricType_updated_30Jan2025.r"))

metlist <- read_csv(paste0(here::here(), "/examples/vegmetlist_bytype.csv")) 
mets <- read_csv(paste0(here::here(), "/examples/Veg_scored_mets.csv")) |> 
  mutate(indvis = if_else(VISIT_NO=='1', 'Yes', 'No'))

prIBI_bytype_test <- prIBI_byType(df = mets, 
                           metList = metlist,
                           idVars = 'UID', 
                           siteVar = 'SITE_ID', 
                           refVar = 'REF_NWCA', 
                           year = NULL, 
                           indexvis = 'indvis', 
                           least = 'L',
                           most = 'M', 
                           nsamp = 1000, 
                           seed = 20160310)

# ProcIBI by number of metrics code
source(paste0(here::here(), "/r/ProcIBI_General_byNumMets.r"))

metlist <- read_csv(paste0(here::here(), "/examples/vegmetlist_bytype.csv"))$METRIC 
mets <- read_csv(paste0(here::here(), "/examples/Veg_scored_mets.csv")) |> 
  mutate(indvis = if_else(VISIT_NO=='1', 'Yes', 'No'))

prIBI_bynum_test <- prIBI_byNumMets(df = mets, 
                           metList = metlist,
                           idVars = 'UID', 
                           siteVar = 'SITE_ID', 
                           refVar = 'REF_NWCA', 
                           year = NULL, 
                           indexvis = 'indvis', 
                           least = 'L',
                           most = 'M', 
                           nummets = c(4, 6, 8),
                           nsamp = 1000, 
                           seed = 20160310)

# Relative scope of impairment
source(paste0(here::here(), "/r/relSOI_test.r"))

soi_test_data <- read_csv(paste0(here::here(), "/examples/Boxplot_comp_test_testdata.csv")) 

soi_test <- relSOI_test(soi_test_data,
                        id_vars = 'UID',
                        ref_var = 'REF_NWCA', 
                        least = 'L', 
                        most = 'M')

