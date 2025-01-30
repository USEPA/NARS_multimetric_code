library(tidyverse)
library(lme4)

source(paste0(here::here(), "/r/range_test.r"))
 # Range test
rg_test_data <- read_csv(paste0(here::here(), "/Testing/Range_test_testdata.csv"))

test_rg <- range_test(rg_test_data, 
                      perc_vars=c('PCTN','XCOV','XABCOV','RFREQ','FREQ','RIMP','IMP','XRCOV'),
                      id_vars='UID', quant_zero=0.75,
                      pass_=0.5, quant_range=1,
                      quant_range_perc=15, quant_range_oth=1/3)

# Save test_rg as expected test results
# write_csv(test_rg, paste0(here::here(), "/Range_test_output.csv"))
# Compare output to expected results 
exp_rg <- read_csv(paste0(here::here(), "/Testing/Range_test_output.csv"))

test_rg_long <- pivot_longer(test_rg, cols = p0:pmid)
exp_rg_long <- pivot_longer(exp_rg, cols = p0:pmid)

# Compare numeric values
inner_join(test_rg_long, exp_rg_long, by = c('METRIC', 'name')) |> 
  filter(value.x!=value.y) # Should be 0

inner_join(test_rg, exp_rg, by = 'METRIC') |> 
  filter(zero_test.x != zero_test.y|
           rg_lim.x!=rg_lim.y|
           RANGE_TEST.x!=RANGE_TEST.y) # Should be 0


# Redundancy test
source(paste0(here::here(), "/r/redund_test.r"))

redund_test_data <- read_csv(paste0(here::here(), "/Testing/Redund_test_testdata.csv"))
test_redund <- redund_test(redund_test_data, id_vars = 'UID', cutoff = 0.70)


# Comparison test
source(paste0(here::here(), "/r/comp_test.r"))

comp_test_data <- read_csv(paste0(here::here(), "/Testing/Boxplot_comp_test_testdata.csv")) |> 
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

sn_test_data <- read_csv(paste0(here::here(), "/Testing/SN_test_testdata.csv"))

test_sn <- sn_test(sn_test_data,
                   id_vars_samp = 'UID',
                   id_vars_site = 'SITE_ID') 

# ProcIBI by type code
source(paste0(here::here(), "/r/ProcIBI_General_byMetricType_updated_30Jan2025.r"))

metlist <- read_csv(paste0(here::here(), "/Testing/vegmetlist_bytype.csv")) 
mets <- read_csv(paste0(here::here(), "/Testing/Veg_scored_mets.csv")) |> 
  mutate(indvis = if_else(VISIT_NO=='1', 'Yes', 'No'))

prIBI_test <- prIBI_byType(df = mets, 
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
