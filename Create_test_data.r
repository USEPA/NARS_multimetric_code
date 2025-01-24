library(tidyverse)
library(tibble)

# Create data for the range test that will trigger each type of fail or pass
range_test_data <- tribble(
  ~Site, ~EPTRICH, ~EPTPIND, ~CHIRPIND, ~CHIRPTAX, ~OLIGRICH, 
  "site1", 10, 50, 10, 25.5, 0,
  "site2", 9, 25, 15, 60, 0, 
  "site3", 5, 20, 12, 90, 0,
  "site4", 5, 20, 11, 75, 0,
  "site5", 1, 25, 13, 60, 0, 
  "site6", 0, 0, 12, 80, 1,
  "site7", 2, 20, 10, 75, 0,
  "site8", 10, 45, 12, 30, 0,
  "site9", 6, 60, 8, 10, 30,
  
)