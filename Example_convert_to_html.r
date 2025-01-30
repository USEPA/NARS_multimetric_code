# First create .Rd files
devtools::document()

# Convert .Rd files to html
tools::Rd2HTML(Rd = paste0(here::here(), "/man/range_test.Rd"), out = paste0(here::here(), "/html/range_test.html"))

tools::Rd2HTML(Rd = paste0(here::here(), "/man/comp_test.Rd"), out = paste0(here::here(), "/html/comp_test.html"))

tools::Rd2HTML(Rd = paste0(here::here(), "/man/sn_test.Rd"), out = paste0(here::here(), "/html/sn_test.html"))

tools::Rd2HTML(Rd = paste0(here::here(), "/man/redund_test.Rd"), out = paste0(here::here(), "/html/redund_test.html"))

tools::Rd2HTML(Rd = paste0(here::here(), "/man/prIBI_byType.Rd"), out = paste0(here::here(), "/html/prIBI_byType.html"))
