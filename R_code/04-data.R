## ----load_washb_data_intro----------------------------------------------------
library(readr)
# read in data via readr::read_csv
dat <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  )
)


## ----skim_washb_data, results="asis", echo=FALSE------------------------------
library(skimr)
# optionally disable sparkline graphs for PDF output
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(dat), format = "latex")
} else {
  skim(dat)
}

