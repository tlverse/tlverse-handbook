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


## ----load_ist_data_intro, message=FALSE, warning=FALSE------------------------
# read in data
ist <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
    "data/ist_sample.csv"
  )
)


## ----skim_ist_data, results="asis", echo=FALSE--------------------------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(ist), format = "latex")
} else {
  skim(ist)
}


## ----load_nhefs_data_intro----------------------------------------------------
# read in data
nhefs_data <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
    "data/NHEFS.csv"
  )
)


## ----skim_nhefs_data, results="asis", echo=FALSE------------------------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(nhefs_data), format = "latex")
} else {
  skim(nhefs_data)
}
