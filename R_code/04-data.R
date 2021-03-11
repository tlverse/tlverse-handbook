## ----load_washb_data_intro----------------------------------------------------
library(readr)
# read in data via readr::read_csv
dat <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  )
)
head(dat)


## ----skim_washb_data, results="asis"------------------------------------------
library(skimr)
# optionally disable sparkline graphs for PDF output
skim_no_sparks <- skim_with(
  numeric = sfl(hist = NULL),
  ts = sfl(line_graph = NULL)
)
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(dat))
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
head(ist)


## ----skim_ist_data, results="asis"--------------------------------------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(ist))
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
head(nhefs_data)


## ----skim_nhefs_data, results="asis"------------------------------------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(nhefs_data))
} else {
  skim(nhefs_data)
}
