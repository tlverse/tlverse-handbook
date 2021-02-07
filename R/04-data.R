## ----load_washb_data_intro, message=FALSE, warning=FALSE----------------------
library(readr)
# read in data via readr::read_csv
dat <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv")
head(dat)


## ----skim_washb_data, message=FALSE, warning=FALSE, results="asis"------------
library(skimr)
# optionally disable sparkline graphs for PDF output
skim_no_sparks <- skim_with(numeric = sfl(hist = NULL),
                            ts = sfl(line_graph = NULL))
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(dat))
} else {
  skim(dat)
}


## ----load_ist_data_intro, message=FALSE, warning=FALSE------------------------
# read in data
ist <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/data/ist_sample.csv")
head(ist)


## ----skim_ist_data, message=FALSE, warning=FALSE, results="asis"--------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(ist))
} else {
  skim(ist)
}


## ----load_nhefs_data_intro, message=FALSE, warning=FALSE----------------------
# read in data
NHEFS <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/data/NHEFS.csv")
head(NHEFS)


## ----skim_nhefs_data, message=FALSE, warning=FALSE, results="asis"------------
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(NHEFS))
} else {
  skim(NHEFS)
}

