## ----load_washb_data_intro, message=FALSE, warning=FALSE----------------------
library(tidyverse)
# read in data
dat <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv")
dat


## ----skim_washb_data, message=FALSE, warning=FALSE----------------------------
library(skimr)
skim(dat)


## ----load_ist_data_intro, message=FALSE, warning=FALSE------------------------
library(tidyverse)
# read in data
ist <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/data/ist_sample.csv")
ist


## ----skim_ist_data, message=FALSE, warning=FALSE------------------------------
skim(ist)


## ----load_nhefs_data_intro, message=FALSE, warning=FALSE----------------------
library(tidyverse)
# read in data
NHEFS <- read_csv("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/data/NHEFS.csv")
NHEFS


## ----skim_nhefs_data, message=FALSE, warning=FALSE----------------------------
skim(NHEFS)

