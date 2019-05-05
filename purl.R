library(here)
library(knitr)
library(tidyverse)

# get list of chapters for which to create .R files
chapters <- setdiff(str_subset(dir(), ".Rmd"),
                    c("index.Rmd", "99-references.Rmd"))

# create .R files
lapply(chapters, function(f) {
  purl(f, output = here("R", str_replace(f, ".Rmd", ".R")))
})
