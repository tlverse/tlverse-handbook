library(here)
library(knitr)
library(stringr)

# get list of chapters for which to create .R files
chapters <- setdiff(
  str_subset(dir(), ".Rmd"),
  c("index.Rmd", "references.Rmd", "tlverse.Rmd")
)

# create .R files
lapply(chapters, function(f) {
  purl(f, output = here("R_code", str_replace(f, ".Rmd", ".R")))
})
