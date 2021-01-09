# Set output options
if (knitr:::is_html_output()) {
  options(width = 80)
}
if (knitr:::is_latex_output()) {
  options(width = 65)
}

options(digits = 7, bookdown.clean_book = TRUE, knitr.kable.NA = "NA")

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  tidy = FALSE,
  fig.align = "center",
  out.width = "\\textwidth"
)

options(dplyr.print_min = 6, dplyr.print_max = 6)

# automatically create a bib database for R packages
knitr::write_bib(c(
  .packages(), "bookdown", "knitr", "rmarkdown"
), "packages.bib")

# same seed across chapters
set.seed(34729)
