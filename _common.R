# same seed across chapters
library(methods)
set.seed(34729)

hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x <- knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x <- strwrap(x, width = n)
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "",
  collapse = TRUE,
  cache = FALSE,
  tidy = FALSE,
  out.width = "80%",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  fig.show = "hold",
  dpi = 600,
  message = FALSE,
  warning = FALSE,
  echo = TRUE,
  linewidth = 80
)

# global options
options(
  digits = 4,
  # scipen = 999,
  dplyr.print_min = 5,
  dplyr.print_max = 5,
  crayon.enabled = FALSE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
  htmltools.dir.version = FALSE, 
  conflicts.policy = FALSE,
  dplyr.summarise.inform = FALSE
)

# overwrite options by output type
if (knitr:::is_html_output()) {
  options(width = 80)
  options(digits = 4)
}
if (knitr:::is_latex_output()) {
  knitr::opts_chunk$set(width = 67)
  options(width = 67)
  options(cli.unicode = TRUE)
  options(digits = 4)
}

# helper for simpler skimr tables
skim_no_sparks <- skimr::skim_with(
  numeric = skimr::sfl(hist = NULL),
  ts = skimr::sfl(line_graph = NULL)
)

# automatically create a bib database for R packages
knitr::write_bib(c(
  .packages(), "bookdown", "knitr", "rmarkdown"
), "packages.bib")

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)

  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background <- trans_rect
  ret$plot.background <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key <- trans_rect

  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 20)
)
