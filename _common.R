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
  digits = 5,
  # scipen = 999,
  dplyr.print_min = 5,
  dplyr.print_max = 5,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
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

# # helper for simpler skimr tables
# skim_no_sparks <- skimr::skim_with(
#   numeric = skimr::sfl(hist = NULL),
#   ts = skimr::sfl(line_graph = NULL)
# )

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
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

## Schematic Data plots
library(viridis)
library(data.table)

# global plot settings
point_alpha = 1
data(schematic_xlim, schematic_ylim, package="tlverse")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbPalette <- cbbPalette[c(2,3,6,7)]

cbPalette <- viridis_pal()(4)[c(1,3,2,4)]

make_plots <- function(plotobj, name, 
                       height=3, width=5, 
                       make_pdf=FALSE, make_png = FALSE){
  print(plotobj)
  
  if(make_pdf){
    pdf(sprintf("figures/%s.pdf",name),height=height,width=width)
    print(plotobj)
    dev.off()
  }
  
  if(make_png){
    png(sprintf("figures/%s.png",name),height=300*height,width=300*width,res=300)
    print(plotobj)
    dev.off()
  }
}

plot_schematic <- function(data, cf_data, cf_mean){
  
  type_guide <- ifelse(length(unique(cf_data$type))>1,"legend",FALSE)
  plotobj <- ggplot(data,aes(x=W,y=Y, color=factor(A)))+
    geom_point(alpha=point_alpha)+
    geom_line(data=cf_data,aes(linetype=type))+theme_bw()+
    xlab("Covariate (W)")+ylab("Outcome (Y)")+
    scale_color_manual("Treatment Level (A)",values=cbPalette)+
    scale_linetype_discrete("Estimation", guide = type_guide)+
    xlim(schematic_xlim)+ylim(schematic_ylim)+
    geom_rug(data=cf_mean,aes(linetype=type), sides="r")
  return(plotobj)
}

