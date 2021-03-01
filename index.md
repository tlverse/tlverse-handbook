---
knit: "bookdown::render_book"
title: "Targeted Learning in R"
subtitle: "Causal Data Science with the tlverse Software Ecosystem"
author: "Mark van der Laan, Jeremy Coyle, Nima Hejazi, Ivana Malenica, Rachael
  Phillips, Alan Hubbard"
date: "March 01, 2021"
documentclass: book
site: bookdown::bookdown_site
bibliography: [book.bib, packages.bib]
biblio-style: apalike
fontsize: '12pt, krantz2'
monofontoptions: "Scale=0.7"
link-citations: yes
colorlinks: yes
lot: yes
lof: yes
always_allow_html: yes
url: 'https\://tlverse.org/tlverse-handbook/'
github-repo: tlverse/tlverse-handbook
twitter-handle: tlverse
graphics: yes
description: "An open source handbook for causal machine learning and data
  science with the Targeted Learning framework using the [`tlverse` software
  ecosystem](https://github.com/tlverse)."
favicon: "img/logos/favicons/favicon.png"
---

# About this book {-}

_Targeted Learning in `R`: Causal Data Science with the `tlverse` Software
Ecosystem_ is an open source, reproducible electronic handbook for applying the
Targeted Learning methodology in practice using the [`tlverse` software
ecosystem](https://github.com/tlverse). This work is currently in an early draft
phase and is available to facilitate input from the community. To view or
contribute to the available content, consider visiting the [GitHub
repository](https://github.com/tlverse/tlverse-handbook).

<img style="float: left; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/tlverse-logo.svg" width="30%" height="30%">
<img style="float: center; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/Rlogo.svg" width="35%" height="35%">
<img style="float: right; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/vdl-logo-transparent.svg" width="30%" height="30%">
<p style="clear: both;">
<br>

## Outline {#outline}

The contents of this handbook are meant to serve as a reference guide for
applied research as well as materials that can be taught in a series of short
courses focused on the applications of Targeted Learning. Each section
introduces a set of distinct causal questions, motivated by a case study,
alongside statistical methodology and software for assessing the causal claim of
interest. The (evolving) set of materials includes

* Motivation: [Why we need a statistical
    revolution](https://senseaboutscienceusa.org/super-learning-and-the-revolution-in-knowledge/)
* The Roadmap and introductory case study: the WASH Beneifits data
* Introduction to the [`tlverse` software
    ecosystem](https://tlverse.org)
* Cross-validation with the [`origami`](https://github.com/tlverse/origami)
    package (_work in progress_)
* Ensemble machine learning with the
    [`sl3`](https://github.com/tlverse/sl3) package
* Targeted learning for causal inference with the
    [`tmle3`](https://github.com/tlverse/tmle3) package
* Optimal treatments regimes and the
    [`tmle3mopttx`](https://github.com/tlverse/tmle3mopttx) package
* Stochastic treatment regimes and the
    [`tmle3shift`](https://github.com/tlverse/tmle3shift) package
* Causal mediation analysis with the
    [`tmle3mediate`](https://github.com/tlverse/tmle3mediate) package
    (_work in progress_)
* _Coda_: [Why we need a statistical
    revolution](https://senseaboutscienceusa.org/super-learning-and-the-revolution-in-knowledge/)

## What this book is not {-}

The focus of this work is __not__ on providing in-depth technical descriptions
of current statistical methodology or recent advancements. Instead, the goal is
to convey key details of state-of-the-art techniques in an manner that is both
clear and complete, without burdening the reader with extraneous information.
We hope that the presentations herein will serve as references for researchers
-- methodologists and domain specialists alike -- that empower them to deploy
the central tools of Targeted Learning in an efficient manner. For technical
details and in-depth descriptions of both classical theory and recent advances
in the field of Targeted Learning, the interested reader is invited to consult
@vdl2011targeted and/or @vdl2018targeted as appropriate. The primary literature
in statistical causal inference, machine learning, and non/semiparametric theory
include many of the most recent advances in Targeted Learning and related areas.

## About the authors {-}

### Mark van der Laan {-}

Mark van der Laan, PhD, is Professor of Biostatistics and Statistics at UC
Berkeley. His research interests include statistical methods in computational
biology, survival analysis, censored data, adaptive designs, targeted maximum
likelihood estimation, causal inference, data-adaptive loss-based learning, and
multiple testing. His research group developed loss-based super learning in
semiparametric models, based on cross-validation, as a generic optimal tool for
the estimation of infinite-dimensional parameters, such as nonparametric density
estimation and prediction with both censored and uncensored data. Building on
this work, his research group developed targeted maximum likelihood estimation
for a target parameter of the data-generating distribution in arbitrary
semiparametric and nonparametric models, as a generic optimal methodology for
statistical and causal inference. Most recently, Mark's group has focused in
part on the development of a centralized, principled set of software tools for
targeted learning, the `tlverse`.

### Jeremy Coyle {-}

Jeremy Coyle, PhD, is a consulting data scientist and statistical programmer,
currently leading the software development effort that has produced the
`tlverse` ecosystem of R packages and related software tools. Jeremy earned his
PhD in Biostatistics from UC Berkeley in 2016, primarily under the supervision
of Alan Hubbard.

### Nima Hejazi {-}

Nima Hejazi is a PhD candidate in biostatistics, working under the collaborative
direction of Mark van der Laan and Alan Hubbard. Nima is affiliated with UC
Berkeley's Center for Computational Biology and NIH Biomedical Big Data training
program, as well as with the Fred Hutchinson Cancer Research Center. Previously,
he earned an MA in Biostatistics and a BA (with majors in Molecular and Cell
Biology, Psychology, and Public Health), both at UC Berkeley.  His research
interests fall at the intersection of causal inference and machine learning,
drawing on ideas from non/semi-parametric estimation in large, flexible
statistical models to develop efficient and robust statistical procedures for
evaluating complex target estimands in observational and randomized studies.
Particular areas of current emphasis include mediation/path analysis,
outcome-dependent sampling designs, targeted loss-based estimation, and vaccine
efficacy trials.  Nima is also passionate about statistical computing and open
source software development for applied statistics.

### Ivana Malenica {-}

Ivana Malenica is a PhD student in biostatistics advised by Mark van der Laan.
Ivana is currently a fellow at the Berkeley Institute for Data Science, after
serving as a NIH Biomedical Big Data and Freeport-McMoRan Genomic Engine fellow.
She earned her Master's in Biostatistics and Bachelor's in Mathematics, and
spent some time at the Translational Genomics Research Institute. Very broadly,
her research interests span non/semi-parametric theory, probability theory,
machine learning, causal inference and high-dimensional statistics. Most of her
current work involves complex dependent settings (dependence through time and
network) and adaptive sequential designs.

### Rachael Phillips {-}

Rachael Phillips is a PhD student in biostatistics, advised by Alan Hubbard and
Mark van der Laan. She has an MA in Biostatistics, BS in Biology, and BA in 
Mathematics. A student of targeted learning and causal inference; her research 
integrates personalized medicine, human-computer interaction, experimental 
design, and regulatory policy. 

### Alan Hubbard {-}

Alan Hubbard is Professor of Biostatistics, former head of the Division of
Biostatistics at UC Berkeley, and head of data analytics core at UC Berkeley's
SuperFund research program. His current research interests include causal
inference, variable importance analysis, statistical machine learning,
estimation of and inference for data-adaptive statistical target parameters, and
targeted minimum loss-based estimation. Research in his group is generally
motivated by applications to problems in computational biology, epidemiology,
and precision medicine.

<!--
# Acknowledgements {-}
-->

## Reproduciblity with the `tlverse` {#repro}

The `tlverse` software ecosystem is a growing collection of packages, several of
which are quite early on in the software lifecycle. The team does its best to
maintain backwards compatibility. Once this work reaches completion, the
specific versions of the `tlverse` packages used will be archived and tagged to
produce it.

This book was written using [bookdown](http://bookdown.org/), and the complete
source is available on [GitHub](https://github.com/tlverse/tlverse-handbook).
This version of the book was built with R version 4.0.2 (2020-06-22),
[pandoc](https://pandoc.org/) version 2.2, and the
following packages:


|package      |version    |source                                |
|:------------|:----------|:-------------------------------------|
|bookdown     |0.21.6     |Github (rstudio/bookdown\@ca0145f)    |
|bslib        |0.2.4.9002 |Github (rstudio/bslib\@aa5a842)       |
|dagitty      |0.3-1      |CRAN (R 4.0.2)                        |
|data.table   |1.14.0     |CRAN (R 4.0.2)                        |
|delayed      |0.3.0      |CRAN (R 4.0.2)                        |
|downlit      |0.2.1      |CRAN (R 4.0.2)                        |
|dplyr        |1.0.4      |CRAN (R 4.0.2)                        |
|forecast     |8.13       |CRAN (R 4.0.2)                        |
|ggdag        |0.2.3      |CRAN (R 4.0.2)                        |
|ggfortify    |0.4.11     |CRAN (R 4.0.2)                        |
|ggplot2      |3.3.3      |CRAN (R 4.0.2)                        |
|kableExtra   |1.3.4      |CRAN (R 4.0.2)                        |
|knitr        |1.31       |CRAN (R 4.0.2)                        |
|mvtnorm      |1.1-1      |CRAN (R 4.0.2)                        |
|origami      |1.0.3      |CRAN (R 4.0.2)                        |
|randomForest |4.6-14     |CRAN (R 4.0.2)                        |
|readr        |1.4.0      |CRAN (R 4.0.2)                        |
|rmarkdown    |2.7.2      |Github (rstudio/rmarkdown\@c0b8584)   |
|skimr        |2.1.2      |CRAN (R 4.0.2)                        |
|sl3          |1.4.3      |Github (tlverse/sl3\@b3f3e92)         |
|stringr      |1.4.0      |CRAN (R 4.0.2)                        |
|tibble       |3.1.0      |CRAN (R 4.0.2)                        |
|tidyr        |1.1.2      |CRAN (R 4.0.2)                        |
|tmle3        |0.1.7      |Github (tlverse/tmle3\@5d86bd4)       |
|tmle3mopttx  |0.1.0      |Github (tlverse/tmle3mopttx\@5ba5f65) |
|tmle3shift   |0.1.9      |Github (tlverse/tmle3shift\@daa0f96)  |

## Learning resources {#learn}

To effectively utilize this handbook, the reader need not be a fully trained
statistician to begin understanding and applying these methods. However, it is
highly recommended for the reader to have an understanding of basic statistical
concepts such as confounding, probability distributions, confidence intervals,
hypothesis tests, and regression. Advanced knowledge of mathematical statistics
may be useful but is not necessary. Familiarity with the `R` programming
language will be essential. We also recommend an understanding of introductory
causal inference.

For learning the `R` programming language we recommend the following (free)
introductory resources:

* [Software Carpentry's _Programming with
    `R`_](http://swcarpentry.github.io/r-novice-inflammation/)
* [Software Carpentry's _`R` for Reproducible Scientific
    Analysis_](http://swcarpentry.github.io/r-novice-gapminder/)
* [Garret Grolemund and Hadley Wickham's _`R` for Data
    Science_](https://r4ds.had.co.nz)

For a general introduction to causal inference, we recommend

* [Miguel A. Hernán and James M. Robins' _Causal Inference_, forthcoming
    2020](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/)
* [Jason A. Roy's _A Crash Course in Causality: Inferring Causal Effects from
  Observational Data_ on
  Coursera](https://www.coursera.org/learn/crash-course-in-causality)

## Setup instructions {#setup}

### R and RStudio

**R** and **RStudio** are separate downloads and installations. R is the
underlying statistical computing environment. RStudio is a graphical integrated
development environment (IDE) that makes using R much easier and more
interactive. You need to install R before you install RStudio.

#### Windows

##### If you already have R and RStudio installed

* Open RStudio, and click on "Help" > "Check for updates". If a new version is
  available, quit RStudio, and download the latest version for RStudio.
* To check which version of R you are using, start RStudio and the first thing
  that appears in the console indicates the version of R you are
  running. Alternatively, you can type `sessionInfo()`, which will also display
  which version of R you are running. Go on the [CRAN
  website](https://cran.r-project.org/bin/windows/base/) and check whether a
  more recent version is available. If so, please download and install it. You
  can [check here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#How-do-I-UNinstall-R_003f)
  for more information on how to remove old versions from your system if you
  wish to do so.

##### If you don't have R and RStudio installed

* Download R from
  the [CRAN website](http://cran.r-project.org/bin/windows/base/release.htm).
* Run the `.exe` file that was just downloaded
* Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
* Under *Installers* select **RStudio x.yy.zzz - Windows
  XP/Vista/7/8** (where x, y, and z represent version numbers)
* Double click the file to install it
* Once it's installed, open RStudio to make sure it works and you don't get any
  error messages.

#### macOS / Mac OS X

##### If you already have R and RStudio installed

* Open RStudio, and click on "Help" > "Check for updates". If a new version is
  available, quit RStudio, and download the latest version for RStudio.
* To check the version of R you are using, start RStudio and the first thing
  that appears on the terminal indicates the version of R you are running.
  Alternatively, you can type `sessionInfo()`, which will also display which
  version of R you are running. Go on the [CRAN
  website](https://cran.r-project.org/bin/macosx/) and check whether a more
  recent version is available. If so, please download and install it.

##### If you don't have R and RStudio installed

* Download R from
  the [CRAN website](http://cran.r-project.org/bin/macosx).
* Select the `.pkg` file for the latest R version
* Double click on the downloaded file to install R
* It is also a good idea to install [XQuartz](https://www.xquartz.org/) (needed
  by some packages)
* Go to the [RStudio download
  page](https://www.rstudio.com/products/rstudio/download/#download)
* Under *Installers* select **RStudio x.yy.zzz - Mac OS X 10.6+ (64-bit)**
  (where x, y, and z represent version numbers)
* Double click the file to install RStudio
* Once it's installed, open RStudio to make sure it works and you don't get any
  error messages.

#### Linux

* Follow the instructions for your distribution
  from [CRAN](https://cloud.r-project.org/bin/linux), they provide information
  to get the most recent version of R for common distributions. For most
  distributions, you could use your package manager (e.g., for Debian/Ubuntu run
  `sudo apt-get install r-base`, and for Fedora `sudo yum install R`), but we
  don't recommend this approach as the versions provided by this are
  usually out of date. In any case, make sure you have at least R 3.3.1.
* Go to the [RStudio download
  page](https://www.rstudio.com/products/rstudio/download/#download)
* Under *Installers* select the version that matches your distribution, and
  install it with your preferred method (e.g., with Debian/Ubuntu `sudo dpkg -i
  rstudio-x.yy.zzz-amd64.deb` at the terminal).
* Once it's installed, open RStudio to make sure it works and you don't get any
  error messages.

These setup instructions are adapted from those written for [Data Carpentry: R
for Data Analysis and Visualization of Ecological
Data](http://www.datacarpentry.org/R-ecology-lesson/).
