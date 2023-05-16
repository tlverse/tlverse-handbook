---
title: "Targeted Learning in R"
subtitle: "Causal Data Science with the tlverse Software Ecosystem"
author: "Mark van der Laan, Jeremy Coyle, Nima Hejazi, Ivana Malenica, Rachael
  Phillips, Alan Hubbard"
date: "May 16, 2023"
documentclass: krantz
lof: yes
fontsize: '12pt, krantz2'
monofontoptions: "Scale=0.7"
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
colorlinks: yes
site: bookdown::bookdown_site
description: "An open source handbook for causal machine learning and data science with the Targeted Learning framework using the [`tlverse` software ecosystem](https://github.com/tlverse)."
favicon: "img/favicons/favicon.png"
github-repo: tlverse/tlverse-handbook
graphics: yes
url: 'https\://tlverse.org/tlverse-handbook/'
twitter-handle: tlverse
---

# Welcome {-}

_Targeted Learning in `R`: Causal Data Science with the `tlverse` Software
Ecosystem_ is an fully reproducible, open source, electronic handbook for
applying Targeted Learning methodology in practice using the software stack
provided by the [`tlverse` ecosystem](https://github.com/tlverse). This work is
a draft phase and is publicly available to solicit input from the community. To
view or contribute, visit the [GitHub
repository](https://github.com/tlverse/tlverse-handbook).

<!--- For HTML Only --->

<img style="float: left; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/tlverse-logo.svg" width="30%" height="30%">
<img style="float: center; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/Rlogo.svg" width="35%" height="35%">
<img style="float: right; margin-right: 1%; margin-bottom: 0.01em"
     src="img/logos/vdl-logo-transparent.svg" width="30%" height="30%">
<p style="clear: both;">
<br>


## Outline {-}

The contents of this handbook are meant to serve as a reference guide for both
applied research and for the teaching of short courses illustrating successful
applications of the Targeted Learning statistical paradigm. Each section
introduces a set of distinct causal inference questions, often motivated by a
case study, alongside statistical methodology and open source software for
assessing the scientific (causal) claim of interest. The set of materials
currently includes

* Motivation: [Why we need a statistical
    revolution](https://senseaboutscienceusa.org/super-learning-and-the-revolution-in-knowledge/)
* The Roadmap and introductory case study: the WASH Benefits Bangladesh dataset
* Introduction to the [`tlverse` software
    ecosystem](https://tlverse.org)
* Cross-validation with the [`origami`](https://github.com/tlverse/origami)
    package
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
* _Coda_: [Why we need a statistical
    revolution](https://senseaboutscienceusa.org/super-learning-and-the-revolution-in-knowledge/)

## What this book is not {-}

This book does __not__ focus on providing in-depth technically sophisticated
descriptions of modern statistical methodology or recent advancements in
Targeted Learning. Instead, the goal is to convey key details of these
state-of-the-art statistical techniques in a manner that is clear, complete, and
intuitive, while simultaneously avoiding the cognitive burden carried by
extraneous details (e.g., mathematically niche theoretical arguments).  Our aim
is for the presentations herein to serve as a coherent reference for researchers
-- applied methodologists and domain specialists alike -- that empower them to
deploy the central statistical tools of Targeted Learning in a manner efficient
for their scientific pursuits.  For a mathematically sophisticated treatment of
some of these topics, inclusive of in-depth technical details, in the field of
Targeted Learning, the interested reader is invited to consult @vdl2011targeted
and @vdl2018targeted, among numerous other works, as appropriate. The primary
literature in causal inference, machine learning, and non/semi-parametric
statistical theory include many of the most recent advances in Targeted Learning
and related areas. For background in causal inference, @hernan2022causal serves
as an introductory modern reference.

## About the authors {-}

### Mark van der Laan {-}

Mark van der Laan is Professor of Biostatistics and of Statistics at UC Berkeley
and co-director of UC Berkeley's [Center for Targeted Machine Learning and
Causal Inference](https://ctml.berkeley.edu/). His research interests include
statistical methods in computational biology, survival analysis, censored data,
adaptive designs, targeted maximum likelihood estimation, causal inference,
data-adaptive loss-based learning, and multiple testing. His research group
developed loss-based super learning in semiparametric models, based on
cross-validation, as a generic optimal tool for the estimation of
infinite-dimensional parameters, such as nonparametric density estimation and
prediction with both censored and uncensored data. Building on this work, his
research group developed targeted maximum likelihood estimation for a target
parameter of the data-generating distribution in arbitrary semiparametric and
nonparametric models, as a generic optimal methodology for statistical and
causal inference. Since mid-2017, Mark's group has focused in part on the
development of a centralized, principled set of software tools for targeted
learning, the `tlverse`.

### Jeremy Coyle {-}

Jeremy Coyle, PhD, is a consulting data scientist and statistical programmer,
currently leading the software development effort that has produced the
`tlverse` ecosystem of R packages and related software tools. Jeremy earned his
PhD in Biostatistics from UC Berkeley in 2017, primarily under the supervision
of Alan Hubbard.

### Nima Hejazi {-}

[Nima Hejazi](https://nimahejazi.org) is Assistant Professor of Biostatistics at
the [Harvard T.H. Chan School of Public
Health](https://www.hsph.harvard.edu/biostatistics/). He obtained his PhD in
Biostatistics at UC Berkeley working with Mark van der Laan and Alan Hubbard,
and held an NSF Mathematical Sciences Postdoctoral Research Fellowship
afterwards. Nima's research interests blend causal inference, machine learning,
non/semi-parametric inference, and computational statistics; areas of recent
emphasis have included nonparametric causal mediation analysis, efficient
estimation under biased sampling designs, and sieve estimation with machine
learning. His methodological work is motivated principally by scientific
collaborations in clinical trials, infectious disease epidemiology, and
computational biology. Nima is also passionate about statistical computing and
open source software design for statistical data science.

### Ivana Malenica {-}

[Ivana Malenica](https://imalenica.github.io/) is a Postdoctoral Researcher in
the [Department of Statistics](https://statistics.fas.harvard.edu/) at Harvard
and a Wojcicki and Troper Data Science Fellow at the [Harvard Data Science
Initiative](https://datascience.harvard.edu/). She obtained her PhD in
Biostatistics at UC Berkeley working with Mark van der Laan, where she was a
Berkeley Institute for Data Science and a NIH Biomedical Big Data Fellow. Her
research interests span non/semi-parametric theory, causal inference and machine
learning, with emphasis on personalized health and dependent settings. Most of
her current work involves causal inference with time and network dependence,
online learning, optimal individualized treatment, reinforcement learning, and
adaptive sequential designs.

### Rachael Phillips {-}

Rachael Phillips is a PhD student in biostatistics, advised by Alan Hubbard and
Mark van der Laan. She has an MA in Biostatistics, BS in Biology, and BA in
Mathematics. As a student of targeted learning, Rachael integrates causal
inference, machine learning, and statistical theory to answer causal questions
with statistical confidence. She is motivated by issues arising in healthcare,
and is especially interested in clinical algorithm frameworks and guidelines.
Related to to this, she is also interested in experimental design;
human-computer interaction; statistical analysis pre-specification, automation,
and reproducibility; and open-source software.

### Alan Hubbard {-}

Alan Hubbard is Professor of Biostatistics at UC Berkeley, current chair of the
Division of Biostatistics of the UC Berkeley School of Public Health, head of
the data analytics core of UC Berkeley's SuperFund research program, and
co-director of UC Berkeley's [Center for Targeted Machine Learning and Causal
Inference](https://ctml.berkeley.edu/). His current research interests include
causal inference, variable importance analysis, statistical machine learning,
estimation of and inference for data-adaptive statistical target parameters, and
targeted minimum loss-based estimation. Research in his group is generally
motivated by applications to problems in computational biology, epidemiology,
and precision medicine.

<!--
# Acknowledgements {-}
-->



## Reproduciblity {-}

The `tlverse` software ecosystem is a growing collection of packages, several of
which are quite early on in the software lifecycle. The team does its best to
maintain backwards compatibility. Once this work reaches completion, the
specific versions of the `tlverse` packages used will be archived and tagged to
produce it.

This book was written using [bookdown](http://bookdown.org/), and the complete
source is available on [GitHub](https://github.com/tlverse/tlverse-handbook).
This version of the book was built with R version 4.3.0 (2023-04-21),
[pandoc](https://pandoc.org/) version 2.19.2, and the
following packages:


|package      |version    |source                                                                  |
|:------------|:----------|:-----------------------------------------------------------------------|
|bookdown     |0.26.3     |Github (rstudio/bookdown\@169c43b6bb95213f2af63a95acd4e977a58a3e1f)     |
|bslib        |0.3.1      |CRAN (R 4.3.0)                                                          |
|dagitty      |0.3-1      |CRAN (R 4.3.0)                                                          |
|data.table   |1.14.2     |CRAN (R 4.3.0)                                                          |
|delayed      |0.3.0      |CRAN (R 4.3.0)                                                          |
|downlit      |0.4.0      |CRAN (R 4.3.0)                                                          |
|dplyr        |1.0.9      |CRAN (R 4.3.0)                                                          |
|forecast     |8.16       |CRAN (R 4.3.0)                                                          |
|future       |1.26.1     |CRAN (R 4.3.0)                                                          |
|ggdag        |0.2.4      |CRAN (R 4.3.0)                                                          |
|ggfortify    |0.4.14     |CRAN (R 4.3.0)                                                          |
|ggplot2      |3.3.6      |CRAN (R 4.3.0)                                                          |
|kableExtra   |1.3.4.9000 |Github (kupietz/kableExtra\@3bf9b21a769c9e6c21c955689bf5f8175dc83350)   |
|knitr        |1.42       |CRAN (R 4.3.0)                                                          |
|mvtnorm      |1.1-3      |CRAN (R 4.3.0)                                                          |
|origami      |1.0.5      |Github (tlverse/origami\@e1b8fe6f5e75fff1d48eed115bb81475c9bd506e)      |
|randomForest |4.7-1.1    |CRAN (R 4.3.0)                                                          |
|readr        |2.1.2      |CRAN (R 4.3.0)                                                          |
|rmarkdown    |2.14       |CRAN (R 4.3.0)                                                          |
|skimr        |2.1.4      |CRAN (R 4.3.0)                                                          |
|sl3          |1.4.5      |Github (tlverse/sl3\@de445c210eefa5aa9dd4c0d1fab8126f0d7c5eeb)          |
|stringr      |1.4.0      |CRAN (R 4.3.0)                                                          |
|tibble       |3.1.7      |CRAN (R 4.3.0)                                                          |
|tidyr        |1.2.0      |CRAN (R 4.3.0)                                                          |
|tmle3        |0.2.0      |Github (tlverse/tmle3\@ed72f8a20e64c914ab25ffe015d865f7a9963d27)        |
|tmle3mediate |0.0.3      |Github (tlverse/tmle3mediate\@70d1151c4adb54d044f355d06d07bcaeb7f8ae07) |
|tmle3mopttx  |1.0.0      |Github (tlverse/tmle3mopttx\@c8c675f051bc5ee6d51fa535fe6dc80791d4d1b7)  |
|tmle3shift   |0.2.0      |Github (tlverse/tmle3shift\@4ed52b50af501a5fa2e6257b568d17fd485d3f42)   |



## Learning resources {-}

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

For a general, modern introduction to causal inference, we recommend

* [Miguel A. Hernán and James M. Robins' _Causal Inference: What If_
    (2022)](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/)
* [Jason A. Roy's _A Crash Course in Causality: Inferring Causal Effects from
  Observational Data_ on
  Coursera](https://www.coursera.org/learn/crash-course-in-causality)

Feel free to [suggest a
resource](https://github.com/tlverse/tlverse-handbook/issues)!

## Want to help? {-}

Any feedback on the book is very welcome. Feel free to [open an
issue](https://github.com/tlverse/tlverse-handbook/issues), or to make a Pull
Request if you spot a typo.
