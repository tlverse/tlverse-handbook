# Welcome to the `tlverse` {#tlverse}

## Learning Objectives {-}

1. Understand the `tlverse` ecosystem conceptually
2. Identify the core components of the `tlverse`
3. Install `tlverse` `R` packages
4. Understand the Targeted Learning roadmap
5. Learn about the WASH Benefits example data

## What is the `tlverse`? {-}

The `tlverse` is a new framework for doing Targeted Learning in R, inspired by
the [`tidyverse` ecosystem](https://tidyverse.org) of R packages.

By analogy to the [`tidyverse`](https://tidyverse.org/):

> The `tidyverse` is an opinionated collection of R packages designed for data
> science. All packages share an underlying design philosophy, grammar, and data
> structures.

So, the [`tlverse`](https://tlverse.org) is

* an opinionated collection of R packages for Targeted Learning
* sharing an underlying philosophy, grammar, and set of data structures

## Anatomy of the `tlverse` {-}

These are the main packages that represent the **core** of the `tlverse`:

* [`sl3`](https://github.com/tlverse/sl3): Modern Super Learning with Pipelines
  * _What?_ A modern object-oriented re-implementation of the Super Learner
    algorithm, employing recently developed paradigms for `R` programming.
  * _Why?_ A design that leverages modern tools for fast computation, is
    forward-looking, and can form one of the cornerstones of the `tlverse`.

* [`tmle3`](https://github.com/tlverse/tmle3): An Engine for Targeted Learning
  * _What?_ A generalized framework that simplifies Targeted Learning by
    identifying and implementing a series of common statistical estimation
    procedures.
  * _Why?_ A common interface and engine that accommodates current algorithmic
    approaches to Targeted Learning and is still flexible enough to remain the
    engine even as new techniques are developed.

In addition to the engines that drive development in the `tlverse`, there are
some supporting packages -- in particular, we have two...

* [`origami`](https://github.com/tlverse/origami): A Generalized Framework for
   Cross-Validation
  * _What?_ A generalized framework for flexible cross-validation
  * _Why?_ Cross-validation is a key part of ensuring error estimates are honest
    and preventing overfitting. It is an essential part of the both the Super
    Learner algorithm and Targeted Learning.

* [`delayed`](https://github.com/tlverse/delayed): Parallelization Framework for
   Dependent Tasks
  * _What?_ A framework for delayed computations (futures) based on task
    dependencies.
  * _Why?_ Efficient allocation of compute resources is essential when deploying
    large-scale, computationally intensive algorithms.

A key principle of the `tlverse` is extensibility. That is, we want to support
new Targeted Learning estimators as they are developed. The model for this is
new estimators are implemented in additional packages using the core packages
above. There are currently two featured examples of this:

* [`tmle3mopttx`](https://github.com/tlverse/tmle3mopttx): Optimal Treatments
  in `tlverse`
  * _What?_ Learn an optimal rule and estimate the mean outcome under the rule
  * _Why?_ Optimal Treatment is a powerful tool in precision healthcare and
    other settings where a one-size-fits-all treatment approach is not
    appropriate.

* [`tmle3shift`](https://github.com/tlverse/tmle3shift): Shift Interventions in
  `tlverse`
  * _What?_ Shift interventions for continuous treatments
  * _Why?_ Not all treatment variables are discrete. Being able to estimate the
    effects of continuous treatment represents a powerful extension of the
    Targeted Learning approach.

## Installation {#installtlverse}

The `tlverse` ecosystem of packages are currently hosted at
https://github.com/tlverse, not yet on [CRAN](https://CRAN.R-project.org/). You
can use the [`usethis` package](https://usethis.r-lib.org/) to install them:


```r
install.packages("devtools")
devtools::install_github("tlverse/tlverse")
```

The `tlverse` depends on a large number of other packages that are also hosted
on GitHub. Because of this, you may see the following error:

```
Error: HTTP error 403.
  API rate limit exceeded for 71.204.135.82. (But here's the good news:
  Authenticated requests get a higher rate limit. Check out the documentation
  for more details.)

  Rate limit remaining: 0/60
  Rate limit reset at: 2019-03-04 19:39:05 UTC

  To increase your GitHub API rate limit
  - Use `usethis::browse_github_pat()` to create a Personal Access Token.
  - Use `usethis::edit_r_environ()` and add the token as `GITHUB_PAT`.
```

This just means that R tried to install too many packages from GitHub in too
short of a window. To fix this, you need to tell R how to use GitHub as your
user (you'll need a GitHub user account). Follow these two steps:

1. Type `usethis::browse_github_pat()` in your R console, which will direct
   you to GitHub's page to create a New Personal Access Token (PAT).
2. Create a PAT simply by clicking "Generate token" at the bottom of the page.
3. Copy your PAT, a long string of lowercase letters and numbers.
4. Type `usethis::edit_r_environ()` in your R console, which will open your
   `.Renviron` file in the source window of RStudio. 
   
    a. If your `.Renviron` file does not pop-up after calling 
       `usethis::edit_r_environ()`; then try inputting 
       `Sys.setenv(GITHUB_PAT = "yourPAT")`, replacing your PAT with inside the 
       quotes. If this does not error, then skip to step 8.
       
5. In your `.Renviron` file, type `GITHUB_PAT=` and then paste your PAT after 
   the equals symbol with no space.
6. In your `.Renviron` file, press the enter key to ensure that your `.Renviron`
   ends with a new line.
7. Save your `.Renviron` file. The example below shows how this syntax should 
   look.
   

```r
GITHUB_PAT=yourPAT

```
8. Restart R. You can restart R via the drop-down menu on RStudio's "Session" 
   tab, which is located at the top of the RStudio interface. You have to 
   restart R for the changes to take effect!
   
After following these steps, you should be able to successfully install the
package which threw the error above.
