# Welcome to the `tlverse` {#tlverse}

Updated: 2022-03-03



## Learning Objectives {-}

This chapter introduces the `tlverse` software ecosystem, including

1. Understanding the `tlverse` ecosystem conceptually.
2. Identifying the core components of the `tlverse`.
3. Installing `tlverse` `R` packages.
4. Understanding the Targeted Learning roadmap.
5. Learning about the WASH Benefits example data.



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

All Targeted Learning methods are targeted maximum likelihood (or minimum
loss-based) estimators (TMLEs). The construction of any Targeted Learning
estimator proceeds through a two-stage process:

1. Flexibly learning particular components of the data-generating distribution
   through macchine learning (e.g., Super Learning), resulting in _initial
   estimates_ of nuisance parameters.
2. Use of a parametric model-based update via maximum likelihood estimation
   (i.e., MLE), incorporating the initial estimates produced by the prior step.

The packages making up the core components of the `tlverse` software ecosystem,
`sl3` and `tmle3`, address the above two goals, respectively. Together, the very
general functionality exposed by both allows one to build specific TMLEs
tailored exactly to a particular estimation problem.

The software packages that make up the **core** of the `tlverse` are

* [`sl3`](https://github.com/tlverse/sl3): Modern Super Machine Learning
  * _What?_ A modern object-oriented re-implementation of the Super Learner
    algorithm, employing recently developed paradigms in `R` programming.
  * _Why?_ A design that leverages modern ideas for faster computation, is
    easily extensible and forward-looking, and forms one of the cornerstones of
    the `tlverse`.

* [`tmle3`](https://github.com/tlverse/tmle3): An Engine for Targeted Learning
  * _What?_ A generalized framework that simplifies Targeted Learning by
    identifying and implementing a series of common statistical estimation
    procedures.
  * _Why?_ A common interface and engine that accommodates current algorithmic
    approaches to Targeted Learning and yet remains a flexible enough engine to
    power the implementation of emerging statistical techniques as they are
    developed.

Beyond these engines that provide the driving force behind the `tlverse`, there
are a few supporting packages that play important roles in the background:

* [`origami`](https://github.com/tlverse/origami): A Generalized Framework for
   Cross-Validation [@coyle2018origami]
  * _What?_ A generalized framework for flexible cross-validation.
  * _Why?_ Cross-validation is a key part of ensuring error estimates are honest
    and in preventing overfitting. It is an essential part of the both the Super
    Learner ensemble modeling algorithm and in the construction of Targeted
    Learning estimators.

* [`delayed`](https://github.com/tlverse/delayed): Parallelization Framework for
   Dependent Tasks
  * _What?_ A framework for delayed computations (i.e., futures) based on task
    dependencies.
  * _Why?_ Efficient allocation of compute resources is essential when deploying
    computationally intensive algorithms at large scale.

A key principle of the `tlverse` is extensibility. That is, the software
ecosystem aims to support the development of new Targeted Learning estimators as
they reaching maturity. To achieve this degree of flexibility, we follow the
model of implementing new classes of estimators, for distinct causal inference
problems, in separate packages, all of which use the core machinery provided by
the `sl3` and `tmle3` packages There are currently three examples:

* [`tmle3mopttx`](https://github.com/tlverse/tmle3mopttx): Optimal Treatments
  in the `tlverse`
  * _What?_ Learn an optimal rule and estimate the mean outcome under the rule.
  * _Why?_ Optimal treatments are a powerful tool in precision healthcare and
    other settings where a one-size-fits-all treatment approach is not
    appropriate.

* [`tmle3shift`](https://github.com/tlverse/tmle3shift): Stochastic Shift
  Interventions in the `tlverse`
  * _What?_ Stochastic shift interventions for continuous-valued treatments.
  * _Why?_ Not all treatment variables are binary or categorical. Estimating the
    total effects of intervening on continuous-valued treatments provides a way
    to probe how an effect changes with shifts in the treatment variable.

* [`tmle3mediate`](https://github.com/tlverse/tmle3mediate): Causal Mediation
  Analysis in the `tlverse`
  * _What?_ Techniques for evaluating the direct and indirect effects of
    treatments through mediating variables.
  * _Why?_ Evaluating the total effect of a treatment does not provide
    information about the pathways through which it may operate. When mediating
    variables have been collected, one can instead evaluate direct and indirect
    effect parameters that speak to the _action mechanism_ of the treatment.

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
  GITHUB_PAT <- yourPAT
  ```

8. Restart R. You can restart R via the drop-down menu on RStudio's "Session"
   tab, which is located at the top of the RStudio interface. You have to
   restart R for the changes to take effect!

After following these steps, you should be able to successfully install the
package which threw the error above.
