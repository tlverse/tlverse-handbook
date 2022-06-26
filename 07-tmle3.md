# The TMLE Framework {#tmle3}

_Jeremy Coyle_

Based on the [`tmle3` `R` package](https://github.com/tlverse/tmle3).

## Learning Objectives {#learn-tmle}

By the end of this chapter, you will be able to

1. Understand why we use TMLE for effect estimation.
2. Use `tmle3` to estimate an Average Treatment Effect (ATE).
3. Understand how to use `tmle3` "Specs" objects.
4. Fit `tmle3` for a custom set of target parameters.
5. Use the delta method to estimate transformations of target parameters.

## Introduction {#tmle-intro}

In the previous chapter on `sl3` we learned how to estimate a regression
function like $\mathbb{E}[Y \mid X]$ from data. That's an important first step
in learning from data, but how can we use this predictive model to estimate
statistical and causal effects?

Going back to [the roadmap for targeted learning](#roadmap), suppose we'd like to
estimate the effect of a treatment variable $A$ on an outcome $Y$. As discussed,
one potential parameter that characterizes that effect is the Average Treatment
Effect (ATE), defined as $\psi_0 = \mathbb{E}_W[\mathbb{E}[Y \mid A=1,W] -
\mathbb{E}[Y \mid A=0,W]]$ and interpreted as the difference in mean outcome
under when treatment $A=1$ and $A=0$, averaging over the distribution of
covariates $W$. We'll illustrate several potential estimators for this
parameter, and motivate the use of the TMLE (targeted maximum likelihood
estimation; targeted minimum loss-based estimation) framework, using the
following example data:

<img src="img/png/schematic_1_truedgd.png" width="80%" style="display: block; margin: auto;" />

The small ticks on the right indicate the mean outcomes (averaging over $W$)
under $A=1$ and $A=0$ respectively, so their difference is the quantity we'd
like to estimate.

While we hope to motivate the application of TMLE in this chapter, we refer the
interested reader to the two Targeted Learning books and associated works for
full technical details.

## Substitution Estimators {#substitution-est}

We can use `sl3` to fit a Super Learner or other regression model to estimate
the outcome regression function $\mathbb{E}_0[Y \mid A,W]$, which we often refer
to as $\overline{Q}_0(A,W)$ and whose estimate we denote $\overline{Q}_n(A,W)$.
To construct an estimate of the ATE $\psi_n$, we need only "plug-in" the
estimates of $\overline{Q}_n(A,W)$, evaluated at the two intervention contrasts,
to the corresponding ATE "plug-in" formula:
$\psi_n = \frac{1}{n}\sum(\overline{Q}_n(1,W)-\overline{Q}_n(0,W))$. This kind
of estimator is called a _plug-in_ or _substitution_ estimator, since accurate
estimates $\psi_n$ of the parameter $\psi_0$ may be obtained by substituting
estimates $\overline{Q}_n(A,W)$ for the relevant regression functions
$\overline{Q}_0(A,W)$ themselves.

Applying `sl3` to estimate the outcome regression in our example, we can see
that the ensemble machine learning predictions fit the data quite well:

<img src="img/png/schematic_2b_sllik.png" width="80%" style="display: block; margin: auto;" />

The solid lines indicate the `sl3` estimate of the regression function, with the
dotted lines indicating the `tmle3` updates [(described below)](#tmle-updates).

While substitution estimators are intuitive, naively using this approach with a
Super Learner estimate of $\overline{Q}_0(A,W)$ has several limitations. First,
Super Learner is selecting learner weights to minimize risk across the entire
regression function, instead of "targeting" the ATE parameter we hope to
estimate, leading to biased estimation. That is, `sl3` is trying to do well on
the full regression curve on the left, instead of focusing on the small ticks on
the right. What's more, the sampling distribution of this approach is not
asymptotically linear, and therefore inference is not possible.

We can see these limitations illustrated in the estimates generated for the
example data:

<img src="img/png/schematic_3_effects.png" width="80%" style="display: block; margin: auto;" />

We see that Super Learner, estimates the true parameter value (indicated by the
dashed vertical line) more accurately than GLM. However, it is still less
accurate than TMLE, and valid inference is not possible. In contrast, TMLE
achieves a less biased estimator and valid inference.

## Targeted Maximum Likelihood Estimation {#tmle}

TMLE takes an initial estimate $\overline{Q}_n(A,W)$ as well as an estimate of
the propensity score $g_n(A \mid W) = \mathbb{P}(A = 1 \mid W)$ and produces an
updated estimate $\overline{Q}^{\star}_n(A,W)$ that is "targeted" to the
parameter of interest. TMLE keeps the benefits of substitution estimators (it is
one), but augments the original, potentially erratic estimates to _correct for
bias_ while also resulting in an _asymptotically linear_ (and thus normally
distributed) estimator that accommodates inference via asymptotically consistent
Wald-style confidence intervals.

### TMLE Updates {#tmle-updates}

There are different types of TMLEs (and, sometimes, multiple for the same set of
target parameters) -- below, we give an example of the algorithm for TML
estimation of the ATE.  $\overline{Q}^{\star}_n(A,W)$ is the TMLE-augmented
estimate $f(\overline{Q}^{\star}_n(A,W)) = f(\overline{Q}_n(A,W)) + \epsilon
\cdot H_n(A,W)$, where $f(\cdot)$ is the appropriate link function (e.g.,
$\text{logit}(x) = \log(x / (1 - x))$), and an estimate $\epsilon_n$ of the
coefficient $\epsilon$ of the "clever covariate" $H_n(A,W)$ is computed. The
form of the covariate $H_n(A,W)$ differs across target parameters; in this case
of the ATE, it is $H_n(A,W) = \frac{A}{g_n(A \mid W)} - \frac{1-A}{1-g_n(A,
W)}$, with $g_n(A,W) = \mathbb{P}(A=1 \mid W)$ being the estimated propensity
score, so the estimator depends both on the initial fit (by `sl3`) of the
outcome regression ($\overline{Q}_n$) and of the propensity score ($g_n$).

There are several robust augmentations that are used across the `tlverse`,
including the use of an additional layer of cross-validation to avoid
over-fitting bias (i.e., CV-TMLE) as well as approaches for more consistently
estimating several parameters simultaneously (e.g., the points on a survival
curve).

### Statistical Inference {#tmle-infer}

Since TMLE yields an **asymptotically linear** estimator, obtaining statistical
inference is very convenient. Each TML estimator has a corresponding
**(efficient) influence function** (often, "EIF", for short) that describes the
asymptotic distribution of the estimator. By using the estimated EIF, Wald-style
inference (asymptotically correct confidence intervals) can be constructed
simply by plugging into the form of the EIF our initial estimates
$\overline{Q}_n$ and $g_n$, then computing the sample standard error.

The following sections describe both a simple and more detailed way of
specifying and estimating a TMLE in the `tlverse`. In designing `tmle3`, we
sought to replicate as closely as possible the very general estimation framework
of TMLE, and so each theoretical object relevant to TMLE is encoded in a
corresponding software object/method. First, we will present the simple
application of `tmle3` to the WASH Benefits example, and then go on to describe
the underlying objects in greater detail.

## Easy-Bake Example: `tmle3` for ATE

We'll illustrate the most basic use of TMLE using the WASH Benefits data
introduced earlier and estimating an average treatment effect.

### Load the Data

We'll use the same WASH Benefits data as the earlier chapters:


```r
library(data.table)
library(dplyr)
library(tmle3)
library(sl3)
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)
```

### Define the variable roles

We'll use the common $W$ (covariates), $A$ (treatment/intervention), $Y$
(outcome) data structure. `tmle3` needs to know what variables in the dataset
correspond to each of these roles. We use a list of character vectors to tell
it. We call this a "Node List" as it corresponds to the nodes in a Directed
Acyclic Graph (DAG), a way of displaying causal relationships between variables.


```r
node_list <- list(
  W = c(
    "month", "aged", "sex", "momage", "momedu",
    "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof", "asset_wardrobe",
    "asset_table", "asset_chair", "asset_khat",
    "asset_chouki", "asset_tv", "asset_refrig",
    "asset_bike", "asset_moto", "asset_sewmach",
    "asset_mobile"
  ),
  A = "tr",
  Y = "whz"
)
```

### Handle Missingness

Currently, missingness in `tmle3` is handled in a fairly simple way:

* Missing covariates are median- (for continuous) or mode- (for discrete)
  imputed, and additional covariates indicating imputation are generated, just
  as described in [the `sl3` chapter](#sl3).
* Missing treatment variables are excluded -- such observations are dropped.
* Missing outcomes are efficiently handled by the automatic calculation (and
  incorporation into estimators) of _inverse probability of censoring weights_
  (IPCW); this is also known as IPCW-TMLE and may be thought of as a joint
  intervention to remove missingness and is analogous to the procedure used with
  classical inverse probability weighted estimators.

These steps are implemented in the `process_missing` function in `tmle3`:


```r
processed <- process_missing(washb_data, node_list)
washb_data <- processed$data
node_list <- processed$node_list
```

### Create a "Spec" Object

`tmle3` is general, and allows most components of the TMLE procedure to be
specified in a modular way. However, most users will not be interested in
manually specifying all of these components. Therefore, `tmle3` implements a
`tmle3_Spec` object that bundles a set of components into a _specification_
("Spec") that, with minimal additional detail, can be run to fit a TMLE.

We'll start with using one of the specs, and then work our way down into the
internals of `tmle3`.


```r
ate_spec <- tmle_ATE(
  treatment_level = "Nutrition + WSH",
  control_level = "Control"
)
```

### Define the learners

Currently, the only other thing a user must define are the `sl3` learners used
to estimate the relevant factors of the likelihood: Q and g.

This takes the form of a list of `sl3` learners, one for each likelihood factor
to be estimated with `sl3`:


```r
# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_rf <- make_learner(Lrnr_ranger)

# define metalearners appropriate to data types
ls_metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_rf),
  metalearner = ls_metalearner
)
sl_A <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_rf),
  metalearner = mn_metalearner
)
learner_list <- list(A = sl_A, Y = sl_Y)
```

Here, we use a Super Learner as defined in the previous chapter. In the future,
we plan to include reasonable defaults learners.

### Fit the TMLE

We now have everything we need to fit the tmle using `tmle3`:


```r
tmle_fit <- tmle3(ate_spec, washb_data, node_list, learner_list)
print(tmle_fit)
A tmle3_Fit that took 1 step(s)
   type                                    param  init_est tmle_est      se
1:  ATE ATE[Y_{A=Nutrition + WSH}-Y_{A=Control}] -0.005233 0.007111 0.05025
      lower  upper psi_transformed lower_transformed upper_transformed
1: -0.09139 0.1056        0.007111          -0.09139            0.1056
```

### Evaluate the Estimates

We can see the summary results by printing the fit object. Alternatively, we
can extra results from the summary by indexing into it:

```r
estimates <- tmle_fit$summary$psi_transformed
print(estimates)
[1] 0.007111
```

## `tmle3` Components

Now that we've successfully used a spec to obtain a TML estimate, let's look
under the hood at the components. The spec has a number of functions that
generate the objects necessary to define and fit a TMLE.

### `tmle3_task`

First is, a `tmle3_Task`, analogous to an `sl3_Task`, containing the data we're
fitting the TMLE to, as well as an NPSEM generated from the `node_list`
defined above, describing the variables and their relationships.


```r
tmle_task <- ate_spec$make_tmle_task(washb_data, node_list)
```


```r
tmle_task$npsem
$W
tmle3_Node: W
	Variables: month, aged, sex, momedu, hfiacat, Nlt18, Ncomp, watmin, elec, floor, walls, roof, asset_wardrobe, asset_table, asset_chair, asset_khat, asset_chouki, asset_tv, asset_refrig, asset_bike, asset_moto, asset_sewmach, asset_mobile, momage, momheight, delta_momage, delta_momheight
	Parents: 

$A
tmle3_Node: A
	Variables: tr
	Parents: W

$Y
tmle3_Node: Y
	Variables: whz
	Parents: A, W
```

### Initial Likelihood

Next, is an object representing the likelihood, factorized according to the
NPSEM described above:


```r
initial_likelihood <- ate_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)
print(initial_likelihood)
W: Lf_emp
A: LF_fit
Y: LF_fit
```

These components of the likelihood indicate how the factors were estimated: the
marginal distribution of $W$ was estimated using NP-MLE, and the conditional
distributions of $A$ and $Y$ were estimated using `sl3` fits (as defined with
the `learner_list`) above.

We can use this in tandem with the `tmle_task` object to obtain likelihood
estimates for each observation:

```r
initial_likelihood$get_likelihoods(tmle_task)
             W      A       Y
   1: 0.000213 0.3302 -0.3550
   2: 0.000213 0.3398 -0.9297
   3: 0.000213 0.3287 -0.8058
   4: 0.000213 0.3247 -0.9373
   5: 0.000213 0.3238 -0.5755
  ---                        
4691: 0.000213 0.2131 -0.5868
4692: 0.000213 0.2130 -0.2243
4693: 0.000213 0.2073 -0.7393
4694: 0.000213 0.2580 -0.9151
4695: 0.000213 0.1821 -1.0360
```

<!-- TODO: make helper to get learners out of fit objects -->

### Targeted Likelihood (updater)

We also need to define a "Targeted Likelihood" object. This is a special type
of likelihood that is able to be updated using an `tmle3_Update` object. This
object defines the update strategy (e.g., submodel, loss function, CV-TMLE or
not).


```r
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
```

When constructing the targeted likelihood, you can specify different update
options. See the documentation for `tmle3_Update` for details of the different
options. For example, you can disable CV-TMLE (the default in `tmle3`) as
follows:


```r
targeted_likelihood_no_cv <-
  Targeted_Likelihood$new(initial_likelihood,
    updater = list(cvtmle = FALSE)
  )
```

### Parameter Mapping

Finally, we need to define the parameters of interest. Here, the spec defines a
single parameter, the ATE. In the next section, we'll see how to add additional
parameters.


```r
tmle_params <- ate_spec$make_params(tmle_task, targeted_likelihood)
print(tmle_params)
[[1]]
Param_ATE: ATE[Y_{A=Nutrition + WSH}-Y_{A=Control}]
```

### Putting it all together

Having used the spec to manually generate all these components, we can now
manually fit a `tmle3`:


```r
tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)
print(tmle_fit_manual)
A tmle3_Fit that took 1 step(s)
   type                                    param  init_est tmle_est      se
1:  ATE ATE[Y_{A=Nutrition + WSH}-Y_{A=Control}] -0.004549  0.01096 0.05046
      lower  upper psi_transformed lower_transformed upper_transformed
1: -0.08794 0.1099         0.01096          -0.08794            0.1099
```

The result is equivalent to fitting using the `tmle3` function as above.

## Fitting `tmle3` with multiple parameters

Above, we fit a `tmle3` with just one parameter. `tmle3` also supports fitting
multiple parameters simultaneously. To illustrate this, we'll use the
`tmle_TSM_all` spec:


```r
tsm_spec <- tmle_TSM_all()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
all_tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
print(all_tsm_params)
[[1]]
Param_TSM: E[Y_{A=Control}]

[[2]]
Param_TSM: E[Y_{A=Handwashing}]

[[3]]
Param_TSM: E[Y_{A=Nutrition}]

[[4]]
Param_TSM: E[Y_{A=Nutrition + WSH}]

[[5]]
Param_TSM: E[Y_{A=Sanitation}]

[[6]]
Param_TSM: E[Y_{A=WSH}]

[[7]]
Param_TSM: E[Y_{A=Water}]
```

This spec generates a Treatment Specific Mean (TSM) for each level of the
exposure variable. Note that we must first generate a new targeted likelihood,
as the old one was targeted to the ATE. However, we can recycle the initial
likelihood we fit above, saving us a super learner step.

### Delta Method

We can also define parameters based on Delta Method Transformations of other
parameters. For instance, we can estimate a ATE using the delta method and two
of the above TSM parameters:


```r
ate_param <- define_param(
  Param_delta, targeted_likelihood,
  delta_param_ATE,
  list(all_tsm_params[[1]], all_tsm_params[[4]])
)
print(ate_param)
Param_delta: E[Y_{A=Nutrition + WSH}] - E[Y_{A=Control}]
```

This can similarly be used to estimate other derived parameters like Relative
Risks, and Population Attributable Risks

### Fit

We can now fit a TMLE simultaneously for all TSM parameters, as well as the
above defined ATE parameter


```r
all_params <- c(all_tsm_params, ate_param)

tmle_fit_multiparam <- fit_tmle3(
  tmle_task, targeted_likelihood, all_params,
  targeted_likelihood$updater
)

print(tmle_fit_multiparam)
A tmle3_Fit that took 1 step(s)
   type                                       param  init_est tmle_est      se
1:  TSM                            E[Y_{A=Control}] -0.592825 -0.62093 0.02978
2:  TSM                        E[Y_{A=Handwashing}] -0.615693 -0.65811 0.04155
3:  TSM                          E[Y_{A=Nutrition}] -0.608009 -0.60779 0.04187
4:  TSM                    E[Y_{A=Nutrition + WSH}] -0.597373 -0.60990 0.04095
5:  TSM                         E[Y_{A=Sanitation}] -0.582596 -0.57852 0.04220
6:  TSM                                E[Y_{A=WSH}] -0.517360 -0.44815 0.04515
7:  TSM                              E[Y_{A=Water}] -0.562570 -0.53743 0.03909
8:  ATE E[Y_{A=Nutrition + WSH}] - E[Y_{A=Control}] -0.004549  0.01103 0.05045
      lower   upper psi_transformed lower_transformed upper_transformed
1: -0.67929 -0.5626        -0.62093          -0.67929           -0.5626
2: -0.73955 -0.5767        -0.65811          -0.73955           -0.5767
3: -0.68986 -0.5257        -0.60779          -0.68986           -0.5257
4: -0.69015 -0.5296        -0.60990          -0.69015           -0.5296
5: -0.66123 -0.4958        -0.57852          -0.66123           -0.4958
6: -0.53664 -0.3597        -0.44815          -0.53664           -0.3597
7: -0.61405 -0.4608        -0.53743          -0.61405           -0.4608
8: -0.08785  0.1099         0.01103          -0.08785            0.1099
```

## Exercises

### Estimation of the ATE with `tmle3` {#tmle3-ex1}

Follow the steps below to estimate an average treatment effect using data from
the Collaborative Perinatal Project (CPP), available in the `sl3` package. To
simplify this example, we define a binary intervention variable, `parity01` --
an indicator of having one or more children before the current child and a
binary outcome, `haz01` -- an indicator of having an above average height for
age.


```r
# load the data set
data(cpp)
cpp <- cpp %>%
  as_tibble() %>%
  dplyr::filter(!is.na(haz)) %>%
  mutate(
    parity01 = as.numeric(parity > 0),
    haz01 = as.numeric(haz > 0)
  )
```
<!--
We're interested in using this simplified data to estimate an Average Treatment
Effect (ATE):
$\Psi(P_0)=\mathbb{E}_0(\mathbb{E}_0[Y \mid A=1,W]-\mathbb{E}_0[Y \mid A=0,W])$

The purely statistical (non-causal) parameter can be interpreted as the average
of the difference in means across the strata for $W$, and only requires the
positivity assumption, that the conditional treatment assignment probabilities
are positive for each possible $w$: $\mathbb{P}_0(A=1 \mid W=w) > 0$ and
$\mathbb{P}_0(A=0 \mid W=w) > 0$ for each possible $w$.

To interpret this parameter as causal, specifically the causal risk difference
$E_0Y_1-E_0Y_0$, then we would also need to make the randomization assumption
stating that $A$ is independent of the counterfactuals $(Y_0,Y_1)$ within
strata of $W$. This assumption might have been included in the original SCM
$\mathcal{M}^F$, but, if one knows there are unmeasured confounders, then the
model $\mathcal{M}^{F\*}$ would be more restrictive by enforcing this "known
to be wrong" randomization assumption. Still, this assumption does not change
the statistical model $\mathcal{M}$, and as a consequence, it does not affect
the estimation problem either. Thus, the theorem's that establish desirable
properties of the TMLE, still hold even when this non-testable randomization
assumption is violated.

We proceed with implementing a targeted minimum loss-based estimator (TMLE),
an efficient substitution estimator which is not only asymptotically
consistent, asymptotically normally distributed, and asymptotically efficient,
but also tailored to have robust finite sample performance.
-->

1. Define the variable roles $(W,A,Y)$ by creating a list of these nodes.
   Include the following baseline covariates in $W$: `apgar1`, `apgar5`,
   `gagebrth`, `mage`, `meducyrs`, `sexn`. Both $A$ and $Y$ are specified
   above. The missingness in the data (specifically, the missingness in the
   columns that are specified in the node list) will need to be taking care of.
   The `process_missing` function can be used to accomplish this, like the
   `washb_data` example above.
2. Define a `tmle3_Spec` object for the ATE, `tmle_ATE()`.
3. Using the same base learning libraries defined above, specify `sl3` base
   learners for estimation of $\overline{Q}_0 = \mathbb{E}_0(Y \mid A, W)$ and
   $g_0 = \mathbb{P}(A = 1 \mid W)$.
4. Define the metalearner like below.


```r
metalearner <- make_learner(
  Lrnr_solnp,
  loss_function = loss_loglik_binomial,
  learner_function = metalearner_logistic_binomial
)
```

5. Define one super learner for estimating $\overline{Q}_0$ and another for
   estimating $g_0$. Use the metalearner above for both super learners.
6. Create a list of the two super learners defined in the step above and call
   this object `learner_list`. The list names should be `A` (defining the super
   learner for estimation of $g_0$) and `Y` (defining the super learner for
   estimation of $\overline{Q}_0$).
7. Fit the TMLE with the `tmle3` function by specifying (1) the `tmle3_Spec`,
   which we defined in Step 2; (2) the data; (3) the list of nodes, which we
   specified in Step 1; and (4) the list of super learners for estimation of
   $g_0$ and $\overline{Q}_0$, which we defined in Step 6. *Note*: Like before,
   you will need to explicitly make a copy of the data (to work around
   `data.table` optimizations), e.g., (`cpp2 <- data.table::copy(cpp)`), then
   use the `cpp2` data going forward.

## Summary

`tmle3` is a general purpose framework for generating TML estimates. The easiest
way to use it is to use a predefined spec, allowing you to just fill in the
blanks for the data, variable roles, and `sl3` learners. However, digging under
the hood allows users to specify a wide range of TMLEs. In the next sections,
we'll see how this framework can be used to estimate advanced parameters such as
optimal treatments and stochastic shift interventions.
