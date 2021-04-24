# Optimal Individualized Treatment Regimes

_Ivana Malenica_

Based on the [`tmle3mopttx` `R` package](https://github.com/tlverse/tmle3mopttx)
by _Ivana Malenica, Jeremy Coyle, and Mark van der Laan_.

Updated: 2021-04-24

## Learning Objectives

1. Differentiate dynamic and optimal dynamic treatment interventions from static
   interventions.
2. Explain the benefits, and challenges, associated with using optimal
   individualized treatment regimes in practice.
3. Contrast the impact of implementing an optimal individualized treatment
   regime in the population with the impact of implementing static and dynamic
   treatment regimes in the population.
4. Estimate causal effects under optimal individualized treatment regimes with
   the `tmle3mopttx` `R` package.
5. Assess the mean under optimal individualized treatment with resource
   constraints.
6. Implement optimal individualized treatment rules based on sub-optimal
   rules, or "simple" rules, and recognize the practical benefit of these rules.
7. Construct "realistic" optimal individualized treatment regimes that respect
   real data and subject-matter knowledge limitations on interventions by
   only considering interventions that are supported by the data.
8. Measure variable importance as defined in terms of the optimal individualized
   treatment interventions.

## Introduction to Optimal Individualized Interventions

Identifying which intervention will be effective for which patient based on
lifestyle, genetic and environmental factors is a common goal in precision
medicine. To put it in context, Abacavir and Tenofovir are commonly prescribed
as part of the antiretroviral therapy to Human Immunodeficiency Virus (HIV)
patients. However, not all individuals benefit from the two medications equally.
In particular, patients with renal dysfunction might further deteriorate if
prescribed Tenofovir, due to the high nephrotoxicity caused by the medication.
While Tenofovir is still highly effective treatment option for HIV patients, in
order to maximize the patient's well-being, it would be beneficial to prescribe
Tenofovir only to individuals with healthy kidney function. Along the same
lines, one might seek to improve retention in HIV care. In a randomized clinical
trial, several interventions show efficacy- including appointment reminders
through text messages, small cash incentives for on time clinic visits, and peer
health workers. Ideally, we want to improve effectiveness by assigning each
patient the intervention they are most likely to benefit from, as well as
improve efficiency by not allocating resources to individuals that do not need
them, or would not benefit from it.

<div class="figure" style="text-align: center">
<img src="img/image/DynamicA_Illustration.png" alt="Dynamic Treatment Regime in a Clinical Setting" width="80%" />
<p class="caption">(\#fig:unnamed-chunk-1)Dynamic Treatment Regime in a Clinical Setting</p>
</div>

One opts to administer the intervention to individuals who will profit from it instead,
instead of assigning treatment on a population level. But how do we know which
intervention works for which patient? This aim motivates a different type of
intervention, as opposed to the static exposures we described in previous chapters. In
particular, in this chapter we learn about dynamic or "individualized"
interventions that tailor the treatment decision based on the collected
covariates. Formally, dynamic treatments represent interventions that at each
treatment-decision stage are allowed to respond to the currently available
treatment and covariate history.

In the statistics community such a treatment strategy is termed an
__individualized treatment regime__ (ITR), and the (counterfactual) population
mean outcome under an ITR is the value of the ITR [@neyman1990; @robins1986;
@pearl2009]. Even more, suppose one wishes to maximize the population mean of an
outcome, where for each individual we have access to some set of measured
covariates. This means, for example, that we can learn for which individual
characteristics assigning treatment increases the probability of a beneficial
outcome. An ITR with the maximal value is referred to as an
optimal ITR or the __optimal individualized treatment__. Consequently, the value
of an optimal ITR is termed the optimal value, or the __mean under the optimal
individualized treatment__.

The problem of estimating the optimal individualized treatment has received much
attention in the statistics literature over the years, especially with the
advancement of precision medicine; see @murphy2003, @robins2004, @laber2012,
@kosorok2012, @moodie2013 and @robins2014 to name a few. However, much of the
early work depends on parametric assumptions. As such, even in a randomized
trial, the statistical inference for the optimal individualized treatment relies
on assumptions that are generally believed to be false, and can lead to biased
results.

In this chapter, we consider estimation of the mean outcome under the optimal
individualized treatment where the candidate rules are restricted to depend only
on user-supplied subset of the baseline covariates. The estimation problem is
addressed in a statistical model for the data distribution that is
nonparametric, and at most places restrictions on the probability of a patient
receiving treatment given covariates (as in a randomized trial). As such, we
don't need to make any assumptions about the relationship of the outcome with
the treatment and covariates, or the relationship between the treatment and
covariates. Further, we provide a Targeted Maximum Likelihood Estimator for the
mean under the optimal individualized treatment that allows us to generate valid
inference for our parameter, without having any parametric assumptions. For a
technical presentation of the algorithm, the interested reader is invited to
further consult @vanderLaanLuedtke15 and @luedtke2016super.

---

## Data Structure and Notation

Suppose we observe $n$ independent and identically distributed observations of
the form $O=(W,A,Y) \sim P_0$. We denote $A$ as categorical treatment, and $Y$
as the final outcome. In particular, we define $A \in \mathcal{A}$ where
$\mathcal{A} \equiv \{a_1, \cdots, a_{n_A} \}$ and $n_A = |\mathcal{A}|$, with
$n_A$ denoting the number of categories (possibly only two, for a binary setup).
Note that we treat $W$ as vector-valued, representing all of our collected
baseline covariates. Therefore, for a single random individual $i$, we have that
their observed data is $O_i$: with corresponding baseline covariates $W_i$,
treatment $A_i$, and final outcome $Y_i$. We say that $O \sim P_0$, or that all
data was drawn from some true probability distribution $P_0$. Let $\mathcal{M}$
denote a statistical model, with $P_0 \in \mathcal{M}$. We emphasize that we
make no assumptions about the distribution of $P_0$, hence $\mathcal{M}$ is a
fully nonparametric model. As previously mentioned, this means that we make no
assumptions on the relationship between variables, but might be able to say
something about the relationship of $A$ and $W$, as is the case of a randomized
trial.  As in previous chapters, we denote $P_n$ as the empirical distribution
which gives each observation weight $1/n$.

We use the nonparametric structural equation model (NPSEM) in order to define
the process that gives rise to the observed (endogenous) and not observed
(exogenous) variables, as described by @pearl2009causality. In particular, we
denote $U=(U_W,U_A,U_Y)$ as the exogenous random variables, and $O=(W,A,Y)$ as
endogenous variables we observe. The joint distribution of exogenous and
endogenous random variables in $\mathcal{M}^F$ (defined by the NPSEM) is
$P_{U,X}$. We can define the relationships between variables with the following
structural equations:
\begin{align}
  W &= f_W(U_W) \\ A &= f_A(W, U_A) \\ Y &= f_Y(A, W, U_Y),
  (\#eq:npsem-mopttx)
\end{align}
where the collection $f=(f_W,f_A,f_Y)$ denotes unspecified functions. Note that
in the case of a randomized trial, we can write the above NPSEM as
\begin{align}
  W &= f_W(U_W) \\ A &= U_A \\ Y &= f_Y(A, W, U_Y),
  (\#eq:npsem-rt-mopttx)
\end{align}
indicating no dependence of treatment on baseline covariates.

The likelihood of the data admits a factorization, implied by the time ordering
of $O$. We denote the density of $O$ as $p_0$, corresponding to the
distribution $P_0$ and dominating measure $\mu$.
\begin{equation}
  p_0(O) = p_{Y,0}(Y \mid A,W) p_{A,0}(A \mid W) p_{W,0}(W) =
    q_{Y,0}(Y \mid A,W) q_{A,0}(A \mid W) q_{W,0}(W),
  (\#eq:likelihood-factorization-mopttx)
\end{equation}
where $p_{Y,0}(Y|A,W)$ is the conditional density of $Y$ given $(A, W)$ with
respect to some dominating measure $\mu_Y$, $p_{A,0}$ is the conditional density
of $A$ given $W$ with respect to dominating measure $\mu_A$, and $p_{W,0}$ is
the density of $W$ with respect to dominating measure $\mu_W$. Consequently, we
define $P_{Y,0}(Y \mid A, W) = Q_{Y,0}(Y \mid A,W)$, $P_{A,0}(A \mid W) =
g_0(A \mid W)$ and $P_{W,0}(W)=Q_{W,0}(W)$ as the corresponding conditional
distribution of $Y$ given $(A,W)$, treatment mechanism $A$ given $W$, and
distribution of baseline covariates. For notational simplicity, we also define
$\bar{Q}_{Y,0}(A,W) \equiv \E_0[Y \mid A,W]$ as the conditional expectation of
$Y$ given $(A,W)$.

Lastly, we define $V$ as a subset of the baseline covariates the optimal
individualized rule depends on, where $V \in W$.  Note that $V$ could be all of
$W$, or an empty set, depending on the subject matter knowledge. In particular,
a researcher might want to consider known effect modifiers available at the time
of treatment decision as possible $V$ covariates. Defining $V$ allows us to
consider possibly sub-optimal rules that are easier to estimate, and thereby
allows for statistical inference for the counterfactual mean outcome under the
sub-optimal rule.

## Defining the Causal Effect of an Optimal Individualized Intervention

Consider dynamic treatment rules $d$ in the set of all possible rules
$\mathcal{D}$. Then, $d$ is a function that takes as input $V$ and outputs a
treatment decision, $V \rightarrow d(V) \in \{a_1, \cdots, a_{n_A} \} \times
\{1\}$. We will use dynamic treatment rules, and the corresponding treatment
decision, to describe an intervention on the treatment mechanism and the
corresponding outcome under a dynamic treatment rule.

As mentioned in the previous section, causal effects are defined in terms of
hypothetical interventions on the NPSEM \@ref(eq:npsem-mopttx). We can define
counterfactuals $Y_{d(V)}$ defined by a modified system in which the equation
for $A$ is replaced by the rule $d(V)$, dependent on covariates $V$. Our
modified system then takes the following form:
\begin{align}
  W &= f_W(U_W) \\ A &= d(V) \\ Y_{d(V)} &= f_Y(d(V), W, U_Y),
  (\#eq:npsem-causal-mopttx)
\end{align}
where the dynamic treatment regime may be viewed as an intervention in which $A$
is set equal to a value based on a hypothetical regime $d(V)$, possibly contrary
to the fact, and $Y_{d(V)}$ is the corresponding outcome under $d(V)$. We
denote the distribution of the counterfactual quantities as $P_{0,d(V)}$.

The goal of any causal analysis motivated by such dynamic interventions is to
estimate a parameter defined as the counterfactual mean of the outcome with
respect to the modified intervention distribution. That is, subject's outcome
if, possibly contrary to the fact, the subject received treatment that would
have been assigned by rule $d(V)$. We can consider different treatment rules,
all in the set $\mathcal{D}$:

1. The true rule, $d_0$, and the corresponding causal parameter
   $\E_{U,X}[Y_{d_0(V)}]$;

2. The estimated rule, $d_n$, and the corresponding causal parameter
   $\E_{U,X}[Y_{d_n(V)}]$.

In this chapter, we will focus on the estimated rule $d_n$, and the
corresponding data-adaptive parameter.

The optimal individualized rule is the rule with the maximal value:
$$d_{opt}(V) \equiv \text{argmax}_{d(V) \in \mathcal{D}}
\E_{P_{U,X}}[Y_{d(V)}]$$.

We note that, in case the problem at hand requires minimizing the mean of an
outcome, our optimal individualized rule will be the rule with the minimal value
instead. Our causal target parameter of interest is the expected outcome under
the estimated optimal individualized rule:

$$\Psi_{d_{n, \text{opt}}(V)}(P_{U,X}) \coloneqq \E_{P_{U,X}}[Y_{d_{n,
\text{opt}}(V)}].$$

### Identification and Statistical Estimand

The optimal individualized rule, as well as the value of a optimal
individualized rule, are causal parameters based on the unobserved
counterfactuals. In order for the causal quantities to be estimated from the
observed data, they need to be identified with statistical parameters. This step
of the roadmap requires me make few assumptions:

1. _Strong ignorability_: $A \indep  Y^{d_n(v)} \mid W$, for all $a \in
   \mathcal{A}$.
2. _Positivity (or overlap)_: $P_0(\min_{a \in \mathcal{A}} g_0(a \mid W) > 0)
   = 1$

Under the above causal assumptions, we can identify $P_{0,d}$ with observed data
using the G-computation formula:

$$P_{0,d_{n, \text{opt}}}(O) = Q_{Y,0}(Y \mid A=d_{n,\text{opt}}(V),W)
g_0(A=d_{n,\text{opt}}(V) \mid W)Q_{W,0}(W).$$
The value of an individualized rule can now be expressed as

$$\E_0[Y_{d_n(V)}] = \E_{0,W}[\bar{Q}_{Y,0}(A=d_n(V),W)],$$

which, under causal assumptions, can is interpreted as the mean outcome if
(possibly contrary to fact), treatment was assigned according to the rule.
Finally, the statistical counterpart to the causal parameter of interest is
defined as
$$\psi_0 = \E_{0,W}[\bar{Q}_{Y,0}(A=d_{n,\text{opt}}(V),W)].$$

Inference for the optimal value has been shown to be difficult at exceptional
laws, defined as probability distributions for which treatment is neither
beneficial nor harmful. Inference is similarly difficult in finite samples if
the treatment effect is very small in all strata, even though valid asymptotic
estimators exist in this setting. With that in mind, we address the estimation
problem under the assumption of non-exceptional laws in effect.

Many methods for learning the optimal rule from data have been developed
[@murphy2003; @robins2004; @laber2012; @kosorok2012; @moodie2013]. In this
chapter, we focus on the methods discussed in @luedtke2016super and
@vanderLaanLuedtke15. Note however, that `tmle3mopttx` also supports the widely
used Q-learning approach, where the optimal individualized rule is based on the
initial estimate of $\bar{Q}_{Y,0}(A,W)$ [@Sutton1998].

We follow the methodology outlined in @luedtke2016super and
@vanderLaanLuedtke15, where we learn the optimal ITR using Super Learner
[@vdl2007super], and estimate its value with cross-validated Targeted Minimum
Loss-based Estimation (CV-TMLE) [@cvtmle2010]. In great generality, we first
need to estimate the true individual treatment regime, $d_0(V)$, which
corresponds to dynamic treatment rule ($d(V)$) that takes a subset of covariates
$V \in W$ and assigns treatment to each individual based on their observed
covariates $v$. With the estimate of the true optimal ITR in hand, we can
estimate its corresponding value.

### Binary treatment

How do we estimate the optimal individualized treatment regime? In the case of a
binary treatment, a key quantity for optimal ITR is the blip function. One can
show that any optimal ITR assigns treatment to individuals falling in strata in
which the stratum specific average treatment effect, the blip function, is
positive and does not assign treatment to individuals for which this quantity is
negative. Therefore for a binary treatment, under causal assumptions, we define
the blip function as:
$$\bar{Q}_0(V) \equiv \E_0[Y_1-Y_0 \mid V] \equiv \E_0[\bar{Q}_{Y,0}(1,W) -
\bar{Q}_{Y,0}(0,W) \mid V],$$
or the average treatment effect within a stratum of $V$. The note that the
optimal individualized rule can now be derived as $d_{n,\text{opt}}(V) =
\mathbb{I}(\bar{Q}_{0}(V) > 0)$.

The package `tmle3mopttx` relies on using the Super Learner to estimate the blip
function, as it easily extends to more general categorical treatment. With that
in mind, the loss function utilized for learning the optimal individualized rule
corresponds to conditional mean type losses. It is however worth mentioning that
@luedtke2016super present three different approaches for learning the optimal
rule. Namely, they focus on:

1. Super Learning the Blip Function,

2. Super Learning the Weighted Classification Problem,

3. Joint Super Learner of the Blip and Weighted Classification Problem.

We refer the interested reader to @luedtke2016super for further reference on
advantages of each approach.

Relying on the Targeted Maximum Likelihood (TML) estimator and the Super Learner
estimate of the blip function, we follow the below steps in order to obtain
value of the ITR:

1. Estimate $\bar{Q}_{Y,0}(A,W)$ and $g_0(A \mid W)$ using `sl3`. We denote such
   estimates as $\bar{Q}_{Y,n}(A,W)$ and $g_n(A \mid W)$.
2. Apply the doubly robust Augmented-Inverse Probability Weighted (A-IPW)
   transform to our outcome, where we define:
   $$D_{\bar{Q}_Y,g,a}(O) \equiv \frac{\mathbb{I}(A=a)}{g(A \mid W)} (Y -
   \bar{Q}_Y(A,W)) + \bar{Q}_Y(A=a,W)$$

Note that under the randomization and positivity assumptions we have that
$\E[D_{\bar{Q}_Y,g,a}(O) \mid V] = \E[Y_a \mid V]$. We emphasize the double
robust nature of the A-IPW transform-consistency of $\E[Y_a \mid V]$ will depend
on correct estimation of either $\bar{Q}_{Y,0}(A,W)$ or $g_0(A \mid W)$. As
such, in a randomized trial, we are guaranteed a consistent estimate of $\E[Y_a
\mid V]$ even if we get $\bar{Q}_{Y,0}(A,W)$ wrong!

Using this transform, we can define the following contrast:
$D_{\bar{Q}_Y,g}(O) = D_{\bar{Q}_Y, g, a=1}(O) - D_{\bar{Q}_Y, g, a=0}(O)$

We estimate the blip function, $\bar{Q}_{0,a}(V)$, by regressing
$D_{\bar{Q}_Y,g}(O)$ on $V$ using the specified `sl3` library of learners and an
appropriate loss function.

3. Our estimated rule corresponds to $\text{argmax}_{a \in \mathcal{A}}
   \bar{Q}_{0,a}(V)$.
4. We obtain inference for the mean outcome under the estimated optimal rule
   using CV-TMLE.

### Categorical treatment

In line with the approach considered for binary treatment, we extend the blip
function to allow for categorical treatment. We denote such blip function
extensions as _pseudo-blips_, which are our new estimation targets in a
categorical setting. We define pseudo-blips as vector-valued entities where the
output for a given $V$ is a vector of length equal to the number of treatment
categories, $n_A$. As such, we define it as:
$$\bar{Q}_0^{pblip}(V) = \{\bar{Q}_{0,a}^{pblip}(V): a \in \mathcal{A} \}$$

We implement three different pseudo-blips in `tmle3mopttx`.

1. _Blip1_ corresponds to choosing a reference category of treatment, and
   defining the blip for all other categories relative to the specified
   reference. Hence we have that:
   $$\bar{Q}_{0,a}^{pblip-ref}(V) \equiv \E_0(Y_a-Y_0 \mid V)$$ where $Y_0$ is
   the specified reference category with $A=0$. Note that, for the case of
   binary treatment, this strategy reduces to the approach described for the
   binary setup.

2. _Blip2_ approach corresponds to defining the blip relative to the average of
   all categories. As such, we can define $\bar{Q}_{0,a}^{pblip-avg}(V)$ as:
   $$\bar{Q}_{0,a}^{pblip-avg}(V) \equiv \E_0(Y_a - \frac{1}{n_A} \sum_{a \in
     \mathcal{A}} Y_a \mid V)$$
   In the case where subject-matter knowledge regarding which reference category
   to use is not available, blip2 might be a viable option.

3. _Blip3_ reflects an extension of Blip2, where the average is now a weighted
   average:
   $$\bar{Q}_{0,a}^{pblip-wavg}(V) \equiv \E_0(Y_a - \frac{1}{n_A} \sum_{a \in
     \mathcal{A}} Y_{a} P(A=a \mid V) \mid V)$$

Just like in the binary case, pseudo-blips are estimated by regressing contrasts
composed using the A-IPW transform on $V$.

### Note on Inference and data-adaptive parameter

In a randomized trial, statistical inference relies on the second-order
difference between the estimator of the optimal individualized treatment and the
optimal individualized treatment itself to be asymptotically negligible. This is
a reasonable condition if we consider rules that depend on a small number of
covariates, or if we are willing to make smoothness assumptions. Alternatively,
we can consider TMLEs and statistical inference for data-adaptive target
parameters defined in terms of an estimate of the optimal individualized
treatment. In particular, instead of trying to estimate the mean under the true
optimal individualized treatment, we aim to estimate the mean under the
estimated optimal individualized treatment. As such, we develop cross-validated
TMLE approach that provides asymptotic inference under minimal conditions for
the mean under the estimate of the optimal individualized treatment. In
particular, considering the data adaptive parameter allows us to avoid
consistency and rate condition for the fitted optimal rule, as required for
asymptotic linearity of the TMLE of the mean under the actual, true optimal
rule. Practically, the estimated (data-adaptive) rule should be preferred, as
this possibly sub-optimal rule is the one implemented in the population.

### Why CV-TMLE?

As discussed in @vanderLaanLuedtke15, CV-TMLE is necessary as the
non-cross-validated TMLE is biased upward for the mean outcome under the rule,
and therefore overly optimistic. More generally however, using CV-TMLE allows us
more freedom in estimation and therefore greater data adaptivity, without
sacrificing inference.

## Interpreting the Causal Effect of an Optimal Individualized Intervention

In summary, the mean outcome under the optimal individualized treatment is a
counterfactual quantity of interest representing what the mean outcome would
have been if everybody, contrary to the fact, received treatment that optimized
their outcome. The optimal individualized treatment regime is a rule that
optimizes the mean outcome under the dynamic treatment, where the candidate
rules are restricted to only respond to a user-supplied subset of the baseline
and intermediate covariates. In essence, our target parameter answers the key
aim of precision medicine: allocating the available treatment by tailoring it to
the individual characteristics of the patient, with the goal of optimizing the
final outcome.

## Evaluating the Causal Effect of an OIT with Binary Treatment {#oit-eval-bin}

Finally, we demonstrate how to evaluate the mean outcome under the optimal
individualized treatment using `tmle3mopptx`. To start, let's load the packages
we'll use and set a seed:


```r
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mopttx)
```

### Simulated Data

First, we load the simulated data. We will start with the more general setup
where the treatment is a binary variable; later in the chapter we will consider
another data-generating distribution where $A$ is categorical. In this example,
our data generating distribution is of the following form:
\begin{align*}
  W &\sim \mathcal{N}(\bf{0},I_{3 \times 3})\\
  \P(A=1 \mid W) &= \frac{1}{1+\exp^{(-0.8*W_1)}}\\
  \P(Y=1 \mid A,W) &= 0.5\text{logit}^{-1}[-5I(A=1)(W_1-0.5)+5I(A=0)(W_1-0.5)] +
     0.5\text{logit}^{-1}(W_2W_3)
\end{align*}


```r
data("data_bin")
```

The above composes our observed data structure $O = (W, A, Y)$. Note that the
mean under the true optimal rule is $\psi=0.578$ for this data generating
distribution.

To formally express this fact using the `tlverse` grammar introduced by the
`tmle3` package, we create a single data object and specify the functional
relationships between the nodes in the _directed acyclic graph_ (DAG) via
_nonparametric structural equation models_ (NPSEMs), reflected in the node list
that we set up:


```r
# organize data and nodes for tmle3
data <- data_bin
node_list <- list(
  W = c("W1", "W2", "W3"),
  A = "A",
  Y = "Y"
)
```

We now have an observed data structure (`data`) and a specification of the role
that each variable in the data set plays as the nodes in a DAG.

### Constructing Optimal Stacked Regressions with `sl3`

To easily incorporate ensemble machine learning into the estimation procedure,
we rely on the facilities provided in the [`sl3` R
package](https://tlverse.org/sl3). Using the framework provided by the [`sl3`
package](https://tlverse.org/sl3), the nuisance parameters of the TML estimator
may be fit with ensemble learning, using the cross-validation framework of the
Super Learner algorithm of @vdl2007super.


```r
# Define sl3 library and metalearners:
lrn_xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_xgboost_300 <- Lrnr_xgboost$new(nrounds = 300)
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm_fast$new()

## Define the Q learner:
Q_learner <- Lrnr_sl$new(
  learners = list(
    lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_300, lrn_mean, lrn_glm
  ),
  metalearner = Lrnr_nnls$new()
)

## Define the g learner:
g_learner <- Lrnr_sl$new(
  learners = list(lrn_xgboost_100, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

## Define the B learner:
b_learner <- Lrnr_sl$new(
  learners = list(
    lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_300, lrn_mean, lrn_glm
  ),
  metalearner = Lrnr_nnls$new()
)
```

As seen above, we generate three different ensemble learners that must be fit,
corresponding to the learners for the outcome regression (Q), propensity score
(g), and the blip function (B). We make the above explicit with respect to
standard notation by bundling the ensemble learners into a list object below:


```r
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)
```

The `learner_list` object above specifies the role that each of the ensemble
learners we've generated is to play in computing initial estimators. Recall that
we need initial estimators of relevant parts of the likelihood in order to
building a TMLE for the parameter of interest. In particular, `learner_list`
makes explicit the fact that our `Y` is used in fitting the outcome regression,
while `A` is used in fitting the treatment mechanism regression, and finally `B`
is used in fitting the blip function.

### Targeted Estimation of the Mean under the Optimal Individualized Interventions Effects

To start, we will initialize a specification for the TMLE of our parameter of
interest simply by calling `tmle3_mopttx_blip_revere`. We specify the argument
`V = c("W1", "W2", "W3")` when initializing the `tmle3_Spec` object in order to
communicate that we're interested in learning a rule dependent on `V`
covariates. Note that we don't have to specify `V`- this will result in a rule
that is not based on any collected covariates. We also need to specify the type
of pseudo-blip we will use in this estimation problem, the list of learners used
to estimate the blip function, whether we want to maximize or minimize the final
outcome, and few other more advanced features including searching for a less
complex rule and realistic interventions.


```r
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  learners = learner_list,
  maximize = TRUE, complex = TRUE, realistic = FALSE
)
```

As seen above, the `tmle3_mopttx_blip_revere` specification object
(like all `tmle3_Spec` objects) does _not_ store the data for our
specific analysis of interest. Later,
we'll see that passing a data object directly to the `tmle3` wrapper function,
alongside the instantiated `tmle_spec`, will serve to construct a `tmle3_Task`
object internally.

We elaborate more on the initialization specifications. In initializing the
specification for the TMLE of our parameter of interest, we have specified the
set of covariates the rule depends on (`V`), the type of pseudo-blip to use
(`type`), and the learners used for estimating the relevant parts of the
likelihood and the blip function. In addition, we need to specify whether we
want to maximize the mean outcome under the rule (`maximize`), and whether we
want to estimate the rule under all the covariates $V$ provided by the user
(`complex`). If `FALSE`, `tmle3mopttx` will instead consider all the possible
rules under a smaller set of covariates including the static rules, and optimize
the mean outcome over all the subsets of $V$. As such, while the user might have
provided a full set of collected covariates as input for $V$, it is possible
that the true rule only depends on a subset of the set provided by the user. In
that case, our returned mean under the optimal individualized rule will be based
on the smaller subset. In addition, we provide an option to search for realistic
optimal individualized interventions via the `realistic` specification. If
`TRUE`, only treatments supported by the data will be considered, therefore
alleviating concerns regarding practical positivity issues. We explore all the
important extensions of `tmle3mopttx` in later sections.


```r
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit
A tmle3_Fit that took 1 step(s)
   type         param init_est tmle_est       se   lower   upper
1:  TSM E[Y_{A=NULL}]  0.42223  0.56606 0.027015 0.51311 0.61901
   psi_transformed lower_transformed upper_transformed
1:         0.56606           0.51311           0.61901
```

We can see that the estimate of $psi_0$ is $0.56$, and that the confidence
interval covers our true mean under the true optimal individualized treatment.

## Evaluating the Causal Effect of an optimal ITR with Categorical Treatment {#oit-eval-cat}

In this section, we consider how to evaluate the mean outcome under the optimal
individualized treatment when $A$ has more than two categories. While the
procedure is analogous to the previously described binary treatment, we now need
to pay attention to the type of blip we define in the estimation stage, as well
as how we construct our learners.

### Simulated Data

First, we load the simulated data. Here, our data generating distribution was
of the following form:
\begin{align*}
  W &\sim \mathcal{N}(\bf{0},I_{4 \times 4})\\
  \P(A=a \mid W) &= \frac{1}{1+\exp^{(-0.8*W_1)}}\\
  \P(Y=1 \mid A,W) = 0.5\text{logit}^{-1}[15I(A=1)(W_1-0.5) -
    3I(A=2)(2W_1+0.5) +
    3I(A=3)(3W_1-0.5)] +\text{logit}^{-1}(W_2W_1)
\end{align*}

We can just load the data available as part of the package as follows:


```r
data("data_cat_realistic")
```

The above composes our observed data structure $O = (W, A, Y)$. Note that the
mean under the true optimal rule is $\psi=0.658$, which is the quantity we aim
to estimate.


```r
# organize data and nodes for tmle3
data <- data_cat_realistic
node_list <- list(
  W = c("W1", "W2", "W3", "W4"),
  A = "A",
  Y = "Y"
)
```

We'll now create new ensemble learners using the
`sl3` learners initialized previously:


```r
## Define the Q learner, which is just a regular learner:
Q_learner <- Lrnr_sl$new(
  learners = list(
    lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_300, lrn_mean, lrn_glm
  ),
  metalearner = Lrnr_nnls$new()
)

# Define the g learner, which is a multinomial learner:
# specify the appropriate loss of the multinomial learner:
mn_metalearner <- make_learner(Lrnr_solnp,
  loss_function = loss_loglik_multinomial,
  learner_function =
    metalearner_linear_multinomial
)
g_learner <- make_learner(
  Lrnr_sl,
  list(lrn_xgboost_100, lrn_xgboost_300, lrn_mean),
  mn_metalearner
)

# Define the Blip learner, which is a multivariate learner:
learners <- list(
  lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_300, lrn_mean, lrn_glm
)
b_learner <- create_mv_learners(learners = learners)
```

As seen above, we generate three different ensemble learners that must be fit,
corresponding to the learners for the outcome regression, propensity score, and
the blip function. Note that we need to estimate $g_0(A \mid W)$ for a
categorical $A$ -- therefore, we use the multinomial Super Learner option
available within the `sl3` package with learners that can address multi-class
classification problems. In order to see which learners can be used to estimate
$g_0(A \mid W)$ in `sl3`, we run the following:


```r
# See which learners support multi-class classification:
sl3_list_learners(c("categorical"))
 [1] "Lrnr_bound"                "Lrnr_caret"               
 [3] "Lrnr_cv_selector"          "Lrnr_glmnet"              
 [5] "Lrnr_grf"                  "Lrnr_gru_keras"           
 [7] "Lrnr_h2o_glm"              "Lrnr_h2o_grid"            
 [9] "Lrnr_independent_binomial" "Lrnr_lightgbm"            
[11] "Lrnr_lstm_keras"           "Lrnr_mean"                
[13] "Lrnr_multivariate"         "Lrnr_nnet"                
[15] "Lrnr_optim"                "Lrnr_polspline"           
[17] "Lrnr_pooled_hazards"       "Lrnr_randomForest"        
[19] "Lrnr_ranger"               "Lrnr_rpart"               
[21] "Lrnr_screener_correlation" "Lrnr_solnp"               
[23] "Lrnr_svm"                  "Lrnr_xgboost"             
```

Note that since the corresponding blip will be vector valued, we will have a
column for each additional level of treatment. As such, we need to create
multivariate learners with the helper function `create_mv_learners` that takes a
list of initialized learners as input.

We make the above explicit with respect to standard notation by bundling the
ensemble learners into a list object below:


```r
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)
```

### Targeted Estimation of the Mean under the Optimal Individualized Interventions Effects


```r
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3", "W4"), type = "blip2",
  learners = learner_list, maximize = TRUE, complex = TRUE,
  realistic = FALSE
)
```


```r
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit
A tmle3_Fit that took 1 step(s)
   type         param init_est tmle_est       se   lower   upper
1:  TSM E[Y_{A=NULL}]  0.54435  0.61851 0.069083 0.48311 0.75391
   psi_transformed lower_transformed upper_transformed
1:         0.61851           0.48311           0.75391
```

We can see that the estimate of $psi_0$ is $0.60$, and that the confidence
interval covers our true mean under the true optimal individualized treatment.

## Extensions to Causal Effect of an OIT

In this section, we consider two extensions to the procedure described for
estimating the value of the OIT. First one considers a setting where the user
might be interested in a grid of possible sub-optimal rules, corresponding to
potentially limited knowledge of potential effect modifiers. The second
extension concerns implementation of a realistic optimal individual
interventions where certain regimes might be preferred, but due to practical or
global positivity restraints are not realistic to implement.

### Simpler Rules

In order to not only consider the most ambitious fully $V$-optimal rule, we
define $S$-optimal rules as the optimal rule that considers all possible subsets
of $V$ covariates, with card($S$) $\leq$ card($V$) and $\emptyset \in S$. This
allows us to consider sub-optimal rules that are easier to estimate and
potentially provide more realistic rules- as such, we allow for statistical
inference for the counterfactual mean outcome under the sub-optimal rule.
Within the `tmle3mopttx` paradigm, we just need to change the `complex`
parameter to `FALSE`:


```r
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = FALSE, realistic = FALSE
)
```


```r
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit
A tmle3_Fit that took 1 step(s)
   type                param init_est tmle_est       se   lower   upper
1:  TSM E[Y_{d(V=W3,W2,W1)}]  0.54841  0.55108 0.062027 0.42951 0.67265
   psi_transformed lower_transformed upper_transformed
1:         0.55108           0.42951           0.67265
```

Therefore even though the user specified all baseline covariates as the basis
for rule estimation, a simpler rule based on only $W_2$ and $W_1$ is sufficient
to maximize the mean under the optimal individualized treatment.

### Realistic Optimal Individual Regimes

In addition to considering less complex rules, `tmle3mopttx` also provides an
option to estimate the mean under the realistic, or implementable, optimal
individualized treatment. It is often the case that assigning particular regime
might have the ability to fully maximize (or minimize) the desired outcome, but
due to global or practical positivity constrains, such treatment can never be
implemented in real life (or is highly unlikely). As such, specifying
`realistic` to `TRUE`, we consider possibly suboptimal treatments that optimize
the outcome in question while being supported by the data.


```r
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = TRUE, realistic = TRUE
)
```


```r
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit
A tmle3_Fit that took 1 step(s)
   type         param init_est tmle_est       se   lower  upper psi_transformed
1:  TSM E[Y_{A=NULL}]  0.55355  0.65922 0.021419 0.61723 0.7012         0.65922
   lower_transformed upper_transformed
1:           0.61723            0.7012

# How many individuals got assigned each treatment?
table(tmle_spec$return_rule)

  2   3 
502 498 
```

### Q-learning

Alternatively, we could estimate the mean under the optimal individualized
treatment using Q-learning. The optimal rule can be learned through fitting the
likelihood, and consequently estimating the optimal rule under this fit of the
likelihood [@Sutton1998; @murphy2003].

Below we outline how to use `tmle3mopttx` package in order to estimate the mean
under the ITR using Q-learning. As demonstrated in the previous sections, we
first need to initialize a specification for the TMLE of our parameter of
interest. As opposed to the previous section however, we will now use
`tmle3_mopttx_Q` instead of `tmle3_mopttx_blip_revere` in order to indicate that
we want to use Q-learning instead of TMLE.


```r
# initialize a tmle specification
tmle_spec_Q <- tmle3_mopttx_Q(maximize = TRUE)

# Define data:
tmle_task <- tmle_spec_Q$make_tmle_task(data, node_list)

# Define likelihood:
initial_likelihood <- tmle_spec_Q$make_initial_likelihood(
  tmle_task,
  learner_list
)

# Estimate the parameter:
Q_learning(tmle_spec_Q, initial_likelihood, tmle_task)[1]
```

## Variable Importance Analysis with OIT

Suppose one wishes to assess the importance of each observed covariate, in
terms of maximizing (or minimizing) the population mean of an outcome under an
optimal individualized treatment regime. In particular, a covariate that
maximizes (or minimizes) the population mean outcome the most under an optimal
individualized treatment out of all other considered covariates under optimal
assignment might be considered _more important_ for the outcome. To put it in
context, perhaps optimal allocation of treatment 1, denoted $A_1$, results in a
larger mean outcome than optimal allocation of another treatment ($A_2$).
Therefore, we would label $A_1$ as having a higher variable importance with
regard to maximizing (minimizing) the mean outcome under the optimal
individualized treatment.

### Simulated Data

In order to run `tmle3mopttx` variable importance measure, we need to consider
covariates to be categorical variables. For illustration purpose, we bin
baseline covariates corresponding to the data-generating distribution
[described previously](#oit-eval):


```r
# bin baseline covariates to 3 categories:
data$W1 <- ifelse(data$W1 < quantile(data$W1)[2], 1,
  ifelse(data$W1 < quantile(data$W1)[3], 2, 3)
)

node_list <- list(
  W = c("W3", "W4", "W2"),
  A = c("W1", "A"),
  Y = "Y"
)
```

Note that our node list now includes $W_1$ as treatments as well! Don't worry,
we will still properly adjust for all baseline covariates.

### Variable Importance using Targeted Estimation of the value of the ITR

In the previous sections we have seen how to obtain a contrast between the mean
under the optimal individualized rule and the mean under the observed outcome
for a single covariate- we are now ready to run the variable importance analysis
for all of our specified covariates. In order to run the variable importance
analysis, we first need to initialize a specification for the TMLE of our
parameter of interest as we have done before. In addition, we need to specify
the data and the corresponding list of nodes, as well as the appropriate
learners for the outcome regression, propensity score, and the blip function.
Finally, we need to specify whether we should adjust for all the other
covariates we are assessing variable importance for. We will adjust for all $W$s
in our analysis, and if `adjust_for_other_A=TRUE`, also for all $A$ covariates
that are not treated as exposure in the variable importance loop.

To start, we will initialize a specification for the TMLE of our parameter of
interest (called a `tmle3_Spec` in the `tlverse` nomenclature) simply by calling
`tmle3_mopttx_vim`. First, we indicate the method used for learning the optimal
individualized treatment by specifying the `method` argument of
`tmle3_mopttx_vim`. If `method="Q"`, then we will be using Q-learning for rule
estimation, and we do not need to specify `V`, `type` and `learners` arguments
in the spec, since they are not important for Q-learning. However, if
`method="SL"`, which corresponds to learning the optimal individualized
treatment using the above outlined methodology, then we need to specify the type
of pseudo-blip we will use in this estimation problem, whether we want to
maximize or minimize the outcome, complex and realistic rules. Finally, for
`method="SL"` we also need to communicate that we're interested in learning a
rule dependent on `V` covariates by specifying the `V` argument. For both
`method="Q"` and `method="SL"`, we need to indicate whether we want to maximize
or minimize the mean under the optimal individualized rule. Finally, we also
need to specify whether the final comparison of the mean under the optimal
individualized rule and the mean under the observed outcome should be on the
multiplicative scale (risk ratio) or linear (similar to average treatment
effect).


```r
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_vim(
  V = c("W2"),
  type = "blip2",
  learners = learner_list,
  contrast = "multiplicative",
  maximize = FALSE,
  method = "SL",
  complex = TRUE,
  realistic = FALSE
)
```


```r
# fit the TML estimator
vim_results <- tmle3_vim(tmle_spec, data, node_list, learner_list,
  adjust_for_other_A = TRUE
)
print(vim_results)
```

The final result of `tmle3_vim` with the `tmle3mopttx` spec is an ordered list
of mean outcomes under the optimal individualized treatment for all categorical
covariates in our dataset.

---

## Exercises

### Real World Data and `tmle3mopttx`

Finally, we cement everything we learned so far with a real data application.

As in the previous sections, we will be using the WASH Benefits data,
corresponding to the effect of water quality, sanitation, hand washing, and
nutritional interventions on child development in rural Bangladesh.

The main aim of the cluster-randomized controlled trial was to assess the
impact of six intervention groups, including:

1. control;

2. hand-washing with soap;

3. improved nutrition through counseling and provision of lipid-based nutrient
   supplements;

4. combined water, sanitation, hand-washing, and nutrition;

5. improved sanitation;

6. combined water, sanitation, and hand-washing;

7. chlorinated drinking water.

We aim to estimate the optimal ITR and the corresponding value under the optimal
ITR for the main intervention in WASH Benefits data.

Our outcome of interest is the weight-for-height Z-score, whereas our primary
treatment is the six intervention groups aimed at improving living conditions.

Questions:

1. Define $V$ as mother's education (`momedu`), current living conditions (`floor`),
   and possession of material items including the refrigerator (`asset_refrig`).
   Why do you think we use these covariates as $V$? Do we want to minimize or
   maximize the outcome?  Which blip type should we use?

2. Load the WASH Benefits data, and define the appropriate nodes for treatment
   and outcome.  Use all the rest of the covariates as $W$ except for
   `momheight` for now. Construct an appropriate `sl3` library for $A$, $Y$ and
   $B$.

3. Based on the $V$ defined in the previous question, estimate the mean under
   the ITR for the main randomized intervention used in the WASH Benefits trial
   with weight-for-height Z-score as the outcome. What's the TMLE value of the
   optimal ITR?  How does it change from the initial estimate? Which
   intervention is the most dominant?  Why do you think that is?

4. Using the same formulation as in questions 1 and 2, estimate the realistic
   optimal ITR and the corresponding value of the realistic ITR. Did the results
   change? Which intervention is the most dominant under realistic rules? Why do
   you think that is?

5. Consider simpler rules for the WASH benefits data example. What set of rules
   are picked?

6. Change the treatment to a binary variable (`asset_sewmach`), and estimate the
   value under the ITR in this setting under a $60\%$ resource constraint.  What
   do the results indicate?

7. Change the treatment once again, now to mother's education (`momedu`), and
   estimate the value under the ITR in this setting. What do the results
   indicate? Can we intervene on such a variable?

### Review of Key Concepts

1. What is the difference between dynamic and optimal individualized regimes?

2. What's the intuition behind using different blip types? Why did we switch
   from `blip1` to `blip2` when considering categorical treatment? What are some
   of the advantages of each?

3. Look back at the results generated in the [section on categorical
   treatments](#oit-eval-cat), and compare then to the mean under the optimal
   individualized treatment in the [section on binary
   treatments](#oit-eval-bin). Why do you think the estimate is higher under the
   less complex rule? How does the set of covariates picked by `tmle3mopttx`
   compare to the baseline covariates the true rule depends on?

4. Compare the distribution of treatments assigned under the true optimal
   individualized treatment and realistic optimal individualized treatment.
   Referring back to the data-generating distribution, why do you think the
   distribution of allocated treatment changed?

5. Using the same simulation, perform a variable importance analysis using
   Q-learning. How do the results change and why?


### Advanced Topics

1. How can we extend the current approach to include exceptional laws?

2. How can we extend the current approach to continuous interventions?

<!--
## Appendix

### Exercise solutions
-->
