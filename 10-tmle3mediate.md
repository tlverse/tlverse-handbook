# Causal Mediation Analysis

_Nima Hejazi_

Based on the [`tmle3mediate` `R`
package](https://github.com/tlverse/tmle3mediate) by _Nima Hejazi, James
Duncan, and David McCoy_.

Updated: 2021-06-30

\begin{VT1}
\VH{Learning Objectives}



1. Appreciate how the presence of a mediating variable in a causal pathway can
   allow direct and indirect effects of the treatment on the outcome to be
   defined.
2. Describe the major differences between direct and indirect causal effects.
3. Differentiate the joint interventions required to define direct and indirect
   effects from the static, dynamic, and stochastic interventions that yield
   the _total_ causal effects previously described.
4. Describe the assumptions needed for identification of the natural direct and
   indirect effects, as well as the limitations of these effect definitions.
5. Estimate the natural direct and indirect effects of a binary treatment using
   the `tmle3mediate` `R` package.
6. Differentiate the population intervention direct and indirect effects of
   stochastic interventions from the natural direct and indirect effects,
   including differences in the assumptions required for their identification.
7. Estimate the population intervention direct effect of a binary treatment
   using the `tmle3mediate` `R` package.

\end{VT1}

## Introduction to Causal Mediation Analysis

A treatment often affects an outcome indirectly, through a particular pathway,
by its effect on _intermediate variables_ (mediators). Causal mediation analysis
concerns the construction and evaluation of these _indirect effects_ and their
complementary _direct effects_. Generally, the indirect effect (IE) of a
treatment on an outcome is the portion of the total effect that is found to work
_through_ mediator variables, while the direct effect often encompasses all
other components of the total effect, including both the effect of the treatment
on the outcome _and_ the effect through all paths not explicitly involving the
mediators). Identifying and quantifying the mechanisms underlying causal effects
is an increasingly desirable endeavor in public health, medicine, and the social
sciences, as such mechanistic knowledge improves understanding of both _why_ and
_how_ treatments may be effective.

While the study of mediation analysis may be traced back quite far, the field
only came into its modern form with the identification and careful study of the
natural direct and indirect effects [@robins1992identifiability;
@pearl2001direct]. The natural direct effect (NDE) and the natural indirect
effect (NIE) are based on a decomposition of the average treatment effect (ATE)
in the presence of mediators [@vanderweele2015explanation]; requisite
theory for the construction of efficient estimators of these quantities only
receiving attention relatively recently [@tchetgen2012semiparametric].

## Data Structure and Notation

Consider $n$ observed units $O_1, \ldots, O_n$, where each observed data random
variable takes the form $O = (W, A, Z, Y)$, for a vector of observed covariates
$W$, a binary or continuous treatment $A$, possibly multivariate mediators $Z$,
and a binary or continuous outcome $Y$. To avoid undue assumptions, we assume
only that $O \sim \mathcal{P} \in \M$ where $\M$ is the nonparametric
statistical model defined as all continuous densities on $O$ with respect to an
arbitrary dominating measure.

We formalize the definition of our counterfactual variables using the following
non-parametric structural equation model (NPSEM):
\begin{align}
  W &= f_W(U_W) \\
  A &= f_A(W, U_A) \\
  Z &= f_Z(W, A, U_Z) \\
  Y &= f_Y(W, A, Z, U_Y).
  (\#eq:npsem-mediate)
\end{align}
This set of equations
represents a mechanistic model generating the observed data $O$; furthermore,
the NPSEM encodes several fundamental assumptions. Firstly, there is an implicit
temporal ordering: $W$ occurs first, depending only on exogenous factors $U_W$;
$A$ happens next, based on both $W$ and exogenous factors $U_A$; then come the
mediators $Z$, which depend on $A$, $W$, and another set of exogenous factors
$U_Z$; and finally appears the outcome $Y$. We assume neither access to the set
of exogenous factors $\{U_W, U_A, U_Z, U_Y\}$ nor knowledge of the forms of the
deterministic generating functions $\{f_W, f_A, f_Z, f_Y\}$. The NPSEM
corresponds to the following DAG:


```r
library(dagitty)
library(ggdag)

# make DAG by specifying dependence structure
dag <- dagitty(
  "dag {
    W -> A
    W -> Z
    W -> Y
    A -> Z
    A -> Y
    Z -> Y
    W -> A -> Y
    W -> A -> Z -> Y
  }"
)
exposures(dag) <- c("A")
outcomes(dag) <- c("Y")
tidy_dag <- tidy_dagitty(dag)

# visualize DAG
ggdag(tidy_dag) +
  theme_dag()
```



\begin{center}\includegraphics[width=0.8\linewidth]{10-tmle3mediate_files/figure-latex/mediation-DAG-1} \end{center}

The likelihood of the data $O$ admits a factorization, wherein, for $p_0^O$,
the density of $O$ with respect to the product measure, the density evaluated
on a particular observation $o$ may be a written
\begin{align}
  p_0^O(x) = &q^O_{0,Y}(y \mid Z = z, A = a, W = w) \cdot \\
    &q^O_{0,Z}(z \mid A = a, W = w) \cdot \\
    &q^O_{0,A}(a \mid W = w) \cdot \\
    &q^O_{0,W}(w),\\
  (\#eq:likelihood-factorization-mediate)
\end{align}
where $q_{0, Y}$ is the conditional density of $Y$ given $(Z, A, W)$, $q_{0, Z}$
is the conditional density of $Z$ given $(A, W)$, $q_{0, A}$ is the conditional
density of $A$ given $W$, and $q_{0, W}$ is the density of $W$. For ease of
notation, we let $\bar{Q}_Y(Z, A, W) = \E[Y \mid Z, A, W]$, $Q_Z(A, W) =
P[Z \mid A, W]$, $g(A \mid W) = \P(A \mid W)$, and $q_W$ the marginal
distribution of $W$.

Finally, note that we have explicitly excluded potential confounders of the
mediator-outcome relationship affected by exposure (i.e., variables affected by
$A$ and affecting both $Z$ and $Y$). Mediation analysis in the presence of such
variables is exceptionally challenging [@avin2005identifiability]; thus, most
efforts to develop definitions of causal direct and indirect effects explicitly
disallowed such a form of confounding. While we will not discuss the matter
here, the interested reader may consult recent advances in the vast literature
on causal mediation analysis, among them @diaz2020nonparametric and
@hejazi2021nonparametric.

## Decomposing the Average Treatment Effect

The natural direct and indirect effects arise from a decomposition of the ATE:
\begin{equation*}
  \E[Y(1) - Y(0)] =
    \underbrace{\E[Y(1, Z(0)) - Y(0, Z(0))]}_{NDE} +
    \underbrace{\E[Y(1, Z(1)) - Y(1, Z(0))]}_{NIE}.
\end{equation*}
In particular, the natural indirect effect (NIE) measures the effect of the
treatment $A \in \{0, 1\}$ on the outcome $Y$ through the mediators $Z$, while
the natural direct effect (NDE) measures the effect of the treatment on the
outcome _through all other paths_. Identification of the natural direct and
indirect effects requires the following non-testable causal assumptions:

1. _Exchangeability_ (randomization): $Y(a, z) \indep (A, Z) \mid W$, further
   implying $\E\{Y(a, z) \mid A=a, W=w, Z=z\} = \E\{Y(a, z) \mid W=w\}$. This
   is a special case of the randomization assumption, extended to observational
   studies with mediators.
2. _Treatment positivity_: For any $a \in \mathcal{A}$ and $w \in
   \mathcal{W}$, $\xi < g(a \mid w) < 1 - \xi$, $\xi > 0$. This mirrors the
   assumption required for static intervention, discussed previously.
3. _Mediator positivity_: For any $z \in \mathcal{Z}$, $a \in \mathcal{A}$, and
   $w \in \mathcal{W}$, $\epsilon < Q(z \mid a, w)$, for $\epsilon > 0$. This
   only requires that the conditional density of the mediators be bounded away
   from zero for all $(z, a, w)$ in their joint support $\mathcal{Z} \times
   \mathcal{A} \times \mathcal{W}$.
4. _Cross-world counterfactual independence_: For all $a \neq a'$, both
   contained in $\mathcal{A}$ and $z \in \mathcal{Z}$, $Y(a', z)$ is independent
   of $Z(a)$, given $W$. That is, the counterfactual outcome under the treatment
   contrast $a' \in \mathcal{A}$ and the counterfactual mediator $Z(a) \in
   \mathcal{Z}$ (under a different contrast $a \in \mathcal{A}$) are
   independent. Note that the counterfactual outcome and mediator are defined
   under differing contrasts, hence the "cross-world" designation.

We note that many attempts have been made to weaken the last assumption, that of
cross-world counterfactual independence, including work by
@petersen2006estimation and @imai2010identification; however, importantly,
@robins2010alternative established that this assumption cannot be satisfied in
randomized experiments. Thus, the natural direct and indirect effects are not
identifiable in randomized experiments, calling into question their utility.
Despite this significant limitation, we will turn to considering estimation of
the statistical functionals corresponding to these effects in observational
studies.

## The Natural Direct Effect

The NDE is defined as
\begin{align*}
  \Psi_{NDE} &= \E[Y(1, Z(0)) - Y(0, Z(0))] \\
  &\overset{\text{rand.}}{=} \sum_w \sum_z
  [\underbrace{\E(Y \mid A = 1, z, w)}_{\bar{Q}_Y(A = 1, z, w)} -
  \underbrace{\E(Y \mid A = 0, z, w)}_{\bar{Q}_Y(A = 0, z, w)}] \times
  \underbrace{p(z \mid A = 0, w)}_{Q_Z(0, w))} \underbrace{p(w)}_{q_W},
\end{align*}
where the likelihood factors $p(z \mid A = 0, w)$ and $p(w)$ (among other
conditional densities) arise from a factorization of the joint likelihood:
\begin{equation*}
  p(w, a, z, y) = \underbrace{p(y \mid w, a, z)}_{Q_Y(A, W, Z)}
  \underbrace{p(z \mid w, a)}_{Q_Z(Z \mid A, W)}
  \underbrace{p(a \mid w)}_{g(A \mid W)}
  \underbrace{p(w)}_{Q_W}.
\end{equation*}

The process of estimating the NDE begins by constructing $\bar{Q}_{Y, n}$, an
estimate of the outcome mechanism $\bar{Q}_Y(Z, A, W) = \E \{Y \mid Z, A,
W\}$ (i.e., the conditional mean of $Y$, given $Z$, $A$, and $W$). With an
estimate of this conditional expectation in hand, predictions of the
counterfactual quantities $\bar{Q}_Y(Z, 1, W)$ (setting $A = 1$) and, likewise,
$\bar{Q}_Y(Z, 0, W)$ (setting $A = 0$) can readily be obtained. We denote the
difference of these counterfactual quantities $\bar{Q}_{\text{diff}}$, i.e.,
$\bar{Q}_{\text{diff}} = \bar{Q}_Y(Z, 1, W) - \bar{Q}_Y(Z, 0, W)$.
$\bar{Q}_{\text{diff}}$ represents the difference in the conditional mean of
$Y$ attributable to changes in $A$ while keeping $Z$ and $W$ at their _natural_
(that is, observed) values.

The estimation procedure treats $\bar{Q}_{\text{diff}}$ itself as a nuisance
parameter, regressing its estimate $\bar{Q}_{\text{diff}, n}$ on $W$, among
control observations only (i.e., those for whom $A = 0$ is observed); the goal
of this step is to remove part of the marginal impact of $Z$ on
$\bar{Q}_{\text{diff}}$, since $W$ is a parent of $Z$. Regressing this
difference on $W$ among the controls recovers the expected
$\bar{Q}_{\text{diff}}$, had all individuals been set to the control condition
$A = 0$. Any residual additive effect of $Z$ on $\bar{Q}_{\text{diff}}$ is
removed during the TML estimation step using the auxiliary (or "clever")
covariate, which accounts for the mediators $Z$. This auxiliary covariate takes
the form
\begin{equation*}
  C_Y(Q_Z, g)(O) = \Bigg\{\frac{\mathbb{I}(A = 1)}{g(1 \mid W)}
  \frac{Q_Z(Z \mid 0, W)}{Q_Z(Z \mid 1, W)} -
  \frac{\mathbb{I}(A = 0)}{g(0 \mid W)} \Bigg\}.
\end{equation*}
Breaking this down, $\frac{\mathbb{I}(A = 1)}{g(1 \mid W)}$ is the inverse
propensity score weight for $A = 1$ and, likewise, $\frac{\mathbb{I}(A = 0)}
{g(0 \mid W)}$ is the inverse propensity score weight for $A = 0$. The middle
term is the ratio of the conditional densities of the mediator under the control
($A = 0$) and treatment ($A = 1$) conditions.

This subtle appearance of a ratio of conditional densities is concerning --
tools to estimate such quantities are sparse in the statistics literature,
unfortunately, and the problem is still more complicated (and computationally
taxing) when $Z$ is high-dimensional. As only the ratio of these conditional
densities is required, a convenient re-parametrization may be achieved, that is,
\begin{equation*}
  \frac{p(A = 0 \mid Z, W) g(0 \mid W)}{p(A = 1 \mid Z, W) g(1 \mid W)}.
\end{equation*}
Going forward, we will denote this re-parameterized conditional probability
$e(A \mid Z, W) := p(A \mid Z, W)$. Similar re-parameterizations have been used
in @zheng2012targeted and @tchetgen2013inverse. This is particularly useful
since this reformulation reduces the problem to one concerning only the
estimation of conditional means, opening the door to the use of a wide range of
machine learning algorithms (e.g., most of those in
[`sl3`](https://github.com/tlverse/sl3)).

Underneath the hood, the counterfactual outcome difference
$\bar{Q}_{\text{diff}}$ and $e(A \mid Z, W)$, the conditional probability of $A$
given $Z$ and $W$, are used in constructing the auxiliary covariate for TML
estimation. These nuisance parameters play an important role in the
bias-correcting update step of the TMLE procedure.

## The Natural Indirect Effect

Derivation and estimation of the NIE is analogous to that of the NDE. The NIE
is the effect of $A$ on $Y$ _only through the mediator(s) $Z$_. This quantity
-- known as the natural indirect effect $\E(Y(Z(1), 1) - \E(Y(Z(0), 1)$ --
corresponds to the difference of the conditional expectation of $Y$ given $A =
1$ and $Z(1)$ (the values the mediator would take under $A = 1$) and the
conditional expectation of $Y$ given $A = 1$ and $Z(0)$ (the values the mediator
would take under $A = 0$).

As with the NDE, the re-parameterization trick can be used to estimate $\E(A
\mid Z, W)$, avoiding estimation of a possibly multivariate conditional density.
However, in this case, the mediated mean outcome difference, denoted
$\Psi_Z(Q)$, is instead estimated as follows
\begin{equation*}
  \Psi_{NIE}(Q) = \E (\Psi_{NIE, Z}(Q)(1, W) - \Psi_{NIE, Z}(Q)(0, W))
\end{equation*}

Here, $\bar{Q}_Y(Z, 1, W)$ (the predicted values for $Y$ given $Z$ and $W$ when
$A = 1$) is regressed on $W$, among the treated units (for whom $A = 1$ is
observed) to obtain the conditional mean $\Psi_{NIE, Z}(Q)(1, W)$. Performing
the same procedure, but now regressing $\bar{Q}_Y(Z, 1, W)$ on $W$ among the
control units (for whom $A = 0$ is observed) yields $\Psi_{NIE,Z}(Q)(0, W)$. The
difference of these two estimates is the NIE and can be thought of as the
additive marginal effect of treatment on the conditional expectation of $Y$
given $W$, $A = 1$, $Z$ through its effects on $Z$. So, in the case of the NIE,
our estimate $\psi_n$ is slightly different, but the same quantity $e(A \mid Z,
W)$ comes into play as the auxiliary covariate.

## The Population Intervention (In)Direct Effects

At times, the natural direct and indirect effects may prove too limiting, as
these effect definitions are based on _static interventions_ (i.e., setting
$A = 0$ or $A = 1$), which may be unrealistic for real-world interventions. In
such cases, one may turn instead to the population intervention direct effect
(PIDE) and the population intervention indirect effect (PIIE), which are based
on decomposing the effect of the population intervention effect (PIE) of
flexible stochastic interventions [@diaz2020causal].

A particular type of stochastic intervention well-suited to working with binary
treatments is the _incremental propensity score intervention_ (IPSI), first
proposed by @kennedy2017nonparametric. Such interventions do not
deterministically set the treatment level of an observed unit to a fixed
quantity (i.e., setting $A = 1$), but instead _alter the odds of receiving the
treatment_ by a fixed amount ($0 \leq \delta \leq \infty$) for each individual.
In particular, this intervention takes the form
\begin{equation*}
  g_{\delta}(1 \mid w) = \frac{\delta g(1 \mid w)}{\delta g(1 \mid w) + 1
  - g(1\mid w)},
\end{equation*}
where the scalar $0 < \delta < \infty$ specifies a _change in the odds of
receiving treatment_. As described by @diaz2020causal, this stochastic
intervention is a special case of exponential tilting, a framework that unifies
post-intervention treatment values that are draws from an altered distribution.

Unlike the natural direct and indirect effects, the conditions required for
identifiability of the population intervention direct and indirect effects are
more lax. Most importantly, these differences involve a (1) treatment positivity
assumption that only requires that the counterfactual treatment be in the
observed support of the treatment $\mathcal{A}$, and (2) no requirement of the
independence any cross-world counterfactuals.

## Decomposing the Population Intervention Effect

We may decompose the population intervention effect (PIE) in terms of the
_population intervention direct effect_ (PIDE) and the _population
intervention indirect effect_ (PIIE):
\begin{equation*}
  \mathbb{E}\{Y(A_\delta)\} - \mathbb{E}Y =
    \overbrace{\mathbb{E}\{Y(A_\delta, Z(A_\delta))
      - Y(A_\delta, Z)\}}^{\text{PIIE}} +
    \overbrace{\mathbb{E}\{Y(A_\delta, Z) - Y(A, Z)\}}^{\text{PIDE}}.
\end{equation*}

This decomposition of the PIE as the sum of the population intervention direct
and indirect effects has an interpretation analogous to the corresponding
standard decomposition of the average treatment effect. In the sequel, we will
compute each of the components of the direct and indirect effects above using
appropriate estimators as follows

* For $\mathbb{E}\{Y(A, Z)\}$, the sample mean $\frac{1}{n}\sum_{i=1}^n Y_i$ is
  consistent;
* for $\mathbb{E}\{Y(A_{\delta}, Z)\}$, a TML estimator for the effect of a
  joint intervention altering the treatment mechanism but not the mediation
  mechanism, based on the proposal in @diaz2020causal; and,
* for $\mathbb{E}\{Y(A_{\delta}, Z_{A_{\delta}})\}$, an efficient estimator for
  the effect of a joint intervention altering both the treatment and mediation
  mechanisms, as proposed in @kennedy2017nonparametric and implemented in the
  [`npcausal` R package](https://github.com/ehkennedy/npcausal).

## Estimating the Effect Decomposition Term

As described by @diaz2020causal, the statistical functional identifying the
decomposition term that appears in both the PIDE and PIIE
$\mathbb{E}\{Y(A_{\delta}, Z)\}$, which corresponds to altering the treatment
mechanism while keeping the mediation mechanism fixed, is
\begin{equation*}
  \theta_0(\delta) = \int m_0(a, z, w) g_{0,\delta}(a \mid w) p_0(z, w)
    d\nu(a, z, w),
\end{equation*}
for which a TML estimator is available. The corresponding _efficient influence
function_ (EIF) with respect to the nonparametric model $\mathcal{M}$ is
$D_{\eta,\delta}(o) = D^Y_{\eta,\delta}(o)
+ D^A_{\eta,\delta}(o) + D^{Z,W}_{\eta,\delta}(o) - \theta(\delta)$.

The TML estimator may be computed basd on the EIF estimating equation and may
incorporate cross-validation [@zheng2011cross; @chernozhukov2018double] to
circumvent possibly restrictive entropy conditions (e.g., Donsker class). The
resultant estimator is
\begin{equation*}
  \hat{\theta}(\delta) = \frac{1}{n} \sum_{i = 1}^n D_{\hat{\eta}_{j(i)},
  \delta}(O_i) = \frac{1}{n} \sum_{i = 1}^n \left\{ D^Y_{\hat{\eta}_{j(i)},
  \delta}(O_i) + D^A_{\hat{\eta}_{j(i)}, \delta}(O_i) +
  D^{Z,W}_{\hat{\eta}_{j(i)}, \delta}(O_i) \right\},
\end{equation*}
which is implemented in `tmle3mediate` (a one-step estimator is also avaialble,
in the [`medshift` R package](https://github.com/nhejazi/medshift)). We
demonstrate the use of `tmle3mediate` to obtain $\mathbb{E}\{Y(A_{\delta}, Z)\}$
via its TML estimator.

## Evaluating the Direct and Indirect Effects

We now turn to estimating the natural direct and indirect effects, as well as
the population intervention direct effect, using the WASH Benefits data,
introduced in earlier chapters. Let's first load the data:


```r
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mediate)

# download data
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)

# make intervention node binary and subsample
washb_data <- washb_data[sample(.N, 600), ]
washb_data[, tr := as.numeric(tr != "Control")]
```

We'll next define the baseline covariates $W$, treatment $A$, mediators $Z$,
and outcome $Y$ nodes of the NPSEM via a "Node List" object:


```r
node_list <- list(
  W = c(
    "momage", "momedu", "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof"
  ),
  A = "tr",
  Z = c("sex", "month", "aged"),
  Y = "whz"
)
```

Here the `node_list` encodes the parents of each node -- for example, $Z$ (the
mediators) have parents $A$ (the treatment) and $W$ (the baseline confounders),
and $Y$ (the outcome) has parents $Z$, $A$, and $W$. We'll also handle any
missingness in the data by invoking `process_missing`:


```r
processed <- process_missing(washb_data, node_list)
washb_data <- processed$data
node_list <- processed$node_list
```

We'll now construct an ensemble learner using a handful of popular machine
learning algorithms:


```r
# SL learners used for continuous data (the nuisance parameter Z)
enet_contin_learner <- Lrnr_glmnet$new(
  alpha = 0.5, family = "gaussian", nfolds = 3
)
lasso_contin_learner <- Lrnr_glmnet$new(
  alpha = 1, family = "gaussian", nfolds = 3
)
fglm_contin_learner <- Lrnr_glm_fast$new(family = gaussian())
mean_learner <- Lrnr_mean$new()
contin_learner_lib <- Stack$new(
  enet_contin_learner, lasso_contin_learner, fglm_contin_learner, mean_learner
)
sl_contin_learner <- Lrnr_sl$new(learners = contin_learner_lib)

# SL learners used for binary data (nuisance parameters G and E in this case)
enet_binary_learner <- Lrnr_glmnet$new(
  alpha = 0.5, family = "binomial", nfolds = 3
)
lasso_binary_learner <- Lrnr_glmnet$new(
  alpha = 1, family = "binomial", nfolds = 3
)
fglm_binary_learner <- Lrnr_glm_fast$new(family = binomial())
binary_learner_lib <- Stack$new(
  enet_binary_learner, lasso_binary_learner, fglm_binary_learner, mean_learner
)
sl_binary_learner <- Lrnr_sl$new(learners = binary_learner_lib)

# create list for treatment and outcome mechanism regressions
learner_list <- list(
  Y = sl_contin_learner,
  A = sl_binary_learner
)
```

## Estimating the Natural Indirect Effect

We demonstrate calculation of the NIE below, starting by instantiating a "Spec"
object that encodes exactly which learners to use for the nuisance parameters
$e(A \mid Z, W)$ and $\Psi_Z$. We then pass our Spec object to the `tmle3`
function, alongside the data, the node list (created above), and a learner list
indicating which machine learning algorithms to use for estimating the nuisance
parameters based on $A$ and $Y$.


```r
tmle_spec_NIE <- tmle_NIE(
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  psi_Z_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE),
  max_iter = 1
)
washb_NIE <- tmle3(
  tmle_spec_NIE, washb_data, node_list, learner_list
)
washb_NIE
A tmle3_Fit that took 1 step(s)
   type                  param  init_est  tmle_est       se     lower    upper
1:  NIE NIE[Y_{A=1} - Y_{A=0}] 0.0022912 0.0026608 0.044295 -0.084156 0.089478
   psi_transformed lower_transformed upper_transformed
1:       0.0026608         -0.084156          0.089478
```

Based on the output, we conclude that the indirect effect of the treatment
through the mediators (sex, month, aged) is
0.00266.

## Estimating the Natural Direct Effect

An analogous procedure applies for estimation of the NDE, only replacing the
Spec object for the NIE with `tmle_spec_NDE` to define learners for the NDE
nuisance parameters:


```r
tmle_spec_NDE <- tmle_NDE(
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  psi_Z_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE),
  max_iter = 1
)
washb_NDE <- tmle3(
  tmle_spec_NDE, washb_data, node_list, learner_list
)
washb_NDE
A tmle3_Fit that took 1 step(s)
   type                  param init_est tmle_est      se   lower   upper
1:  NDE NDE[Y_{A=1} - Y_{A=0}] 0.012983 0.012983 0.10285 -0.1886 0.21457
   psi_transformed lower_transformed upper_transformed
1:        0.012983           -0.1886           0.21457
```

From this, we can draw the conclusion that the direct effect of the treatment
(through all paths not involving the mediators (sex, month, aged)) is
0.01298. Note that, together, the estimates of
the natural direct and indirect effects approximately recover the _average
treatment effect_, that is, based on these estimates of the NDE and NIE, the
ATE is roughly
0.01564.

## Estimating the Population Intervention Direct Effect

As previously noted, the assumptions underlying the natural direct and indirect
effects may be challenging to justify; moreover, the effect definitions
themselves depend on the application of a static intervention to the treatment,
sharply limiting their flexibility. When considering binary treatments,
incremental propensity score shifts provide an alternative class of flexible,
stochastic interventions. We'll now consider estimating the PIDE with an IPSI
that modulates the odds of receiving treatment by $\delta = 3$.  Such an
intervention may be interpreted (hypothetically) as the effect of a design that
encourages study participants to opt in to receiving the treatment, thus
increasing their relative odds of receiving said treatment. To exemplify our
approach, we postulate a motivational intervention that _triples the odds_
(i.e., $\delta = 3$) of receiving the treatment for each individual:


```r
# set the IPSI multiplicative shift
delta_ipsi <- 3

# instantiate tmle3 spec for stochastic mediation
tmle_spec_pie_decomp <- tmle_medshift(
  delta = delta_ipsi,
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  phi_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE)
)

# compute the TML estimate
washb_pie_decomp <- tmle3(
  tmle_spec_pie_decomp, washb_data, node_list, learner_list
)
washb_pie_decomp

# get the PIDE
washb_pie_decomp$summary$tmle_est - mean(washb_data[, get(node_list$Y)])
```

Recall that, based on the decomposition outlined previously, the PIDE may be
denoted $\beta_{\text{PIDE}}(\delta) = \theta_0(\delta) - \mathbb{E}Y$. Thus,
an estimator of the PIDE, $\hat{\beta}_{\text{PIDE}}(\delta)$ may be expressed
as a composition of estimators of its constituent parameters:
\begin{equation*}
  \hat{\beta}_{\text{PIDE}}({\delta}) = \hat{\theta}(\delta) -
  \frac{1}{n} \sum_{i = 1}^n Y_i.
\end{equation*}

<!--
Based on the above, we may construct an estimator of the PIDE using the already
estimated decomposition term and the empirical (marginal) mean of the outcome.
Note that this is a straightforward application of the delta method and could
equivalently be performed using the functionality exposed in the [`tmle3`
package](https://github.com/tlverse/tmle3).
-->

<!--

```r
tmle_task <- tmle_spec_pie_decomp$make_tmle_task(
  weight_behavior_complete, node_list
)
initial_likelihood <- tmle_spec_pie_decomp$make_initial_likelihood(
  tmle_task, learner_list
)
```
-->

<!--
## Exercises

### Review of Key Concepts

1. Examine the WASH Benefits dataset and choose a different set of potential
   mediators of the effect of the treatment on weight-for-height Z-score. Using
   this newly chosen set of mediators (or single mediator), estimate the
   natural direct and indirect effects. Provide an interpretation of these
   estimates.

2. Additivity of the natural (in)direct effects and the ATE.

3. Incremental propensity score interventions.

### The Ideas in Action

1. TODO

2. TODO

3. TODO
-->

<!--
## Appendix

### Exercise solutions
-->
