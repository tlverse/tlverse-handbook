## ---- fig.cap="Animation of how a counterfactual outcome changes as the natural treatment distribution is subjected to a simple stochastic intervention", echo=FALSE, eval=TRUE, out.width='60%'----
knitr::include_graphics(path = "img/gif/shift_animation.gif")


## ----setup-shift, message=FALSE, warning=FALSE--------------------------------
library(data.table)
library(haldensify)
library(sl3)
library(tmle3)
library(tmle3shift)


## ----sl3_lrnrs-Qfit-shift, message=FALSE, warning=FALSE-----------------------
# learners used for conditional mean of the outcome
mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new()
rf_lrnr <- Lrnr_ranger$new()
hal_lrnr <- Lrnr_hal9001$new()

# SL for the outcome regression
sl_reg_lrnr <- Lrnr_sl$new(
  learners = list(mean_lrnr, fglm_lrnr, rf_lrnr, hal_lrnr),
  metalearner = Lrnr_nnls$new()
)


## ----sl3_density_lrnrs_search-shift, message=FALSE, warning=FALSE-------------
sl3_list_learners("density")


## ----sl3_lrnrs-gfit-shift, message=FALSE, warning=FALSE-----------------------
# learners used for conditional densities (i.e., generalized propensity score)
haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = 5, grid_type = "equal_mass",
  lambda_seq = exp(seq(-1, -10, length = 200))
)
# semiparametric density estimator based on homoscedastic errors (HOSE)
hose_hal_lrnr <- make_learner(Lrnr_density_semiparametric,
  mean_learner = hal_lrnr
)
# semiparametric density estimator based on heteroscedastic errors (HESE)
hese_rf_glm_lrnr <- make_learner(Lrnr_density_semiparametric,
  mean_learner = rf_lrnr,
  var_learner = fglm_lrnr
)

# SL for the conditional treatment density
sl_dens_lrnr <- Lrnr_sl$new(
  learners = list(haldensify_lrnr, hose_hal_lrnr, hese_rf_glm_lrnr),
  metalearner = Lrnr_solnp_density$new()
)


## ----learner-list-shift, message=FALSE, warning=FALSE-------------------------
learner_list <- list(Y = sl_reg_lrnr, A = sl_dens_lrnr)


## ----sim_data, message=FALSE, warning=FALSE-----------------------------------
# simulate simple data for tmle-shift sketch
n_obs <- 500 # number of observations
tx_mult <- 2 # multiplier for the effect of W = 1 on the treatment

## baseline covariates -- simple, binary
W <- replicate(2, rbinom(n_obs, 1, 0.5))

## create treatment based on baseline W
A <- rnorm(n_obs, mean = tx_mult * W, sd = 1)

## create outcome as a linear function of A, W + white noise
Y <- rbinom(n_obs, 1, prob = plogis(A + W))

# organize data and nodes for tmle3
data <- data.table(W, A, Y)
setnames(data, c("W1", "W2", "A", "Y"))
node_list <- list(W = c("W1", "W2"), A = "A", Y = "Y")
head(data)


## ----spec_init-shift, message=FALSE, warning=FALSE----------------------------
# initialize a tmle specification
tmle_spec <- tmle_shift(
  shift_val = 0.5,
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)


## ----fit_tmle-shift, message=FALSE, warning=FALSE, cache=FALSE----------------
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)
tmle_fit


## ----vim_spec_init, message=FALSE, warning=FALSE------------------------------
# what's the grid of shifts we wish to consider?
delta_grid <- seq(-1, 1, 1)

# initialize a tmle specification
tmle_spec <- tmle_vimshift_delta(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)


## ----fit_tmle_wrapper_vimshift, message=FALSE, warning=FALSE, cache=FALSE-----
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)
tmle_fit


## ----msm_fit, message=FALSE, warning=FALSE------------------------------------
tmle_fit$summary[4:5, ]


## ----vim_targeted_msm_fit, message=FALSE, warning=FALSE, cache=FALSE----------
# initialize a tmle specification
tmle_msm_spec <- tmle_vimshift_msm(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)

# fit the TML estimator and examine the results
tmle_msm_fit <- tmle3(tmle_msm_spec, data, node_list, learner_list)
tmle_msm_fit


## ----load_washb_data_shift, message=FALSE, warning=FALSE, cache=FALSE---------
washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv",
                    stringsAsFactors = TRUE)
washb_data <- washb_data[!is.na(momage), lapply(.SD, as.numeric)]
head(washb_data, 3)


## ----washb_data_npsem_shift, message=FALSE, warning=FALSE, cache=FALSE--------
node_list <- list(
  W = names(washb_data)[!(names(washb_data) %in%
    c("whz", "momage"))],
  A = "momage", Y = "whz"
)


## ----vim_spec_init_washb_shift, message=FALSE, warning=FALSE------------------
# initialize a tmle specification for the variable importance parameter
washb_vim_spec <- tmle_vimshift_delta(
  shift_grid = c(-2, 2),
  max_shifted_ratio = 2
)


## ----sl3_lrnrs_gfit_washb_shift, message=FALSE, warning=FALSE-----------------
# we need to turn on cross-validation for the HOSE learner
cv_hose_hal_lrnr <- Lrnr_cv$new(
  learner = hose_hal_lrnr,
  full_fit = TRUE
)

# modify learner list, using existing SL for Q fit
learner_list <- list(Y = sl_reg_lrnr, A = cv_hose_hal_lrnr)


## ----fit_tmle_wrapper_washb_shift, message=FALSE, warning=FALSE, eval=FALSE----
## washb_tmle_fit <- tmle3(washb_vim_spec, washb_data, node_list, learner_list)
## washb_tmle_fit

