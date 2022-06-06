## ---- fig.cap="Dynamic Treatment Regime in a Clinical Setting", results="asis", echo=FALSE----
knitr::include_graphics(path = "img/png/DynamicA_Illustration.png")


## ----setup-mopttx-------------------------------------------------------------
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mopttx)
library(devtools)

set.seed(111)


## ----load sim_bin_data--------------------------------------------------------
data("data_bin")


## ----data_nodes2-mopttx-------------------------------------------------------
# organize data and nodes for tmle3
data <- data_bin
node_list <- list(
  W = c("W1", "W2", "W3"),
  A = "A",
  Y = "Y"
)


## ----mopttx_sl3_lrnrs2--------------------------------------------------------
# Define sl3 library and metalearners:
lrn_xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)

lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm_fast$new()
lrn_lasso <- Lrnr_glmnet$new()

## Define the Q learner:
Q_learner <- Lrnr_sl$new(
  learners = list(lrn_lasso, lrn_mean, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

## Define the g learner:
g_learner <- Lrnr_sl$new(
  learners = list(lrn_lasso, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

## Define the B learner:
b_learner <- Lrnr_sl$new(
  learners = list(lrn_lasso,lrn_mean, lrn_glm),
  metalearner = Lrnr_nnls$new()
)


## ----mopttx_make_lrnr_list----------------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----mopttx_spec_init_complex-------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  learners = learner_list,
  maximize = TRUE, complex = TRUE,
  realistic = FALSE, resource = 1
)


## ----mopttx_fit_tmle_auto_blip_revere_complex---------------------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----mopttx_spec_init_complex_resource----------------------------------------
# initialize a tmle specification
tmle_spec_resource <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  learners = learner_list,
  maximize = TRUE, complex = TRUE,
  realistic = FALSE, resource = 0.90
)


## ----mopttx_fit_tmle_auto_blip_revere_complex_resource, eval=T----------------
# fit the TML estimator
fit_resource <- tmle3(tmle_spec_resource, data, node_list, learner_list)
fit_resource


## ----mopttx_compare_resource--------------------------------------------------
# Number of individuals getting treatment (no resource constraint):
table(tmle_spec$return_rule)

# Number of individuals getting treatment (resource constraint):
table(tmle_spec_resource$return_rule)


## ----mopttx_spec_init_complex_V_empty-----------------------------------------
# initialize a tmle specification
tmle_spec_V_empty <- tmle3_mopttx_blip_revere(
  type = "blip1",
  learners = learner_list,
  maximize = TRUE, complex = TRUE,
  realistic = FALSE, resource = 0.90
)


## ----mopttx_fit_tmle_auto_blip_revere_complex_V_empty, eval=T-----------------
# fit the TML estimator
fit_V_empty <- tmle3(tmle_spec_V_empty, data, node_list, learner_list)
fit_V_empty


## ----load sim_cat_data--------------------------------------------------------
data("data_cat_realistic")


## ----data_nodes-mopttx--------------------------------------------------------
# organize data and nodes for tmle3
data <- data_cat_realistic
node_list <- list(
  W = c("W1", "W2", "W3", "W4"),
  A = "A",
  Y = "Y"
)


## ----data_cats-mopttx---------------------------------------------------------
# organize data and nodes for tmle3
table(data$A)


## ----sl3_lrnrs-mopttx---------------------------------------------------------
# Initialize few of the learners:
lrn_xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm_fast$new()

## Define the Q learner, which is just a regular learner:
Q_learner <- Lrnr_sl$new(
  learners = list(lrn_xgboost_100, lrn_mean, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

# Define the g learner, which is a multinomial learner:
# specify the appropriate loss of the multinomial learner:
mn_metalearner <- make_learner(Lrnr_solnp,
  eval_function = loss_loglik_multinomial,
  learner_function = metalearner_linear_multinomial
)
g_learner <- make_learner(Lrnr_sl, list(lrn_xgboost_100, lrn_xgboost_500, lrn_mean), mn_metalearner)

# Define the Blip learner, which is a multivariate learner:
learners <- list(lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_500, lrn_mean, lrn_glm)
b_learner <- create_mv_learners(learners = learners)


## ----cat_learners-------------------------------------------------------------
# See which learners support multi-class classification:
sl3_list_learners(c("categorical"))


## ----make_lrnr_list-mopttx----------------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----spec_init----------------------------------------------------------------
# initialize a tmle specification
tmle_spec_cat <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3", "W4"), type = "blip2",
  learners = learner_list, maximize = TRUE, complex = TRUE,
  realistic = FALSE
)


## ----fit_tmle_auto------------------------------------------------------------
# fit the TML estimator
fit_cat <- tmle3(tmle_spec_cat, data, node_list, learner_list)
fit_cat

# How many individuals got assigned each treatment?
table(tmle_spec_cat$return_rule)


## ----mopttx_spec_init_noncomplex----------------------------------------------
# initialize a tmle specification
tmle_spec_cat_simple <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = FALSE, realistic = FALSE
)


## ----mopttx_fit_tmle_auto_blip_revere_noncomplex------------------------------
# fit the TML estimator
fit_cat_simple <- tmle3(tmle_spec_cat_simple, data, node_list, learner_list)
fit_cat_simple


## ----mopttx_spec_init_realistic-----------------------------------------------
# initialize a tmle specification
tmle_spec_cat_realistic <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = TRUE, realistic = TRUE
)


## ----mopttx_fit_tmle_auto_blip_revere_realistic-------------------------------
# fit the TML estimator
fit_cat_realistic <- tmle3(tmle_spec_cat_realistic, data, node_list, learner_list)
fit_cat_realistic

# How many individuals got assigned each treatment?
table(tmle_spec_cat_realistic$return_rule)


## ----data_nodes-add-missigness-mopttx-----------------------------------------
data_missing <- data_cat_realistic

#Add some random missingless:
rr <- sample(nrow(data_missing), 100, replace = FALSE)
data_missing[rr,"Y"]<-NA

summary(data_missing$Y)


## ----sl3_lrnrs-add-mopttx-----------------------------------------------------
delta_learner <- Lrnr_sl$new(
  learners = list(lrn_mean, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, delta_Y=delta_learner)


## ----spec_init_missingness----------------------------------------------------
# initialize a tmle specification
tmle_spec_cat_miss <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3", "W4"), type = "blip2",
  learners = learner_list, maximize = TRUE, complex = TRUE,
  realistic = FALSE
)


## ----fit_tmle_auto2, eval=T---------------------------------------------------
# fit the TML estimator
fit_cat_miss <- tmle3(tmle_spec_cat_miss, data_missing, node_list, learner_list)
fit_cat_miss


## ----spec_init_Qlearning2, eval=FALSE-----------------------------------------
## # initialize a tmle specification
## tmle_spec_Q <- tmle3_mopttx_Q(maximize = TRUE)
## 
## # Define data:
## tmle_task <- tmle_spec_Q$make_tmle_task(data, node_list)
## 
## # Define likelihood:
## initial_likelihood <- tmle_spec_Q$make_initial_likelihood(
##   tmle_task,
##   learner_list
## )
## 
## # Estimate the parameter:
## Q_learning(tmle_spec_Q, initial_likelihood, tmle_task)[1]


## ----data_vim-nodes-mopttx----------------------------------------------------
# bin baseline covariates to 3 categories:
data$W1<-ifelse(data$W1<quantile(data$W1)[2],1,ifelse(data$W1<quantile(data$W1)[3],2,3))

node_list <- list(
  W = c("W3", "W4", "W2"),
  A = c("W1", "A"),
  Y = "Y"
)


## ----mopttx_spec_init_vim-----------------------------------------------------
# initialize a tmle specification
tmle_spec_vim <- tmle3_mopttx_vim(
  V=c("W2"),
  type = "blip2",
  learners = learner_list,
  maximize = FALSE,
  method = "SL",
  complex = TRUE,
  realistic = FALSE
)


## ----mopttx_fit_tmle_auto_vim, eval=TRUE--------------------------------------
# fit the TML estimator
vim_results <- tmle3_vim(tmle_spec_vim, data, node_list, learner_list,
  adjust_for_other_A = TRUE
)

print(vim_results)

