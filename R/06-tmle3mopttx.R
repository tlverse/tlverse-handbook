## ---- fig.cap="Illustration of a Dynamic Treatment Regime in a Clinical Setting", echo=FALSE, eval=TRUE, out.width='60%'----
knitr::include_graphics(path = "img/image/DynamicA_Illustration.png")


## ----setup-mopttx, message=FALSE, warning=FALSE--------------------------
library(here)
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mopttx)
library(devtools)
set.seed(111)


## ----load sim_cat_data---------------------------------------------------
data("data_cat")


## ----data_nodes-mopttx---------------------------------------------------
# organize data and nodes for tmle3
data <- data_cat
node_list <- list(W = c("W1", "W2", "W3", "W4"), A = "A", Y = "Y")


## ----sl3_lrnrs-mopttx----------------------------------------------------
# Initialize some of the learners.
# Here we use xgboost with various parameters: glm, HAL and the mean
xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
lrn1 <- Lrnr_mean$new()
lrn2 <- Lrnr_glm_fast$new()
lrn3 <- Lrnr_hal9001$new()

# Define the Q learner, which is just a regular learner:
Q_learner <- Lrnr_sl$new(
  learners = list(
    xgboost_50, xgboost_100, xgboost_500,
    lrn1, lrn2
  ),
  metalearner = Lrnr_nnls$new()
)

# Define the g learner, which is a multinomial learner:
glib <- list(
  rf = make_learner(Lrnr_randomForest),
  xgb = make_learner(Lrnr_xgboost),
  glmnet = make_learner(Lrnr_glmnet),
  multinom_gf = make_learner(
    Lrnr_independent_binomial,
    make_learner(Lrnr_glm_fast)
  ),
  mean = make_learner(Lrnr_mean)
)

# Specify the appropriate loss of the multinomial learner:
mn_metalearner <- make_learner(Lrnr_solnp,
  loss_function = loss_loglik_multinomial,
  learner_function =
    metalearner_linear_multinomial
)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)

# Define the Blip learner, which is a multivariate learner:
learners <- list(xgboost_50, xgboost_100, xgboost_500, lrn1, lrn2)
b_learner <- create_mv_learners(learners = learners)


## ----cat_learners--------------------------------------------------------
# See which learners support multi-class classification:
sl3_list_learners(c("categorical"))


## ----make_lrnr_list-mopttx-----------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----spec_init-----------------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3", "W4"), type =
    "blip2", b_learner = learner_list$B,
  maximize = TRUE, complex = TRUE
)


## ----mopttx_fit_tmle_auto_blip2, eval=T----------------------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----load sim_bin_data---------------------------------------------------
data("data_bin")


## ----data_nodes2-mopttx--------------------------------------------------
# organize data and nodes for tmle3
data <- data_bin
node_list <- list(W = c("W1", "W2", "W3"), A = "A", Y = "Y")


## ----mopttx_sl3_lrnrs2---------------------------------------------------
# Define sl3 library and metalearners:
xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
lrn1 <- Lrnr_mean$new()
lrn2 <- Lrnr_glm_fast$new()
lrn3 <- Lrnr_hal9001$new()

Q_learner <- Lrnr_sl$new(
  learners = list(
    xgboost_50, xgboost_100, xgboost_500,
    lrn1, lrn2
  ),
  metalearner = Lrnr_nnls$new()
)

g_learner <- Lrnr_sl$new(
  learners = list(xgboost_100, lrn2),
  metalearner = Lrnr_nnls$new()
)

b_learner <- Lrnr_sl$new(
  learners = list(
    xgboost_50, xgboost_100, xgboost_500,
    lrn1, lrn2
  ),
  metalearner = Lrnr_nnls$new()
)


## ----mopttx_make_lrnr_list-----------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----mopttx_spec_init_complex--------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  b_learner = learner_list$B,
  maximize = TRUE,
  complex = TRUE
)


## ----mopttx_fit_tmle_auto_blip_revere_complex, eval=T--------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----spec_init_Qlearning2------------------------------------------------
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
Q_learning(tmle_spec_Q, initial_likelihood, tmle_task)


## ----mopttx_spec_init_noncomplex-----------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  b_learner = learner_list$B, maximize =
    TRUE, complex = FALSE
)


## ----mopttx_fit_tmle_auto_blip_revere_noncomplex, eval=T-----------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----mopttx_sl3_lrnrs3---------------------------------------------------
# Define sl3 library and metalearners:
qlib <- make_learner_stack("Lrnr_mean", "Lrnr_glm_fast")
glib <- make_learner_stack("Lrnr_mean", "Lrnr_glmnet", "Lrnr_xgboost")
blib <- make_learner_stack("Lrnr_mean", "Lrnr_glm_fast")

metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(Lrnr_solnp,
  loss_function = loss_loglik_multinomial,
  learner_function =
    metalearner_linear_multinomial
)

Q_learner <- make_learner(Lrnr_sl, qlib, metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
b_learner <- make_learner(Lrnr_sl, blib, metalearner)


## ----mopttx_make_lrnr_list3----------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----mopttx_spec_init_vim------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_vim(
  V = c("W1", "W2", "W3"), type = "blip1",
  b_learner = learner_list$B, contrast =
    "multiplicative", maximize = FALSE, method =
    "SL"
)


## ----mopttx_fit_tmle_auto_vim, eval=FALSE--------------------------------
## # fit the TML estimator
## vim_results <- tmle3_vim(tmle_spec, data, node_list, learner_list,
##   adjust_for_other_A = FALSE
## )
## vim_results


## ----load_washb_data_mopttx, message=FALSE, warning=FALSE, cache=FALSE----
washb_data <- fread(here("data", "washb_data.csv"), stringsAsFactors = TRUE)
washb_data <- washb_data[!is.na(momage), lapply(.SD, as.numeric)]
head(washb_data, 3)


## ----washb_data_npsem_mopttx, message=FALSE, warning=FALSE, cache=FALSE----
node_list <- list(
  W = names(washb_data)[!(names(washb_data) %in%
    c("whz", "tr"))],
  A = "tr",
  Y = "whz"
)


## ----sl3_lrnrs_washb_mopttx----------------------------------------------
xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
glm_fast <- Lrnr_glm_fast$new()
lrn_mean <- Lrnr_mean$new()

# Define the Q learner, which is just a regular learner:
Q_learner <- Lrnr_sl$new(
  learners = list(xgboost_100, xgboost_500, lrn_mean),
  metalearner = Lrnr_nnls$new()
)

# Define the g learner, which is a multinomial learner:
glib <- list(
  xgb = make_learner(Lrnr_xgboost),
  mean = make_learner(Lrnr_mean)
)

# Specify the appropriate loss of the multinomial learner:
mn_metalearner <- make_learner(Lrnr_solnp,
  loss_function =
    loss_loglik_multinomial, learner_function =
    metalearner_linear_multinomial
)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)

# Define the Blip learner, which is a multivariate learner:
learners <- list(xgboost_100, xgboost_500, lrn_mean)
b_learner <- create_mv_learners(learners = learners)

# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----spec_init_washb_mopttx----------------------------------------------
table(washb_data$momedu)
table(washb_data$floor)
table(washb_data$asset_refrig)

# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("momedu", "floor", "asset_refrig"),
  type = "blip2",
  b_learner = learner_list$B,
  maximize = TRUE, complex = TRUE
)


## ----fit_tmle_auto_washb_mopttx, eval=T----------------------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data = washb_data, node_list, learner_list)
fit

