## ----mediation-DAG------------------------------------------------------------
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


## ----tmle3mediate-load-data---------------------------------------------------
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


## ----tmle3mediate-node-list---------------------------------------------------
node_list <- list(
  W = c(
    "momage", "momedu", "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof"
  ),
  A = "tr",
  Z = c("sex", "month", "aged"),
  Y = "whz"
)


## ----tmle3mediate-process_missing---------------------------------------------
processed <- process_missing(washb_data, node_list)
washb_data <- processed$data
node_list <- processed$node_list


## ----tmle3mediate-sl-learners-------------------------------------------------
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


## ----tmle3mediate-NIE---------------------------------------------------------
tmle_spec_NIE <- tmle_NIE(
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  psi_Z_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE),
  max_iter = 1
)
washb_NIE <- tmle3(
  tmle_spec_NIE, washb_data, node_list, learner_list
)
washb_NIE


## ----tmle3mediate-NDE---------------------------------------------------------
tmle_spec_NDE <- tmle_NDE(
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  psi_Z_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE),
  max_iter = 1
)
washb_NDE <- tmle3(
  tmle_spec_NDE, washb_data, node_list, learner_list
)
washb_NDE


## ----tmle3mediate-pide-decomp, eval=FALSE-------------------------------------
## # set the IPSI multiplicative shift
## delta_ipsi <- 3
##
## # instantiate tmle3 spec for stochastic mediation
## tmle_spec_pie_decomp <- tmle_medshift(
##   delta = delta_ipsi,
##   e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
##   phi_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE)
## )
##
## # compute the TML estimate
## washb_pie_decomp <- tmle3(
##   tmle_spec_pie_decomp, washb_data, node_list, learner_list
## )
## washb_pie_decomp
##
## # get the PIDE
## washb_pie_decomp$summary$tmle_est - mean(washb_data[, get(node_list$Y)])


## ----pide_delta, message=FALSE, warning=FALSE, eval=FALSE---------------------
## tmle_task <- tmle_spec_pie_decomp$make_tmle_task(
##   weight_behavior_complete, node_list
## )
## initial_likelihood <- tmle_spec_pie_decomp$make_initial_likelihood(
##   tmle_task, learner_list
## )
