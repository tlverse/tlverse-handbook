## ----setup, message=FALSE, warning=FALSE---------------------------------
library(here)
library(tidyverse)
library(data.table)
library(sl3)
library(SuperLearner)
library(origami)
library(knitr)
set.seed(7194)

# load data set and take a peek
washb_data <- fread(here("data", "washb_data.csv"), stringsAsFactors = TRUE)
head(washb_data, 3) %>%
  kable(format = "markdown", digits = 3)


## ----task----------------------------------------------------------------
# specify the outcome and covariates
outcome <- "whz"
covars <- colnames(washb_data)[-which(names(washb_data) == outcome)]

# create the sl3 task
washb_task <- make_sl3_Task(
  data = washb_data,
  covariates = covars,
  outcome = outcome
)

# examine it
washb_task


## ----list-properties-----------------------------------------------------
sl3_list_properties()


## ----list-learners-------------------------------------------------------
sl3_list_learners("continuous")


## ----baselearners--------------------------------------------------------
# choose base learners
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_ranger <- make_learner(Lrnr_ranger)
lrnr_glmnet <- make_learner(Lrnr_glmnet)


## ----stack---------------------------------------------------------------
stack <- make_learner(
  Stack,
  lrnr_glm, lrnr_mean, lrnr_ranger, lrnr_glmnet
)


## ----metalearner---------------------------------------------------------
metalearner <- make_learner(Lrnr_nnls)


## ----make-sl-------------------------------------------------------------
sl <- make_learner(Lrnr_sl,
  learners = stack,
  metalearner = metalearner
)


## ----slbasic-------------------------------------------------------------
sl_fit <- sl$train(washb_task)


## ----slbasic-summary-----------------------------------------------------
sl_preds <- sl_fit$predict()
head(sl_preds)
sl_fit$print() %>%
  kable(format = "markdown", digits = 3)


## ----extra-lrnr----------------------------------------------------------
lrnr_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
lrnr_ranger1k <- make_learner(Lrnr_ranger, num.trees = 1000)
lrnr_gam <- Lrnr_pkg_SuperLearner$new("SL.gam")
lrnr_bayesglm <- Lrnr_pkg_SuperLearner$new("SL.bayesglm")


## ----new-stack-----------------------------------------------------------
new_stack <- make_learner(
  Stack,
  lrnr_glm, lrnr_mean, lrnr_ranger, lrnr_glmnet, lrnr_ranger1k, lrnr_ranger100,
  lrnr_gam, lrnr_bayesglm
)


## ----screeners-----------------------------------------------------------
screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")


## ----screeners-pipe------------------------------------------------------
cor_pipeline <- make_learner(Pipeline, screen_cor, new_stack)


## ----screeners-stack-----------------------------------------------------
fancy_stack <- make_learner(Stack, cor_pipeline, new_stack)
dt <- delayed_learner_train(fancy_stack, washb_task)
plot(dt, color = FALSE, height = "300px")


## ----sl-fancy, eval=FALSE------------------------------------------------
## # run sl and predict on WASH
## sl_fancy <- Lrnr_sl$new(learners = fancy_stack, metalearner = metalearner)
## sl_fancy_feast <- sl_fancy$train(washb_task)
## sl_fancy_feast$print() %>%
##   kable(format = "markdown", digits = 3)


## ----CVsl----------------------------------------------------------------
washb_task_new <- make_sl3_Task(
  data = washb_data,
  covariates = covars,
  outcome = outcome,
  folds = make_folds(washb_data, fold_fun = folds_vfold, V = 2)
)
CVsl_fancy <- CV_lrnr_sl(sl_fit_fancy, washb_task_new, loss_squared_error)
CVsl_fancy %>%
  kable(format = "markdown", digits = 3)


## ----varimp--------------------------------------------------------------
washb_varimp <- varimp(sl_fit_fancy, loss_squared_error)
washb_varimp %>%
  kable(format = "markdown", digits = 3)


## ----ex-setup------------------------------------------------------------
# load the data set
db_data <-
 url("https://raw.githubusercontent.com/benkeser/sllecture/master/chspred.csv")
chspred <- read_csv(file = db_data, col_names = TRUE)
# take a quick peek
head(chspred, 3) %>%
  kable(format = "markdown", digits = 3)


## ----ex-key, eval=FALSE--------------------------------------------------
## chspred_task <- make_sl3_Task(
##   data = chspred,
##   covariates = head(colnames(chspred), -1),
##   outcome = "mi"
## )
## 
## glm_learner <- Lrnr_glm$new()
## lasso_learner <- Lrnr_glmnet$new(alpha = 1)
## ridge_learner <- Lrnr_glmnet$new(alpha = 0)
## enet_learner <- Lrnr_glmnet$new(alpha = 0.5)
## curated_glm_learner <- Lrnr_glm_fast$new(formula = "mi ~ smoke + beta + waist")
## mean_learner <- Lrnr_mean$new() # That is one mean learner!
## glm_fast_learner <- Lrnr_glm_fast$new()
## ranger_learner <- Lrnr_ranger$new()
## svm_learner <- Lrnr_svm$new()
## xgb_learner <- Lrnr_xgboost$new()
## 
## screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")
## glm_pipeline <- make_learner(Pipeline, screen_cor, glm_learner)
## ranger_pipeline <- make_learner(Pipeline, screen_cor, ranger_learner)
## 
## stack <- make_learner(
##   Stack,
##   glm_pipeline, ranger_pipeline, glm_learner,
##   lasso_learner, ridge_learner, enet_learner,
##   curated_glm_learner, mean_learner, glm_fast_learner,
##   ranger_learner, svm_learner, xgb_learner
## )
## 
## metalearner <- make_learner(Lrnr_nnls)
## 
## sl <- Lrnr_sl$new(
##   learners = stack,
##   metalearner = metalearner
## )
## sl_fit <- sl$train(task)
## sl_fit$print() %>%
##   kable(format = "markdown", digits = 3)

