## ----cv_fig, fig.show="hold", echo = FALSE------------------------------------
knitr::include_graphics("img/misc/SLKaiserNew.pdf")


## ----setup--------------------------------------------------------------------
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(SuperLearner)
library(origami)
library(sl3)
library(knitr)
library(kableExtra)

# load data set and take a peek
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)
head(washb_data) %>%
  kable() %>%
  kableExtra::kable_styling(fixed_thead = T) %>%
  scroll_box(width = "100%", height = "300px")


## ----task, warning=TRUE-------------------------------------------------------
# specify the outcome and covariates
outcome <- "whz"
covars <- colnames(washb_data)[-which(names(washb_data) == outcome)]

# create the sl3 task
washb_task <- make_sl3_Task(
  data = washb_data,
  covariates = covars,
  outcome = outcome
)


## ----task-examine-------------------------------------------------------------
washb_task


## ----task-folds-examine-------------------------------------------------------
length(washb_task$folds) # how many folds?

head(washb_task$folds[[1]]$training_set) # row indexes for fold 1 training
head(washb_task$folds[[1]]$validation_set) # row indexes for fold 1 validation

any(
  washb_task$folds[[1]]$training_set %in%
    washb_task$folds[[1]]$validation_set
)


## ----list-properties----------------------------------------------------------
sl3_list_properties()


## ----list-learners------------------------------------------------------------
sl3_list_learners("continuous")


## ----baselearners-------------------------------------------------------------
# choose base learners
lrn_glm <- make_learner(Lrnr_glm)
lrn_mean <- Lrnr_mean$new()


## ----extra-lrnr-awesome-------------------------------------------------------
lrn_lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_enet.5 <- make_learner(Lrnr_glmnet, alpha = 0.5)

lrn_polspline <- Lrnr_polspline$new()

lrn_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)

lrn_hal_faster <- Lrnr_hal9001$new(max_degree = 2, reduce_basis = 0.05)

xgb_fast <- Lrnr_xgboost$new() # default with nrounds = 20 is pretty fast
xgb_50 <- Lrnr_xgboost$new(nrounds = 50)


## ----interaction-learner------------------------------------------------------
interactions <- list(c("elec", "tr"), c("tr", "hfiacat"))
# main terms as well as the interactions above will be included
lrn_interaction <- make_learner(Lrnr_define_interactions, interactions)


## ----interaction-pipe---------------------------------------------------------
# we already instantiated a linear model learner above, no need to do it again
lrn_glm_interaction <- make_learner(Pipeline, lrn_interaction, lrn_glm)
lrn_glm_interaction


## ----extra-lrnr-woah----------------------------------------------------------
lrn_bayesglm <- Lrnr_pkg_SuperLearner$new("SL.bayesglm")


## ----extra-lrnr-mindblown-svm, eval = FALSE-----------------------------------
## # I like to crock pot my SLs
## grid_params <- list(
##   cost = c(0.01, 0.1, 1, 10, 100, 1000),
##   gamma = c(0.001, 0.01, 0.1, 1),
##   kernel = c("polynomial", "radial", "sigmoid"),
##   degree = c(1, 2, 3)
## )
## grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
## svm_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
##   do.call(Lrnr_svm$new, as.list(tuning_params))
## })

## ----extra-lrnr-mindblown-xgboost---------------------------------------------
grid_params <- list(
  max_depth = c(2, 4, 6),
  eta = c(0.001, 0.1, 0.3),
  nrounds = 100
)
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
grid

xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})
xgb_learners


## ----carotene, eval=FALSE-----------------------------------------------------
## # Unlike xgboost, I have no idea how to tune a neural net or BART machine, so
## # I let caret take the reins
## lrnr_caret_nnet <- make_learner(Lrnr_caret, algorithm = "nnet")
## lrnr_caret_bartMachine <- make_learner(Lrnr_caret,
##   algorithm = "bartMachine",
##   method = "boot", metric = "Accuracy",
##   tuneLength = 10
## )


## ----stack--------------------------------------------------------------------
stack <- make_learner(
  Stack, lrn_glm, lrn_polspline, lrn_enet.5, lrn_ridge, lrn_lasso, xgb_50
)
stack


## ----alt-stack----------------------------------------------------------------
# named vector of learners first
learners <- c(lrn_glm, lrn_polspline, lrn_enet.5, lrn_ridge, lrn_lasso, xgb_50)
names(learners) <- c(
  "glm", "polspline", "enet.5", "ridge", "lasso", "xgboost50"
)
# next make the stack
stack <- make_learner(Stack, learners)
# now the names are pretty
stack


## ----alt-stack-cv-------------------------------------------------------------
cv_stack <- Lrnr_cv$new(stack)
cv_stack


## ----screener-----------------------------------------------------------------
miniforest <- Lrnr_ranger$new(
  num.trees = 20, write.forest = FALSE,
  importance = "impurity_corrected"
)

# learner must already be instantiated, we did this when we created miniforest
screen_rf <- Lrnr_screener_importance$new(learner = miniforest, num_screen = 5)
screen_rf

# which covariates are selected on the full data?
screen_rf$train(washb_task)


## ----screener-augment---------------------------------------------------------
keepme <- c("aged", "momage")
# screener must already be instantiated, we did this when we created screen_rf
screen_augment_rf <- Lrnr_screener_augment$new(
  screener = screen_rf, default_covariates = keepme
)
screen_augment_rf


## ----screener-coefs-----------------------------------------------------------
# we already instantiated a lasso learner above, no need to do it again
screen_lasso <- Lrnr_screener_coefs$new(learner = lrn_lasso, threshold = 0)
screen_lasso


## ----screener-pipe------------------------------------------------------------
screen_rf_pipe <- make_learner(Pipeline, screen_rf, stack)
screen_lasso_pipe <- make_learner(Pipeline, screen_lasso, stack)


## ----screeners-stack----------------------------------------------------------
# pretty names again
learners2 <- c(learners, screen_rf_pipe, screen_lasso_pipe)
names(learners2) <- c(names(learners), "randomforest_screen", "lasso_screen")

fancy_stack <- make_learner(Stack, learners2)
fancy_stack


## ----make-sl------------------------------------------------------------------
sl <- make_learner(Lrnr_sl, learners = fancy_stack)


## ----make-sl-discrete---------------------------------------------------------
discrete_sl_metalrn <- Lrnr_cv_selector$new()
discrete_sl <- Lrnr_sl$new(
  learners = fancy_stack,
  metalearner = discrete_sl_metalrn
)


## ----make-sl-plot-------------------------------------------------------------
dt_sl <- delayed_learner_train(sl, washb_task)
plot(dt_sl, color = FALSE, height = "400px", width = "90%")


## ----sl-----------------------------------------------------------------------
set.seed(4197)
sl_fit <- sl$train(washb_task)


## ----sl-predictions-----------------------------------------------------------
# we did it! now we have SL predictions
sl_preds <- sl_fit$predict()
head(sl_preds)


## ---- plot-predvobs-woohoo, eval=FALSE----------------------------------------
##
## # df_plot <- data.frame(Observed = washb_data[["whz"]], Predicted = sl_preds,
## #                        count = seq(1:nrow(washb_data))
##
## # df_plot_melted <- melt(df_plot, id.vars = "count",
## #                         measure.vars = c("Observed", "Predicted"))
##
## # ggplot(df_plot_melted, aes(value, count, color = variable)) + geom_point()


## ---- sl-summary--------------------------------------------------------------
sl_fit_summary <- sl_fit$print()


## ----CVsl---------------------------------------------------------------------
washb_task_new <- make_sl3_Task(
  data = washb_data,
  covariates = covars,
  outcome = outcome,
  folds = origami::make_folds(washb_data, fold_fun = folds_vfold, V = 2)
)
CVsl <- CV_lrnr_sl(
  lrnr_sl = sl_fit, task = washb_task_new, loss_fun = loss_squared_error
)
CVsl %>%
  kable(digits = 4) %>%
  kableExtra::kable_styling(fixed_thead = TRUE) %>%
  scroll_box(width = "100%", height = "300px")


## ----varimp-------------------------------------------------------------------
washb_varimp <- importance(sl_fit, loss = loss_squared_error, type = "permute")
washb_varimp %>%
  kable(digits = 4) %>%
  kableExtra::kable_styling(fixed_thead = TRUE) %>%
  scroll_box(width = "100%", height = "300px")


## ----varimp-plot--------------------------------------------------------------
# plot variable importance
importance_plot(
  washb_varimp,
  main = "sl3 Variable Importance for WASH Benefits Example Data"
)


## ----ex-setup-----------------------------------------------------------------
# load the data set
db_data <- url(
  paste0(
    "https://raw.githubusercontent.com/benkeser/sllecture/master/",
    "chspred.csv"
  )
)
chspred <- read_csv(file = db_data, col_names = TRUE)

# take a quick peek
head(chspred) %>%
  kable(digits = 4) %>%
  kableExtra::kable_styling(fixed_thead = TRUE) %>%
  scroll_box(width = "100%", height = "300px")


## ----ex-setup2----------------------------------------------------------------
ist_data <- paste0(
  "https://raw.githubusercontent.com/tlverse/",
  "tlverse-handbook/master/data/ist_sample.csv"
) %>% fread()

# number 3 help
ist_task_CVsl <- make_sl3_Task(
  data = ist_data,
  outcome = "DRSISC",
  covariates = colnames(ist_data)[-which(names(ist_data) == "DRSISC")],
  drop_missing_outcome = TRUE,
  folds = origami::make_folds(
    n = sum(!is.na(ist_data$DRSISC)),
    fold_fun = folds_vfold,
    V = 5
  )
)


## ----ex-key, eval=FALSE-------------------------------------------------------
## db_data <- url(
##   "https://raw.githubusercontent.com/benkeser/sllecture/master/chspred.csv"
## )
## chspred <- read_csv(file = db_data, col_names = TRUE)
##
## # make task
## chspred_task <- make_sl3_Task(
##   data = chspred,
##   covariates = head(colnames(chspred), -1),
##   outcome = "mi"
## )
##
## # make learners
## glm_learner <- Lrnr_glm$new()
## lasso_learner <- Lrnr_glmnet$new(alpha = 1)
## ridge_learner <- Lrnr_glmnet$new(alpha = 0)
## enet_learner <- Lrnr_glmnet$new(alpha = 0.5)
## # curated_glm_learner uses formula = "mi ~ smoke + beta + waist"
## curated_glm_learner <- Lrnr_glm_fast$new(covariates = c("smoke, beta, waist"))
## mean_learner <- Lrnr_mean$new() # That is one mean learner!
## glm_fast_learner <- Lrnr_glm_fast$new()
## ranger_learner <- Lrnr_ranger$new()
## svm_learner <- Lrnr_svm$new()
## xgb_learner <- Lrnr_xgboost$new()
##
## # screening
## screen_cor <- make_learner(Lrnr_screener_correlation)
## glm_pipeline <- make_learner(Pipeline, screen_cor, glm_learner)
##
## # stack learners together
## stack <- make_learner(
##   Stack,
##   glm_pipeline, glm_learner,
##   lasso_learner, ridge_learner, enet_learner,
##   curated_glm_learner, mean_learner, glm_fast_learner,
##   ranger_learner, svm_learner, xgb_learner
## )
##
## # make and train SL
## sl <- Lrnr_sl$new(
##   learners = stack
## )
## sl_fit <- sl$train(chspred_task)
## sl_fit$print()
##
## CVsl <- CV_lrnr_sl(sl_fit, chspred_task, loss_loglik_binomial)
## CVsl
##
## varimp <- importance(sl_fit, type = "permute")
## varimp %>%
##   importance_plot(
##     main = "sl3 Variable Importance for Myocardial Infarction Prediction"
##   )


## ----ex2-key, eval=FALSE------------------------------------------------------
## library(ROCR) # for AUC calculation
##
## ist_data <- paste0(
##   "https://raw.githubusercontent.com/tlverse/",
##   "tlverse-handbook/master/data/ist_sample.csv"
## ) %>% fread()
##
## # stack
## ist_task <- make_sl3_Task(
##   data = ist_data,
##   outcome = "DRSISC",
##   covariates = colnames(ist_data)[-which(names(ist_data) == "DRSISC")],
##   drop_missing_outcome = TRUE
## )
##
## # learner library
## lrn_glm <- Lrnr_glm$new()
## lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
## lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
## lrn_enet <- Lrnr_glmnet$new(alpha = 0.5)
## lrn_mean <- Lrnr_mean$new()
## lrn_ranger <- Lrnr_ranger$new()
## lrn_svm <- Lrnr_svm$new()
## # xgboost grid
## grid_params <- list(
##   max_depth = c(2, 5, 8),
##   eta = c(0.01, 0.15, 0.3)
## )
## grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
## params_default <- list(nthread = getOption("sl.cores.learners", 1))
## xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
##   do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))
## })
## learners <- unlist(list(
##   xgb_learners, lrn_ridge, lrn_mean, lrn_lasso,
##   lrn_glm, lrn_enet, lrn_ranger, lrn_svm
## ),
## recursive = TRUE
## )
##
## # SL
## sl <- Lrnr_sl$new(learners)
## sl_fit <- sl$train(ist_task)
##
## # AUC
## preds <- sl_fit$predict()
## obs <- c(na.omit(ist_data$DRSISC))
## AUC <- performance(prediction(sl_preds, obs), measure = "auc")@y.values[[1]]
## plot(performance(prediction(sl_preds, obs), "tpr", "fpr"))
##
## # CVsl
## ist_task_CVsl <- make_sl3_Task(
##   data = ist_data,
##   outcome = "DRSISC",
##   covariates = colnames(ist_data)[-which(names(ist_data) == "DRSISC")],
##   drop_missing_outcome = TRUE,
##   folds = origami::make_folds(
##     n = sum(!is.na(ist_data$DRSISC)),
##     fold_fun = folds_vfold,
##     V = 5
##   )
## )
## CVsl <- CV_lrnr_sl(sl_fit, ist_task_CVsl, loss_loglik_binomial)
## CVsl
##
## # sl3 variable importance plot
## ist_varimp <- importance(sl_fit, type = "permute")
## ist_varimp %>%
##   importance_plot(
##     main = "Variable Importance for Predicting Recurrent Ischemic Stroke"
##   )


## ----ex3-key, eval=FALSE------------------------------------------------------
## # TODO
