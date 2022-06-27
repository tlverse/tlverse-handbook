## ----setup-handbook-utils-noecho, echo = FALSE--------------------------------
library(knitr)
library(kableExtra)
library(data.table)


## ----load-data----------------------------------------------------------------
library(data.table)
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)


## ----show-data-normal-noeval, eval = FALSE------------------------------------
## head(washb_data)

## ----show-data-handbook, echo = FALSE-----------------------------------------
if (knitr::is_latex_output()) {
  head(washb_data) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(washb_data) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----install-sl3, eval = FALSE------------------------------------------------
## library(devtools)
## install_github("tlverse/sl3@devel")


## ----load-sl3-----------------------------------------------------------------
library(sl3)


## ----task---------------------------------------------------------------------
# create the task (i.e., use washb_data to predict outcome using covariates)
task <- make_sl3_Task(
  data = washb_data,
  outcome = "whz",
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)

# let's examine the task
task


## ----list-properties----------------------------------------------------------
sl3_list_properties()


## ----list-learners------------------------------------------------------------
sl3_list_learners(properties = "continuous")


## ----learners-----------------------------------------------------------------
lrn_glm <- Lrnr_glm$new()
lrn_mean <- Lrnr_mean$new()


## ----more-learners------------------------------------------------------------
# penalized regressions:
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)


## ----more-learners-np---------------------------------------------------------
# spline regressions:
lrn_polspline <- Lrnr_polspline$new()
lrn_earth <- Lrnr_earth$new()

# fast highly adaptive lasso (HAL) implementation
lrn_hal <- Lrnr_hal9001$new(max_degree = 2, num_knots = c(3,2), nfolds = 5)

# tree-based methods
lrn_ranger <- Lrnr_ranger$new()
lrn_xgb <- Lrnr_xgboost$new()


## ----more-learners-final------------------------------------------------------
lrn_gam <- Lrnr_gam$new()
lrn_bayesglm <- Lrnr_bayesglm$new()


## ----stack--------------------------------------------------------------------
stack <- Stack$new(
  lrn_glm, lrn_mean, lrn_ridge, lrn_lasso, lrn_polspline, lrn_earth, lrn_hal, 
  lrn_ranger, lrn_xgb, lrn_gam, lrn_bayesglm
)
stack


## ----make-sl------------------------------------------------------------------
sl <- Lrnr_sl$new(learners = stack, metalearner = Lrnr_nnls$new())


## ----train-sl-----------------------------------------------------------------
start_time <- proc.time() # start time

set.seed(4197)
sl_fit <- sl$train(task = task)

runtime_sl_fit <- proc.time() - start_time # end time - start time = run time
runtime_sl_fit


## ----sl-predictions-----------------------------------------------------------
sl_preds <- sl_fit$predict(task = task)
head(sl_preds)


## ----glm-predictions----------------------------------------------------------
glm_preds <- sl_fit$learner_fits$Lrnr_glm_TRUE$predict(task = task)
head(glm_preds)


## ----glm-predictions-fullfit--------------------------------------------------
# we can also access the candidate learner full fits directly and obtain
# the same "full fit" candidate predictions from there 
# (we split this into two lines to avoid overflow)
stack_full_fits <- sl_fit$fit_object$full_fit$learner_fits$Stack$learner_fits
glm_preds_full_fit <- stack_full_fits$Lrnr_glm_TRUE$predict(task)

# check that they are identical
identical(glm_preds, glm_preds_full_fit)


## ----predvobs-----------------------------------------------------------------
# table of observed and predicted outcome values and arrange by observed values
df_plot <- data.table(
  Obs = washb_data[["whz"]], SL_Pred = sl_preds, GLM_Pred = glm_preds,
  Mean_Pred = sl_fit$learner_fits$Lrnr_mean$predict(task)
)
df_plot <- df_plot[order(df_plot$Obs), ] 

## ----predvobs-head, eval = FALSE----------------------------------------------
## head(df_plot)

## ----predvobs-head-handbook, echo = FALSE-------------------------------------
if (knitr::is_latex_output()) {
  head(df_plot) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(df_plot) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}

## ----predobs-plot, fig.asp = .55, fig.cap = "Observed and predicted values for weight-for-height z-score (whz)"----
# melt the table so we can plot observed and predicted values
df_plot$id <- seq(1:nrow(df_plot))
df_plot_melted <- melt(
  df_plot, id.vars = "id",
  measure.vars = c("Obs", "SL_Pred", "GLM_Pred", "Mean_Pred")
)

library(ggplot2)
ggplot(df_plot_melted, aes(id, value, color = variable)) + 
  geom_point(size = 0.1) + 
  labs(x = "Subjects (ordered by increasing whz)", 
       y = "whz") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 1)))


## ----cv-predictions-----------------------------------------------------------
# one way to obtain the CV predictions for the candidate learners
cv_preds_option1 <- sl_fit$fit_object$cv_fit$predict_fold(
  task = task, fold_number = "validation"
)
# another way to obtain the CV predictions for the candidate learners
cv_preds_option2 <- sl_fit$fit_object$cv_fit$predict(task = task)

# we can check that they are identical
identical(cv_preds_option1, cv_preds_option2)


## ----cv-predictions-head, eval = FALSE----------------------------------------
## head(cv_preds_option1)

## ----cv-predictions-head-handbook, echo = FALSE-------------------------------
if (knitr::is_latex_output()) {
  head(cv_preds_option1) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(cv_preds_option1) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----glm-predict-fold---------------------------------------------------------
full_fit_preds <- sl_fit$fit_object$cv_fit$predict_fold(
  task = task, fold_number = "full"
)
glm_full_fit_preds <- full_fit_preds$Lrnr_glm_TRUE

# check that they are identical
identical(glm_preds, glm_full_fit_preds)


## ----cv-predictions-long------------------------------------------------------
##### CV predictions "by hand" #####
# for each fold, i, we obtain validation set predictions:
cv_preds_list <- lapply(seq_along(task$folds), function(i){
  
  # get validation dataset for fold i:
  v_data <- task$data[task$folds[[i]]$validation_set, ]
  
  # get observed outcomes in fold i's validation dataset:
  v_outcomes <- v_data[["whz"]]

  # make task (for prediction) using fold i's validation dataset as data, 
  # and keeping all else the same:
  v_task <- make_sl3_Task(covariates = task$nodes$covariates, data = v_data)
  
  # get predicted outcomes for fold i's validation dataset, using candidates 
  # trained to fold i's training dataset
  v_preds <- sl_fit$fit_object$cv_fit$predict_fold(
    task = v_task, fold_number = i
  )
  # note: v_preds is a matrix of candidate learner predictions, where the 
  # number of rows is the number of observations in fold i's validation dataset 
  # and the number of columns is the number of candidate learners (excluding 
  # any that might have failed)
  
  # an identical way to get v_preds, which is used when we calculate the 
  # cv risk by hand in a later part of this chapter:
  # v_preds <- sl_fit$fit_object$cv_fit$fit_object$fold_fits[[i]]$predict(
  #   task = v_task
  # )
  
  # we will also return the row indices for fold i's validation set, so we 
  # can later reorder the CV predictions and make sure they are equal to what 
  # we obtained above
  return(list("v_preds" = v_preds, "v_index" = task$folds[[i]]$validation_set))
})

# extract the validation set predictions across all folds
cv_preds_byhand <- do.call(rbind, lapply(cv_preds_list, "[[", "v_preds"))

# extract the indices of validation set observations across all folds
# then reorder cv_preds_byhand to correspond to the ordering in the data
row_index_in_data <- unlist(lapply(cv_preds_list, "[[", "v_index"))
cv_preds_byhand_ordered <- cv_preds_byhand[order(row_index_in_data), ]
# now we can check that they are identical
identical(cv_preds_option1, cv_preds_byhand_ordered)


## ----predictions-new-task, eval = FALSE---------------------------------------
## # we do not evaluate this code chunk, as `washb_data_new` does not exist
## prediction_task <- make_sl3_Task(
##   data = washb_data_new, # assuming we have some new data for predictions
##   covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu",
##                  "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec",
##                  "floor", "walls", "roof", "asset_wardrobe", "asset_table",
##                  "asset_chair", "asset_khat", "asset_chouki", "asset_tv",
##                  "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach",
##                  "asset_mobile")
## )
## sl_preds_new_task <- sl_fit$predict(task = prediction_task)


## ----cf-predictions-static----------------------------------------------------
### 1. Copy data
tr_intervention_data <- data.table::copy(washb_data) 

### 2. Define intervention in copied dataset
tr_intervention <- rep("Nutrition + WSH", nrow(washb_data))
# NOTE: When we intervene on a categorical variable (such as "tr"), we need to 
#       define the intervention as a categorical variable (ie a factor).
#       Also, even though not all levels of the factor will be represented in 
#       the intervention, we still need this factor to reflect all of the 
#       levels that are present in the observed data
tr_levels <- levels(washb_data[["tr"]])
tr_levels
tr_intervention <- factor(tr_intervention, levels = tr_levels)
tr_intervention_data[,"tr" := tr_intervention, ]

### 3. Create a new sl3_Task
# note that we do not need to specify the outcome in this new task since we are 
# only using it to obtain predictions
tr_intervention_task <- make_sl3_Task(
  data = tr_intervention_data, 
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)
### 4. Get predicted values under intervention of interest
# SL predictions of what "whz" would have been had everyone received "tr" 
# equal to "Nutrition + WSH"
counterfactual_pred <- sl_fit$predict(tr_intervention_task)


## ----cf-predictions-dynamic---------------------------------------------------
dynamic_tr_intervention_data <- data.table::copy(washb_data) 

dynamic_tr_intervention <- ifelse(
  washb_data[["asset_refrig"]] == 1, "Nutrition + WSH", "WSH"
)
dynamic_tr_intervention <- factor(dynamic_tr_intervention, levels = tr_levels)
dynamic_tr_intervention_data[,"tr" := dynamic_tr_intervention, ]

dynamic_tr_intervention_task <- make_sl3_Task(
  data = dynamic_tr_intervention_data, 
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)
### 4. Get predicted values under intervention of interest
# SL predictions of what "whz" would have been had every subject received "tr" 
# equal to "Nutrition + WSH" if they had a fridge and "WSH" if they didn't have 
# a fridge
counterfactual_pred <- sl_fit$predict(dynamic_tr_intervention_task)


## ----sl-coefs-simple----------------------------------------------------------
round(sl_fit$coefficients, 3)


## ----metalearner-fit----------------------------------------------------------
metalrnr_fit <- sl_fit$fit_object$cv_meta_fit$fit_object
round(metalrnr_fit$coefficients, 3)


## ----sl-summary---------------------------------------------------------------
cv_risk_table <- sl_fit$cv_risk(eval_fun = loss_squared_error)

## ----cv-risk-summary, eval = FALSE--------------------------------------------
## cv_risk_table[,c(1:3)]

## ----cv-risk-summary-handbook, echo = FALSE-----------------------------------
if (knitr::is_latex_output()) {
  cv_risk_table[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  cv_risk_table[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----cv-risk-byhand-----------------------------------------------------------
##### CV risk "by hand" #####
# for each fold, i, we obtain predictive performance/risk for each candidate:
cv_risks_list <- lapply(seq_along(task$folds), function(i){
  
  # get validation dataset for fold i:
  v_data <- task$data[task$folds[[i]]$validation_set, ]
  
  # get observed outcomes in fold i's validation dataset:
  v_outcomes <- v_data[["whz"]]

  # make task (for prediction) using fold i's validation dataset as data, 
  # and keeping all else the same:
  v_task <- make_sl3_Task(covariates = task$nodes$covariates, data = v_data)
  
  # get predicted outcomes for fold i's validation dataset, using candidates 
  # trained to fold i's training dataset
  v_preds <- sl_fit$fit_object$cv_fit$fit_object$fold_fits[[i]]$predict(v_task)
  # note: v_preds is a matrix of candidate learner predictions, where the 
  # number of rows is the number of observations in fold i's validation dataset 
  # and the number of columns is the number of candidate learners (excluding 
  # any that might have failed)
  
  # calculate predictive performance for fold i for each candidate
  eval_function <- loss_squared_error # valid for estimation of conditional mean
  v_losses <- apply(v_preds, 2, eval_function, v_outcomes)
  cv_risks <- colMeans(v_losses)
  return(cv_risks)
})
# average the predictive performance across all folds for each candidate
cv_risks_byhand <- colMeans(do.call(rbind, cv_risks_list))
cv_risk_table_byhand <- data.table(
  learner = names(cv_risks_byhand), MSE = cv_risks_byhand
)
# check that the CV risks are identical when calculated by hand and function
# (ignoring small differences by rounding to the fourth decimal place)
identical(
  round(cv_risk_table_byhand$MSE,4), round(as.numeric(cv_risk_table$MSE),4)
)


## ----sl-summary-plot, eval = F------------------------------------------------
## 
## # Column "se" in the CV risk table is the standard error across all losses for
## # a learner, i.e., se = sd(loss)/sqrt(n), where loss is an n length vector of
## # validation set predictions across all folds, and n is the number of
## # validation set observations across all folds. We can use this to
## cv_risk_table[, "lower" := MSE - qnorm(.975)*se]
## cv_risk_table[, "upper" := MSE + qnorm(.975)*se]
## 
## ggplot(cv_risk_table,
##        aes_string(x = "learner", y = "MSE", ymin = "lower", ymax = "upper")) +
##   geom_pointrange() +
##   coord_flip() +
##   ylab("V-fold CV Risk Estimate") +
##   xlab("Learner")


## ----cvsl, eval = FALSE-------------------------------------------------------
## start_time <- proc.time()
## 
## set.seed(569)
## cv_sl_fit <- cv_sl(lrnr_sl = sl_fit, task = task, eval_fun = loss_squared_error)
## 
## runtime_cv_sl_fit <- proc.time() - start_time
## runtime_cv_sl_fit


## ----cvsl-save, eval = FALSE, echo = FALSE------------------------------------
## library(here)
## save(cv_sl_fit, file=here("data", "fit_objects", "cv_sl_fit.Rdata"), compress=T)
## save(runtime_cv_sl_fit, file=here("data", "fit_objects", "runtime_cv_sl_fit.Rdata"))


## ----cvsl-load, eval = TRUE, echo = FALSE-------------------------------------
library(here)
load(here("data", "fit_objects", "cv_sl_fit.Rdata"))
load(here("data", "fit_objects", "runtime_cv_sl_fit.Rdata"))
runtime_cv_sl_fit


## ----cvsl-risk-summary, eval = FALSE------------------------------------------
## cv_sl_fit$cv_risk[,c(1:3)]


## ----cvsl-risk-summary-handbook, echo = FALSE---------------------------------
if (knitr::is_latex_output()) {
  cv_sl_fit$cv_risk[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  cv_sl_fit$cv_risk[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----cvsl-risk-summary-coefs, eval = FALSE------------------------------------
## round(cv_sl_fit$coef, 3)

## ----cvsl-risk-summary-coefs-handbook, echo = FALSE---------------------------
if (knitr::is_latex_output()) {
  round(cv_sl_fit$coef, 3) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  round(cv_sl_fit$coef, 3) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----sl-revere-risk-----------------------------------------------------------
cv_risk_w_sl_revere <- sl_fit$cv_risk(
  eval_fun = loss_squared_error, get_sl_revere_risk = TRUE
)


## ----sl-revere-risk-summary, eval = FALSE-------------------------------------
## cv_risk_w_sl_revere[,c(1:3)]


## ----sl-revere-risk-handbook, echo = FALSE------------------------------------
if (knitr::is_latex_output()) {
  cv_risk_w_sl_revere[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  cv_risk_w_sl_revere[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----sl-revere-risk-byhand----------------------------------------------------
##### revere-based risk "by hand" #####
# for each fold, i, we obtain predictive performance/risk for the SL
sl_revere_risk_list <- lapply(seq_along(task$folds), function(i){
  # get validation dataset for fold i:
  v_data <- task$data[task$folds[[i]]$validation_set, ]
  
  # get observed outcomes in fold i's validation dataset:
  v_outcomes <- v_data[["whz"]]
  
  # make task (for prediction) using fold i's validation dataset as data, 
  # and keeping all else the same:
  v_task <- make_sl3_Task(
    covariates = task$nodes$covariates, data = v_data
  )
  
  # get predicted outcomes for fold i's validation dataset, using candidates 
  # trained to fold i's training dataset
  v_preds <- sl_fit$fit_object$cv_fit$fit_object$fold_fits[[i]]$predict(v_task)

  # make a metalevel task (for prediction with sl):
  v_meta_task <- make_sl3_Task(
    covariates = sl_fit$fit_object$cv_meta_task$nodes$covariates,
    data = v_preds
  )
  
  # get predicted outcomes for fold i's metalevel dataset, using the fitted
  # metalearner, cv_meta_fit 
  sl_revere_v_preds <- sl_fit$fit_object$cv_meta_fit$predict(task=v_meta_task)
  # note: cv_meta_fit was trained on the metalevel dataset, which contains the
  # candidates' cv predictions and validation dataset outcomes across ALL folds, 
  # so cv_meta_fit has already seen fold i's validation dataset outcomes.
  
  # calculate predictive performance for fold i for the SL
  eval_function <- loss_squared_error # valid for estimation of conditional mean
  # note: by evaluating the predictive performance of the SL using outcomes 
  # that were already seen by the metalearner, this is not a cross-validated 
  # measure of predictive performance for the SL. 
  sl_revere_v_loss <- eval_function(
    pred = sl_revere_v_preds, observed = v_outcomes
  )
  sl_revere_v_risk <- mean(sl_revere_v_loss)
  return(sl_revere_v_risk)
})
# average the predictive performance across all folds for the SL
sl_revere_risk_byhand <- mean(unlist(sl_revere_risk_list))
sl_revere_risk_byhand

# check that our calculation by hand equals what is output in cv_risk_table_revere
sl_revere_risk <- as.numeric(cv_risk_w_sl_revere[learner=="SuperLearner","MSE"])
sl_revere_risk


## ----make-dSL-----------------------------------------------------------------
cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
dSL <- Lrnr_sl$new(learners = stack, metalearner = cv_selector)


## ----fit-dSL------------------------------------------------------------------
set.seed(4197)
dSL_fit <- dSL$train(task)


## ----summarize-dSL-coefs------------------------------------------------------
round(dSL_fit$coefficients, 3)


## ----summarize-dSL-cv-risk----------------------------------------------------
dSL_cv_risk_table <- dSL_fit$cv_risk(eval_fun = loss_squared_error)

## ----summarize-dSL-cv-risk-tbl, eval = FALSE----------------------------------
## dSL_cv_risk_table[,c(1:3)]

## ----summarize-dSL-cv-risk-tbl-handbook, echo = FALSE-------------------------
if (knitr::is_latex_output()) {
  dSL_cv_risk_table[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  dSL_cv_risk_table[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----verify-dSL-preds---------------------------------------------------------
dSL_pred <- dSL_fit$predict(task)
earth_pred <- dSL_fit$learner_fits$Lrnr_earth_2_3_backward_0_1_0_0$predict(task)
identical(dSL_pred, earth_pred)


## ----recall-eSL---------------------------------------------------------------
# in the section 2 we defined Lrnr_sl as
# sl <- Lrnr_sl$new(learners = stack, metalearner = Lrnr_nnls$new())


## ----rename-eSL---------------------------------------------------------------
# let's rename it to clarify that this is an eSL that uses NNLS as meta-learner
eSL_metaNNLS <- sl


## ----eSL-in-stack-------------------------------------------------------------
stack_with_eSL <- Stack$new(stack, eSL_metaNNLS)


## ----eSL-in-dSL---------------------------------------------------------------
cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
dSL <- Lrnr_sl$new(learners = stack_with_eSL, metalearner = cv_selector)


## ----make-sl-discrete-multi-esl-----------------------------------------------
# instantiate more eSLs
eSL_metaNNLSconvex <- Lrnr_sl$new(
  learners = stack, metalearner = Lrnr_nnls$new(convex = TRUE)
)
eSL_metaLasso <- Lrnr_sl$new(learners = stack, metalearner = lrn_lasso)
eSL_metaEarth <- Lrnr_sl$new(learners = stack, metalearner = lrn_earth)
eSL_metaRanger <- Lrnr_sl$new(learners = stack, metalearner = lrn_ranger)
eSL_metaHAL <- Lrnr_sl$new(learners = stack, metalearner = lrn_hal)
# adding the eSLs to the stack that defined them
stack_with_eSLs <- Stack$new(
  stack, eSL_metaNNLS, eSL_metaNNLSconvex, eSL_metaLasso, eSL_metaEarth, 
  eSL_metaRanger, eSL_metaHAL
)
# specify dSL
dSL <- Lrnr_sl$new(learners = stack_with_eSLs, metalearner = cv_selector)


## ----fit-sl-parallel----------------------------------------------------------
# let's load the future package and set n-1 cores for parallel processing
library(future)
ncores <- availableCores()-1
ncores
plan(multicore, workers = ncores)
# now, let's re-train sl in parallel for demonstrative purposes
# we will also set a stopwatch so we can see how long this takes
start_time <- proc.time()

set.seed(4197)
sl_fit_parallel <- sl$train(task)

runtime_sl_fit_parallel <- proc.time() - start_time
runtime_sl_fit_parallel


## ----task-with-warning, warning=TRUE------------------------------------------
# create the task (i.e., use washb_data to predict outcome using covariates)
task <- make_sl3_Task(
  data = washb_data,
  outcome = "whz",
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)


## ----which-data-missing-------------------------------------------------------
# which columns have missing values, and how many observations are missing?
colSums(is.na(washb_data))


## ----data-missing-------------------------------------------------------------
some_rows_with_missingness <- which(!complete.cases(washb_data))[31:33]
# note: we chose 31:33 because missingness in momage & momheight is there
washb_data[some_rows_with_missingness, c("momage", "momheight")]


## ----task-data-imputed--------------------------------------------------------
task$data[some_rows_with_missingness,
          c("momage", "momheight", "delta_momage", "delta_momheight")]
colSums(is.na(task$data))


## ----kitty--------------------------------------------------------------------
cats <- c("calico", "tabby", "cow", "ragdoll", "mancoon", "dwarf", "calico")
cats <- factor(cats)
cats_onehot <- factor_to_indicators(cats)
cats_onehot


## ----show-X, eval = FALSE-----------------------------------------------------
## head(task$X)


## ----show-X-handbook, echo = FALSE--------------------------------------------
if (knitr::is_latex_output()) {
  head(task$X) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(task$X) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----stack-names--------------------------------------------------------------
stack


## ----name-glm, eval = FALSE---------------------------------------------------
## lrn_glm <- Lrnr_glm$new(name = "GLM")


## ----stack-pretty-------------------------------------------------------------
learners_pretty_names <- c(
  "GLM" = lrn_glm, "Mean" = lrn_mean, "Ridge" = lrn_ridge, 
  "Lasso" = lrn_lasso, "Polspline" = lrn_polspline, "Earth" = lrn_earth, 
  "HAL" = lrn_hal, "RF" = lrn_ranger, "XGBoost" = lrn_xgb, "GAM" = lrn_gam, 
  "BayesGLM" = lrn_bayesglm
)
stack_pretty_names <- Stack$new(learners_pretty_names)
stack_pretty_names


## ----lrnr-grid-diy------------------------------------------------------------
grid_params <- list(
  max_depth = c(3, 5, 8),
  eta = c(0.001, 0.1, 0.3),
  nrounds = 100
)
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})
xgb_learners


## ----lrnr-grid-diy-names------------------------------------------------------
names(xgb_learners) <- c(
  "XGBoost_depth3_eta.001", "XGBoost_depth5_eta.001", "XGBoost_depth8_eta.001", 
  "XGBoost_depth3_eta.1", "XGBoost_depth5_eta.1", "XGBoost_depth8_eta.1", 
  "XGBoost_depth3_eta.3", "XGBoost_depth5_eta.3", "XGBoost_depth8_eta.3"
)


## ----lrnr-grid-caret, eval = FALSE--------------------------------------------
## lrnr_nnet_autotune <- Lrnr_caret$new(method = "nnet", name = "NNET_autotune")


## ----interaction-learner------------------------------------------------------
lrnr_glm_interaction <- Lrnr_glm$new(formula = "~.^2")


## ----screener-properties------------------------------------------------------
sl3_list_learners(properties = "importance")


## ----screener-importance------------------------------------------------------
ranger_with_importance <- Lrnr_ranger$new(importance = "impurity_corrected")
RFscreen_top10 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 10
)
RFscreen_top10_glm <- Pipeline$new(RFscreen_top10, lrn_glm)


## ----screener-importance-stack------------------------------------------------
RFscreen_top10_stack <- Pipeline$new(RFscreen_top10, stack)


## ----screener-coefs-----------------------------------------------------------
lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso, threshold = 0)
lasso_screen_glm <- Pipeline$new(lasso_screen, lrn_glm)


## ----screener-stack-----------------------------------------------------------
lasso_screen_stack <- Pipeline$new(lasso_screen, stack)


## ----screener-corr------------------------------------------------------------
# select top 10 most correlated covariates
corRank_screen <- Lrnr_screener_correlation$new(
  type = "rank", num_screen = 10
)
corRank_screen_stack <- Pipeline$new(corRank_screen, stack)

# select covariates with correlation p-value below 0.05, and a minimum of 3
corP_screen <- Lrnr_screener_correlation$new(
  type = "threshold", pvalue_threshold = 0.05, min_screen = 3
)
corP_screen_stack <- Pipeline$new(corP_screen, stack)


## ----screener-augment---------------------------------------------------------
keepme <- c("aged", "momage")
# using corRank_screen as an example, but any instantiated screener can be 
# supplied as screener.
corRank_screen_augmented <- Lrnr_screener_augment$new(
  screener = corRank_screen, default_covariates = keepme
)
corRank_screen_augmented_glm <- Pipeline$new(corRank_screen_augmented, lrn_glm)


## ----screeners-stack----------------------------------------------------------
screeners_stack <- Stack$new(stack, corP_screen_stack, corRank_screen_stack, 
                             lasso_screen_stack, RFscreen_top10_stack)


## ----varimp-------------------------------------------------------------------
assets <- c("asset_wardrobe", "asset_table", "asset_chair", "asset_khat",
            "asset_chouki", "asset_tv", "asset_refrig", "asset_bike", 
            "asset_moto", "asset_sewmach", "asset_mobile", "Nlt18", "Ncomp", 
            "watmin", "elec", "floor", "walls", "roof")
set.seed(983)
washb_varimp <- importance(
  fit = sl_fit, eval_fun = loss_squared_error, type = "permute", 
  covariate_groups = list("assets" = assets)
)


## ----varimp-print, eval = FALSE-----------------------------------------------
## washb_varimp

## ----varimp-print-handbook, echo = FALSE--------------------------------------
if (knitr::is_latex_output()) {
  washb_varimp %>%
    kable(digits = 4, format = "latex")
} else if (knitr::is_html_output()) {
  washb_varimp %>%
    kable(digits = 4) %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----varimp-plot, fig.asp = .62, fig.cap = "sl3 variable importance for predicting weight-for-height z-score with WASH Benefits example dataset"----
# plot variable importance
importance_plot(x = washb_varimp)


## ----cde-using-locscale, eval = FALSE-----------------------------------------
## # semiparametric density estimator with homoscedastic errors (HOSE)
## hose_hal_lrnr <- Lrnr_density_semiparametric$new(
##   mean_learner = Lrnr_hal9001$new()
## )
## # semiparametric density estimator with heteroscedastic errors (HESE)
## hese_rf_glm_lrnr <- Lrnr_density_semiparametric$new(
##   mean_learner = Lrnr_ranger$new()
##   var_learner = Lrnr_glm$new()
## )
## 
## # SL for the conditional treatment density
## sl_dens_lrnr <- Lrnr_sl$new(
##   learners = list(hose_hal_lrnr, hese_rf_glm_lrnr),
##   metalearner = Lrnr_solnp_density$new()
## )


## ----cde-using-pooledhaz, eval = FALSE----------------------------------------
## # learners used for conditional densities for (g_n)
## haldensify_lrnr <- Lrnr_haldensify$new(
##   n_bins = c(5, 10)
## )


## ----ex-setup-----------------------------------------------------------------
# load the data set
library(readr)
db_data <- url(
  paste0(
    "https://raw.githubusercontent.com/benkeser/sllecture/master/",
    "chspred.csv"
  )
)
chspred <- read_csv(file = db_data, col_names = TRUE)


## ----ex-head, eval = FALSE----------------------------------------------------
## head(chspred)

## ----ex-head-handbook, echo = FALSE-------------------------------------------
if (knitr::is_latex_output()) {
  head(chspred) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(chspred) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----ex-key, eval=FALSE-------------------------------------------------------
## db_data <- url(
##   "https://raw.githubusercontent.com/benkeser/sllecture/master/chspred.csv"
## )
## chspred <- read_csv(file = db_data, col_names = TRUE)
## data.table::setDT(chspred)
## 
## # make task
## chspred_task <- make_sl3_Task(
##   data = chspred,
##   covariates = colnames(chspred)[-1],
##   outcome = "mi"
## )
## 
## # make learners
## glm_learner <- Lrnr_glm$new()
## lasso_learner <- Lrnr_glmnet$new(alpha = 1)
## ridge_learner <- Lrnr_glmnet$new(alpha = 0)
## enet_learner <- Lrnr_glmnet$new(alpha = 0.5)
## # curated_glm_learner uses formula = "mi ~ smoke + beta"
## curated_glm_learner <- Lrnr_glm_fast$new(covariates = c("smoke", "beta"))
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
## sl_fit$cv_risk(loss_squared_error)

