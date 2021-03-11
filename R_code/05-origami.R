## ----setup--------------------------------------------------------------------
library(data.table)
library(origami)
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

washb_data <- washb_data[1:30, ]
head(washb_data) %>%
  kable() %>%
  kableExtra::kable_styling(fixed_thead = TRUE) %>%
  scroll_box(width = "100%", height = "300px")


## ----resubstitution-----------------------------------------------------------
folds_resubstitution(nrow(washb_data))


## ----loo----------------------------------------------------------------------
folds <- folds_loo(nrow(washb_data))
folds[[1]]
folds[[2]]


## ----cv-----------------------------------------------------------------------
folds <- folds_vfold(nrow(washb_data), V = 2)
folds[[1]]
folds[[2]]


## ----montecarlo---------------------------------------------------------------
folds <- folds_montecarlo(nrow(washb_data), V = 2, pvalidation = 0.2)
folds[[1]]
folds[[2]]


## ----bootstrap----------------------------------------------------------------
folds <- folds_bootstrap(nrow(washb_data), V = 2)
folds[[1]]
folds[[2]]


## ----plot_airpass-------------------------------------------------------------
library(ggfortify)

data(AirPassengers)
AP <- AirPassengers

autoplot(AP) +
  labs(
    x = "Date",
    y = "Passenger numbers (1000's)",
    title = "Air Passengers from 1949 to 1961"
  )

t <- length(AP)


## ---- fig.cap="Rolling origin CV", results="asis", echo=FALSE-----------------
knitr::include_graphics(path = "img/image/rolling_origin.png")


## ----rolling_origin-----------------------------------------------------------
folds <- folds_rolling_origin(
  t,
  first_window = 50, validation_size = 10, gap = 5, batch = 20
)
folds[[1]]
folds[[2]]


## ---- fig.cap="Rolling window CV", results="asis", echo=FALSE-----------------
knitr::include_graphics(path = "img/image/rolling_window.png")


## ----rolling_window-----------------------------------------------------------
folds <- folds_rolling_window(
  t,
  window_size = 50, validation_size = 10, gap = 5, batch = 20
)
folds[[1]]
folds[[2]]


## ---- fig.cap="Rolling origin V-fold CV", results="asis", echo=FALSE----------
knitr::include_graphics(path = "img/image/rolling_origin_v_fold.png")


## ---- fig.cap="Rolling window V-fold CV", results="asis", echo=FALSE----------
knitr::include_graphics(path = "img/image/rolling_window_v_fold.png")


## ----setup_ex-----------------------------------------------------------------
library(stringr)
library(dplyr)
library(tidyr)

# load data set and take a peek
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)

# Remove missing data, then pick just the first 500 rows
washb_data <- washb_data %>%
  drop_na() %>%
  slice(1:500)

outcome <- "whz"
covars <- colnames(washb_data)[-which(names(washb_data) == outcome)]

head(washb_data) %>%
  kable() %>%
  kableExtra::kable_styling(fixed_thead = TRUE) %>%
  scroll_box(width = "100%", height = "300px")


## ----covariates---------------------------------------------------------------
outcome
covars


## ----linear_mod---------------------------------------------------------------
lm_mod <- lm(whz ~ ., data = washb_data)
summary(lm_mod)


## ----get_naive_error----------------------------------------------------------
(err <- mean(resid(lm_mod)^2))


## ----define_fun_cv_lm---------------------------------------------------------
cv_lm <- function(fold, data, reg_form) {
  # get name and index of outcome variable from regression formula
  out_var <- as.character(unlist(str_split(reg_form, " "))[1])
  out_var_ind <- as.numeric(which(colnames(data) == out_var))

  # split up data into training and validation sets
  train_data <- training(data)
  valid_data <- validation(data)

  # fit linear model on training set and predict on validation set
  mod <- lm(as.formula(reg_form), data = train_data)
  preds <- predict(mod, newdata = valid_data)
  valid_data <- as.data.frame(valid_data)

  # capture results to be returned as output
  out <- list(
    coef = data.frame(t(coef(mod))),
    SE = (preds - valid_data[, out_var_ind])^2
  )
  return(out)
}


## ----cv_lm_resub--------------------------------------------------------------
# re-substitution estimate
resub <- make_folds(washb_data, fold_fun = folds_resubstitution)[[1]]
resub_results <- cv_lm(fold = resub, data = washb_data, reg_form = "whz ~ .")
mean(resub_results$SE, na.rm = TRUE)


## ----cv_lm_cross_valdate------------------------------------------------------
# cross-validated estimate
folds <- make_folds(washb_data)
cvlm_results <- cross_validate(
  cv_fun = cv_lm, folds = folds, data = washb_data, reg_form = "whz ~ .",
  use_future = FALSE
)
mean(cvlm_results$SE, na.rm = TRUE)


## ----cv_fun_randomForest------------------------------------------------------
# make sure to load the package!
library(randomForest)

cv_rf <- function(fold, data, reg_form) {
  # get name and index of outcome variable from regression formula
  out_var <- as.character(unlist(str_split(reg_form, " "))[1])
  out_var_ind <- as.numeric(which(colnames(data) == out_var))

  # define training and validation sets based on input object of class "folds"
  train_data <- training(data)
  valid_data <- validation(data)

  # fit Random Forest regression on training set and predict on holdout set
  mod <- randomForest(formula = as.formula(reg_form), data = train_data)
  preds <- predict(mod, newdata = valid_data)
  valid_data <- as.data.frame(valid_data)

  # define output object to be returned as list (for flexibility)
  out <- list(
    coef = data.frame(mod$coefs),
    SE = ((preds - valid_data[, out_var_ind])^2)
  )
  return(out)
}


## ----cv_fun_randomForest_run--------------------------------------------------
# now, let's cross-validate...
folds <- make_folds(washb_data)
cvrf_results <- cross_validate(
  cv_fun = cv_rf, folds = folds, data = washb_data, reg_form = "whz ~ .",
  use_future = FALSE
)
mean(cvrf_results$SE)


## ----load_airpass-------------------------------------------------------------
data(AirPassengers)
print(AirPassengers)


## ----folds_airpass------------------------------------------------------------
folds <- make_folds(AirPassengers,
  fold_fun = folds_rolling_origin,
  first_window = 36, validation_size = 24, batch = 10
)

# How many folds where generated?
length(folds)

# Examine the first 2 folds.
folds[[1]]
folds[[2]]


## ----fit_airpass--------------------------------------------------------------
# make sure to load the package!
library(forecast)

# function to calculate cross-validated squared error
cv_forecasts <- function(fold, data) {
  # Get training and validation data
  train_data <- training(data)
  valid_data <- validation(data)
  valid_size <- length(valid_data)

  train_ts <- ts(log10(train_data), frequency = 12)

  # First arima model
  arima_fit <- arima(train_ts, c(0, 1, 1),
    seasonal = list(
      order = c(0, 1, 1),
      period = 12
    )
  )
  raw_arima_pred <- predict(arima_fit, n.ahead = valid_size)
  arima_pred <- 10^raw_arima_pred$pred
  arima_MSE <- mean((arima_pred - valid_data)^2)

  # Second arima model
  arima_fit2 <- arima(train_ts, c(5, 1, 1),
    seasonal = list(
      order = c(0, 1, 1),
      period = 12
    )
  )
  raw_arima_pred2 <- predict(arima_fit2, n.ahead = valid_size)
  arima_pred2 <- 10^raw_arima_pred2$pred
  arima_MSE2 <- mean((arima_pred2 - valid_data)^2)

  out <- list(mse = data.frame(
    fold = fold_index(),
    arima = arima_MSE, arima2 = arima_MSE2
  ))
  return(out)
}

mses <- cross_validate(
  cv_fun = cv_forecasts, folds = folds, data = AirPassengers,
  use_future = FALSE
)
mses$mse
colMeans(mses$mse[, c("arima", "arima2")])
