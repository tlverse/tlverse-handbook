## -----------------------------------------------------------------------------
library(tmle3)
library(sl3)


## -----------------------------------------------------------------------------
vet_data <- read.csv("https://raw.githubusercontent.com/tlverse/deming2019-workshop/master/data/veteran.csv")
vet_data$trt <- vet_data$trt - 1
# make fewer times for illustration
vet_data$time <- ceiling(vet_data$time / 20) 
k_grid <- 1:max(vet_data$time)
head(vet_data)


## -----------------------------------------------------------------------------
var_types <- list(T_tilde = Variable_Type$new("continuous"),
                  t = Variable_Type$new("continuous"),
                  Delta = Variable_Type$new("binomial"))
survival_spec <- tmle_survival(treatment_level = 1, control_level = 0,
                                 target_times = intersect(1:10, k_grid),
                                 variable_types = var_types)
node_list <- list(W = c("celltype", "karno", "diagtime", "age", "prior"),
                  A = "trt", T_tilde = "time", Delta = "status", id ="X")

long_data_tuple <- survival_spec$transform_data(vet_data, node_list)
df_long <- long_data_tuple$long_data
long_node_list <- long_data_tuple$long_node_list


## -----------------------------------------------------------------------------
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(A = sl_A, N = sl_A, A_c = sl_A)


## -----------------------------------------------------------------------------
tmle_task <- survival_spec$make_tmle_task(df_long, long_node_list)
initial_likelihood <- survival_spec$make_initial_likelihood(tmle_task,
                                                            learner_list)


## -----------------------------------------------------------------------------
up <- tmle3_Update_survival$new(
    maxit = 3e1, 
    cvtmle = TRUE,
    convergence_type = "scaled_var",
    delta_epsilon = 1e-2,
    fit_method = "l2",
    use_best = TRUE,
    verbose=FALSE
  )

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                               updater = up)
tmle_params <- survival_spec$make_params(tmle_task, targeted_likelihood)
tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)


## -----------------------------------------------------------------------------
print(tmle_fit_manual)

