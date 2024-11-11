# Simulation example of CUSUM chart for observations following a piecewise constant baseline excess hazard model with proportional alternative
# First: Simulate observations between 2000 and 2010 from the model above used for fitting the semi-parametric model 
# Second: Simulate observations from the out-of-control alternative with rho = 0.75 and run the CUSUM chart based on the output from the semi-parametric model
rm(list = ls())
library(relsurv)
source("simple_examples_scripts/read_population_table.R")

start_year_train <- 2000
end_year_train <- 2010

new_pop_data_male <- pop.data.male[pop.data.male$Year >= 2010, ] # Only look at rows relevant for the monitoring time period
new_pop_data_male$Age[new_pop_data_male$Age == "110+"] <- 110
new_pop_data_male$Age <- as.integer(new_pop_data_male$Age)
new_pop_data_male <- as.matrix(new_pop_data_male)

new_pop_data_female <- pop.data.female[pop.data.female$Year >= 2010, ] 
new_pop_data_female$Age[new_pop_data_female$Age == "110+"] <- 110
new_pop_data_female$Age <- as.integer(new_pop_data_female$Age)
new_pop_data_female <- as.matrix(new_pop_data_female)

set.seed(1)

# Simulate arrival times as a Poisson process 

# From the 10 year survival piecewise model on real data --> on average 3660 observations per year
# Therefore choose lambda_arrival = 3750

arrival_times_vec_train <- cumsum(rexp(100000, rate = 3750))
arrival_times_vec_train <- arrival_times_vec_train[arrival_times_vec_train <= 10]
range(arrival_times_vec_train)

max_follow_up_vec_train <- end_year_train - (start_year_train + arrival_times_vec_train)
max_follow_up_vec_train[max_follow_up_vec_train >= 10] <- 10 # none of these here, but this is more general if the monitoring time period is longer than the number of years of interest when it comes to survival

n_obs_train <- length(arrival_times_vec_train)

# Simulate other relevant variables like age, gender etc. 

age_sim_vec_train <- rnorm(n_obs_train, mean = 75, sd = 10)
age_sim_vec_train[age_sim_vec_train <= 50] <- 50
age_sim_vec_train[age_sim_vec_train >= 105] <- 105

gender_sim_vec_train <- factor(sample(1:2, n_obs_train, replace = T)) # 1 = male, 2 = female like in relsurv.

seer_sim_vec_train <- sample(c("Distant", "Localised", "Regional", "Unknown"), n_obs_train, replace = T, prob = c(0.2, 0.2, 0.55, 0.05))
seer_sim_vec_train <- factor(seer_sim_vec_train, levels = c("Distant", "Localised", "Regional", "Unknown"))

icd_sim_vec_train <- factor(sample(0:3, n_obs_train, replace = T, prob = c(0.27, 0.44, 0.28, 0.01)))

surgery_sim_vec_train <- factor(sample(0:2, n_obs_train, replace = T, prob = c(0.8275, 0.1715, 0.001)))

morphology_type_train <- factor(sample(c("adeno-ca", "muc-adeno-ca"), n_obs_train, replace = T, prob = c(0.9, 0.1)))

# Collect all the simulated variables as a data frame

x_data_train <- data.frame(KJOENN = gender_sim_vec_train, SEER_STADIUM = seer_sim_vec_train, ICD7.indicator = icd_sim_vec_train, surgery.type = surgery_sim_vec_train, morphology_type = morphology_type_train)

x_matrix_train <- model.matrix( ~., data = x_data_train)[, -1] # need a covariate matrix for the CUSUM calculation 

# Prepare quantities for the excess hazard

beta_vec <- c(
  0.005, # Kjoenn = 2 (Female)
  -3, # SEER_STADIUM = Localised
  -1.75, # SEER_STADIUM = Regional
  -1, # SEER_STADIUM = Unknown
  0.5, # ICD7.indicator = 1
  0.2, # ICD7.indicator = 2
  0.3, # ICD7.indicator = 3
  1.5, # surgery.type = 1
  2.5, # surgery.type = 2
  -0.05 # morphology type = muc-adeno-ca
)

linear_predictor_vec_train <- x_matrix_train %*% beta_vec

follow_up_interval_partition <- c(seq(0, 5, by = 1), 10) # 6 bands, yearly bands for the first 5 years and a larger band between 5 and 10 years of follow up

chi_vec <- c(-1.4, -1.6, -1.8, -2, -2.1, -3)

baseline_excess_vec <- exp(chi_vec) # baseline excess hazard for a time point in a given band

# Simulate population survival times and population censoring indicator

u_pop_vec_train <- runif(n_obs_train)

matrix_time_pop_sim_train <- CUSUMrelsurv::vec_pop_sim(
  age_vec = age_sim_vec_train,
  gender_vec = gender_sim_vec_train,
  start_year = start_year_train,
  end_year = end_year_train,
  table_end_year = 2020,
  arrival_time_vec = arrival_times_vec_train,
  pop_matrix_male = new_pop_data_male,
  pop_matrix_female = new_pop_data_female,
  u = u_pop_vec_train
)

time_pop_sim_train <- matrix_time_pop_sim_train[, 1]
delta_p_vec_train <- matrix_time_pop_sim_train[, 2]

u_excess_vec_train <- runif(n_obs_train)

vec_time_excess_sim_train <- CUSUMrelsurv::vec_excess_sim_piecewise(
  partition_t_vec = follow_up_interval_partition,
  baseline_vec = baseline_excess_vec,
  max_follow_up_vec = max_follow_up_vec_train,
  u_vec = u_excess_vec_train,
  linear_predictor_vec = linear_predictor_vec_train
)

interim_censoring_rate <- 0.000275
time_censoring_sim_train <- cbind(rexp(n_obs_train, rate = interim_censoring_rate), max_follow_up_vec_train)
time_censoring_sim_train <- apply(time_censoring_sim_train, 1, FUN = min)

time_observed_sim_train <- pmin(time_pop_sim_train, vec_time_excess_sim_train, time_censoring_sim_train)
delta_i_vec_train <- pmax(delta_p_vec_train, as.numeric(vec_time_excess_sim_train < time_pop_sim_train)) * as.numeric(time_censoring_sim_train > time_observed_sim_train)
censoring_indices_train <- which(delta_i_vec_train == 0 & time_observed_sim_train < max_follow_up_vec_train) # consider interim censoring, not censoring due to no event at the end of follow up
estimated_censoring_rate_train <- length(censoring_indices_train) / sum(time_observed_sim_train) # MLE of exponential distributed event times with censoring data

time_observed_sim_train <- round(time_observed_sim_train * 365.241) / 365.241 # relsurv works faster with days, multiply back byb 365.241 when fitting the model using relsurv

data_train <- cbind(x_data_train, diagnosis_year = arrival_times_vec_train + start_year_train, age = age_sim_vec_train, time = time_observed_sim_train, status = delta_i_vec_train)

# Fit the semi-parametric model with the EM-algorithm based on this dataset 

em_excess_model <-
  rsadd(
    Surv(time * 365.241, status) ~ as.factor(KJOENN) +
      as.factor(SEER_STADIUM)
    + as.factor(ICD7.indicator) + as.factor(surgery.type) +
      as.factor(morphology_type),
    data = data_train,
    ratetable = nortab,
    method = "EM",
    bwin = 25,
    rmap = list(
      age = age * 365.241,
      sex = as.integer(KJOENN),
      year = (diagnosis_year - 1970) * 365.241
    )
  )

summary(em_excess_model)

estimated_beta_vec <- em_excess_model$coefficients

# Need to create some objects required by the CUSUM calculations from the output of the semiparametric model 

# First: Baseline hazard fit

estimated_baseline_values <- as.data.frame(cbind(em_excess_model$lambda0 * 365.241, em_excess_model$times / 365.241)) # everything in year scale 
colnames(estimated_baseline_values) <- c("baseline_haz_values", "time_in_years")

baseline_haz_fit <- smooth.spline(estimated_baseline_values$time_in_years, estimated_baseline_values$baseline_haz_values, df = 25) # also need the baseline hazard at any time point other than the times to event from data
# This is done by using smooth splines on the values obtained at the times to event in the training data. 

baseline_haz_vec_general_func <- function(t){ # t in years - create a wrapper function that returns estimated baseline hazard for any given time point t
  baseline_value <- predict(baseline_haz_fit, t)$y
  final_baseline_value <- ifelse(baseline_value >= 0, baseline_value, 0) # non-negative baseline
  final_baseline_value
}

# Second: Cumulative baseline hazard fit

cumulative_baseline_values <- as.data.frame(cbind(em_excess_model$Lambda0, em_excess_model$times / 365.241)) # everything in year scale 
colnames(cumulative_baseline_values) <- c("cumulative_baseline_haz_values", "time_in_years")

cumulative_baseline_haz_fit <- smooth.spline(cumulative_baseline_values$time_in_years, cumulative_baseline_values$cumulative_baseline_haz_values, df = 25) # similar as above

cumulative_baseline_haz_vec_general_func <- function(t){ # t in years - create a wrapper function that returns estimated cumulative baseline hazard for any given time point t
  cumulative_baseline_value <- ifelse(t > min(cumulative_baseline_values$time_in_years),
                                      predict(cumulative_baseline_haz_fit, t)$y,
                                      predict(cumulative_baseline_haz_fit, min(cumulative_baseline_values$time_in_years))$y / min(cumulative_baseline_values$time_in_years) * t) # linear interpolation
  cumulative_baseline_value
}

predict_cumulative_baseline_haz_vec_general_func <- function(t){ # t in years 
  predict(cumulative_baseline_haz_fit, t)$y
}

min_time_data_cumulative_baseline_haz_fit <- min(cumulative_baseline_values$time_in_years)

# Now simulate test data in the period 2010-2020 that will be used to calculate the CUSUM chart.
# To illustrate the method, the test observations experience the out-of-control after eta = 5 years of monitoring (i.e. 2015).
# The proportional constant rho is set to 0.5 for illustration purposes. 

set.seed(2)

start_year_test <- 2010
end_year_test <- 2020

rho <- 0.5
eta <- 5

arrival_times_vec_test <- cumsum(rexp(100000, rate = 3750))
arrival_times_vec_test <- arrival_times_vec_test[arrival_times_vec_test <= 10]
range(arrival_times_vec_test)

max_follow_up_vec_test <- end_year_test - (start_year_test + arrival_times_vec_test)
max_follow_up_vec_test[max_follow_up_vec_test >= 10] <- 10 # none of these here, but this is more general if the monitoring time period is longer than the number of years of interest when it comes to survival

n_obs_test <- length(arrival_times_vec_test)

# Simulate other relevant variables like age, gender etc. 

age_sim_vec_test <- rnorm(n_obs_test, mean = 75, sd = 10)
age_sim_vec_test[age_sim_vec_test <= 50] <- 50
age_sim_vec_test[age_sim_vec_test >= 105] <- 105

gender_sim_vec_test <- factor(sample(1:2, n_obs_test, replace = T)) # 1 = male, 2 = female like in relsurv.

seer_sim_vec_test <- sample(c("Distant", "Localised", "Regional", "Unknown"), n_obs_test, replace = T, prob = c(0.2, 0.2, 0.55, 0.05))
seer_sim_vec_test <- factor(seer_sim_vec_test, levels = c("Distant", "Localised", "Regional", "Unknown"))

icd_sim_vec_test <- factor(sample(0:3, n_obs_test, replace = T, prob = c(0.27, 0.44, 0.28, 0.01)))

surgery_sim_vec_test <- factor(sample(0:2, n_obs_test, replace = T, prob = c(0.8275, 0.1715, 0.001)))

morphology_type_test <- factor(sample(c("adeno-ca", "muc-adeno-ca"), n_obs_test, replace = T, prob = c(0.9, 0.1)))

# Collect all the simulated variables as a data frame

x_data_test <- data.frame(KJOENN = gender_sim_vec_test, SEER_STADIUM = seer_sim_vec_test, ICD7.indicator = icd_sim_vec_test, surgery.type = surgery_sim_vec_test, morphology_type = morphology_type_test)

x_matrix_test <- model.matrix( ~., data = x_data_test)[, -1] # need a covariate matrix for the CUSUM calculation 

# Prepare quantities for the excess hazard

linear_predictor_vec_test <- x_matrix_test %*% beta_vec # still simulate using true parameter values

# Simulate population survival times and population censoring indicator

u_pop_vec_test <- runif(n_obs_test)

matrix_time_pop_sim_test <- CUSUMrelsurv::vec_pop_sim(
  age_vec = age_sim_vec_test,
  gender_vec = gender_sim_vec_test,
  start_year = start_year_test,
  end_year = end_year_test,
  table_end_year = 2020,
  arrival_time_vec = arrival_times_vec_test,
  pop_matrix_male = new_pop_data_male,
  pop_matrix_female = new_pop_data_female,
  u = u_pop_vec_test
)

time_pop_sim_test <- matrix_time_pop_sim_test[, 1]
delta_p_vec_test <- matrix_time_pop_sim_test[, 2]

u_excess_vec_test <- runif(n_obs_test)

vec_time_excess_sim_test <- CUSUMrelsurv::vec_excess_sim_prop_piecewise_in_and_out(
  partition_t_vec = follow_up_interval_partition,
  baseline_vec = baseline_excess_vec,
  max_follow_up_vec = max_follow_up_vec_test,
  u_vec = u_excess_vec_test,
  arrival_time_vec = arrival_times_vec_test,
  linear_predictor_vec = linear_predictor_vec_test,
  rho = rho, 
  eta = eta
)

interim_censoring_rate <- 0.000275
time_censoring_sim_test <- cbind(rexp(n_obs_test, rate = interim_censoring_rate), max_follow_up_vec_test)
time_censoring_sim_test <- apply(time_censoring_sim_test, 1, FUN = min)

time_observed_sim_test <- pmin(time_pop_sim_test, vec_time_excess_sim_test, time_censoring_sim_test)
delta_i_vec_test <- pmax(delta_p_vec_test, as.numeric(vec_time_excess_sim_test < time_pop_sim_test)) * as.numeric(time_censoring_sim_test > time_observed_sim_test)
censoring_indices_test <- which(delta_i_vec_test == 0 & time_observed_sim_test < max_follow_up_vec_test) # consider interim censoring, not censoring due to no event at the end of follow up
estimated_censoring_rate_test <- length(censoring_indices_test) / sum(time_observed_sim_test) # MLE of exponential distributed event times with censoring data

time_observed_sim_test <- round(time_observed_sim_test * 365.241) / 365.241 # relsurv works faster with days, multiply back byb 365.241 when fitting the model using relsurv

data_test <- cbind(x_data_test, diagnosis_year = arrival_times_vec_test + start_year_test, age = age_sim_vec_test, time = time_observed_sim_test, status = delta_i_vec_test)

# Run CUSUM chart

monitoring_t_grid <- seq(from = 0, to = 10, by = 0.01)

cusum_chart <- CUSUMrelsurv::cusum_prop_r_t_EM_smoothing_spline(
  cumulative_baseline_haz_vec_general_func = cumulative_baseline_haz_vec_general_func,
  baseline_haz_vec_general_func = baseline_haz_vec_general_func,
  start_year = start_year_test,
  age_vec = age_sim_vec_test,
  gender_vec = gender_sim_vec_test,
  x_matrix = x_matrix_test,
  time_obs_vec = time_observed_sim_test,
  arrival_time_vec = arrival_times_vec_test,
  delta_i_vec = delta_i_vec_test,
  beta_vec = estimated_beta_vec, # Use estimated beta coefficients here
  rho = rho,
  t_grid = monitoring_t_grid,
  pop_data_male = new_pop_data_male,
  pop_data_female = new_pop_data_female,
  end_year_table = 2020
)

plot(monitoring_t_grid, cusum_chart[, 1], type = "l")
