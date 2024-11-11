# Simulation example of CUSUM chart for a piecewise constant baseline excess hazard model with proportional alternative
rm(list = ls())
source("simple_examples_scripts/read_population_table.R")

start_year <- 2010
end_year <- 2020

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

arrival_times_vec <- cumsum(rexp(100000, rate = 3750))
arrival_times_vec <- arrival_times_vec[arrival_times_vec <= 10]
range(arrival_times_vec)

max_follow_up_vec <- end_year - (start_year + arrival_times_vec)
max_follow_up_vec[max_follow_up_vec >= 10] <- 10 # none of these here, but this is more general if the monitoring time period is longer than the number of years of interest when it comes to survival

n_obs <- length(arrival_times_vec)

# Simulate other relevant variables like age, gender etc. 

age_sim_vec <- rnorm(n_obs, mean = 75, sd = 10)
age_sim_vec[age_sim_vec <= 50] <- 50
age_sim_vec[age_sim_vec >= 105] <- 105

gender_sim_vec <- factor(sample(1:2, n_obs, replace = T)) # 1 = male, 2 = female like in relsurv.

seer_sim_vec <- sample(c("Distant", "Localised", "Regional", "Unknown"), n_obs, replace = T, prob = c(0.2, 0.2, 0.55, 0.05))
seer_sim_vec <- factor(seer_sim_vec, levels = c("Distant", "Localised", "Regional", "Unknown"))

icd_sim_vec <- factor(sample(0:3, n_obs, replace = T, prob = c(0.27, 0.44, 0.28, 0.01)))

surgery_sim_vec <- factor(sample(0:2, n_obs, replace = T, prob = c(0.8275, 0.1715, 0.001)))

morphology_type <- factor(sample(c("adeno-ca", "muc-adeno-ca"), n_obs, replace = T, prob = c(0.9, 0.1)))

# Collect all the simulated variables as a data frame

x_data <- data.frame(KJOENN = gender_sim_vec, SEER_STADIUM = seer_sim_vec, ICD7.indicator = icd_sim_vec, surgery.type = surgery_sim_vec, morphology_type = morphology_type)

x_matrix <- model.matrix( ~., data = x_data)[, -1] # need a covariate matrix for the CUSUM calculation 

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

linear_predictor_vec <- x_matrix %*% beta_vec

follow_up_interval_partition <- c(seq(0, 5, by = 1), 10) # 6 bands, yearly bands for the first 5 years and a larger band between 5 and 10 years of follow up

chi_vec <- c(-1.4, -1.6, -1.8, -2, -2.1, -3)

baseline_excess_vec <- exp(chi_vec) # baseline excess hazard for a time point in a given band

# Simulate population survival times and population censoring indicator

u_pop_vec <- runif(n_obs)

matrix_time_pop_sim <- CUSUMrelsurv::vec_pop_sim(
  age_vec = age_sim_vec,
  gender_vec = gender_sim_vec,
  start_year = start_year,
  end_year = end_year,
  table_end_year = 2020,
  arrival_time_vec = arrival_times_vec,
  pop_matrix_male = new_pop_data_male,
  pop_matrix_female = new_pop_data_female,
  u = u_pop_vec
)

time_pop_sim <- matrix_time_pop_sim[, 1]
delta_p_vec <- matrix_time_pop_sim[, 2]

# First 5 year - in control (from 2010 - 2015), from 2015 and beyond we have out of control

rho <- 0.75 # the hazard is reduced by 25% when out of control in comparison to in control period
eta <- 5

u_excess_vec <- runif(n_obs)

vec_time_excess_sim <- CUSUMrelsurv::vec_excess_sim_prop_piecewise_in_and_out(
  partition_t_vec = follow_up_interval_partition,
  baseline_vec = baseline_excess_vec,
  max_follow_up_vec = max_follow_up_vec,
  u_vec = u_excess_vec,
  arrival_time_vec = arrival_times_vec,
  linear_predictor_vec = linear_predictor_vec,
  rho = rho,
  eta = eta
)

interim_censoring_rate <- 0.000275
time_censoring_sim <- cbind(rexp(n_obs, rate = interim_censoring_rate), max_follow_up_vec)
time_censoring_sim <- apply(time_censoring_sim, 1, FUN = min)

time_observed_sim <- pmin(time_pop_sim, vec_time_excess_sim, time_censoring_sim)
delta_i_vec <- pmax(delta_p_vec, as.numeric(vec_time_excess_sim < time_pop_sim)) * as.numeric(time_censoring_sim > time_observed_sim)

# Run CUSUM chart

monitoring_t_grid <- seq(from = 0, to = 10, by = 0.01)

cusum_chart <- CUSUMrelsurv::cusum_prop_r_t_piecewise(
  partition_t_vec = follow_up_interval_partition,
  baseline_vec = baseline_excess_vec,
  start_year = start_year,
  age_vec = age_sim_vec,
  gender_vec = gender_sim_vec,
  x_matrix = x_matrix,
  time_obs_vec = time_observed_sim,
  arrival_time_vec = arrival_times_vec,
  delta_i_vec = delta_i_vec,
  beta_vec = beta_vec,
  rho = rho,
  t_grid = monitoring_t_grid,
  pop_data_male = new_pop_data_male,
  pop_data_female = new_pop_data_female,
  end_year_table = 2020
)

plot(monitoring_t_grid, cusum_chart[, 1], type = "l")
