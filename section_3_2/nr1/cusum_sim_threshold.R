# Simulation example of CUSUM piecewise inspired by the 10 year survival piecewise model based on the last 15 years of observations
rm(list = ls())
source("simple_examples_scripts/read_population_table.R")

new_pop_data_male <- pop.data.male[pop.data.male$Year >= 2010, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_male$Age[new_pop_data_male$Age == "110+"] <- 110
new_pop_data_male$Age <- as.integer(new_pop_data_male$Age)
new_pop_data_male <- as.matrix(new_pop_data_male)

new_pop_data_female <- pop.data.female[pop.data.female$Year >= 2010, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_female$Age[new_pop_data_female$Age == "110+"] <- 110
new_pop_data_female$Age <- as.integer(new_pop_data_female$Age)
new_pop_data_female <- as.matrix(new_pop_data_female)

# We simulate the threshold by choosing the value such that the false probability after 5 years is 0.05. Need therefore to simulate the distribution of max cusum value when in control. 

# Function to simulate patients/observations

sim_obs_func <- function() {
  
  # Simulate arrival times as a Poisson process 
  
  start_year <- 2010
  
  false_prob_t_grid <- seq(from = 0, to = 10, by = 0.01) # want to find the threshold based on the probability of false alarm during a 10 year period
  
  num_of_years <- false_prob_t_grid[length(false_prob_t_grid)]
  
  end_year <- start_year + num_of_years
  
  arrival_times_vec <- cumsum(rexp(5000, rate = 250)) # set a slightly larger arrival rate for this example - around 250 observations per year
  arrival_times_vec <- arrival_times_vec[arrival_times_vec <= num_of_years]
  range(arrival_times_vec)
  
  max_follow_up_vec <- num_of_years - arrival_times_vec # largest possible follow up time for an observation
  
  n_obs <- length(arrival_times_vec)
  
  # Simulate other relevant variables like age, gender etc. based on the following distributions
  
  age_sim_vec <- rnorm(n_obs, mean = 75, sd = 10)
  age_sim_vec[age_sim_vec <= 50] <- 50
  age_sim_vec[age_sim_vec >= 105] <- 105
  
  gender_sim_vec <- factor(sample(1:2, n_obs, replace = T)) # 1 = male, 2 = female like in relsurv.
  
  seer_sim_vec <- sample(c("Distant", "Localised", "Regional", "Unknown"), n_obs, replace = T, prob = c(0.2, 0.2, 0.55, 0.05))
  seer_sim_vec <- factor(seer_sim_vec, levels = c("Distant", "Localised", "Regional", "Unknown"))
  
  icd_sim_vec <- factor(sample(0:3, n_obs, replace = T, prob = c(0.27, 0.44, 0.28, 0.01)), levels = 0:3)
  
  surgery_sim_vec <- factor(sample(0:2, n_obs, replace = T, prob = c(0.8275, 0.1715, 0.001)), levels = 0:2)
  
  morphology_type <- factor(sample(c("adeno-ca", "muc-adeno-ca"), n_obs, replace = T, prob = c(0.9, 0.1)), levels = c("adeno-ca", "muc-adeno-ca"))
  
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
  
  baseline_excess_vec <- exp(chi_vec) # baseline excess hazard for a time point in a given band - chi above is defined to coincide with e.g. excess hazard models with piecewise constant baseline from relsurv
  
  interim_censoring_rate <- 0.000275 # interim censoring rate
  
  list( # need these arguments to run the sim_predefined_obs_func_cusum_threshold_piecewise_function_parallel - function.
    start_year = start_year,
    end_year = end_year,
    age_vec = age_sim_vec,
    gender_vec = gender_sim_vec,
    partition_t_vec = follow_up_interval_partition,
    baseline_vec = baseline_excess_vec,
    max_follow_up_vec = max_follow_up_vec,
    arrival_time_vec = arrival_times_vec,
    interim_censoring_rate = interim_censoring_rate,
    beta_vec = beta_vec,
    x_matrix = x_matrix,
    linear_predictor_vec = linear_predictor_vec,
    t_grid = false_prob_t_grid,
    n_obs = n_obs
  )
  
}

system.time(
  max_cusum_matrix <- CUSUMrelsurv::sim_predefined_obs_func_cusum_prop_threshold_piecewise_function_parallel(
    sim_obs_func = sim_obs_func,
    n_iterations = 10000,
    n_cores = 10, 
    rho_vec = c(0.8, 0.9, 1.1, 1.2),
    pop_data_male = new_pop_data_male,
    pop_data_female = new_pop_data_female,
    end_year_table = 2020,
    random_state = 1
  )
)

saveRDS(max_cusum_matrix, "section_3_2/nr1/sim_cusum_in_control_samples_for_threshold.RDS")
