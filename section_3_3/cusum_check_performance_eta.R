# Simulation example of CUSUM piecewise inspired by the 10 year survival piecewise model based on the last 15 years of observations
rm(list = ls())
source("simple_examples_scripts/read_population_table.R")
library(doParallel)

start_year <- 2010
end_year <- 2020

new_pop_data_male <- pop.data.male[pop.data.male$Year >= 2010, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_male$Age[new_pop_data_male$Age == "110+"] <- 110
new_pop_data_male$Age <- as.integer(new_pop_data_male$Age)
new_pop_data_male <- as.matrix(new_pop_data_male)

new_pop_data_female <- pop.data.female[pop.data.female$Year >= 2010, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_female$Age[new_pop_data_female$Age == "110+"] <- 110
new_pop_data_female$Age <- as.integer(new_pop_data_female$Age)
new_pop_data_female <- as.matrix(new_pop_data_female)

# Function to simulate patients/observations

sim_obs_func <- function() {
  
  # Simulate arrival times as a Poisson process 
  
  start_year <- 2010
  
  end_year <- 2020
  
  num_of_years <- end_year - start_year
  
  t_grid <- seq(from = 0, to = num_of_years, by = 0.01)
  
  arrival_times_vec <- cumsum(rexp(50000, rate = 3700)) # set a slightly larger arrival rate for this example - around 250 observations per year
  arrival_times_vec <- arrival_times_vec[arrival_times_vec <= num_of_years]
  range(arrival_times_vec)
  
  max_follow_up_vec <- end_year - (start_year + arrival_times_vec)
  max_follow_up_vec[max_follow_up_vec >= num_of_years] <- num_of_years # none of these here, but this is more general if the monitoring time period is longer than the number of years of interest when it comes to survival
  
  n_obs <- length(arrival_times_vec)
  
  # Simulate other relevant variables like age, gender etc. 
  
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
  
  baseline_excess_vec <- exp(chi_vec) # baseline excess hazard for a time point in a given band
  
  interim_censoring_rate <- 0.000275
  
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
    t_grid = t_grid,
    n_obs = n_obs
  )
  
}

n_cores <- parallel::detectCores() - 2

init_cluster <- parallel::makeCluster(n_cores)

doParallel::registerDoParallel(cl = init_cluster)

n_iterations <- 1000
gamma_vec <- c(-0.002, 0.002, 0.005) # test out for different values of gamma
num_of_gamma_values <- length(gamma_vec)

doRNG::registerDoRNG(seed = 42)

system.time(
  
  max_cusum_matrix <- foreach::foreach(
    i = 1:n_iterations,
    .combine = "rbind"
  ) %dopar% {
    
    store_max_cusum <- numeric(num_of_gamma_values)
    
    # Simulate population survival times and population censoring indicator
    
    sim_obs_run <- sim_obs_func()
    
    u_pop_vec <- runif(sim_obs_run$n_obs)
    
    matrix_time_pop_sim <- CUSUMrelsurv::vec_pop_sim(
      age_vec = sim_obs_run$age_vec,
      gender_vec = sim_obs_run$gender_vec,
      start_year = sim_obs_run$start_year,
      end_year = sim_obs_run$end_year,
      table_end_year = 2020,
      arrival_time_vec = sim_obs_run$arrival_time_vec,
      pop_matrix_male = new_pop_data_male,
      pop_matrix_female = new_pop_data_female,
      u = u_pop_vec
    )
    
    time_pop_sim <- matrix_time_pop_sim[, 1]
    delta_p_vec <- matrix_time_pop_sim[, 2]
    
    # First 5 year - in control (from 2010 - 2015), from 2015 and beyond we have out of control
    
    eta <- 5
    
    for(j in 1:num_of_gamma_values) {
      
      u_excess_vec <- runif(sim_obs_run$n_obs)
      
      vec_time_excess_sim <- CUSUMrelsurv::vec_excess_sim_add_excess_non_neg_piecewise_in_and_out(
        partition_t_vec = sim_obs_run$partition_t_vec,
        baseline_vec = sim_obs_run$baseline_vec,
        max_follow_up_vec = sim_obs_run$max_follow_up_vec,
        u_vec = u_excess_vec,
        arrival_time_vec = sim_obs_run$arrival_time_vec,
        linear_predictor_vec = sim_obs_run$linear_predictor_vec,
        gamma = gamma_vec[j],
        eta = eta
      )
      
      interim_censoring_rate <- 0.000275
      time_censoring_sim <- cbind(rexp(sim_obs_run$n_obs, rate = interim_censoring_rate), sim_obs_run$max_follow_up_vec)
      time_censoring_sim <- apply(time_censoring_sim, 1, FUN = min)
      
      time_observed_sim <- pmin(time_pop_sim, vec_time_excess_sim, time_censoring_sim)
      delta_i_vec <- pmax(delta_p_vec, as.numeric(vec_time_excess_sim < time_pop_sim)) * as.numeric(time_censoring_sim > time_observed_sim)
      
      # Run CUSUM chart
      
      cusum_chart <- CUSUMrelsurv::cusum_add_exact_r_t_piecewise(
        partition_t_vec = sim_obs_run$partition_t_vec,
        baseline_vec = sim_obs_run$baseline_vec,
        start_year = sim_obs_run$start_year,
        age_vec = sim_obs_run$age_vec,
        gender_vec = sim_obs_run$gender_vec,
        x_matrix = sim_obs_run$x_matrix,
        time_obs_vec = time_observed_sim,
        arrival_time_vec = sim_obs_run$arrival_time_vec,
        delta_i_vec = delta_i_vec,
        beta_vec = sim_obs_run$beta_vec,
        gamma = gamma_vec[j],
        t_grid = sim_obs_run$t_grid,
        pop_data_male = new_pop_data_male,
        pop_data_female = new_pop_data_female,
        end_year_table = 2020
      )
      
      store_max_cusum[j] <- max(cusum_chart[, 1])
      
    }
    
    store_max_cusum
    
  }
  
)

parallel::stopCluster(init_cluster)

colnames(max_cusum_matrix) <- paste("gamma_", gamma_vec, sep = "")

saveRDS(max_cusum_matrix, file = "section_3_3/sim_cusum_check_performance_eta.RDS")

# If the simulations has been performed

max_cusum_matrix <- readRDS("section_3_3/sim_cusum_check_performance_eta.RDS")

max_cusum_matrix_threshold <- readRDS("section_3_3/sim_cusum_in_control_samples_for_threshold.RDS")

false_prob_vec <- c(0.05, 0.01)

final_cusum_threshold_matrix <- apply(max_cusum_matrix_threshold, 2, function(col) quantile(col, 1 - false_prob_vec)) # consider 5% (first row) and 1% (second row) false probability 

cusum_check_combinations <- expand.grid(false_prob = false_prob_vec, gamma = gamma_vec)
cusum_check_combinations_indices <- expand.grid(false_prob_index = 1:length(false_prob_vec), gamma_index = 1:length(gamma_vec))

cusum_check_combinations <- cusum_check_combinations[order(cusum_check_combinations_indices$false_prob_index), ]
cusum_check_combinations_indices <- cusum_check_combinations_indices[order(cusum_check_combinations_indices$false_prob_index), ]

cusum_signal_df <- data.frame(
  eta = rep("eta = 5", nrow(cusum_check_combinations)),
  cusum_check_combinations,
  signal_ratio = apply(cusum_check_combinations_indices, 1, function(indices) mean(max_cusum_matrix[, indices[2]] >= final_cusum_threshold_matrix[indices[1], indices[2]])) * 100
)

saveRDS(cusum_signal_df, "section_3_3/cusum_signal_df_eta_normal.RDS")
