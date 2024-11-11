# Simulation of in control performance

# We specify a true baseline excess hazard, in this case based on a Weibull baseline hazard.
# Then, we simulate a dataset from the underlying distribution from above. 
# Based on this dataset, a piecewise excess baseline hazard model is estimated, which of course does not correspond to the true shape of the baseline excess hazard that the data is generated from.
# Using the results from the estimated piecewise excess baseline hazard model, simulations are done to obtain a threshold leading to 5% false signals in a 5-year period when in control. 
# Finally, this threshold is then used to see if the ratio of false signal is similar to the specified value from above when we run the chart using the estimated piecewise model, but the observations are generated using the true Weibull now.
library(doParallel)
library(ggplot2)

rm(list = ls())
source("simple_examples_scripts/read_population_table.R")

# Function to simulate excess times based on a Weibull baseline hazard

vec_excess_sim_weibull <- function(lambda, phi, beta_vec, x_matrix, u_vec = NULL) { # lambda = shape, phi = scale
  
  n_elements <- nrow(x_matrix)
  if (is.null(u_vec)) {
    u_vec <- runif(n_elements)
  } 
  inner_prod_vec <- x_matrix %*% beta_vec
  phi_vec <- phi * exp(inner_prod_vec)
  t_vec <- (-log(1 - u_vec) / phi_vec)^(1 / lambda) # inverse transform
  
  t_vec
  
} 

base_start_year <- 2000
base_end_year <- 2010

new_pop_data_male <- pop.data.male[pop.data.male$Year >= 2000, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_male$Age[new_pop_data_male$Age == "110+"] <- 110
new_pop_data_male$Age <- as.integer(new_pop_data_male$Age)
new_pop_data_male <- as.matrix(new_pop_data_male)

new_pop_data_female <- pop.data.female[pop.data.female$Year >= 2000, ] # Only look at rows relevant for the monitoring time period? Filter the rest out?
new_pop_data_female$Age[new_pop_data_female$Age == "110+"] <- 110
new_pop_data_female$Age <- as.integer(new_pop_data_female$Age)
new_pop_data_female <- as.matrix(new_pop_data_female)

base_beta_vec <- c(
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

follow_up_interval_partition <- c(seq(0, 5, by = 1), 10)

n_sims <- 1000

proportion_signals_in_control_estimated_vs_true <- numeric(n_sims)
store_threshold <- numeric(n_sims)

false_prob_t_grid <- seq(from = 0, to = 5, by = 0.01) # want to find the threshold based on the probability of false alarm during a 5 year period

set.seed(1)

seed_values <- sample(1:100000, n_sims)

system.time(
  for (i in 1:n_sims) {
    
    set.seed(seed_values[i])
    
    # Simulate a dataset using Weibull excess baseline hazard from the period 2000 to 2010 - fit a model based on this dataset and use it as basis for simulation of thresholds
    
    base_arrival_times_vec <- cumsum(rexp(50000, rate = 3750))
    base_arrival_times_vec <- base_arrival_times_vec[base_arrival_times_vec <= 10]
    range(base_arrival_times_vec)
    
    base_max_follow_up_vec <- base_end_year - (base_start_year + base_arrival_times_vec)
    base_max_follow_up_vec[base_max_follow_up_vec >= 10] <- 10 # none of these here, but this is more general if the monitoring time period is longer than the number of years of interest when it comes to survival
    
    base_n_obs <- length(base_arrival_times_vec)
    
    # Simulate other relevant variables like age, gender etc. 
    
    base_age_sim_vec <- rnorm(base_n_obs, mean = 75, sd = 10)
    base_age_sim_vec[base_age_sim_vec <= 50] <- 50
    base_age_sim_vec[base_age_sim_vec >= 105] <- 105
    
    base_gender_sim_vec <- factor(sample(1:2, base_n_obs, replace = T)) # 1 = male, 2 = female like in relsurv.
    
    base_seer_sim_vec <- sample(c("Distant", "Localised", "Regional", "Unknown"), base_n_obs, replace = T, prob = c(0.2, 0.2, 0.55, 0.05))
    base_seer_sim_vec <- factor(base_seer_sim_vec, levels = c("Distant", "Localised", "Regional", "Unknown"))
    
    base_icd_sim_vec <- factor(sample(0:3, base_n_obs, replace = T, prob = c(0.27, 0.44, 0.28, 0.01)))
    
    base_surgery_sim_vec <- factor(sample(0:2, base_n_obs, replace = T, prob = c(0.8275, 0.1715, 0.001)))
    
    base_morphology_type <- factor(sample(c("adeno-ca", "muc-adeno-ca"), base_n_obs, replace = T, prob = c(0.9, 0.1)))
    
    # Collect all the simulated variables as a data frame
    
    base_x_data <- data.frame(KJOENN = base_gender_sim_vec, SEER_STADIUM = base_seer_sim_vec, ICD7.indicator = base_icd_sim_vec, surgery.type = base_surgery_sim_vec, morphology_type = base_morphology_type)
    
    base_x_matrix <- model.matrix( ~., data = base_x_data)[, -1] # need a covariate matrix for the CUSUM calculation 
    
    # Prepare quantities for the excess hazard
    
    # Simulate population survival times and population censoring indicator
    
    base_u_pop_vec <- runif(base_n_obs)
    
    base_matrix_time_pop_sim <- CUSUMrelsurv::vec_pop_sim(
      age_vec = base_age_sim_vec,
      gender_vec = base_gender_sim_vec,
      start_year = base_start_year,
      end_year = base_end_year,
      table_end_year = 2020,
      arrival_time_vec = base_arrival_times_vec,
      pop_matrix_male = new_pop_data_male,
      pop_matrix_female = new_pop_data_female,
      u = base_u_pop_vec
    )
    
    base_time_pop_sim <- base_matrix_time_pop_sim[, 1]
    base_delta_p_vec <- base_matrix_time_pop_sim[, 2]
    
    # Simulate excess times from the true model
    
    base_u_excess_vec <- runif(base_n_obs)
    
    base_vec_time_excess_sim <- vec_excess_sim_weibull(
      lambda = 0.65,
      phi = 0.25,
      beta_vec = base_beta_vec,
      x_matrix = base_x_matrix,
      u_vec = base_u_excess_vec
    )
    
    base_interim_censoring_rate <- 0.000275
    base_time_censoring_sim <- cbind(rexp(base_n_obs, rate = base_interim_censoring_rate), base_max_follow_up_vec)
    base_time_censoring_sim <- apply(base_time_censoring_sim, 1, FUN = min)
    
    # base_time_observed_sim <- round(pmin(base_time_pop_sim, base_vec_time_excess_sim, base_time_censoring_sim) * 365.241) / 365.241
    base_time_observed_sim <- pmin(base_time_pop_sim, base_vec_time_excess_sim, base_time_censoring_sim)
    base_delta_i_vec <- pmax(base_delta_p_vec, as.numeric(base_vec_time_excess_sim < base_time_pop_sim)) * as.numeric(base_time_censoring_sim > base_time_observed_sim)
    base_censoring_indices <- which(base_delta_i_vec == 0 & base_time_observed_sim < base_max_follow_up_vec) # consider interim censoring, not censoring due to no event at the end of follow up
    base_estimated_censoring_rate <- length(base_censoring_indices) / sum(base_time_observed_sim) # MLE of exponential distributed event times with censoring data
    
    base_time_observed_sim <- round(base_time_observed_sim * 365.241) / 365.241
    
    base_data <- cbind(base_x_data, diagnosis_year = base_arrival_times_vec + base_start_year, age = base_age_sim_vec, time = base_time_observed_sim, status = base_delta_i_vec)
    
    # Fit a piecewise excess baseline model based on this dataset 
    
    base_piecewise_excess_model <-
      rsadd(
        Surv(time * 365.241, status) ~ as.factor(KJOENN) +
          as.factor(SEER_STADIUM)
        + as.factor(ICD7.indicator) + as.factor(surgery.type) +
          as.factor(morphology_type),
        data = base_data,
        ratetable = nortab,
        method = "glm.poi",
        int = follow_up_interval_partition,
        rmap = list(
          age = age * 365.241,
          sex = as.integer(KJOENN),
          year = (diagnosis_year - 1970) * 365.241
        )
      )
    
    base_estimated_beta_vec <- base_piecewise_excess_model$coefficients[1:length(base_beta_vec)]
    
    base_estimated_baseline_vec <- exp(base_piecewise_excess_model$coefficients[-(1:length(base_beta_vec))])
    
    # Bootstrap the base simulated baseline data from above to calculate the threshold
    
    max_cusum_matrix <- CUSUMrelsurv::sim_bootstrap_cusum_prop_threshold_piecewise_function_parallel(
      partition_t_vec = follow_up_interval_partition,
      baseline_vec = base_estimated_baseline_vec,
      n_iterations = 1000,
      n_cores = 8, 
      start_year = 2010,
      estimated_rate_censoring = base_estimated_censoring_rate,
      age_vec = base_age_sim_vec,
      gender_vec = base_gender_sim_vec,
      x_matrix = base_x_matrix,
      arrival_time_vec = base_arrival_times_vec,
      beta_vec = base_estimated_beta_vec,
      rho_vec = c(0.9),
      t_grid = false_prob_t_grid,
      pop_data_male = new_pop_data_male,
      pop_data_female = new_pop_data_female,
      end_year_table = 2020,
      random_state = 1
    )
    
    threshold <- quantile(max_cusum_matrix, probs = c(0.95)) # threshold to get 5% false probability under in control using estimated model as the true underlying generating excess hazard
    store_threshold[i] <- threshold
    
    # Now simulate under true model with Weibull baseline and run the CUSUM chart based on the estimated quantities from the estimated piecewise baseline excess hazard model
    
    init_cluster <- parallel::makeCluster(8)
    
    doParallel::registerDoParallel(cl = init_cluster)
    
    print(foreach::getDoParRegistered())
    
    rho_vec <- 0.9
    
    n_iterations <- 1000
    
    n_elements_rho_vec <- length(rho_vec)
    
    doRNG::registerDoRNG(seed = 42)
    
    store_matrix <- foreach::foreach(
      i = 1:n_iterations,
      .combine = "rbind"
    ) %dopar% {
      
      store_vec <- numeric(n_elements_rho_vec)
      
      start_year <- 2010
      
      num_of_years <- false_prob_t_grid[length(false_prob_t_grid)]
      
      end_year <- start_year + num_of_years
      
      arrival_times_vec <- cumsum(rexp(37500, rate = 3750)) 
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
      
      # Simulate excess times from the true model with Weibull baseline
      
      u_excess_vec <- runif(n_obs)
      
      vec_time_excess_sim <- vec_excess_sim_weibull(
        lambda = 0.65,
        phi = 0.25,
        beta_vec = base_beta_vec,
        x_matrix = x_matrix,
        u_vec = u_excess_vec
      )
      
      interim_censoring_rate <- 0.000275
      time_censoring_sim <- cbind(rexp(n_obs, rate = interim_censoring_rate), max_follow_up_vec)
      time_censoring_sim <- apply(time_censoring_sim, 1, FUN = min)
      
      time_observed_sim <- pmin(time_pop_sim, vec_time_excess_sim, time_censoring_sim)
      delta_i_vec <- pmax(delta_p_vec, as.numeric(vec_time_excess_sim < time_pop_sim)) * as.numeric(time_censoring_sim > time_observed_sim)
      
      for (j in 1:n_elements_rho_vec) {
        cusum_chart <- CUSUMrelsurv::cusum_prop_r_t_piecewise(
          partition_t_vec = follow_up_interval_partition,
          baseline_vec = base_estimated_baseline_vec,
          start_year = start_year,
          age_vec = age_sim_vec,
          gender_vec = gender_sim_vec,
          x_matrix = x_matrix,
          time_obs_vec = time_observed_sim,
          arrival_time_vec = arrival_times_vec,
          delta_i_vec = delta_i_vec,
          beta_vec = base_estimated_beta_vec,
          rho = rho_vec[j],
          t_grid = false_prob_t_grid,
          pop_data_male = new_pop_data_male,
          pop_data_female = new_pop_data_female,
          end_year_table = 2020
        )
        
        store_vec[j] <- max(cusum_chart[, 1])
      }
      
      store_vec
      
    }
    
    parallel::stopCluster(init_cluster)
    
    proportion_signals_in_control_estimated_vs_true[i] <- mean(store_matrix >= threshold) # check the proportion of signals given the estimated threshold when in control
    
    print(i)
    
  }
  
)

list_weibull_in_control_piecewise_bootstrap <- list(threshold_values = store_threshold, proportion_in_control_signal = proportion_signals_in_control_estimated_vs_true)
saveRDS(list_weibull_in_control_piecewise_bootstrap, "section_3_5/weibull_baseline/parametric_piecewise_constant_model/list_weibull_in_control_piecewise_bootstrap_baseline_period.RDS")

# If the simulations have been performed

list_weibull_in_control_piecewise_bootstrap <- readRDS("section_3_5/weibull_baseline/parametric_piecewise_constant_model/list_weibull_in_control_piecewise_bootstrap_baseline_period.RDS")

df_weibull_in_control_piecewise_bootstrap <- data.frame(proportion_in_control_signal = list_weibull_in_control_piecewise_bootstrap$proportion_in_control_signal)

hist(df_weibull_in_control_piecewise_bootstrap$proportion_in_control_signal, breaks = seq(from = 0, to = 0.4, by = 0.01), probability = T)

ggplot(df_weibull_in_control_piecewise_bootstrap, aes(x = proportion_in_control_signal, y = after_stat(density))) + 
  geom_histogram(binwidth = 0.01, fill="#69b3a2", color="#e9ecef", alpha=1, center = 0.005) + 
  geom_vline(aes(xintercept = mean(proportion_in_control_signal)), color = "red", linewidth = 1) + 
  geom_vline(aes(xintercept = median(proportion_in_control_signal)), color = "blue", linewidth = 1) +
  labs(
    title = "Histogram of signal ratio under in-control (piecewise baseline excess hazard model and bootstrap)",
    x = "Signal ratio",
    y = "Density"
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5)
  )
