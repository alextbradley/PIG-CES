# normalize_simulation_data.R
#
# Normalize the simulation input and output data according to specified mean and standard deviation
#
# inputs:
# 
#   model_input               : (n_simulations x n_input) array of model input parameters
# 
#   input_normalization_mean  : (1 x n_input) array, mean of model input parameter normalization 
#
#   input_normalization_sd    : (1 x n_input) array, sd of model input parameter normalization 
#
#   model_output              : (n_simulations x n_output) array of model outputs
# 
#   output_normalization_mean : (1 x n_output) array, mean of model output normalization 
#
#   output_normalization_sd   : (1 x n_output) array, sd of model output normalization
#
# outputs:
#   
#   normalized_model_input    : (n_simulations x n_input) array of normalized model input parameters
#
#   normalized_model_output   : (n_simulations x n_output) array of normalized model output 
#
#   normalized_observations   : (1 x n_output) array of normalized observations

### normalize inputs
input_normalization_mean_array <- matrix(rep(input_normalization_mean, n_members*n_iterations), nrow = n_members*n_iterations, byrow = TRUE)
input_normalization_sd_array   <- matrix(rep(input_normalization_sd, n_members*n_iterations), nrow = n_members*n_iterations, byrow = TRUE) #arrays with mean and sd repeated row wise
normalized_model_input <- (model_input - input_normalization_mean_array)/input_normalization_sd_array

### normalize outputs
output_normalization_mean_array <- matrix(rep(output_normalization_mean, n_members*n_iterations), nrow = n_members*n_iterations, byrow = TRUE)
output_normalization_sd_array   <- matrix(rep(output_normalization_sd, n_members*n_iterations), nrow = n_members*n_iterations, byrow = TRUE) #arrays with mean and sd repeated row wise
normalized_model_output <- (model_output - output_normalization_mean_array)/output_normalization_sd_array

### normalize observations
normalized_observations <- (observations - output_normalization_mean)/output_normalization_sd

