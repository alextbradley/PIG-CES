##### Script to emulate model output and sample using MCMC


#### libraries
library(RobustGaSP)
library(ggplot2)
library(RColorBrewer)
library(dplyr) #for filtering data frames
library(MASS)

#### include the functions
source("./emulate_sample_functions.R")

#### specify model results
realization <- 30
iterations  <- 1:5
members     <- 1:20
verbose     <- 0

#### emulator parameters
training_iterations <- 4:5 #which iterations to train on
method              <- 'post_mode'
nugget_est          <- F
kernel_type         <- 'matern_5_2'
max_eval            <- 100
alpha               <- NA  #alpha value for exponential kernels 

#### mcmc parameters
N_steps <- 50000

#### specify plots
leave_one_out_validation_plot <- 1
mcmc_traceplot                <- 1
mcmc_histogram                <- 1

#### prior info 
dimensional_prior_mean <- c(1.0, 1.0, 1.0, 0.0,0.0 ,200.0, 5.0);
dimensional_prior_covariance <- diag(c(0.3, 0.3, 0.3, 1.2, 200.0, 100.0, 2.5)) #give large for agnositicism
input_headers <- c("weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",
                   "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration")

#### obs info
actual_dimensional <- 0 #dimensional value of actual (set to be relative to the truth)
actual_dimensional_error_cov <- 1

#### run the thing
emulate_sample(realization,
               iterations,
               members, 
               training_iterations, 
               method, 
               nugget_est, 
               kernel_type, 
               max_eval, 
               alpha, 
               leave_one_out_validation_plot,
               mcmc_traceplot,
               mcmc_histogram,
               N_steps,
               dimensional_prior_mean,
               dimensional_prior_covariance,
               input_headers,
               actual_dimensional,
               actual_dimensional_error_cov)


emulate_sample <- function(realization,
                            iterations,
                            members, 
                            training_iterations, 
                            method, 
                            nugget_est, 
                            kernel_type, 
                            max_eval, 
                            alpha, 
                            leave_one_out_validation_plot,
                            mcmc_traceplot,
                            mcmc_histogram,
                            N_steps,
                            dimensional_prior_mean,
                            dimensional_prior_covariance,
                            input_headers,
                            actual_dimensional,
                           actual_dimensional_error_cov){
    
    #####
    ##### get the model results and restrict to only grounding line volume and to training data
    #####
    data <- get_data(realization,iterations,members, verbose)
    indices <- which(data$meta_data[,1] %in% training_iterations)
    model_input <- data$model_input
    model_input <- model_input[indices,]
    model_output <- data$model_output
    model_output <- model_output[indices,3]
    model_output <- matrix(model_output, nrow = length(model_output), ncol = 1 ) #make into a matrix
    meta_data <- data$meta_data
    meta_data <- meta_data[indices,]
    

    #####
    ##### normalize the data
    ##### 
    normalize_z_score <- function(x) {
      (x - mean(x)) / sd(x)
    }
    model_input_sd <- apply(model_input, 2, sd)
    model_input_mean <- apply(model_input, 2, mean)
    model_output_sd <- apply(model_output, 2, sd)
    model_output_mean <- apply(model_output, 2, mean)
    normalized_model_output <- apply(model_output, 2, normalize_z_score)
    normalized_model_input  <- apply(model_input, 2, normalize_z_score)
    normalized_actual <- (actual_dimensional_output - model_output_mean)/model_output_sd
    normalized_error_cov <- actual_dimensional_error_cov/model_output_sd
    
    

    #####
    ##### perform the leave one out validation and make plots
    #####
    if (leave_one_out_validation_plot){

        leave_one_out_data <- leave_one_out_validation(normalized_model_input, normalized_model_output, method, nuggest_est, kernel_type, max_eval,alpha)

        #print the RMSE
        rmse_fn <- function(actual, predicted) {
            sqrt(mean((actual - predicted)^2))
            }
        rmse <- rmse_fn(leave_one_out_data$y, leave_one_out_data$x)
        
        print(paste("normalised RMSE:", rmse))
        print(paste("actual RMSE:", (rmse * sd(model_output,3))))

        P1 <- plot_modelled_vs_emulated(leave_one_out_data)
        P1 <- P1 + xlim(-2, 2) + ylim(-2,2)
        print(P1)

        P2 <- plot_normalized_errors(leave_one_out_data)
        P2 <- P2 + xlim(-4, 4) 
        print(P2)

    }
    
    ##### 
    ##### create the emulator
    #####
    model <- rgasp(design = normalized_model_input, response = normalized_model_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)
    
    #####
    ##### initialize the mcmc
    #####
    indices <- which(meta_data[,1] == max(meta_data[,1])); #indices corresponding to final iteration
    final_iteration_parameters <- normalized_model_input[indices,]
    naive_posterior_mean <- apply(final_iteration_parameters, 2, mean)
    naive_posterior_sd <- apply(final_iteration_parameters, 2, sd)
    theta0 <- apply(final_iteration_parameters, 2, mean) 
    
    #compute C, which is based on the final iteration of the EKS
    C <- array(0, dim = c(7, 7))
    for (i in members) {
      err <- final_iteration_parameters[i,] - mean(final_iteration_parameters[i,])
      C <- C + matrix(kronecker(err, err), 7, 7)/length(members)
    }
    
    ##### do the mcmc
    normalized_prior_mean  <- (dimensional_prior_mean - model_input_mean)/model_input_sd;
    normalized_prior_covariance  <- dimensional_prior_covariance/model_input_sd
    mcmc_data <- run_mcmc(theta0, N_steps, C, model, normalized_prior_mean, normalized_prior_covariance, normalized_actual, normalized_error_cov)
  
    ##### get the dimensional data
    dimensional_mcmc_data <- mcmc_data
    dimensional_naive_posterior_mean <- naive_posterior_mean
    dimensional_naive_posterior_sd <- naive_posterior_sd
    
    for (i in 1:7){
      dimensional_mcmc_data[,i] <- (mcmc_data[,i]*model_input_sd[i] + model_input_mean[i])
      dimensional_naive_posterior_mean[i] <- naive_posterior_mean[i] + model_input_mean[i]
      dimensional_naive_posterior_sd[i] <- naive_posterior_sd[i]*model_input_sd[i]
      
    }
    
    #### make some plots
    if (mcmc_traceplot){
      P <- make_mcmc_traceplot(dimensional_mcmc_data)
    }
      
    if (mcmc_histogram){
      dimensional_prior_sd = diag(dimensional_prior_covariance)
      P <- make_mcmc_histogram(dimensional_mcmc_data,dimensional_prior_mean, dimensional_prior_sd,dimensional_naive_posterior_mean,dimensional_naive_posterior_sd, input_headers)
      
    }
      
    return(list(mcmc_data = mcmc_data, dimensional_mcmc_data = dimensional_mcmc_data))
}

