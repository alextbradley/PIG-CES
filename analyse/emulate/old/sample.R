###########################################################
########## Run the MCMC on the normalized data ############
###########################################################

#
# get initial guess from final iteration
#
indices <- which(meta_data[,1] == max(meta_data[,1])); #indices corresponding to final iteration
dimensional_final_iteration_parameters <- dimensional_model_input[indices,]
normalized_final_iteration_parameters <- normalized_model_input[indices,]
normalized_theta0 <- apply(normalized_final_iteration_parameters, 2, mean) 
normalized_naive_posterior_mean <- apply(normalized_final_iteration_parameters, 2, mean)
normalized_naive_posterior_sd <- apply(normalized_final_iteration_parameters, 2, sd)

#
# update the prior from earlier?
#
#dimensional_prior_mean <- apply(dimensional_final_iteration_parameters, 2, mean) 
#dimensional_prior_sd   <- apply(dimensional_final_iteration_parameters, 2, sd) 
#normalized_prior_mean <- (dimensional_prior_mean - dimensional_model_input_mean)/dimensional_model_input_sd 
#normalized_prior_sd   <- dimensional_prior_sd/dimensional_model_input_sd

#
#compute martix C, which is based on the final iteration of the EKS and set the update size
#
C <- array(0, dim = c(7, 7))
for (i in members) {
  err <- normalized_final_iteration_parameters[i,] - mean(normalized_final_iteration_parameters[i,])
  C <- C + matrix(kronecker(err, err), 7, 7)/length(members)
}

# 
# intialize mcmc
#
normalized_prior_covariance <- diag(normalized_prior_sd)
inv_prior_covariance <- solve(normalized_prior_covariance)
normalized_theta_values <- matrix(0, nrow = N_steps, ncol = length(normalized_theta0))
normalized_theta_values[1,] <- normalized_theta0
normalized_theta_current <- normalized_theta0

#shift the first entry away
#update_size <- mvrnorm(1, c(0,0,0,0,0,0,0), C) #draw from multivariable random normal with zero mean
#normalized_theta_current <- normalized_theta_current + 2*update_size

#
# perform steps
#
for (i in 1:(N_steps-1)) {
  
  #determine the trial theta
  update_size <- mvrnorm(1, c(0,0,0,0,0,0,0), C) #draw from multivariable random normal with zero mean
  normalized_theta_trial <- normalized_theta_current + update_size
  
  #compute likelihoods
  uncertainty_trial <- predict(model,t(normalized_theta_trial))$upper95 - predict(model,t(normalized_theta_trial))$lower95
  uncertainty_current <- predict(model,t(normalized_theta_current))$upper95 - predict(model,t(normalized_theta_current))$lower95
  
  likelihood_term_trial   = exp(-(normalized_actual - predict(model,t(normalized_theta_trial))$mean )^2 / 2 /(normalized_error_cov+uncertainty_trial)^2 - log(uncertainty_trial)/2)
  likelihood_term_current = exp(-(normalized_actual - predict(model,t(normalized_theta_current))$mean )^2 / 2 /(normalized_error_cov+uncertainty_current)^2 - log(uncertainty_current)/2)
  
  #compute prior terms
  normalized_prior_diff_trial   <-  matrix(normalized_theta_trial - normalized_prior_mean, 1, 7) 
  normalized_prior_diff_current <-  matrix(normalized_theta_current - normalized_prior_mean, 1, 7) 
  
  prior_term_trial    <- exp(- normalized_prior_diff_trial    %*% inv_prior_covariance %*% t(normalized_prior_diff_trial) / 2);
  prior_term_current  <- exp(- normalized_prior_diff_current  %*% inv_prior_covariance %*% t(normalized_prior_diff_current) / 2);
  #prior_term_trial <- 1
  #prior_term_current <- 1
  
  #likelihood ratio
  a <- min(c(1, (likelihood_term_trial*prior_term_trial / (likelihood_term_current*prior_term_current)) ));    
  
  random_number <- runif(1, min = 0, max = 1)
  if (random_number < a){ #accept
    #  cat('accepted update at step ', i)
    normalized_theta_current <- normalized_theta_trial
  }
  
  # store the new value
  # print(normalized_theta_current)
  #  Sys.sleep(0.2)
  normalized_theta_values[i,] <- normalized_theta_current
  
  #progress tracker
  if (i%%1000 == 0)
  {
    print(i)
  } 
}

#
# post process
#
normalized_mcmc_data <- as.data.frame(normalized_theta_values) #make into dataframe
normalized_mcmc_data$Index <- 1:dim(normalized_theta_values)[1] #set an index

#rescale to recover dimensional data
dimensional_mcmc_data <- normalized_mcmc_data
dimensional_naive_posterior_mean <- normalized_naive_posterior_mean
dimensional_naive_posterior_sd <- normalized_naive_posterior_sd

for (i in 1:7){
  dimensional_mcmc_data[,i] <- (normalized_mcmc_data[,i]*dimensional_model_input_sd[i] + dimensional_model_input_mean[i])
  dimensional_naive_posterior_mean[i] <- dimensional_naive_posterior_mean[i]*dimensional_model_input_sd[i] + dimensional_model_input_mean[i]
  dimensional_naive_posterior_sd[i] <- dimensional_naive_posterior_sd[i]*dimensional_model_input_sd[i] 
}

###########################################################
##################### Make MCMC plots #####################
###########################################################

#
# make mcmc traceplot
#
if (mcmc_traceplot){
  plot_list <- list()
  ylims <- matrix(c(0,0,0,-1, -200, 0,0,2,2.5,2, 2, 1000,600, 15), nrow = 7, ncol = 2)
  for (i in 1:7) {
    p <- ggplot(dimensional_mcmc_data, aes_string(x = "Index", y = paste0("V", i))) +
      geom_line() +
      theme_minimal() + 
      theme(plot.title = element_blank(),axis.title.x = element_blank()) + # Remove plot titles
      ylim(ylims[i,]) # Set y-axis limits
    
    plot_list[[i]] <- p
  }
  p <- grid.arrange(grobs = plot_list, ncol = 1)
} #end traceplot flag

#
# make distribution plots
#
if (mcmc_histogram){
  plot_list <- list()
  nbars <- 200; #number of histogram bars
  
  for (i in 1:7) {
    
    #get the prior and put in a data frame
    xx <- seq(from = min(dimensional_mcmc_data[,i]), to = max(dimensional_mcmc_data[,i]), length.out = 100)
    yy <- 1/sqrt(2*pi*dimensional_prior_sd[i]^2) * exp(-(xx - dimensional_prior_mean[i])^2 / 2 /dimensional_prior_sd[i]^2)
    data_prior <- data.frame(
      x = xx,
      y = yy
    )
    
    data_naive_posterior <- data.frame(
      x = xx,
      y =  1/sqrt(2*pi*dimensional_naive_posterior_sd[i]^2) * exp(-(xx - dimensional_naive_posterior_mean[i])^2 / 2 /dimensional_naive_posterior_sd[i]^2)
    )
    
    
    #make plot
    xlims <- matrix(c(0.5, 0.5, 0.5, 0, -200, 0,0,
                      1.25, 2.0, 2.0,  1, 600,  500, 8), nrow = 7, ncol = 2)
    
    p <- ggplot(dimensional_mcmc_data, aes_string(x = colnames(dimensional_mcmc_data)[i])) +
      geom_histogram(aes(y = after_stat(density)), binwidth = (max( dimensional_mcmc_data[,i]) - min( dimensional_mcmc_data[,i]))/nbars, fill = "blue", alpha = 0.7) +
      geom_density(alpha = 1, color = "black", linewidth =0.75) +
      #geom_line(aes(x = xx, y = yy), color = "black", linetype = "dashed") +
      theme_minimal() +
      xlim(xlims[i,])+
      ggtitle(input_headers[i]) 
    p <- p + geom_line(data = data_prior, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.75) 
    p <- p + geom_line(data = data_naive_posterior, aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 0.75)
    
    plot_list[[i]] <- p
  } #end loop over variables
  grid.arrange(grobs = plot_list, ncol = 3)
}

###########################################################
## plot emulator error as a function of model parameters ##
##########################################################
# plot the emulated model output as a function of input parameter for each of the parameters, at the naive posterior mean
np <- 1000 #number of x points
dimensionless_prior_mean <- normalized_prior_mean
dimensionless_prior_sd <- normalized_prior_sd
normalized_xx <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
dimensional_xx <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_mean <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_u95 <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_l95 <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))


for (i in 1:7){
  normalized_xx[,i] <- seq(from = (dimensionless_prior_mean[i] - 3*dimensionless_prior_sd[i]), to = (dimensionless_prior_mean[i] + 3*dimensionless_prior_sd[i]) , length.out = np)
  dimensional_xx[,i] <- (normalized_xx[,i]*dimensional_model_input_sd[i] + dimensional_model_input_mean[i])
  test_values <- matrix(rep(normalized_naive_posterior_mean, times = np), nrow = np, ncol = 7, byrow = TRUE) #test values at naive posterior mean everywhere
  test_values[,i] <- normalized_xx[,i] #replace the ith column with normalized xx values
  
  #test the model on these values
  tmp <- predict(model,test_values)
  
  normalized_yy_mean[,i] <- tmp$mean
  normalized_yy_u95[,i] <- tmp$upper95
  normalized_yy_l95[,i] <- tmp$lower95
}

normalized_xx <- data.frame(normalized_xx)
dimensional_xx <- data.frame(dimensional_xx)
normalized_yy_mean <- data.frame(normalized_yy_mean)
normalized_yy_l95 <- data.frame(normalized_yy_l95)
normalized_yy_u95 <- data.frame(normalized_yy_u95)

plot_list <- list()

for (i in 1:7) {
  # Create a data frame for the current column
  df <- data.frame(
    x = dimensional_xx[[i]],
    y = normalized_yy_mean[[i]] - normalized_actual,
    ymin = normalized_yy_l95[[i]] - normalized_actual,
    ymax = normalized_yy_u95[[i]] - normalized_actual
  )
  
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = dimensional_naive_posterior_mean[i], linetype = "dashed") +
    
    ylim(-.1, .1) +
    xlim(xlims[i,]) + 
    labs(x = input_headers[i], y = "emu error") +
    theme_minimal()
  
  # Add the plot to the list
  plot_list[[i]] <- p
}

grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)


if (output_mcmc){
  dimensional_mcmc_data_noindex <- dimensional_mcmc_data[,1:7];
  file_name <- paste0("mcmc_output_realization", padded_realization, ".csv")
  colnames(dimensional_mcmc_data_noindex) <- input_headers
  
  write.csv(dimensional_mcmc_data_noindex, file = file_name, row.names = FALSE)
  
}
#return(list(dimensional_mcmc_data = dimensional_mcmc_data, model = model))

#}

#emulate_sample_data <- emulate_sample(realization,iterations,members,verbose,
#                                       training_iterations,method,nugget_est,kernel_type,max_eval,alpha,
#                                       N_steps,
#                                       leave_one_out_validation_plot,mcmc_traceplot,mcmc_histogram,
#                                       dimensional_prior_mean,dimensional_prior_sd,input_headers,
#                                       dimensional_observation,dimensional_error_cov)