# run_mcmc.R
#
# Run the mcmc.
#
# inputs: init_sample (dimensional), final_iteration_parameters, n_steps, n_burn, fac
# outputs: 
# plots:

nvar <- length(init_sample)
fac  <- 10   #how big is the factor relating model and obs error? 0 implies no model error

#normalize input stuff
normalized_init_sample                <- (init_sample - input_normalization_mean)/input_normalization_sd
normalized_final_iteration_parameters <- (final_iteration_parameters - matrix(rep(input_normalization_mean, times = dim(final_iteration_parameters)[1]), nrow = dim(final_iteration_parameters)[1], ncol = nvar, byrow = TRUE))/matrix(rep(input_normalization_sd, times = dim(final_iteration_parameters)[1]), nrow = dim(final_iteration_parameters)[1], ncol = nvar, byrow = TRUE)
normalized_prior_mean                 <- (prior_mean -  input_normalization_mean)/input_normalization_sd
normalized_prior_sd                   <- prior_sd/input_normalization_sd
normalized_prior_cov                  <- normalized_prior_sd^2
normalized_prior_covariance           <- diag(normalized_prior_sd)
inv_normalized_prior_covariance       <- diag(1 / normalized_prior_cov) #inverse of the prior covariance matrix, NB only works bc assume no correlation between inputs

#normalized output stuff
normalized_observations <- (observations - output_normalization_mean)/output_normalization_sd
normalized_observation_error <- observations_noise/output_normalization_sd


#
# get the step size
#
C <- array(0, dim = c(nvar, nvar))
for (i in n_members) {
  err <- normalized_final_iteration_parameters[i,] - mean(normalized_final_iteration_parameters[i,])
  C <- C + matrix(kronecker(err, err), nvar, nvar)/n_members
}

#
# seed the distribution
#
normalized_current <- normalized_init_sample
normalized_samples <- matrix(0, nrow = n_steps, ncol = nvar)

accept <- c()
ice_loss <- c()


#
# loop over steps
#
for (i in 1:n_steps){
  
  # get the step size
  step_size <- mvrnorm(1, mu = rep(0, nvar), Sigma = C)
  
  # update the trial value
  normalized_trial <- normalized_current + step_size
  
  # get the emulator predictions for each
  model.predict<-predict(model, array(unname(normalized_trial), dim = c(1, length(normalized_trial)))) #we have to unname because the emulator is trained on unnamed data
  normalized_trial_emulator_value <- model.predict$mean
  normalized_trial_emulator_sd    <- model.predict$mean
  
  
  model.predict<-predict(model, array(unname(normalized_current), dim = c(1, length(normalized_current)))) #we have to unname because the emulator is trained on unnamed data
  normalized_current_emulator_value <- model.predict$mean
  normalized_current_emulator_sd   <- model.predict$mean
  
  
  #work out the terms in the log-likelihood (equation 2.14 in Cleary et al.)
  phi_gp_trial <- 1/2 * (normalized_observations - normalized_trial_emulator_value)^2 / (normalized_observation_error^2 + fac*normalized_observation_error^2 + normalized_trial_emulator_sd^2) + 
                1/2 * log(normalized_observation_error^2 + fac*normalized_observation_error^2 + normalized_trial_emulator_sd^2)
  
  phi_gp_current <- 1/2 * (normalized_observations - normalized_current_emulator_value)^2 / (normalized_observation_error^2 + fac*normalized_observation_error^2 + normalized_current_emulator_sd^2) + 
                1/2 * log(normalized_observation_error^2 + fac*normalized_observation_error^2 + normalized_current_emulator_sd^2)
  
  
  #work out the prior terms in the log-likelihood (second term in the round brackets in (2.15) of Cleary et al.)
  prior_term_trial   <- 1/2 * array((normalized_trial   - normalized_prior_mean), dim = c(1,nvar)) %*% inv_normalized_prior_covariance %*% array((normalized_trial   - normalized_prior_mean), dim = c(nvar,1))
  prior_term_current <- 1/2 * array((normalized_current - normalized_prior_mean), dim = c(1,nvar)) %*% inv_normalized_prior_covariance %*% array((normalized_current - normalized_prior_mean), dim = c(nvar,1))
  
  # compute the acceptance probability (equation 2.15 in Cleary et al.)
  a = min(1, exp(-(phi_gp_trial + prior_term_trial) + (phi_gp_current + prior_term_current)))
  
  # generate random number and accept/reject
  u <- runif(1)
  if (u < a) {
    normalized_current <- normalized_trial
    accept <- c(accept, 1)
  } else{
    accept <- c(accept, 0)
  }
  
  normalized_samples[i,] <- normalized_current
  
  #store the model output
  model.predict<-predict(model, array(unname(normalized_current), dim = c(1, length(normalized_current)))) #we have to unname because the emulator is trained on unnamed data
  ice_loss <- c(ice_loss, model.predict$mean)
  
  
  #progress tracker
  if (i%%1000 == 0)
  {
    cat("Completed", i, "steps of the MCMC\n")
    
  } 
  
}


# print the acceptance ratio to console
print(paste('Acceptance ratio is', sum(accept)/length(accept)))

#remove burn-in period
normalized_posterior_samples <- normalized_samples[n_burn:dim(normalized_samples)[1],]

#store in a data frame
normalized_posterior_samples <- as.data.frame(normalized_posterior_samples) #make into dataframe
normalized_posterior_samples$Index <- 1:dim(normalized_posterior_samples)[1] #set an index

#rescale to recover dimensional data
posterior_samples <- normalized_posterior_samples
for (i in 1:nvar){
  posterior_samples[,i] <- (normalized_posterior_samples[,i]*input_normalization_sd[i] + input_normalization_mean[i])
}
colnames(posterior_samples) <- c(input_colnames, "Index")

# make a plot showing the trace of the mcmc
plot_list <- list()

ylims <- matrix(c(0,  0,-1, -200, 0. ,0,
                  2,2, 2, 1000,600, 15), nrow = nvar, ncol = 2) #ylimits of plot

for (i in 1:nvar) {
  p <- ggplot(posterior_samples, aes_string(x = "Index", y = input_colnames[i])) +
    geom_line() +
    theme_minimal() + 
    theme(plot.title = element_blank(),axis.title.x = element_blank()) + # Remove plot titles
    ylim(ylims[i,]) # Set y-axis limits
  
  plot_list[[i]] <- p
}
#p <- grid.arrange(grobs = plot_list, ncol = 1)
#plot_list

# make a plot showing the posterior histograms alongside prior, kde, and naive posterior
plot_list <- list()
nbars <- 100; #number of histogram bars

#naive posterior
naive_posterior_sd   = apply(final_iteration_parameters, 2, sd)
naive_posterior_mean = apply(final_iteration_parameters, 2, mean)

for (i in 1:nvar) {
  
  #get the prior and put in a data frame
  #xx <- seq(from = min(posterior_samples[,i]), to = max(posterior_samples[,i]), length.out = 100)
  xx <- seq(from = prior_mean[i] - 2 *prior_sd[i], to =  prior_mean[i] + 2 *prior_sd[i], length.out = 100)
  
  yy <- 1/sqrt(2*pi*prior_sd[i]^2) * exp(-(xx - prior_mean[i])^2 / 2 /prior_sd[i]^2)
  data_prior <- data.frame(
    x = xx,
    y = yy
  )

  
  data_naive_posterior <- data.frame(
    x = xx,
    y =  1/sqrt(2*pi*naive_posterior_sd[i]^2) * exp(-(xx - naive_posterior_mean[i])^2 / 2 /naive_posterior_sd[i]^2)
  )
  
  
  #make plot
  xlims <- matrix(c((prior_mean - 2 *prior_sd), (prior_mean + 2*prior_sd)), nrow = nvar, ncol = 2)
  
  p <- ggplot(posterior_samples, aes_string(x = colnames(posterior_samples)[i])) +
    geom_histogram(aes(y = after_stat(density)), binwidth = (max( posterior_samples[,i]) - min( posterior_samples[,i]))/nbars, fill = "blue", alpha = 0.5) +
    geom_density(alpha = 1, color = "black", linewidth =0.75) +
    #geom_line(aes(x = xx, y = yy), color = "black", linetype = "dashed") +
    theme_minimal() +
    xlim(xlims[i,])+
   # ylim(c(0, 3/sqrt(2*pi*prior_sd[i]^2))) + 
    ggtitle(input_colnames[i]) 
  
  p <- p + geom_line(data = data_prior, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.75) 
  p <- p + geom_line(data = data_naive_posterior, aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 0.75)
  
  plot_list[[i]] <- p
} #end loop over variables
grid.arrange(grobs = plot_list, ncol = 3)
