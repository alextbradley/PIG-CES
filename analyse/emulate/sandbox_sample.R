library(MASS)

model_input_sd <- apply(model_input, 2, sd)
model_input_mean <- apply(model_input, 2, mean)
model_output_sd <- apply(model_output, 2, sd)
model_output_mean <- apply(model_output, 2, mean)

#get the initial guess as the mean of the final iteration

indices <- which(meta_data[,1] == max(meta_data[,1])); #indices corresponding to final iteration
final_iteration_parameters <- normalized_model_input[indices,]
theta0 <- apply(final_iteration_parameters, 2, mean) #initial guess is mean of parameters in final iteration
dimensional_theta0 <- model_input_mean + model_input_sd * theta0 #dimensional version of the model input (we only compute this for curiosity)


dimensional_prior_mean <- c(1.0, 1.0, 1.0, 0.0, 0.0, 200.0, 5.0) #dimnesional prior means
prior_mean <- (dimensional_prior_mean - model_input_mean)/model_input_sd

theta0_output <- predict(model_2015grv,t(theta0)) #model output at the mean

# in the languange of Cleary et al (2020):
# J: number of ensemble members
# N: number of iterations
# M: size of the subset of model pairs (M <= JN, Cleary et al suggest using M = J using final iteration only)
# p: number of input parameters (here, 7)

#compute C, which is based on the final iteration of the EKS
C <- array(0, dim = c(7, 7))

for (i in members) {
  err <- final_iteration_parameters[i,] - mean(final_iteration_parameters[i,])
  C <- C + matrix(kronecker(err, err), 7, 7)/length(members)
}

# set up corresponding multivariable normal
mu <- c(0,0,0,0,0,0,0)
update <- mvrnorm(1, mu, C)

theta_trial <- theta0 + update
theta_current <- theta0
theta_values <- c(theta0)

actual_dimensional_output <- 0 #dimensional value of actual (set to be relative to the truth)
normalized_actual <- (actual_dimensional_output - model_output_mean[3])/model_output_sd[3]

actual_dimensional_error_cov <- 1
normalized_error_cov <- actual_dimensional_error_cov/model_output_sd[3]

likelihood_trial = exp(-(normalized_actual - predict(model_2015grv,t(theta_new_trial))$mean )^2 / 2 /normalized_error_cov^2)
likelihood_prev = exp(-(normalized_actual - predict(model_2015grv,t(theta_current))$mean )^2 / 2 /normalized_error_cov^2)

# define prior stuff
#input_headers <- c("weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",
#                   "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration") #just to remember order

dimensional_prior_mean <- c(1.0, 1.0, 1.0, 0.0,0.0 ,200.0, 5.0);
dimensional_prior_covariance <- diag(c(0.3, 0.3, 0.3, 1.2, 200.0, 100.0, 2.5))
normalized_prior_mean <- (dimensional_prior_mean - model_input_mean)/model_input_sd;
normalized_prior_covariance  <- dimensional_prior_covariance/model_input_sd
inv_normalized_prior_covariance <- solve(normalized_prior_covariance)

prior_diff_trial <-  matrix(theta_new_trial - normalized_prior_mean, 1, 7) 
prior_diff_prev  <-  matrix(theta_current - normalized_prior_mean, 1, 7) 

prior_term_trial <- exp(- prior_diff_trial %*% inv_normalized_prior_covariance %*% t(prior_diff_trial) / 2);
prior_term_prev  <- exp(- prior_diff_prev  %*% inv_normalized_prior_covariance %*% t(prior_diff_prev) / 2);

a <- min(c(1, (likelihood_trial*prior_term_trial / (likelihood_prev*prior_term_prev)) ));

# compute a random number and check whether to accept/reject
random_number <- runif(1, min = 0, max = 1)

if (random_number < a){ #accept
  theta_current <- theta_trial
}

# store the new value
theta_values <- rbind(theta_values, theta_current)

# loop over lots of values
N_steps = 50000

for (i in 1:N_steps) {

  #determine the trial theta_new
  update <- mvrnorm(1, mu, C)
  theta_trial <- theta_current + update
  
  #compute emulator likelihoods
  likelihood_term_trial = exp(-(normalized_actual - predict(model_2015grv,t(theta_trial))$mean )^2 / 2 /normalized_error_cov^2)
  likelihood_term_current = exp(-(normalized_actual - predict(model_2015grv,t(theta_current))$mean )^2 / 2 /normalized_error_cov^2)
  
  #compute prior terms
  prior_diff_trial    <-  matrix(theta_trial - normalized_prior_mean, 1, 7) 
  prior_diff_current  <-  matrix(theta_current - normalized_prior_mean, 1, 7) 
  prior_term_trial    <- exp(- prior_diff_trial    %*% inv_normalized_prior_covariance %*% t(prior_diff_trial) / 2);
  prior_term_current  <- exp(- prior_diff_current  %*% inv_normalized_prior_covariance %*% t(prior_diff_current) / 2);
  
  #update ratio
  a <- min(c(1, (likelihood_term_trial*prior_term_trial / (likelihood_term_current*prior_term_current)) ));
  
  #accept/reject
  if (random_number < a){ #accept
  theta_current <- theta_trial
  }
  
  # store the new value
  theta_values <- rbind(theta_values, theta_current)
  
  #progress tracker
  if (i%%1000 == 0)
  {
    print(i)
  }
  
}

# plot the results
#data_df <- as.data.frame(theta_values*model_input_sd + model_input_mean) #rescaled data
data_df <- as.data.frame(theta_values) #rescaled data
data_df$Index <- 1:dim(theta_values)[1]


# plot histograms
#plot_list <- list()
#for (i in 1:7) {
#  p <- ggplot(data_df, aes_string(x = colnames(data_df)[i])) +
#    geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
#    theme_minimal() +
#    ggtitle(input_headers[i]) + 
#    xlim(-2.5, 2.5) # Set y-axis limits
  
#  plot_list[[i]] <- p
#}
#grid.arrange(grobs = plot_list, ncol = 3)

# make dimensional again
dimensional_data_df = data_df
for (i in 1:7){
  dimensional_data_df[,i] <- (data_df[,i]*model_input_sd[i] + model_input_mean[i])
}

# plot dimensional histogram
if (1){
plot_list <- list()
nbars < 30; #number of histogram bars

#priors
mu <- c(1,1,1,0,0,200,5)
sigma <- c(0.3, 0.3, 0.3, 1.2, 200, 100, 2.5)

#naive posteriors
mu_naivepost <- apply(model_input[indices,], 2, mean)
sigma_naivepost <- apply(model_input[indices,], 2, sd)

for (i in 1:7) {
  
  #get the prior and put in a data frame
  xx <- seq(from = min( dimensional_data_df[,i]), to = max( dimensional_data_df[,i]), length.out = 100)
  yy <- 1/sqrt(2*pi*sigma[i]^2) * exp(-(xx - mu[i])^2 / 2 /sigma[i]^2)
  data_prior <- data.frame(
    x = xx,
    y = yy
  )
  
  data_naivepost <- data.frame(
    x = xx,
    y =  1/sqrt(2*pi*sigma_naivepost[i]^2) * exp(-(xx - mu_naivepost[i])^2 / 2 /sigma_naivepost[i]^2)
  )
  
  #add the naive posterior based on final iteration
  
  
  #make plot
  p <- ggplot(dimensional_data_df, aes_string(x = colnames(data_df)[i])) +
    geom_histogram(aes(y = after_stat(density)), binwidth = (max( dimensional_data_df[,i]) - min( dimensional_data_df[,i]))/nbars, fill = "blue", alpha = 0.7) +
    geom_density(alpha = 1, color = "black", linewidth =0.75) +
    #geom_line(aes(x = xx, y = yy), color = "black", linetype = "dashed") +
    theme_minimal() +
    ggtitle(input_headers[i]) 
  p <- p + geom_line(data = data_prior, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.75) 
  p <- p + geom_line(data = data_naivepost, aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 0.75)
  
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 3)
}



plot_list <- list()
ylims <- matrix(c(0,0,0,-1, -200, 0,0,2,2.5,2, 2, 1000,600, 15), nrow = 7, ncol = 2)
for (i in 1:7) {
  p <- ggplot(dimensional_data_df, aes_string(x = "Index", y = paste0("V", i))) +
    geom_line() +
    theme_minimal() + 
    theme(plot.title = element_blank(),axis.title.x = element_blank()) + # Remove plot titles
    ylim(ylims[i,]) # Set y-axis limits
  
  plot_list[[i]] <- p
}
P <- grid.arrange(grobs = plot_list, ncol = 1)

#


##### get a kernel density estimate
kde <- density(dimensional_data_df[,5])
#plot(kde, main="Kernel Density Estimate", xlab="Value", ylab="Density")


