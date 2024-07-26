##### Script to emulate model output and sample using MCMC

###########################################################
##################### Preliminaries #######################
###########################################################

rm(list = ls())

#
# libraries
#
library(RobustGaSP)
library(ggplot2)
library(RColorBrewer)
library(dplyr) #for filtering data frames
library(MASS)
library(gridExtra)


#
# make sure you run from the right place
#
setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

#
# specify model results to use
#
realization <- 26
iterations  <- 1:5
members     <- 1:20
verbose     <- 0

#
# emulator parameters
#
training_iterations <- 1:5 #which iterations to train on
method              <- 'post_mode'
nugget_est          <- F
kernel_type         <- 'matern_5_2'
max_eval            <- 100
alpha               <- NA  #alpha value for exponential kernels 

#
# mcmc parameters
#
N_steps <- 100000

#
# specify plots
#
leave_one_out_validation_plot <- 0
mcmc_traceplot                <- 1
mcmc_histogram                <- 1

#
# prior info 
#
dimensional_prior_mean <- c(1.0, 1.0, 1.0, 0.0,0.0 ,200.0, 5.0);
dimensional_prior_sd   <- c(0.3, 0.3, 0.3, 1.2, 200.0, 100.0, 2.5)
input_headers <- c("weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",
                   "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration")
#print("!!!!! warning: prior is updated in the code to equal the final iteration mean and sd !!!! ")

#
# obs info
#
dimensional_observation <- 0 #dimensional value of actual (set to be relative to the truth)
dimensional_error_cov   <- 1 #error covariance 

#emulate_sample <- function(realization,iterations,members,verbose,
#                           training_iterations,method,nugget_est,kernel_type,max_eval,alpha,
#                           N_steps,
#                           leave_one_out_validation_plot,mcmc_traceplot,mcmc_histogram,
#                           dimensional_prior_mean,dimensional_prior_sd,input_headers,
#                           dimensional_observation,dimensional_error_cov){

                           
###########################################################
################ get the model results ####################
###########################################################

#
# set up stuff 
#
count <- 1 
n_output <- 3 #3 output variables
n_input  <- 7 #7 input variables
n_meta   <- 2 #number of meta data (iteration, member)
n_iterations <- length(iterations)
n_members <- length(members)
dimensional_model_output <- array(NA, dim = c(n_iterations*n_members,n_output))
dimensional_model_input  <- array(NA, dim = c(n_iterations*n_members,n_input))
meta_data    <- array(NA, dim = c(n_iterations*n_members,n_meta))

# 
# loop over iterations
#
for (iter in iterations) {

#pad zeros to realization and iteration
padded_realization <- sprintf("%03d", realization)
padded_iteration  <- sprintf("%03d", iter)

#get the params file
params_file_path <- paste0("../../model-inputs-and-outputs/realization", padded_realization,
                            "/iteration", padded_iteration, "/params.csv")
params <- as.matrix(read.csv(params_file_path, header = TRUE))

for (mem in members) {
    
    # pad zeros to member
    padded_member     <- sprintf("%03d", mem)
    
    # get the file path
    file_path <- paste0("../../model-inputs-and-outputs/realization", padded_realization,
                        "/iteration", padded_iteration, "/member", padded_member, "/outputs_cts.csv")
    
    # Check if the file exists and output 1 if it does
    if (verbose)
    if (file.exists(file_path)) {
        cat("Found the file at", file_path, "\n")
        
    } else {
        cat("Did not find the file at", file_path, "\n")
    }
    
    # Load the CSV file
    data <- read.csv(file_path,header = FALSE)
    
    # Add the model input and output to the arrays
    dimensional_model_output[count,] <- as.matrix(data)
    dimensional_model_input[count,]  <- params[mem, ]
    meta_data[count,1]   <- iter
    meta_data[count,2]   <- mem
    
    count <- count + 1
    
} #end loop over members   
} #end loop over iterations

#
# restrict to only the grounded volume
#
indices <- which(meta_data[,1] %in% training_iterations)
dimensional_model_input <- dimensional_model_input[indices,]
dimensional_model_output <- dimensional_model_output[indices,3]
dimensional_model_output <- matrix(dimensional_model_output, nrow = length(dimensional_model_output), ncol = 1 ) #make into a matrix
meta_data <- meta_data[indices,]

#
# normalize the data
#
normalize_z_score <- function(x) {
    (x - mean(x)) / sd(x)
}
#get sd and mean of input and output for use later
dimensional_model_input_sd <- apply(dimensional_model_input, 2, sd)
dimensional_model_input_mean <- apply(dimensional_model_input, 2, mean)
dimensional_model_output_sd <- apply(dimensional_model_output, 2, sd)
dimensional_model_output_mean <- apply(dimensional_model_output, 2, mean)

normalized_model_output <- apply(dimensional_model_output, 2, normalize_z_score)
normalized_model_input  <- apply(dimensional_model_input, 2, normalize_z_score)

normalized_prior_mean <- (dimensional_prior_mean - dimensional_model_input_mean)/dimensional_model_input_sd 
normalized_prior_sd   <- dimensional_prior_sd/dimensional_model_input_sd

normalized_actual <- (dimensional_observation - dimensional_model_output_mean)/dimensional_model_output_sd
normalized_error_cov <- dimensional_error_cov/dimensional_model_output_sd

###########################################################
################ Leave one out validation #################
###########################################################
if (leave_one_out_validation_plot){
#
#initialize storage (LOO = leave one out)
#
LOO_pred_mean = rep(NaN, length(normalized_model_output))
LOO_pred_u95  = rep(NaN, length(normalized_model_output))
LOO_pred_l95  = rep(NaN, length(normalized_model_output))
LOO_pred_sd   = rep(NaN, length(normalized_model_output))
  
# 
# loop over output members
#
for (i in 1:length(normalized_model_output)) {
    
    # Get a dataset which doesn't include this point
    all_except_i <- setdiff(1:length(normalized_model_output), i)  # Indices except i
    
    # training and testing data
    training_input <- normalized_model_input[all_except_i,]
    testing_input  <- t(matrix(normalized_model_input[i,]))
    training_output <- normalized_model_output[all_except_i]
    testing_output  <- normalized_model_output[i]
    
    #train the model on this data and evaluate on the left out point. Condition on whether exponential (needs extra argument) or not
    if (kernel_type == "pow_exp") {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval,alpha=alpha)
    } else {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)
    }
    model.predict<-predict(model,testing_input)
    
    #get the prediction and the upper and lower bounds
    LOO_pred_mean[i] = model.predict$mean
    LOO_pred_l95[i]  = model.predict$lower95
    LOO_pred_u95[i]  = model.predict$upper95
    LOO_pred_sd[i]   = model.predict$sd
    
} #end loop over output points
  
#
#create data frame
#
LOO_data <- data.frame(
x = normalized_model_output ,  #actual output values
y = LOO_pred_mean,      # model output
y_lower = LOO_pred_l95,# upper error bars (nb length of error bar, not the end of it for plotting)
y_upper = LOO_pred_u95, #lower error bars
sd = LOO_pred_sd) #standard deviation in the training points

#
# make the plots
# 

# filter the data into two data frames, depending on whether y is in the range (y_lower, y_upper) or not
LOO_data_in_range <- LOO_data %>%
filter(x > y_lower & x < y_upper)

LOO_data_out_of_range <- LOO_data %>%
filter(!(x > y_lower & x < y_upper))

inrange_color = rgb(42/255, 103/255, 131/255)
outrange_color = rgb(228/255,128/255,111/255)
  
#rmse function
rmse_fn <- function(actual, predicted) {
sqrt(mean((actual - predicted)^2))
}
rmse <- rmse_fn(LOO_data$y, LOO_data$x)

#make the plot of modelled versus emulated
P <- ggplot() + 
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  #add the 1-1 line

    geom_point(data = LOO_data_in_range, aes(x = x, y = y), color = inrange_color) +
    geom_errorbar(data = LOO_data_in_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = inrange_color, width = 0.05) +
    
    geom_point(data = LOO_data_out_of_range, aes(x = x, y = y), color = outrange_color) +
    geom_errorbar(data = LOO_data_out_of_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = outrange_color, width = 0.05) +
    
    geom_point(aes(x = normalized_actual, y = normalized_actual), color = "red", size = 3) +

    labs(title = paste("leave on out validation. RMSE:", rmse),
         x = "normalized model output",
         y = "normalized emulated output") +
    
    theme_minimal() 

P <- P + xlim(-2, 2) + ylim(-2,2)
print(P)

# make the normalized error plot 
P <- ggplot(rbind(data.frame(value = (LOO_data_in_range$x - LOO_data_in_range$y)/(LOO_data_in_range$sd), group = "inside 95% CI"),data.frame(value = (LOO_data_out_of_range$x - LOO_data_out_of_range$y)/(LOO_data_out_of_range$sd), group = "outside 95% CI")), aes(x = value, fill = group)) +
    geom_histogram(position = "identity", alpha = 1.0, binwidth = 0.25 ) +
    scale_fill_manual(values = c("inside 95% CI" = inrange_color, "outside 95% CI" = outrange_color)) +
    #labs(title = "Overlapping Histograms", x = "Value", y = "Frequency", fill = "Group") +
    theme_minimal()

P <- P + xlim(-4, 4) 
print(P)

} #end leave_one_out_validation_plot flag

###########################################################
################# Create the emulator #####################
###########################################################

model <- rgasp(design = normalized_model_input, 
            response = normalized_model_output, 
            nugget.est = nugget_est, 
            method=method, 
            kernel_type= kernel_type,
            max_eval = max_eval)

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
###########################################################
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

#return(list(dimensional_mcmc_data = dimensional_mcmc_data, model = model))

#}

#emulate_sample_data <- emulate_sample(realization,iterations,members,verbose,
#                                       training_iterations,method,nugget_est,kernel_type,max_eval,alpha,
#                                       N_steps,
#                                       leave_one_out_validation_plot,mcmc_traceplot,mcmc_histogram,
#                                       dimensional_prior_mean,dimensional_prior_sd,input_headers,
#                                       dimensional_observation,dimensional_error_cov)