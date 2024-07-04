
#####
##### make mcmc traceplots
#####
make_mcmc_traceplot <- function(mcmc_data){
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
  print(P)
}

make_mcmc_histogram <- function(mcmc_data,prior_mean,prior_sd,naive_posterior_mean,naive_posterior_sd, input_headers){
  plot_list <- list()
  nbars < 30; #number of histogram bars
  
  
  for (i in 1:7) {
    
    #get the prior and put in a data frame
    xx <- seq(from = min(mcmc_data[,i]), to = max(mcmc_data[,i]), length.out = 100)
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
    p <- ggplot(mcmc_data, aes_string(x = colnames(mcmc_data)[i])) +
      geom_histogram(aes(y = after_stat(density)), binwidth = (max( mcmc_data[,i]) - min( mcmc_data[,i]))/nbars, fill = "blue", alpha = 0.7) +
      geom_density(alpha = 1, color = "black", linewidth =0.75) +
      #geom_line(aes(x = xx, y = yy), color = "black", linetype = "dashed") +
      theme_minimal() +
      ggtitle(input_headers[i]) 
    #p <- p + geom_line(data = data_prior, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.75) 
    #p <- p + geom_line(data = data_naivepost, aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 0.75)
    
    plot_list[[i]] <- p
  } #end loop over variables
  grid.arrange(grobs = plot_list, ncol = 3)
}


#####
##### run mcmc function
#####
run_mcmc <- function(theta0, N_steps,C, model, prior_mean, prior_covariance, actual, error_cov){
  theta_values <- matrix(0, nrow = N_steps, ncol = length(theta0))
  #theta_values[1,] <- theta0
  update_size <- mvrnorm(1, c(0,0,0,0,0,0,0), C) 
  theta_current <- theta0 + update_size #seed the current value
  inv_prior_covariance <- solve(prior_covariance)
  
  for (i in 1:(N_steps-1)) {
    
    #determine the trial theta
    update_size <- mvrnorm(1, c(0,0,0,0,0,0,0), C) #draw from multivariable random normal with zero mean
    theta_trial <- theta_current + update_size
    
    #compute emulator likelihoods
    likelihood_term_trial = exp(-(actual - predict(model,t(theta_trial))$mean )^2 / 2 /error_cov^2)
    likelihood_term_current = exp(-(actual - predict(model,t(theta_current))$mean )^2 / 2 /error_cov^2)
    
    #compute prior terms
    prior_diff_trial    <-  matrix(theta_trial - prior_mean, 1, 7) 
    prior_diff_current  <-  matrix(theta_current - prior_mean, 1, 7) 
    prior_term_trial    <- exp(- prior_diff_trial    %*% inv_prior_covariance %*% t(prior_diff_trial) / 2);
    prior_term_current  <- exp(- prior_diff_current  %*% inv_prior_covariance %*% t(prior_diff_current) / 2);
    
    #likelihood ratio
    a <- min(c(1, (likelihood_trial*prior_term_trial / (likelihood_prev*prior_term_prev)) ));
    
    #perform update and store
    random_number <- runif(1, min = 0, max = 1)
    if (random_number < a){ #accept
      theta_current <- theta_trial
    }
    theta_values[i,] <- theta_current
    
    #progress tracker
    if (i%%1000 == 0)
    {
      print(i)
    } 
  } #end loop over values
  
  data_df <- as.data.frame(theta_values) #make into dataframe
  data_df$Index <- 1:dim(theta_values)[1] #set an index
  return(data_df)
}

#####
##### get data function
#####
get_data <- function(realization, iterations, members, verbose){
  count <- 1 
  n_output <- 3 #3 output variables
  n_input  <- 7 #7 input variables
  n_meta   <- 2 #number of meta data (iteration, member)
  model_output <- array(NA, dim = c(n_iterations*n_members,n_output))
  model_input  <- array(NA, dim = c(n_iterations*n_members,n_input))
  meta_data    <- array(NA, dim = c(n_iterations*n_members,n_meta))
  
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
      model_output[count,] <- as.matrix(data)
      model_input[count,]  <- params[mem, ]
      meta_data[count,1]   <- iter
      meta_data[count,2]   <- mem
      
      count <- count + 1
      
    } #end loop over members   
  } #end loop over iterations
  
  return(list(meta_data = meta_data, model_input = model_input, model_output = model_output))
} #end function 

#####
##### plot modelled vs emulated 
#####
plot_modelled_vs_emulated <- function(data) {
  
  # filter the data into two data frames, depending on whether y is in the range (y_lower, y_upper) or not
  data_in_range <- data %>%
    filter(x > y_lower & x < y_upper)
  
  data_out_of_range <- data %>%
    filter(!(x > y_lower & x < y_upper))
  
  inrange_color = rgb(42/255, 103/255, 131/255)
  outrange_color = rgb(228/255,128/255,111/255)
  
  rmse_fn <- function(actual, predicted) {
    sqrt(mean((actual - predicted)^2))
  }
  rmse <- rmse_fn(data$y, data$x)
  
  P <- ggplot() + 
    geom_point(data = data_in_range, aes(x = x, y = y), color = inrange_color) +
    geom_errorbar(data = data_in_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = inrange_color, width = 0.05) +
    
    geom_point(data = data_out_of_range, aes(x = x, y = y), color = outrange_color) +
    geom_errorbar(data = data_out_of_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = outrange_color, width = 0.05) +
    
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  #add the 1-1 line
    
    labs(title = paste("leave on out validation. RMSE:", rmse),
         x = "model output",
         y = "emulated output") +
    
    theme_minimal() 
  
  return(P)
} #plotting script for modelled vs emulated


##### 
##### plot normalized errors
#####
plot_normalized_errors <- function(data){
  
  # filter the data into two data frames, depending on whether y is in the range (y_lower, y_upper) or not
  data_in_range <- data %>%
    filter(x > y_lower & x < y_upper)
  
  data_out_of_range <- data %>%
    filter(!(x > y_lower & x < y_upper))
  
  inrange_color = rgb(42/255, 103/255, 131/255)
  outrange_color = rgb(228/255,128/255,111/255)
  
  data1 <- data.frame(value = (data_in_range$x - data_in_range$y)/(data_in_range$sd), group = "inside 95% CI")
  data2 <- data.frame(value = (data_out_of_range$x - data_out_of_range$y)/(data_out_of_range$sd), group = "outside 95% CI")
  data <- rbind(data1, data2)
  
  P <- ggplot(data, aes(x = value, fill = group)) +
    geom_histogram(position = "identity", alpha = 1.0, binwidth = 0.25 ) +
    scale_fill_manual(values = c("inside 95% CI" = inrange_color, "outside 95% CI" = outrange_color)) +
    #labs(title = "Overlapping Histograms", x = "Value", y = "Frequency", fill = "Group") +
    theme_minimal()
  
  return(P)
  
}

##### 
##### leave one out validation
#####
leave_one_out_validation <- function(input, output, method, nugget_est, kernel_type, max_eval, alpha){
  
  #initialize storage
  pred_mean = rep(NaN, length(output))
  pred_u95  = rep(NaN, length(output))
  pred_l95  = rep(NaN, length(output))
  pred_sd   = rep(NaN, length(output))
  
  for (i in 1:length(output)) {
    
    # Get a dataset which doesn't include this point
    all_except_i <- setdiff(1:length(output), i)  # Indices except i
    
    # training and testing data
    training_input <- input[all_except_i,]
    testing_input  <- t(matrix(input[i,]))
    training_output <- output[all_except_i]
    testing_output  <- output[i]
    
    #train the model on this data and evaluate on the left out point. Condition on whether exponential (needs extra argument) or not
    if (kernel_type == "pow_exp") {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval,alpha=alpha)
    } else {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)
    }
    model.predict<-predict(model,testing_input)
    
    #get the prediction and the upper and lower bounds
    pred_mean[i] = model.predict$mean
    pred_l95[i]  = model.predict$lower95
    pred_u95[i]  = model.predict$upper95
    pred_sd[i]   = model.predict$sd
    
  } #end loop over output points
  
  #create data frame
  data <- data.frame(
    x = output ,  #actual output values
    y = pred_mean,      # model output
    y_lower = pred_l95,# upper error bars (nb length of error bar, not the end of it for plotting)
    y_upper = pred_u95, #lower error bars
    sd = pred_sd) #standard deviation in the training points 
  
  
  return(data)
} 