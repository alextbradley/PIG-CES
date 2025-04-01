# loocv.R
#
# Perform the leave one out cross validation on normalized model output. 
# Produces one plot for each of the model inputs. Assumes same emulator structure for each.
# 
# inputs: normalized_model_output, normalized_model_input, training_iterations, method, nugget_est, kernel_type, max_eval, alpha, normalization mean and sd for inputs and outputs
#
#   
# outputs: LOO_data
# plots made: plot of modelled vs emulated 

#
# initialized
for (iO in 1:length(observations)){
  res_normalized_model_output = normalized_model_output[,iO] #restrict to this particular output
  
  LOO_normalized_pred_mean = rep(NaN, length(res_normalized_model_output))
  LOO_normalized_pred_u95  = rep(NaN, length(res_normalized_model_output))
  LOO_normalized_pred_l95  = rep(NaN, length(res_normalized_model_output))
  LOO_normalized_pred_sd   = rep(NaN, length(res_normalized_model_output))
  
  # 
  # loop over output members
  #
  for (i in 1:length(res_normalized_model_output)) {
    
    # Get a dataset which doesn't include this point
    all_except_i <- setdiff(1:length(res_normalized_model_output), i)  # Indices except i
    
    # training and testing data
    training_input <- normalized_model_input[all_except_i,]
    testing_input  <- t(matrix(normalized_model_input[i,]))
    training_output <- res_normalized_model_output[all_except_i]
    testing_output  <- res_normalized_model_output[i]
    
    #train the model on this data and evaluate on the left out point. Condition on whether exponential (needs extra argument) or not
    if (kernel_type == "pow_exp") {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval,alpha=alpha)
    } else {
      model <- rgasp(design = training_input, response = training_output, nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)
    }
    model.predict<-predict(model,testing_input)
    
    #get the prediction and the upper and lower bounds
    LOO_normalized_pred_mean[i] = model.predict$mean
    LOO_normalized_pred_l95[i]  = model.predict$lower95
    LOO_normalized_pred_u95[i]  = model.predict$upper95
    LOO_normalized_pred_sd[i]   = model.predict$sd
    
  } #end loop over output points
  
  #
  #create data frame (regular and normalized)
  #
  LOO_normalized_data <- data.frame(
    x = as.vector(res_normalized_model_output) ,    #actual output values
    y = LOO_normalized_pred_mean,             # model output
    y_lower = LOO_normalized_pred_l95,        # upper error bars (nb length of error bar, not the end of it for plotting)
    y_upper = LOO_normalized_pred_u95,        #lower error bars
    sd = LOO_normalized_pred_sd) #standard deviation in the training points
  
  LOO_data <- data.frame(
    x = as.vector(model_output[,iO]) ,  #actual output values
    y = LOO_normalized_pred_mean*output_normalization_sd[iO] + output_normalization_mean[iO],      # model output
    y_lower = LOO_normalized_pred_l95*output_normalization_sd[iO] + output_normalization_mean[iO], # upper error bars (nb length of error bar, not the end of it for plotting)
    y_upper = LOO_normalized_pred_u95*output_normalization_sd[iO] + output_normalization_mean[iO], #lower error bars
    sd = LOO_normalized_pred_sd*output_normalization_sd[iO]) #standard deviation in the training points
  
  #
  # make plot
  #
  
  # filter the data into two data frames, depending on whether y is in the range (y_lower, y_upper) or not
  LOO_data_in_range <- LOO_data %>%
    filter(x > y_lower & x < y_upper)
  
  LOO_data_out_of_range <- LOO_data %>%
    filter(!(x > y_lower & x < y_upper))
  
  inrange_color = rgb(42/255, 103/255, 131/255)
  outrange_color = rgb(228/255,128/255,111/255)
  
  #get rmse, coverage, and r^2
  rmse_fn <- function(actual, predicted) {
    sqrt(mean((actual - predicted)^2))
  }
  rmse <- rmse_fn(LOO_data$y, LOO_data$x)
  coverage <- dim(LOO_data_in_range)[1]/dim(LOO_data)[1]
  model <- lm(LOO_data$y ~ LOO_data$x)
  r_squared <- summary(model)$r.squared
  
  #make plot 
  P <- ggplot() + 
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  #add the 1-1 line
    
    geom_point(data = LOO_data_in_range, aes(x = x, y = y), color = inrange_color) +
    geom_errorbar(data = LOO_data_in_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = inrange_color, width = 0.05) +
    
    geom_point(data = LOO_data_out_of_range, aes(x = x, y = y), color = outrange_color) +
    geom_errorbar(data = LOO_data_out_of_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = outrange_color, width = 0.05) +
    
    #geom_point(aes(x = array(observations[i0]), y = array(observations[i0])), color = "red", size = 4) +
    
    ggtitle(colnames(model_output)[iO]) +
    
    labs(
      x = "model error",
      y = "emulator error") +
    
    theme_minimal()  + 
    theme(
      axis.title.x = element_text(size = 10), # Change X-axis label font size
      axis.title.y = element_text(size = 10)  # Change Y-axis label font size
    )
  
  P <- P + geom_text(aes(x = -30, y = 40), label = sprintf("RMSE:%.3f \ncoverage: %.3f \n R^2: %.3f", rmse, coverage,r_squared), color = "black", size = 4)
  P <- P + xlim(-40, 40) + ylim(-60,60)
  show(P)
  
  if (iO == 1) {
    LOO_data_gl_1930 <- LOO_data}
  
  if (iO == 2) {
    LOO_data_gl_2015 <- LOO_data}
  
  if (iO == 3) {
    LOO_data_grv_2015 <- LOO_data}
  
  
  
}
