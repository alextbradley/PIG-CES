#### libraries
library(RobustGaSP)
library(ggplot2)
library(RColorBrewer)
library(dplyr) #for filtering data frames
  setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

verbose = 0;  #set to one for lots of model output

#### get model outputs
realization <- 30;
iterations  <- 1:5;
members     <- 1:20;

#initialize structure to store outputs
n_output <- 3 #3 output variables
n_input  <- 7 #7 input variables
n_meta   <- 2 #number of meta data (iteration, member)

n_iterations <- length(iterations)
n_members <- length(members)
model_output <- array(NA, dim = c(n_iterations*n_members,n_output))
model_input  <- array(NA, dim = c(n_iterations*n_members,n_input))
meta_data    <- array(NA, dim = c(n_iterations*n_members,n_meta))
count <- 1 
input_headers <- c("weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",
           "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration")
output_headers <- c("dimensionless_gl_error_1930", "dimensionless_gl_error_2015", "dimensionless_grv_error_2015")

#get the data
for (iter in iterations) {
  #pad zeros to member
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
    
  }
}

##### normalise the data
normalize_z_score <- function(x) {
  (x - mean(x)) / sd(x)
} # Z-Score Normalization function
normalized_model_output <- apply(model_output, 2, normalize_z_score)
normalized_model_input  <- apply(model_input, 2, normalize_z_score)

##### make scatter plots of parameter values and grounded volume error as a function of iteration
plot_input_vs_output <- function(input, output, iteration) {
 # if (ncol(input) != 7 || ncol(output) < 1) {
#    stop("Input must have 7 columns and output must have at least 3 columns")
 # }
  
  #compute the colours
  palette <- brewer.pal(max(iteration), "YlOrRd")
  colors <- unlist(sapply(iteration, function(x) palette[x]))

  
  par(mfrow = c(3, 3))  # Set up a 3x3 plotting area
  
  for (i in 1:7) {
    plot(input[, i], output, xlab = input_headers[i], ylab = "normalized output",ylim = c(-1, 1), pch = 21, bg =colors,col=NA)
    abline(h = 0, col = "black", lwd = 1 ,lty= 2)  # Add horizontal line at y = 0
  }
  
  # Reset plotting area to default
  #par(mfrow = c(1, 1))
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(2, 2, 2, 2), mgp = c(2, 0.5, 0))  # Set up a 3x3 plotting area with reduced margins
  
}
#plot_input_vs_output(model_input, normalized_model_output[,3],meta_data[,1])


#### do the model emulation via leave one out
# build the model and do the leave one out validation, and plot the results
method = 'post_mode'
nuggest_est = F
kernel_type = 'matern_5_2'
max_eval = 100
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

#do the leave one out validation
method = 'post_mode'
nuggest_est = F
max_eval = 100
alpha = NA
kernel_type = "matern_5_2"
leave_one_out_data <- leave_one_out_validation(normalized_model_input, normalized_model_output[,3], method, nuggest_est, kernel_type, max_eval,alpha)

#print the RMSE
rmse_fn <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}
rmse <- rmse_fn(leave_one_out_data$y, leave_one_out_data$x)
print(paste("normalised RMSE:", rmse))
print(paste("actual RMSE:", (rmse * sd(model_output,3))))

##### plot the modelled versus emulated
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

P1 <- plot_modelled_vs_emulated(leave_one_out_data)
P1 <- P1 + xlim(-1, 1) + ylim(-1,1)
P1

##### plot the normalised error histogram
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
P2 <- plot_normalized_errors(leave_one_out_data)
P2 <- P2 + xlim(-4, 4) 
P2


##### repeat for different types of kernel
#method = 'post_mode'
#nuggest_est = F
#max_eval = 100
#kernel_list <- list("matern_5_2", "matern_3_2", "pow_exp", "pow_exp", "pow_exp")
#alphas <- c(NA,NA,1.9, 1.0, 0.1)

#count <- 1
#plot_list <- list()
#rmses <- list()
#for (kernel in kernel_list) {
#  leave_one_out_data <- leave_one_out_validation(normalized_model_input, normalized_model_output[,3], method, nuggest_est, kernel, max_eval, alphas[count])
#  P1 <- plot_modelled_vs_emulated(leave_one_out_data)
#  P1 <- P1 + xlim(-1, 1) + ylim(-1,1)
#  
#  rmse <- rmse_fn(leave_one_out_data$y, leave_one_out_data$x)
#  
#  rmses[[count]] <- rmse
#  plot_list[[count]] <- P1
#  
#  count <- count + 1
#}

##### emulate the other variables
method = 'post_mode'
nuggest_est = N
max_eval = 100
alpha = 0.1
kernel_type = "pow_exp"
kernel_type = "matern_3_2"

leave_one_out_data_1930gl <- leave_one_out_validation(normalized_model_input, normalized_model_output[,1], method, nuggest_est, kernel_type, max_eval,alpha)
P1 <- plot_modelled_vs_emulated(leave_one_out_data_1930gl)
P1 <- P1 + xlim(-1, 1) + ylim(-2,2)
P1

##### make the emulators
method = 'post_mode'
nugget_est = F
max_eval = 100
alpha = NA
kernel_type = "matern_3_2"
model_2015grv <- rgasp(design = normalized_model_input, response = normalized_model_output[,3], nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)


##### what if we only emulate the final iteration?
indices <- which(meta_data[,1] >= max(meta_data[,1])-1); #indices corresponding to final iteration or final two?
leave_one_out_data_final_iter <- leave_one_out_validation(normalized_model_input[indices, ], normalized_model_output[indices,3], method, nuggest_est, kernel_type, max_eval,alpha)
P1 <- plot_modelled_vs_emulated(leave_one_out_data_final_iter)
P1 <- P1 + xlim(-0.2, 0.2) + ylim(-0.2,0.2)
P1

#make model based in this
model_2015grv <- rgasp(design = normalized_model_input[indices,], response = normalized_model_output[indices,3], nugget.est = nugget_est, method=method, kernel_type= kernel_type,max_eval = max_eval)


