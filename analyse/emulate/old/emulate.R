##### Script to emulate model output and sample using MCMC

###########################################################
##################### Preliminaries #######################
###########################################################

rm(list = ls())

#
# libraries
#
include("shared.R")


#
# make sure you run from the right place
#
setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

#
# specify model results to use
#
realization <- 30
iterations  <- 1:5
members     <- 1:20
verbose     <- 0

#
# emulator parameters
#
training_iterations <- 1:5 #which iterations to train on
method              <- 'post_mode'
nugget_est          <- F
kernel_type         <- 'matern_3_2'
max_eval            <- 100
alpha               <- NA  #alpha value for exponential kernels 

#
# mcmc parameters
#
N_steps <- 100000

#
# specify plots and output
#
leave_one_out_validation_plot <- 1
mcmc_traceplot                <- 1
mcmc_histogram                <- 1

output_mcmc <-0 

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

coverage <- dim(LOO_data_in_range)[1]/dim(LOO_data)[1]

model <- lm(LOO_data$y ~ LOO_data$x)
r_squared <- summary(model)$r.squared



#make the plot of modelled versus emulated
P <- ggplot() + 
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  #add the 1-1 line

    geom_point(data = LOO_data_in_range, aes(x = x, y = y), color = inrange_color) +
    geom_errorbar(data = LOO_data_in_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = inrange_color, width = 0.05) +
    
    geom_point(data = LOO_data_out_of_range, aes(x = x, y = y), color = outrange_color) +
    geom_errorbar(data = LOO_data_out_of_range, aes(x = x, ymin = y_lower, ymax = y_upper), color = outrange_color, width = 0.05) +
    
    geom_point(aes(x = normalized_actual, y = normalized_actual), color = "red", size = 3) +

  
    labs(title = sprintf("RMSE:%.3f, coverage: %.3f, R^2: %.3f", rmse, coverage,r_squared),
         x = "normalized model output",
         y = "normalized emulated output") +
    
    theme_minimal() 

P <- P + xlim(-2, 2) + ylim(-4,4)
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

