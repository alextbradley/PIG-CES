library(RobustGaSP)
library(ggplot2)
setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

######## data and preprocessing
data <- read.csv("../../realization002/model_simulations.csv")
inputs  <- data[,1:7] #input features
outputs <- data[,8:9] #output features

# normalize the inputs
normalize_z_score <- function(x) {
  (x - mean(x)) / sd(x)
} # Z-Score Normalization function
normalized_matrix_z_score <- apply(data, 2, normalize_z_score)
normalized_input <- normalized_matrix_z_score[,1:7] #normalized input

# normalize the outputs
observations <- read.csv("../../observations/truth.csv", header = FALSE);
observations_repmat <- data.frame(matrix(rep(observations, each = 30), nrow = 30, byrow = FALSE));
colnames(observations_repmat) <- colnames(outputs)

outputs$gl_position_1930  <- as.numeric(outputs$gl_position_1930)
outputs$gl_position_2015  <- as.numeric(outputs$gl_position_2015)
observations_repmat$gl_position_1930  <- as.numeric(observations_repmat$gl_position_1930)
observations_repmat$gl_position_2015  <- as.numeric(observations_repmat$gl_position_2015)
shifted_output <- outputs - observations_repmat

std_devs <- sapply(shifted_output, sd)
normalized_output <- sweep(shifted_output, 2, std_devs, FUN = "/")
#normalized_output <- shifted_output/10000 #std not a good metric bc of outliers?

####### leave one out cross validation for independently emulated estimates
pred_mean_1930 = rep(NaN, length(normalized_output[,1]))
pred_u95_1930  = rep(NaN, length(normalized_output[,1]))
pred_l95_1930  = rep(NaN, length(normalized_output[,1]))
is_within      = rep(1, length(normalized_output[,1]))

for (i in 1:length(normalized_output[,1])) {
  
  # Get a dataset which doesn't include this point
  all_except_i <- setdiff(1:length(normalized_output[,1]), i)  # Indices except i
  
  # training and testing data
  training_input <- normalized_input[all_except_i,]
  testing_input  <- t(matrix(normalized_input[i,]))
  training_output <- normalized_output[all_except_i,1]
  testing_output  <- normalized_output[i,1]
  
  #train the model on this data and evaluate on the left out point
  model <- rgasp(design = training_input, response = training_output, nugget.est = T, method = 'post_mode', kernel_type= 'matern_5_2',max_eval = 100)
  model.predict<-predict(model,testing_input)
  
  #get the prediction and the upper and lower bounds
  pred_mean_1930[i] <- model.predict$mean
  pred_l95_1930[i]  <- model.predict$lower95
  pred_u95_1930[i]  <- model.predict$upper95
  
  #store whether this prediction has bounds intersecting
  if (model.predict$lower95 > normalized_output[i, 1] || model.predict$upper95 < normalized_output[i, 1]) {
    is_within[i] <- 0
  }
  

  
}

#print the rmse
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}
rmse <- rmse(normalized_output[,1] , pred_mean_1930)
print(paste("RMSE for 1930:", rmse))

#print the coverage proportion
prop_emulator_1930 <- length(which((pred_l95_1930<=normalized_output[,1])
                              &(pred_u95_1930>=normalized_output[,1])))/length(normalized_output[, 1])

print(paste("proportion emulated within 95% CI for 1930:", prop_emulator_1930))


data <- data.frame(
  x = normalized_output[,1] ,  #actual output values
  y = pred_mean_1930,          # model output
  y_lower = pred_l95_1930 ,   # upper error bars
  y_upper = pred_u95_1930)    #lower error bars


P <- plot_modelled_vs_emulated(data, 1930, is_within)
print(P)

######### repeat for 2015 #################################
pred_mean_2015 = rep(NaN, length(normalized_output[,2]))
pred_u95_2015  = rep(NaN, length(normalized_output[,2]))
pred_l95_2015  = rep(NaN, length(normalized_output[,2]))
is_within_2015      = rep(1, length(normalized_output[,2]))

for (i in 1:length(normalized_output[,2])) {
  
  # Get a dataset which doesn't include this point
  all_except_i <- setdiff(1:length(normalized_output[,2]), i)  # Indices except i
  
  # training and testing data
  training_input <- normalized_input[all_except_i,]
  testing_input  <- t(matrix(normalized_input[i,]))
  training_output <- normalized_output[all_except_i,2]
  testing_output  <- normalized_output[i,2]
  
  #train the model on this data and evaluate on the left out point
  model <- rgasp(design = training_input, response = training_output, nugget.est = T, method = 'post_mode', kernel_type= 'matern_5_2',max_eval = 100)
  model.predict<-predict(model,testing_input)
  
  #get the prediction and the upper and lower bounds
  pred_mean_2015[i] <- model.predict$mean
  pred_l95_2015[i]  <- model.predict$lower95
  pred_u95_2015[i]  <- model.predict$upper95
  
  #store whether this prediction has bounds intersecting
  if (model.predict$lower95 > normalized_output[i, 2] || model.predict$upper95 < normalized_output[i, 2]) {
    is_within_2015[i] <- 0
  }
  
}

#print the rmse
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}
rmse <- rmse(normalized_output[,2] , pred_mean_2015)
print(paste("RMSE for 2015:", rmse))

#print the coverage proportion
prop_emulator_2015 <- length(which((pred_l95_2015<=normalized_output[,2])
                              &(pred_u95_2015>=normalized_output[,2])))/length(normalized_output[, 2])
print(paste("proportion emulated within 95% CI for 2015:", prop_emulator_2015))


data_2015 <- data.frame(
  x = normalized_output[,2] ,  #actual output values
  y = pred_mean_2015,          # model output
  y_lower = pred_l95_2015 ,   # upper error bars
  y_upper = pred_u95_2015)    #lower error bars


P <- plot_modelled_vs_emulated(data_2015, 2015, is_within_2015)
print(P)




####### functions
plot_modelled_vs_emulated <- function(data,year, is_within) {
  P <- ggplot(data, aes(x = x, y = y, color = factor(is_within))) +
    geom_point() +# Scatter plot points
    geom_abline(slope = 1, intercept = 0, color = "black") +  #add the 1-1 line
    geom_errorbar(aes(ymin = y_lower, ymax = y_upper), width = 0.1) +  # Error bars
    scale_color_manual(values = c("red", "blue"))
    labs(title = paste("leave one out validation for", as.character(year) , sep = " "),
         x = "model output",
         y = "emulated output") +
    theme_minimal() #+                                      # Minimal theme
  #  xlim(c(0, 20))  # Set x-axis limits
  
  return(P)
}

