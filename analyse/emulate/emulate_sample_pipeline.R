# run the emulate-sample pipeline for a single realization of forcing
rm(list = ls())
setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

source("shared.R")

################################################################################
### set run info
################################################################################
realization <- 30
iterations  <- 1:5
members     <- 1:20

################################################################################
### load the observations
################################################################################
observations <- read.csv("../../observations/truth.csv", header = FALSE) 
#observations of:
# - error in the grounding line position in 1930, in number of grid cells (3km in size) 
# - error in the grounding line position in 2015, in number of grid cells (3km in size) 
# - error in the grounded volume in 2015 / 1e12 (m^3)
observations_noise <- read.csv("../../observations/noise.csv", header = FALSE) 


################################################################################
### get the simulation data
################################################################################
source("get_simulation_data.R")


################################################################################
### restrict to only the grounded volume output
################################################################################

model_output <- model_output[,3]
model_output <- matrix(model_output, nrow = length(model_output), ncol = 1 ) #make into a matrix
colnames(model_output) <- c("grv_error_2015")
observations <- observations[,3]
observations_noise <- observations_noise[,3]

################################################################################
### normalize the input and output data
################################################################################

prior_mean <- c(1.0, 1.0, 1.0, 0.0,0.0 ,200.0, 5.0);
prior_sd   <- c(0.3, 0.3, 0.3, 1.2, 200.0, 100.0, 2.5)

input_normalization_mean <- prior_mean
input_normalization_sd   <- prior_sd
#input_normalization_mean <- colMeans(model_input)
#input_normalization_sd   <- apply(model_input, 2, sd)

output_normalization_mean <- colMeans(model_output)
output_normalization_sd   <- apply(model_output, 2, sd)

source("normalize_simulation_data.R")

################################################################################
### perform the loocv
################################################################################

training_iterations <- 1:5 #which iterations to train on
method              <- 'post_mode'
nugget_est          <- F
kernel_type         <- 'matern_3_2'
max_eval            <- 100
alpha               <- NA  #alpha value for exponential kernels 
 
source("loocv.R")

################################################################################
### create the emulator
################################################################################
model <- rgasp(design = normalized_model_input, 
               response = normalized_model_output, 
               nugget.est = nugget_est, 
               method=method, 
               kernel_type= kernel_type,
               max_eval = max_eval)

################################################################################
### make main effects plots
################################################################################

# use the prior mean as the nominal values
nominal_values <- prior_mean #values along which to take the plots

# or use the mean of the final iteration params
final_iteration_parameters <- model_input[which(meta_data[,1] == max(meta_data[,1])),]
nominal_values <- colMeans(final_iteration_parameters)

source("meff.R")

################################################################################
### run the mcmc
################################################################################

final_iteration_parameters <- model_input[which(meta_data[,1] == max(meta_data[,1])),]
init_sample                <- colMeans(final_iteration_parameters)

n_steps <- 200000
n_burn  <- 1000

source("run_mcmc.R")



