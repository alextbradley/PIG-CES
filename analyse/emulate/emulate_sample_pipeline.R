# run the emulate-sample pipeline for a single realization of forcing
rm(list = ls())
setwd("/Users/user/Dropbox/BAS_Postdoc/Projects/AttributionRealWorld/manual-EKI/analyse/emulate/")

source("shared.R")

################################################################################
### set run info
################################################################################
write_output <- 1
for (realization in 36){
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
### split the output into parts
################################################################################

model_output <- matrix(model_output, nrow = length(model_output[,1]), ncol = 3 ) #make into a matrix
colnames(model_output) <- c("gl_error_1930","gl_error_2015","grv_error_2015")
#observations <- observations[,3]
#observations_noise <- observations_noise[,3]

# adjust the observation noise? Make it 1% of the observed mass?
grv_2015 <- read.csv("../../observations/truth_actual.csv", header = FALSE)/1e12 #everything is normalized by 1e12
grv_2015 <- grv_2015[,3]
observations_noise[,3] <- 0.01 * grv_2015 # 1 percent error

################################################################################
### remove the ungrounded_weertman_c_prefactor
################################################################################
model_input <- model_input[,c(1,3,4,5,6,7)]
input_colnames <- input_colnames[c(1,3,4,5,6,7)]

################################################################################
### normalize the input and output data
################################################################################

#prior_mean <- c(1.0, 1.0, 1.0, 0.0,0.0 ,200.0, 5.0);
#prior_sd   <- c(0.3, 0.3, 0.3, 1.2, 200.0, 100.0, 2.5) #with the ungrounded_weertman_c_prefactor

prior_mean <- c(1.0, 1.0, 0.0, 0.0,   200.0, 5.0);
prior_sd   <- c(0.3, 0.3, 1.2, 200.0, 100.0, 2.5) #with the ungrounded_weertman_c_prefactor

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

if (write_output){ #write the output of the LOOCV
  file_name_gl1930 <- paste0("./mcmc_output/realization",padded_realization, "_LOOCV_gl1930.csv")
  file_name_gl2015 <- paste0("./mcmc_output/realization",padded_realization, "_LOOCV_gl2015.csv")
  file_name_grv2015 <- paste0("./mcmc_output/realization",padded_realization, "_LOOCV_grv2015.csv")

  write.csv(LOO_data_gl_1930, file = file_name_gl1930, row.names = FALSE)
  write.csv(LOO_data_gl_2015, file = file_name_gl2015, row.names = FALSE)
  write.csv(LOO_data_grv_2015, file = file_name_grv2015, row.names = FALSE)
  
}

################################################################################
### create the emulators
################################################################################
#1930 grounding line
model_gl1930 <- rgasp(design = normalized_model_input, 
                      response = normalized_model_output[,1], 
                      nugget.est = nugget_est, 
                      method=method, 
                      kernel_type= kernel_type,
                      max_eval = max_eval)


#2015 grounding line
model_gl2015 <- rgasp(design = normalized_model_input, 
                       response = normalized_model_output[,2], 
                       nugget.est = nugget_est, 
                       method=method, 
                       kernel_type= kernel_type,
                       max_eval = max_eval)


#2015 grounded volume
model_grv2015 <- rgasp(design = normalized_model_input, 
               response = normalized_model_output[,3], 
               nugget.est = nugget_est, 
               method=method, 
               kernel_type= kernel_type,
               max_eval = max_eval)

################################################################################
### make main effects plots (f)
########################################################################### #####

# use the prior mean as the nominal values
nominal_values <- prior_mean #values along which to take the plots

# or use the mean of the final iteration params
final_iteration_parameters <- model_input[which(meta_data[,1] == max(meta_data[,1])),]
nominal_values <- colMeans(final_iteration_parameters)

#source("meff.R")

################################################################################
### run the mcmc
################################################################################

final_iteration_parameters <- model_input[which(meta_data[,1] == max(meta_data[,1])),]
init_sample                <- colMeans(final_iteration_parameters)

n_steps <- 51000
n_burn  <- 1000

#n_steps <- 30000
#n_burn  <- 1000

source("run_mcmc.R")


################################################################################
### output the samples
################################################################################
posterior_samples_noindex <- posterior_samples[,1:6]
file_name <- paste0("./mcmc_output/mcmc_output_realization", padded_realization, ".csv")
colnames(posterior_samples_noindex) <- input_colnames
  
if (write_output){
  write.csv(posterior_samples_noindex, file = file_name, row.names = FALSE)
  }
################################################################################
### sample from the posterior with and without a trend 
################################################################################
#source("sample_w_wo_trend.R")
}