# get_simulation_data.R
# 
# Return the simulation data for a specified realization, iterations and members.
#
# inputs:
#
# realization : scalar
#               Realization number
#
# iterations  : (n_iterations x 1) array
#               list of the iteration numbers used. n_iterations: number of 
#               iterations used
#
# members     : (n_members x 1) array
#               list of the member numbers used. n_members: number of members 
#               used
#
#
# outputs (currently not a function):
#
# model_input   : (n_iterations * n_members x n_input) array
#                 array of model input. Each row contains n_input (=7) model 
#                 inputs (see repository README for more info):
#                 - weertman_c_prefactor: prefactor on the sliding law
#                 - ungrounded_weertmanC_prefactor: << redundant -- ignore >>
#                 - glen_a_ref_prefactor: prefactor on the ice viscosity
#                 - melt_rate_prefactor_exponent: prefactor exponent on the melt rate
#                 - per_century_trend: trend in the forcing
#                 - bump_amplitude: amplitude in m of the 1940s anomaly
#                 - bump_duration: duration of the 1940s anomaly
#
# model_output  : (n_iterations * n_members x n_output) array
#                 array of model output. Each row contains n_output (= 3) entries:
#                 - error in the grounding line position in 1930, in number of grid cells (3km in size)
#                 - error in the grounding line position in 2015, in number of grid cells (3km in size)
#                 - error in the grounded volume in 2015 / 1e12 (m^3)
#
# meta_data     : (n_iterations * n_members x 2) array
#                 simulation meta data, containing the iteration and member no
#             

################################################################################
################################ Preliminaries #################################
################################################################################

count        <- 1 #initialize counter for the number of simulations
n_output     <- 3 #number of output variables
n_input      <- 7 #number of input variables
n_meta       <- 2 #number of meta data (iteration, member)

n_iterations <- length(iterations)
n_members    <- length(members)

# initilize storage of simulation inputs, outputs and metadata
model_output <- array(NA, dim = c(n_iterations*n_members,n_output))
model_input  <- array(NA, dim = c(n_iterations*n_members,n_input))
meta_data    <- array(NA, dim = c(n_iterations*n_members,n_meta))

# column names
output_colnames <- c("gl_error_1930",  #error in the grounding line position in 1930, in number of grid cells (3km in size)
                     "gl_error_2015",  #error in the grounding line position in 2015, in number of grid cells (3km in size)
                     "grv_error_2015") #error in the grounded volume in 2015 / 1e12 (m^3)

input_colnames <- c("weertman_c_prefactor", 
                   "ungrounded_weertmanC_prefactor" ,
                   "glen_a_ref_prefactor",
                   "melt_rate_prefactor_exponent",
                   "per_century_trend",
                   "bump_amplitude",
                   "bump_duration")

################################################################################
############################ Loop over simulations #############################
################################################################################
for (iter in iterations) #loop over the iterations
  {
  
  #pad zeros to realization and iteration for file names
  padded_realization <- sprintf("%03d", realization)
  padded_iteration  <- sprintf("%03d", iter)
  
  #get the file containing the model input parameters
  params_file_path <- paste0("../../model-inputs-and-outputs/realization", padded_realization,
                             "/iteration", padded_iteration, "/params.csv")
  params <- as.matrix(read.csv(params_file_path, header = TRUE))
  
  for (mem in members) #loop over each of the members
    {
    
    # pad zeros to member for file names
    padded_member <- sprintf("%03d", mem)
    
    # get the file path
    file_path <- paste0("../../model-inputs-and-outputs/realization", padded_realization,
                        "/iteration", padded_iteration, "/member", padded_member, "/outputs_cts.csv")
    
    # Check if the file exists and output accordingly
    #if (file.exists(file_path)) {
    #  cat("Found the file at", file_path, "\n")
    #  
    #} else {
    #  cat("Did not find the file at", file_path, "\n")
    #}
    
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

colnames(model_input)  <- input_colnames
colnames(model_output) <- output_colnames
colnames(meta_data) <- c("iteration", "member")


