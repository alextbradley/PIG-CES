# for playing with the CES stuff
using JLD2
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample
using Printf
using Plots
const EKP = CalibrateEmulateSample.EnsembleKalmanProcesses
const PD = EKP.ParameterDistributions
const CES = CalibrateEmulateSample


    # check we're running from the directory
if pwd() != dirname(@__FILE__)
    println("changing directory to ", dirname(@__FILE__))
    cd(dirname(@__FILE__))
    return
end

#load the EKI
N_iterations = 5
realization  = 26
padded_iteration = @sprintf("%03d", N_iterations+1); #need to ensure you've done the final update
padded_realization = @sprintf("%03d", realization);
eki_path =  "../../model-inputs-and-outputs/realization" * padded_realization * "/iteration" * padded_iteration * "/eki.jld2"
@load eki_path eki param_dict prior

# get the error covariance from the obs_mean
obs_path = "../../observations/truth.jld2"
@load obs_path Γ

input_output_pairs = CES.Utilities.get_training_points(eki, N_iterations)
unconstrained_inputs = CES.Utilities.get_inputs(input_output_pairs)
inputs = Emulators.transform_unconstrained_to_constrained(prior, unconstrained_inputs)
outputs = CES.Utilities.get_outputs(input_output_pairs)

####### emulation 
gppackage = Emulators.GPJL()
gauss_proc = Emulators.GaussianProcess(gppackage, noise_learn = false)
emulator_gp = Emulator(gauss_proc, input_output_pairs, normalize_inputs = true,  obs_noise_cov = Γ)
optimize_hyperparameters!(emulator_gp)


###### reduced data
eki_copy = deepcopy(eki)

# loop over the iterations, and for each store only the last row
outputs_2015_grv_error = ones(1,N_iterations*eki_copy.N_ens) #this will hold only the last rows

count = 1
for i = 1:N_iterations
    outputs_2015_grv_error[:,count:(count+eki_copy.N_ens  - 1)] = eki.g[i].stored_data[end,:]
    count = count + eki_copy.N_ens
end

outputs_2015_grv_error_data_container = EnsembleKalmanProcesses.DataContainers.DataContainer(outputs_2015_grv_error)
input_output_pairs_onlygrv = EnsembleKalmanProcesses.DataContainers.PairedDataContainer(EnsembleKalmanProcesses.DataContainers.DataContainer(inputs),outputs_2015_grv_error_data_container)

# make new error covariance? Maybe ok because will be same format?

#try random features?
input_dim = 7
output_dim = 1
# Select number of features
n_features = 60
nugget = 1e-9
kernel_structure = NonseparableKernel(LowRankFactor(2, nugget))
optimizer_options = Dict(
    "n_ensemble" => 50,
    "cov_sample_multiplier" => 10,
    "scheduler" => EKP.DataMisfitController(on_terminate = "continue"),
    "n_iteration" => 50,
    "verbose" => true,
)
random_features = VectorRandomFeatureInterface(
    n_features,
    input_dim,
    output_dim,
    kernel_structure = kernel_structure,
    optimizer_options = optimizer_options,
)
emulator_random_features =
    Emulator(random_features, input_output_pairs_onlygrv, normalize_inputs = true, obs_noise_cov = Γ, decorrelate = false)
optimize_hyperparameters!(emulator_random_features)