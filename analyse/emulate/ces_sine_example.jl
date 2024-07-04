#running thru the CES example 
using LinearAlgebra, Random
using JLD2
using Statistics, Distributions
using Plots
# CES
using CalibrateEmulateSample
using CalibrateEmulateSample.Emulators
const CES = CalibrateEmulateSample
const EKP = CalibrateEmulateSample.EnsembleKalmanProcesses
const PD = EKP.ParameterDistributions


# Define x-axis
dt = 0.01
trange = 0:dt:(2 * pi + dt)

function model(amplitude, vert_shift)
    # Set phi, the random phase
    phi = 2 * pi * rand()
    return amplitude * sin.(trange .+ phi) .+ vert_shift
end

amplitude_true = 3.0
vert_shift_true = 7.0
# Our input parameters are 2d and our outputs are 2d
theta_true = [amplitude_true, vert_shift_true]
dim_params = 2
# Generate the "true" signal for these parameters
signal_true = model(amplitude_true, vert_shift_true)

y1_true = maximum(signal_true) - minimum(signal_true)
y2_true = mean(signal_true)

dim_output = 2
Γ = 0.2 * I
white_noise = MvNormal(zeros(dim_output), Γ)
y_obs = [y1_true, y2_true] .+ rand(white_noise)
println("Observations:", y_obs)

function G(theta)
    amplitude, vert_shift = theta
    sincurve = model(amplitude, vert_shift)
    return [maximum(sincurve) - minimum(sincurve), mean(sincurve)]
end

##### Calibrate 
prior_u1 = PD.constrained_gaussian("amplitude", 2, 1, 0, Inf)
prior_u2 = PD.constrained_gaussian("vert_shift", 0, 5, -Inf, Inf)
prior = PD.combine_distributions([prior_u1, prior_u2])
# Plot priors
p = plot(prior, fill = :lightgray)

#set up eki
N_ensemble = 10
N_iterations = 5
initial_ensemble = EKP.construct_initial_ensemble(prior, N_ensemble)
ensemble_kalman_process = EKP.EnsembleKalmanProcess(initial_ensemble, y_obs, Γ, EKP.Inversion())

for i in 1:N_iterations
    params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
    G_ens = hcat([G(params_i[:, i]) for i in 1:N_ensemble]...)
    EKP.update_ensemble!(ensemble_kalman_process, G_ens)
end

final_ensemble = EKP.get_ϕ_final(prior, ensemble_kalman_process)

# Check that the ensemble mean is close to the theta_true
println("Ensemble mean: ", mean(final_ensemble, dims=2))   # [3.05, 6.37]
println("True parameters: ", theta_true)   # [3.0, 7.0]


input_output_pairs = CES.Utilities.get_training_points(ensemble_kalman_process, N_iterations)
unconstrained_inputs = CES.Utilities.get_inputs(input_output_pairs)
inputs = Emulators.transform_unconstrained_to_constrained(prior, unconstrained_inputs)
outputs = CES.Utilities.get_outputs(input_output_pairs)

gppackage = Emulators.GPJL()
gauss_proc = Emulators.GaussianProcess(gppackage, noise_learn = false)

# Build emulator with data
emulator_gp = Emulator(gauss_proc, input_output_pairs, normalize_inputs = true,  obs_noise_cov = Γ)
optimize_hyperparameters!(emulator_gp)