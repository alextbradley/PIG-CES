# Code and Data for "Detection and attribution of anthropogenic climate change in retreat of the Pine Island Glacier"

This repository contains code (to run simulations and produce figures) and data (from model simulations) for analysis contained in "Detection and attribution of anthropogenic climate change in retreat of the Pine Island Glacier" by Bradley et al. 

Note that full simulation data is too large to be hosted here, and can be obtained on request from the corresponding author (alex.bradley@kcl.ac.uk). Simulation output for the three observational constraints (1930s grounding line position, 2015 grounding line position, and 2015 grounded ice volume) are hosted here.


The folder `model-inputs-and-outputs/realizationXXX/iterationYYY/memberZZZ` folder contains csv files `outputs.csv` and `outputs_cts.csv`, which contain the 3x1 model output for the the XXXth realization, YYYth iteration and ZZZth member of this iteration. These files are identical except that `outputs.csv` has the output when the grounding line position is measured in a discrete way (with the grounding line taken as the first not-fully grounded cell) and `outputs_cts.csv` has the output when the grounding line is measured as the weighted average of not-fully grounded cells, as outlined in the paper. Note that these outputs are in dimensionless form, they are:  
     * the number of grid cells away from the observations for the 1930 grounding line position  
     * the number of grid cells away from the observations for the 2015s grounding line position  
     * the difference between the modelled and observed grounded volume divided by 10^12 
     
This folder also contains a matlab file "output_trajectory.mat", which contains:  
     * temporal evolution of the grounding line position (discrete), called `gl_pos_discrete`  
     * temporal evolution of the grounding line position (continuous), called `gl_pos_cts`  
     * temporal evolution of the grounded volume, called `grv`  
     * corresponding time trajectory, in years after 1750, called `t`  
This files are generated using the scripts in the "update" folder, in combination with full model output (hosted elsewhere).

The 'model-inputs-and-outputs' folder structure is initialized using the 'initialize.jl' script:
```julia
include("initialize.jl")
initialize_EKI.(realizations, n_ensemble)
```
where `realizations` is an array of realization numbers and `n_enemble` is an integer number of ensemble members per iteration. You only have to run this once. 

For the EKI, we build upon the EnsembleKalmanProcesses library (https://github.com/CliMA/EnsembleKalmanProcesses.jl). The script `initialize.jl` script initializes the EKI. With model outputs computed, the EKI is updated with `update_EKI.jl`. The state of the Ensemble Kalman Process for the XXXth realization of forcing and YYYth iteration is stored as `model-inputs-and-outputs/realizationXXX/iterationYYY/eki.jld2`.

The `analyse/emulate` folder contains the full pipeline for the emulation and sampling steps of the CES, see `analyse/emulate/emulate_sample_pipeline.R`. 

The `observations` folder contains information on the observations of ice sheet retreat. This folder contains the following files used in the CES:  
     * `truth.csv`: a vector of zeros, which is the value of the model output if it perfectly matches observations (the CES is run is dimensionless form, so we use this file as observations)  
     * `truth.jld2`: as in `truth.csv` but in julia native format  
     * `truth_actual.csv`: actual values of observed values, where the grounding line position is measured relative to a polar stereographic grid, and grounded volume is in m^3  
     * `truth_actual_cts.csv`: as in `truth_actual.csv`, but with the grounding line position as measured using the continuous method  
     * `noise.csv`: dimensionless observation noise  
     * `noise_actual`: actual observation noise (note that this file says that the noise in the grounded volume is 10^12, but we override this to be 1% of observed grounded volume, see line 46 of `emulate_sample_pipeline.R`.  
     






- Priors.toml stores the information on the priors, which are the same for each realization of forcing.

- observations/truth.jld2 contains the observations, used to update the EKI:
     Γ = 1.0 * 3*10.0^3 * I #noisy observation of grounding line position, where I is a 2 x 2identity matrix
     y is a 2 x 1 noisy observation of grounding line position

- results are stored in /realizationXXX/iterationYYY/memberZZZ, with the eki stored in /realizationXXX/iterationYYY/eki.jld2
