include("shared.jl")

"""
    Run this script to generate the folder structure and initialize the eki for different realizations of forcing. 
    In theory, this only needs to be run once! Must be in main directory to do so. 

"""
function initialize_EKI(realization, n_ensemble)

    # check we're running from the directory
    if pwd() != dirname(@__FILE__)
        println("Please run from main directory")
        return
    end

    # check that this folder does not already exist. If so, terminal and do nothing
    padded_realization = @sprintf("%03d", realization);
    folder_path = "./model-inputs-and-outputs/realization" * string(padded_realization)
    iter_path   = folder_path * "/iteration001"
    if isdir(folder_path)
        println("Folder exists: ", folder_path, ", terminating")
        return  # Exit the function early
    else
        mkdir(folder_path)
        mkdir(iter_path)
    end

    # read in the priors
    toml_path = "./priors.toml";
    param_dict = TOML.parsefile(toml_path)
    names = ["weertman_c_prefactor", "ungrounded_weertmanC_prefactor", "glen_a_ref_prefactor",  "melt_rate_prefactor_exponent", "per_century_trend","bump_amplitude","bump_duration"] #must match the priors names
    prior_vec = [get_parameter_distribution(param_dict, n) for n in names]
    prior = combine_distributions(prior_vec)

    #get the initial ensemble
    rng_seed = realization #use a different random seed for each realization of forcing? Or is it easier to use the same?
    rng_ekp = Random.MersenneTwister(rng_seed)
    initial_ensemble = EKP.construct_initial_ensemble(rng_ekp, prior, n_ensemble)
    
    # set up the folder structure
    for member = 1:n_ensemble
        member_padded = @sprintf("%03d", member);
        mkdir(iter_path * "/member" * member_padded)

        #write a tmp file so that git knows about folder structure
        path = iter_path*"/member" * member_padded*"/tmp.txt"
        open(path, "w") do file
            write(file, "temp file\n")
        end
    end

    # save the parameter values in txt format
    df = DataFrame(initial_ensemble', names)
    csv_path = iter_path * "/params.csv"
    CSV.write(csv_path, df)

    # load in the data
    @load "./observations/truth.jld2" y Γ

    # generate the EKI 
    eki = EKP.EnsembleKalmanProcess(initial_ensemble, vec(y), Γ, Inversion(); rng = rng_ekp)

    # save the EKI
    eki_path = iter_path * "/eki.jld2"
    @save eki_path eki param_dict prior

    return nothing
end

#realization = 1:20
#n_ensemble  = 10
#initialize_EKI.(realization, n_ensemble)