include("shared.jl")


function update_EKI(realization,completed_iteration)
    """
    update the EKI and prepare the parameters csv for the next iteration
    """

    padded_realization = @sprintf("%03d", realization);
    padded_completed_iteration = @sprintf("%03d", completed_iteration);
    eki_path =  "./realization" * padded_realization * "/iteration" * padded_completed_iteration * "/eki.jld2"

    @load eki_path eki param_dict prior
    N_ensemble = eki.N_ens
    dim_output = size(eki.obs_mean)[1]

    #loop over the model runs, check the files exist for them all and, if so, pull the results
    G_ens = zeros(dim_output, N_ensemble)
    for member in 1:N_ensemble
        padded_member = @sprintf("%03d", member);
        output_path =  "./realization" * padded_realization * "/iteration" * padded_completed_iteration * "/member" * padded_member *"/outputs.csv"
        if isfile(output_path)
        model_output = Matrix(CSV.read(output_path, DataFrame, header = false))'
        G_ens[:, member] = model_output
        else
        println("The file '$output_path' does not exist.")
        return
        end
    end

    # perform the update
    EKP.update_ensemble!(eki, G_ens)

    # generate the folder structure
    padded_next_iteration = @sprintf("%03d", completed_iteration+1);
    iter_path = "./realization" * padded_realization * "/iteration" * padded_next_iteration 
    mkdir(iter_path)

    for member = 1:N_ensemble
        member_padded = @sprintf("%03d", member);
        mkdir(iter_path * "/member" * member_padded)
    end


    # generate the next params file
    parameters_next_iteration = eki.u[completed_iteration + 1].stored_data

    # save the parameter values in csv format
    colnames = ["weertman_c_prefactor", "ungrounded_weertmanC_prefactor", "glen_a_ref_prefactor",  "melt_rate_prefactor_exponent", "per_century_trend","bump_amplitude","bump_duration"] #must match the priors names
    df = DataFrame(parameters_next_iteration', colnames)
    csv_path = iter_path * "/params.csv"
    CSV.write(csv_path, df)

    # save the EKI
    eki_path = iter_path * "/eki.jld2"
    @save eki_path eki param_dict prior    

end