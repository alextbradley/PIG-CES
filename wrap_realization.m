% Wrap up all of the outputs for a given realization into a single csv
% file. This code produces an output file names "model_simulations.csv",
% which is saved in the appropriate realization folder, i.e. at
% "./realizationXXX/model_simulations.csv"

realization = "002";
realization_folder = strcat('realization', realization);
realization_dir = dir(realization_folder);

outputs = [];
parameters = [];

for i = 1:length(realization_dir) %for each folder in the directory
    iteration_folder = realization_dir(i).name;
    prefix = "iteration";
    if strncmp(iteration_folder, prefix, length(prefix)) %filter only for iterations by their prefix
        iteration_dir = dir(fullfile(".", realization_folder, iteration_folder));

        for j = 1:length(iteration_dir)

            %load the parameter file
            params_iteration = readmatrix(fullfile(".", realization_folder, iteration_folder, "params.csv"));
            
            member_folder = iteration_dir(j).name;

            
            if strncmp(member_folder, "member", length("member")) %filter only for members based on their prefix
                output_path = fullfile(".", realization_folder, iteration_folder, member_folder, "outputs.csv"); 
                
                %append this output and the corresponding parameters
                if exist("output_path")
                    outputs = [outputs; csvread(output_path)];                  
                    member_number = str2num(member_folder(end-2:end));
                    parameters = [parameters; params_iteration(member_number,:)];
                end %end file exists flag
            end %end member prefix flag           
        end %end loop over iteration directory
    end
end

%glue them together
data = [parameters, outputs];
data = array2table(data, 'VariableNames', {'weertman_c_prefactor', 'ungrounded_weertmanC_prefactor', 'glen_a_ref_prefactor', 'melt_rate_prefactor_exponent', 'per_century_trend', 'bump_amplitude', 'bump_duration', 'gl_position_1930', 'gl_position_2015'});

%save the output
save_path = fullfile(realization_folder, "model_simulations.csv");
writetable(data, save_path);