% output a set of params.csv files in the correct format to be used by the
% ensembler

realizations = 11:15;
iterations   = 1;

params = [];
for ir = 1:length(realizations)
    for ii = 1:length(iterations)
        %pad the inputs to get the file
        padded_realization = sprintf('%03d', realizations(ir));
        padded_iteration   = sprintf('%03d', iterations(ii));

        %find params file
        fpath = strcat("./", "realization", padded_realization, "/iteration", padded_iteration , "/params.csv");
        
        if isfile(fpath)
            arr = table2array(readtable(fpath));

            %put the random seed on the end, which is the realization
            sz = size(arr);
            for is = 1:sz(1)
                arr(is,8) = str2num(padded_realization);
   
            end
            params = [params;arr];
        else
            disp(strcat('did not find a params file at', fpath, ", skipping..."))
        end

    end
end

% loop over the entries of params and print
for i= 1:length(params)
    fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
    fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',params(i,2) )
    fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,3) )
    fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,4) )
    fprintf('        per_century_trend: %.4f \n',params(i,5) )
    fprintf('        bump_amplitude: %.4f \n',params(i,6) )
    fprintf('        bump_duration: %.4f \n',params(i,7) )
    fprintf('        random_seed: %02d \n',params(i,8) )
end
