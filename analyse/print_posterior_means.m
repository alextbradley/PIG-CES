% output the posterior mean parameters for use in the ensembles, including
% with (1) no trend in forcing and (2) with no bump

realizations = ["021","022", "023", "024" ,"025", "026", "027", "028", "029", "030", "031", "032", "033", "038", "039", "040"];


params = zeros(length(realizations),8);
params(:,2) = 1.0; %fix the ungrounded weertman c
pars = [1,3,4,5,6,7]; %tell the code to skip 2nd entry (the ungrounded weertman c)

for ir = 1:length(realizations)

    %load data
    fname = strcat("emulate/mcmc_output/mcmc_output_realization_lowN", realizations(ir), ".csv");
    A = readmatrix(fname);

    for i = 1:6

        %kde = fitdist(A(:,i),'kernel'); 
        %params(ir, pars(i)) = mean(kde); %to use mean of the kde, not max

        [density, x_values] = ksdensity(A(:,i));
        [max_density, idx] = max(density);
        params(ir, pars(i)) = x_values(idx);

    end
    params(ir, 8) = str2num(realizations(ir));
    
end


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

%% repeat this with zero trend
for i= 1:length(params)
    fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
    fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',params(i,2) )
    fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,3) )
    fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,4) )
    fprintf('        per_century_trend: %.4f \n',0 )
    fprintf('        bump_amplitude: %.4f \n',params(i,6) )
    fprintf('        bump_duration: %.4f \n',params(i,7) )
    fprintf('        random_seed: %02d \n',params(i,8) )
end

%% repeat this with zero bump
for i= 1:length(params)
    fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
    fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',params(i,2) )
    fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,3) )
    fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,4) )
    fprintf('        per_century_trend: %.4f \n',params(i,5) )
    fprintf('        bump_amplitude: %.4f \n', 0 )
    fprintf('        bump_duration: %.4f \n',params(i,7) )
    fprintf('        random_seed: %02d \n',params(i,8) )
end

%% repeat this with zero bump or trend
for i= 1:length(params)
    fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
    fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',params(i,2) )
    fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,3) )
    fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,4) )
    fprintf('        per_century_trend: %.4f \n',0 )
    fprintf('        bump_amplitude: %.4f \n', 0 )
    fprintf('        bump_duration: %.4f \n',params(i,7) )
    fprintf('        random_seed: %02d \n',params(i,8) )
end