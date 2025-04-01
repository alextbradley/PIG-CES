% Output the latin hypercube info for the sytem
save_as_csv = 0; 
fname       = 'lhc.csv';

% Number of samples
numSamples = 100;

% Number of variables
numVariables = 6;

% Priors from EKI
prior_means = [1,1, 0, 200, 5, 0]; 
prior_stds  = [0.3, 0.3, 1.2, 100.0, 2.5, 200];

variableRanges = nan(6,2);
for i = 1:6
    variableRanges(i,:) = prior_means(i) + 2*[-prior_stds(i), prior_stds(i)];
end


rng(0) %set random seed

% Generate Latin Hypercube Sample

lhsSample = lhsdesign(numSamples, numVariables, 'criterion', 'maximin');

% Rescale the LHS samples to match variable ranges
for i = 1:numVariables
    lhsSample(:, i) = lhsSample(:, i) * (variableRanges(i, 2) - variableRanges(i, 1)) + variableRanges(i, 1);
end

% Display the Latin Hypercube Sample
%disp('Latin Hypercube Sample:');
%disp(lhsSample);
params = lhsSample;
%%
for i= 1:length(params)
    fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
    fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',1.0 )
    fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,2) )
    fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,3) )
    fprintf('        per_century_trend: %.4f \n',params(i,4) )
    fprintf('        bump_amplitude: %.4f \n',params(i,5) )
    fprintf('        bump_duration: %.4f \n',params(i,6) )
    fprintf('        random_seed: %02d \n',21 )
end
