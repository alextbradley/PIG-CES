%Generate forward run samples

realizations = ["011", "012", "013", "014", "015", "016", "017", "018", "019",...
                "020", "021", "022", "023", "024", "025", "026", "027", "028",...
                "029", "030", "031", "032", "033", "036", "038", "039", "040"];


trend2 = 100;
trend2 = '1960_to_2015'; %continue the trend from 1960 to 2015 



rng(10)
numSamples = 5; %number of samples from each realization

params = nan(numSamples*length(realizations),8); %six model parameters, random seed and trend in forcing
count = 1;
for ir = 1:length(realizations)

    %load the mcmc output
    fname = strcat("./emulate/mcmc_output/mcmc_output_realization", realizations(ir), ".csv");
    mcmc_output = readmatrix(fname);

    %loop over each of the mcmc outputs
    for i = 1:6
        %generate kde
        kde = fitdist(mcmc_output(:,i), 'Kernel', 'Kernel', 'normal'); % specify 'normal' or other kernel type
    
        %sample from the kde
        samples = random(kde, numSamples, 1);

        %store 
        params(count:(count + numSamples-1),i) = samples;

    end

    params(count:(count + numSamples-1), 7) = str2num(realizations(ir));

    %work out the trend in forcing
    if strcmp(trend2, '1960_to_2015')
        %load the realization of forcing
        rf_data = load(strcat('../model-inputs-and-outputs/realization',realizations(ir),'/realization.mat'));

        %compute the trend linear trend from 1960 (210) to 2015 (265)
        [~,idx1] = min(abs(rf_data.time - 210)); %index of start
        [~,idx2] = min(abs(rf_data.time - 265)); %index of end
        vv = rf_data.pycnocline_center; 
        vv = vv(idx1:idx2); %pycnocline centres over 1960 to 2015
        tt = rf_data.time;
        tt = tt(idx1:idx2);

        %get the linear trend 
        p = polyfit(tt,vv,1); %p(1) is the natural trend per year
        natural_trend = p(1)*100;

        %store the trend as the trend in the sample minus the natural trend
        anthro_trend = params(count:(count + numSamples-1), 4);
        params(count:(count + numSamples-1), 8) = anthro_trend - natural_trend; 

    else
       params(count:(count + numSamples-1), 8) = trend2;

    end

    count = count + numSamples;


end

% %% output results
% for i= 1:length(params)
% fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
% fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',1.0 )
% fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,2) )
% fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,3) )
% fprintf('        per_century_trend: %.4f \n',params(i,4) )
% fprintf('        bump_amplitude: %.4f \n',params(i,5) )
% fprintf('        bump_duration: %.4f \n',params(i,6) )
% fprintf('        random_seed: %02d \n',params(i,7) )
% fprintf('        per_century_trend2: %.4f \n',params(i,4) )
% end


%% output results
for i= 1:length(params)
fprintf('      - weertman_c_prefactor: %.4f \n',params(i,1) )
fprintf('        ungrounded_weertmanC_prefactor: %.4f \n',1.0 )
fprintf('        glen_a_ref_prefactor: %.4f \n',params(i,2) )
fprintf('        melt_rate_prefactor_exponent: %.4f \n',params(i,3) )
fprintf('        per_century_trend: %.4f \n',0.0 ) %fprintf('        per_century_trend: %.4f \n',params(i,4) )
fprintf('        bump_amplitude: %.4f \n',0.0 )
fprintf('        bump_duration: %.4f \n',params(i,6) )
fprintf('        random_seed: %02d \n',params(i,7) )
fprintf('        per_century_trend2: %.4f \n',params(i,4) )
end
warning("above is no bump and no trend")