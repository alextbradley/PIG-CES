% attribution stuff

addpath('../functions/');
realizations = ["021","022", "023", "024" ,"025", "026", "027", "028", "029", ...
    "030", "031", "032", "033", "038", "039", "040"];

%% T-test on the ensemble mean trends

mean_trends = nan(1,length(realizations));

for ir = 1:length(realizations)

    %load data
    fname = strcat("emulate/mcmc_output/mcmc_output_realization", realizations(ir), ".csv");
    A = readmatrix(fname);

    kde = fitdist(A(:,4),'kernel'); %4th one is the trend
    mean_trends(ir) = mean(kde);
end

%do the t-test
[h,p,ci,stats] = ttest(mean_trends);

fprintf("Reject null hypothesis (trends not significantly different from zero) at the %.3g level \n", p);


%% Plot the posterior climate 
climates = nan(3001,length(realizations));
central_estimates = nan(length(realizations), 6);
for ir = 1:length(realizations)

    %get the mcmc output
    fname = strcat("emulate/mcmc_output/mcmc_output_realization", realizations(ir), ".csv");
    A = readmatrix(fname);


    for i = 1:6
        %construct the kde
        kde = fitdist(A(:,i),'kernel');
        central_estimates(ir,i) = mean(kde); 
    end

    %get the random forcing anomaly 
    fname = strcat('../model-inputs-and-outputs/realization', realizations(ir), "/realization.mat");
    data = load(fname);
    t = data.time;
    t = t + 1750;

    pc_random = data.pycnocline_center;

    %add the central components
    trend = central_estimates(ir,4);
    bump_amplitude = central_estimates(ir,5);
    bump_duration  = central_estimates(ir,6);
    
    pc_random = pc_random + trend/100* (t-1960) .* (t > 1960) +...
                 bump_amplitude*exp(-(t - 1945).^2 /2 /bump_duration^2);


    climates(:,ir) = pc_random;
end

figure(2);clf; hold on;
for ir = 1:length(realizations)
    p = plot(t, climates(:,ir));
    p.Color(4) = 0.3;  
end

plot(t, mean(climates, 2), 'k', 'linewidth',1.5)
plot([min(t), max(t)], [-500,-500], 'k--', 'linewidth', 1.5)


%fill in the std
%xf = [t; flip(t)];
%yf = [mean(climates, 2) - (std(climates')');flip(mean(climates, 2) + (std(climates')'))];
%fill(xf, yf, 'k','FaceAlpha',0.2, );

%% Plot histories of melting 


for it = 1:length(t)
    

end
