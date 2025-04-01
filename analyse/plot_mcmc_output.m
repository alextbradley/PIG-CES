% Plot the MCMC output for a given realization

realization = "011";


% load the MCMC output
fname = strcat("./emulate/mcmc_output/mcmc_output_realization", realization, ".csv");
mcmc_output = readmatrix(fname);
colnames = ["weertman_c_prefactor", "glen_a_ref_prefactor", "melt_rate_prefactor_exponent", "per_century_trend", "bump_amplitude", "bump_duration"];

%
% make a plot of the MCMC output
%
make_plot = 0; 

kdes = struct;
figure(1);clf;
for i = 1:6
    ax(i) = subplot(2,3,i);

    histogram(ax(i), mcmc_output(:,i), 50, 'Normalization','pdf');
    ax(i).XLabel.String= colnames(i);
    ax(i).XLabel.Interpreter = 'none';
    ax(i).FontSize = 12;
    hold(ax(i), 'on')


    kde = fitdist(mcmc_output(:,i), 'Kernel', 'Kernel', 'normal'); % specify 'normal' or other kernel type
    kdes(i).kde = kde;

    % add the distribution of the kde for a sanity check
    numSamples = 500;
    samples = random(kde, numSamples, 1);
    xx = ax(i).XLim;
    xx = linspace(xx(1), xx(2), 100);
    y_values = pdf(kde, xx);
    plot(xx, y_values, 'k', 'LineWidth', 2); % KDE
    histogram(samples, 'Normalization', 'pdf', 'FaceAlpha',0.2); % sampled data histogram

end

% add the naive posterior
params = readmatrix(strcat("../model-inputs-and-outputs/realization", realization, "/iteration005/params.csv"));
params = params( :,[1,3:7]); %remove unused ungrounded weertman c
naive_posterior_mean = mean(params,1);
naive_posterior_std  = std(params,1);
expfn = @(x,mu,sigma) 1/sqrt(2*pi*sigma^2)*exp(-(x-mu).^2 / 2/ sigma^2);

for i = 1:6
    xx = ax(i).XLim;
    xx = linspace(xx(1), xx(2), 100);
    hold(ax(i),'on')
    plot(ax(i), xx, expfn(xx,naive_posterior_mean(i), naive_posterior_std(i)),'k--', 'linewidth', 1.5);
end









