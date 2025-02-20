% create a plot of the posterior distributions of the parameters by
% averaging over kdes of individual realizations

addpath('../functions/');
realizations = ["021","022", "023", "024" ,"025", "026", "027", "028", "029", "030", "031", "032", "033", "038", "039", "040"];

input_headers = ["weertman_c_prefactor" ,"glen_a_ref_prefactor",...
    "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration"];

figure(1); clf;
vars = [1,2,3,4,5,6]; %skip the ungrounded weertman c

lims = [0.5, 2;
    0.5, 2;
    -2,  2;
    -400, 400;
    -400, 400;
    0, 10];  %limits of the kde


np = 1000;
xx = zeros(np,6);
dx = diff(xx);
dx = dx(1,:);
alphaval = 0.4;
colmap = parula(length(realizations));

for i = 1:6
    ax(i) = subplot(2,3,i);
    hold(ax(i), "on");
    box(ax(i), 'on');
    xx(:,i) = linspace(lims(i,1),lims(i,2), np);
    ax(i).XLabel.String = input_headers(i);
    ax(i).XLabel.Interpreter = 'none';

end

ss = nan(length(realizations), 6, np);
kdes = struct;
kde_mean = nan(length(realizations), 6);
kde_std  = nan(length(realizations), 6);

for ir = 1:length(realizations)
    fname = strcat("./emulate/mcmc_output/mcmc_output_realization", realizations(ir), ".csv");
    A = readmatrix(fname);

    for i = 1:6

        pts = xx(:,i);
        [f,xi] = ksdensity(A(:,vars(i)),pts);
        ss(ir, i, :) = f;

        kde = fitdist(A(:,i),'kernel');
        y = pdf(kde,pts);
        kdes(ir,i).kde = kde;

        kde_mean(ir,i) = mean(kde);
        kde_std(ir,i) = std(kde);
       

        p = plot(ax(i), xi, f, 'linewidth', 1.5);
        p.Color = [colmap(ir,:), alphaval];

    end
end


%% add ensemble means
f_mean = squeeze(mean(ss,1));
f_sd   = squeeze(std(ss,1));

for i = 1:6
    xf = [xx(:,i); flip(xx(:,i))];
    yf = [f_mean(i,:) + f_sd(i,:), flip(f_mean(i,:) - f_sd(i,:))];
    fill(ax(i),xf, yf, 'k', 'LineStyle','none', 'FaceAlpha',0.1);
    plot(ax(i), xx(:,i), f_mean(i,:), 'k', 'linewidth', 1.5);
end



%% add priors
prior_mean = [1.0, 1.0, 0.0,0.0 ,200.0, 5.0];
prior_sd   = [0.3, 0.3, 1.2, 200.0, 100.0, 2.5];
expfn = @(x,mu, sigma) (1/sqrt(2*pi*sigma^2)*exp(-(x-mu).^2 ./ 2 ./sigma^2));

for i = 1:6
    plot(ax(i), xx(:,i), expfn(xx(:,i),prior_mean(i), prior_sd(i)), 'k--', 'linewidth', 1.5);
    ax(i).YLim = [0, 5*max(expfn(xx(:,i),prior_mean(i), prior_sd(i)))];
end

%% what if every realization is gaussian..?
gaussian_posterior_mean = mean(kde_mean, 1);
gaussian_posterior_std  = std(kde_mean, 1);

for i = 1:6
    plot(ax(i), xx(:,i), expfn(xx(:,i),gaussian_posterior_mean(i), gaussian_posterior_std(i)), 'r--', 'linewidth', 1.5)
end


%% plot the cdfs%
figure(2); clf;

dx = diff(xx);
dx = dx(1,:);
for i = 1:6
    axs(i)=subplot(2,3,i); hold on; box on;
    for ir = 1:length(realizations)
        p = plot(axs(i), xx(:,i), cumsum(squeeze(ss(ir, i, :)))*dx(i), 'linewidth', 1.5);
        p.Color = [colmap(ir,:), alphaval];
    end
    plot(axs(i), xx(:,i), cumsum(f_mean(i,:))*dx(i), 'k', 'linewidth', 1.5);
    plot(axs(i), xx(:,i), cumsum(expfn(xx(:,i),prior_mean(i), prior_sd(i)))/max(cumsum(expfn(xx(:,i),prior_mean(i), prior_sd(i)))), 'k--', 'linewidth', 1.5);
    plot(axs(i), xx(:,i), cumsum(expfn(xx(:,i),gaussian_posterior_mean(i), gaussian_posterior_std(i)))/max(cumsum(expfn(xx(:,i),gaussian_posterior_mean(i), gaussian_posterior_std(i)))), 'r--', 'linewidth', 1.5);
    
    axs(i).XLabel.String = input_headers(i);
    axs(i).XLabel.Interpreter = 'none';
end



