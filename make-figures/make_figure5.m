% Make figure 5 of the manuscript, showing spread of posterior
% distributions for climate and model parameters.
%
% 26/02/25 ATB alex.bradley@kcl.ac.uk, MIT licence.
%
%% Preliminaries
realizations = 1:14;


param_names  = ["basal sliding prefactor",...
    "ice viscosity prefactor",...
    "melt rate exponent prefactor",...
    "Anthropogenic trend since 1960s",...
    "1940s event magnitude",...
    "1940s event duration"];

% priors
prior_mean = [1.0, 1.0, 0.0, 0.0, 200.0, 5.0];
prior_std  = [0.3, 0.3, 1.2, 200.0,100, 2.5];

% setup figure
fig = figure(1); clf;
for i = 1:6
    ax(i) = subplot(2,3,i);
    hold(ax(i), 'on');
    box(ax(i), 'on');
    ax(i).FontSize = 14;
    ax(i).YTick = [realizations, 15,16]; %add for prior and mean posterior
    ax(i).YLim = [0.5, 16.5];
    ax(i).YTickLabels{15} = 'combined';
    ax(i).YTickLabels{16} = 'prior';
    ax(i).XLim = [prior_mean(i) - 2.1*prior_std(i), prior_mean(i) + 2.1*prior_std(i)];

    ax(i).XLabel.String = param_names(i);
end

barcol = [255,127,80]/255;



%% Get the MCMC data

%setup storage for the intervals
int_05 = nan(6,length(realizations));
int_25 = nan(6,length(realizations));
int_50 = nan(6,length(realizations));
int_75 = nan(6,length(realizations));
int_95 = nan(6,length(realizations));

all_samples = [];
for ir = 1:length(realizations)
    padded_realization = sprintf("%03d", realizations(ir));

    mcmc_data = readmatrix(strcat('../mcmc_output/mcmc_output_realization', padded_realization, '.csv'));
    all_samples = [all_samples;mcmc_data];

    for iP = 1:6 %loop over parameters
        Q = quantile(mcmc_data(:,iP),[0.05 0.25 0.5 0.75 0.95]);
        int_05(iP,ir) = Q(1);
        int_25(iP,ir) = Q(2);
        int_50(iP,ir) = Q(3);
        int_75(iP,ir) = Q(4);
        int_95(iP,ir) = Q(5);
    end


end %end loop over realizations

Q_all = nan(6, 5);
for iP = 1:6
    Q_all(iP,:) = quantile(all_samples(:,iP),[0.05 0.25 0.5 0.75 0.95]);
end



%% Make panels

tlw = 3;
wlw = 6; %wide line width
for iP = 1:6
    for ir = 1:length(realizations)
        % 95% interval
        plot(ax(iP), [int_05(iP,ir), int_95(iP, ir)], [realizations(ir), realizations(ir)], 'color',barcol,  'LineWidth',tlw);

        % IQR
        plot(ax(iP), [int_25(iP,ir), int_75(iP, ir)], [realizations(ir), realizations(ir)], 'color',barcol,  'LineWidth',wlw);


        %central estimate
        plot(ax(iP), int_50(iP,ir), realizations(ir), 'ko', 'MarkerFaceColor','w',  'MarkerSize',8, 'linewidth',3);

    end

    % priors: 95% interval
    plot(ax(iP),[(prior_mean(iP) - 2*prior_std(iP)), (prior_mean(iP) + 2*prior_std(iP))], [16,16], 'LineWidth',tlw, 'Color',[26,64,98]/255);


    % priors: IQR
    plot(ax(iP),[(prior_mean(iP) - 1.35*prior_std(iP)), (prior_mean(iP) + 1.35*prior_std(iP))], [16,16], 'LineWidth',wlw, 'Color',[26,64,98]/255);

    %prior central 
    plot(ax(iP),prior_mean(iP),16,'ko', 'MarkerFaceColor','w',  'MarkerSize',8, 'linewidth',3);

    % combined 95% interval
    plot(ax(iP), [Q_all(iP,1), Q_all(iP, end)], [15, 15],'color',[125, 180, 223]/255,  'LineWidth',tlw);

    % combined IQR
    plot(ax(iP), [Q_all(iP,2), int_75(iP, ir)], [15, 15], 'color',[125, 180, 223]/255,  'LineWidth',wlw);


    %combined central estimate
    plot(ax(iP), Q_all(iP,3), 15, 'ko', 'MarkerFaceColor','w',  'MarkerSize',8, 'linewidth',3);


end





%exportgraphics(fig, "figures/raw/figure6.pdf", 'ContentType','vector');


%% attribution stuff
[h,p] = ttest(int_50(4,:));

all_trend_sample = (all_samples(:,4));
frac_negative = sum(all_trend_sample<0)/length(all_trend_sample);