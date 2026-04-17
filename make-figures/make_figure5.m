% Make figure 5 of the manuscript, showing spread of posterior
% distributions for climate and model parameters.
%
% 26/02/25 ATB alex.bradley@kcl.ac.uk, MIT licence.
%
%% Preliminaries
realizations = 1:14;


param_names  = ["sliding prefactor",...
    "viscosity prefactor",...
    "melt prefactor",...
    "per century trend",...
    "1940s event amplitude",...
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
    grid on
    ax(i).XLabel.String = param_names(i);
end

combcol = [255,127,80]/255;
barcol = [125, 180, 223]/255;

fig = gcf;
fig.Position(3:4) = [1300, 660];

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
    plot(ax(iP), [Q_all(iP,1), Q_all(iP, end)], [15, 15],'color',combcol,  'LineWidth',tlw);

    % combined IQR
    plot(ax(iP), [Q_all(iP,2), Q_all(iP, end-1)], [15, 15], 'color',combcol,  'LineWidth',wlw);


    %combined central estimate
    plot(ax(iP), Q_all(iP,3), 15, 'ko', 'MarkerFaceColor','w',  'MarkerSize',8, 'linewidth',3);


end





%exportgraphics(fig, "figures/raw/figure6.pdf", 'ContentType','vector');


%% attribution stuff
[h,p] = ttest(int_50(4,:));

all_trend_sample = (all_samples(:,4));
frac_negative = sum(all_trend_sample<0)/length(all_trend_sample);

%% Make the posterior forcing for trend
fig2 = figure(2); clf; hold on;
fig2.Position(3:4) = [600, 250]; box on
ax(3) = gca;
% add the prior and stuff from figure 2 c
realizations = 1:14;


for ir = 1:length(realizations)
    fpath = strcat("../model-inputs-and-outputs/realization",  sprintf('%03d', realizations(ir)), "/realization.mat");
    rr = load(fpath);

    pc = rr.pycnocline_center;
    time = rr.time; 

    %restrict a bit for each of plotting
    pc = pc(1:10:end);
    time = time(1:10:end);
    time = time + 1750;

    if ir == 1
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5,'Marker','none');
    else        
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5, "HandleVisibility","off", 'Marker', 'none');

    end

    pp.Color(4) = 0.1; %set the alpha

    if ir == 1
        ens_mean = pc;
    else
        ens_mean = ens_mean + pc;
    end

end

ens_mean = ens_mean/length(realizations);
plot(ax(3), time, ens_mean, 'k', 'linewidth', 2);
ax(3).YLim = [-650, -250];
ax(3).YTick = -600:100:-200;
ax(3).YLabel.String = 'pycnocline depth (m)';
ax(3).FontSize = 14;
ax(3).XLabel.String = 'year';
ax(3).YLabel.String = 'pycnocline depth (m)';

% add the prior shading
prior95 = (time > 1960).*(time-1960)*2*200/100;
prior75 = (time > 1960).*(time-1960)*1.73*200/100;
prior50 = (time > 1960).*(time-1960)*0*200/100;
prior25 = (time > 1960).*(time-1960)*-1.73*200/100;
prior05 = (time > 1960).*(time-1960)*-2*200/100;

xf05 = [time; flip(time)];
yf05 = -500 + [prior05; flip(prior95)];
yfiqr = -500 + [prior25; flip(prior75)];
fill(xf05, yf05, [26,64,98]/255, 'FaceAlpha', 0.2, 'Linestyle', 'none');
%fill(xf05, yfiqr, [26,64,98]/255, 'FaceAlpha', 0.3, 'Linestyle', 'none');


% compute the components for the trend
trend95 = (time > 1960).*(time - 1960)*Q_all(4,end)/100; %95% interval
trend75 = (time > 1960).*(time - 1960)*Q_all(4,end-1)/100; %75% interval
trend50 = (time > 1960).*(time - 1960)*Q_all(4,3)/100; %75% interval
trend25 = (time > 1960).*(time - 1960)*Q_all(4,2)/100; %75% interval
trend05 = (time > 1960).*(time - 1960)*Q_all(4,1)/100; %75% interval


yf05 = -500 + [trend05; flip(trend95)];
yfiqr = -500 + [trend25; flip(trend75)];
fill(xf05, yf05, combcol, 'FaceAlpha', 0.2, 'Linestyle', 'none');
%fill(xf05, yfiqr, barcol, 'FaceAlpha', 0.3, 'Linestyle', 'none');
plot(ax(3), time,-500+trend50,'Color',combcol, 'LineWidth',2)
plot(ax(3), time,-500+prior50,'Color',[26,64,98]/255, 'LineWidth',2)

%% Make the posterior forcing for bump
fig2 = figure(3); clf; hold on;
fig2.Position(3:4) = [600, 250]; box on
ax(3) = gca;
% add the prior and stuff from figure 2 c
realizations = 1:14;


for ir = 1:length(realizations)
    fpath = strcat("../model-inputs-and-outputs/realization",  sprintf('%03d', realizations(ir)), "/realization.mat");
    rr = load(fpath);

    pc = rr.pycnocline_center;
    time = rr.time; 

    %restrict a bit for each of plotting
    pc = pc(1:10:end);
    time = time(1:10:end);
    time = time + 1750;

    if ir == 1
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5,'Marker','none');
    else        
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5, "HandleVisibility","off", 'Marker', 'none');

    end

    pp.Color(4) = 0.1; %set the alpha

    if ir == 1
        ens_mean = pc;
    else
        ens_mean = ens_mean + pc;
    end

end

ens_mean = ens_mean/length(realizations);
plot(ax(3), time, ens_mean, 'k', 'linewidth', 2);
ax(3).YLim = [-650, -250];
ax(3).YTick = -600:100:-200;
ax(3).YLabel.String = 'pycnocline depth (m)';
ax(3).FontSize = 14;
ax(3).XLabel.String = 'year';
ax(3).YLabel.String = 'pycnocline depth (m)';

% add the prior shading
prior95 =  (200 + 2*100) * exp(-(time - 1945).^2 / 2 / (5.0 + 2*2.5)^2);
prior75 =  (200 + 1.73*100) * exp(-(time - 1945).^2 / 2 / (5.0 + 1.73*2.5)^2);
prior50 =  200 * exp(-(time - 1945).^2 / 2 / 5.0^2);
prior25 =  (200 - 1.73*100) * exp(-(time - 1945).^2 / 2 / (5.0 - 1.73*2.5)^2);
prior05 =  max((200 - 2*100) * exp(-(time - 1945).^2 / 2 / 5.0^2),0);



xf05 = [time; flip(time)];
yf05 = -500 + [prior05; flip(prior95)];
yfiqr = -500 + [prior25; flip(prior75)];
fill(xf05, yf05, [26,64,98]/255, 'FaceAlpha', 0.2, 'Linestyle', 'none');
%fill(xf05, yfiqr, [26,64,98]/255, 'FaceAlpha', 0.3, 'Linestyle', 'none');


% compute the components for the trend
trend95 =  Q_all(5,end) * exp(-(time - 1945).^2 / 2 / (5.0 + 2*2.5)^2);
trend75 =  Q_all(5,end-1) * exp(-(time - 1945).^2 / 2 / (5.0 + 1.73*2.5)^2);
trend50 =  Q_all(5,end-2) * exp(-(time - 1945).^2 / 2 / 5.0^2);
trend25 =  Q_all(5,end-3) * exp(-(time - 1945).^2 / 2 / (5.0 - 1.73*2.5)^2);
trend05 =  max(Q_all(5,end-4) * exp(-(time - 1945).^2 / 2 / 5.0^2),0);

yf05 = -500 + [trend05; flip(trend95)];
yfiqr = -500 + [trend25; flip(trend75)];
fill(xf05, yf05, combcol, 'FaceAlpha', 0.2, 'Linestyle', 'none');
%fill(xf05, yfiqr, barcol, 'FaceAlpha', 0.3, 'Linestyle', 'none');
plot(ax(3), time,-500+trend50,'Color',combcol, 'LineWidth',2)
plot(ax(3), time,-500+prior50,'Color',[26,64,98]/255, 'LineWidth',2)