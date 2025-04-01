% Make figure 6 of the manuscript:  reconstructed grounding line retreat
% for each of the realizations of forcing.

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license. 

%% Preliminaries
realizations = 1:14;
samples      = 1:10; 
outpath = "../model-inputs-and-outputs/forward-runs/default/";
init_time = 1750; 
plot_colour = [0, 47, 167]/255;
fig =figure(1); clf; fig.Position(3:4) = [1400, 560];

%% Pull data
ss = zeros(length(realizations), length(samples), 276);
for ir = 1:length(realizations)
    for is = 1:length(samples)
        fpath = strcat(outpath, "realization_", sprintf('%03d', realizations(ir)),"/sample_", sprintf('%03d',samples(is)), "/output_trajectory.mat");
        data = load(fpath);
        gl_retreat = data.gl_pos_cts;
        gl_retreat = (gl_retreat - gl_retreat(1))/1e3;% turn into retreat in km 
        ss(ir,is,:) = gl_retreat;

    end
end
%% get the obs
observations = readmatrix("../observations/truth_actual.csv");
observations = observations(1:2); %remove the gl constraint
observations = (observations - observations(1))/1e3; % turn into retreat in km 

obs_noise    = readmatrix("../observations/noise_actual.csv");
obs_noise    = obs_noise(1:2)/1e3;

obs_time     = readmatrix("../observations/truth_times.csv");
obs_time     = obs_time + init_time;
%% Make the plot
tt = tiledlayout(3,5);
time = data.t;
time = time + init_time;


for ir = 1:length(realizations)
    ens_mean_for_this_realization = squeeze(mean(ss(ir,:,:),2));
    ens_std_for_this_realization = squeeze(std(ss(ir,:,:),0,2));

    xf = [time; flip(time)];
    yf = [(ens_mean_for_this_realization - ens_std_for_this_realization);flip(ens_mean_for_this_realization + ens_std_for_this_realization) ];

    nexttile

    %plot the ens mean and spread
    fill(xf, yf, plot_colour, 'FaceAlpha',0.1, 'LineStyle','none');
    hold on
    plot(time, ens_mean_for_this_realization, 'color', plot_colour, 'LineWidth',2)

    %plot the observations
    for i = 1:2
        plot(obs_time(i), observations(i), 'ko', 'markersize', 10, 'markerfacecolor', 'k');
        plot(obs_time(i)*[1,1], observations(i) + [-obs_noise(i), obs_noise(i)], 'k', 'linewidth', 2);
    end

 


    ax = gca;
    ax.YLim = [0,50];
    ax.XLim = [1750,2050];

    ax.XLabel.String = 'time';
    ax.YLabel.String = 'GL retreat (km)';
    ax.FontSize = 12; 
    ax.XTick = 1800:100:2300;
   
end


