% Make figure 8 of the manuscript: (a) reconstructed grounding line
% positions for the four cases: default, no trend, no 1940s ENSO, and no
% trend nor 1940s ENSO. Timeframe: 1750-2200. (b) 

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.

%% Preliminaries
addpath("../functions")

realizations = 1:14;
samples      = 1:10;
outpath   = "../model-inputs-and-outputs/forward-runs/";
scenarios = ["default", "nobump", "notrend", "notrend_nobump"];

colmap  = linspecer(length(scenarios));

init_time = 1750;
fig =figure(1); clf; fig.Position(3:4) =[1150, 460];
fs = 16;
for i = 1:2
    ax(i) = subplot(1,2,i);
    hold(ax(i), 'on')
    box(ax(i), 'on')
    ax(i).FontSize = fs;
end

%% Pull data
ss_gl  = zeros(length(scenarios), length(realizations), length(samples), 276);
for iSc = 1:length(scenarios)
    for ir = 1:length(realizations)
        for is = 1:length(samples)
            fpath = strcat(outpath,scenarios(iSc), "/realization_", sprintf('%03d', realizations(ir)),"/sample_", sprintf('%03d',samples(is)), "/output_trajectory.mat");
            data = load(fpath);
            gl_retreat = data.gl_pos_cts;
            gl_retreat = (gl_retreat - gl_retreat(1))/1e3;% turn into retreat in km
            ss_gl(iSc, ir,is,:) = gl_retreat;       

        end
    end
end

%average over realizations (3rd entry) and members (4th entry) 
ens_mean_gl  = squeeze(mean(mean(ss_gl,3),2));  
ens_std_gl   = squeeze(std(std(ss_gl, 0, 2), 0, 3)); 

time = data.t;
time = time + init_time;

%% get the obs
observations = readmatrix("../observations/truth_actual.csv");
observations_gl = observations(1:2); %remove the gl constraint
observations_gl = (observations_gl - observations_gl(1))/1e3; % turn into retreat in km

obs_noise     = readmatrix("../observations/noise_actual.csv");
obs_noise_gl  = obs_noise(1:2)/1e3;

obs_time     = readmatrix("../observations/truth_times.csv");
obs_time     = obs_time + init_time;


%% Make gl plot

for iSc = 1:length(scenarios)
    xf = [time; flip(time)];
    yf = [(ens_mean_gl(iSc, :) - ens_std_gl(iSc, :)),flip(ens_mean_gl(iSc, :) + ens_std_gl(iSc, :)) ]';
    fill(ax(1), xf, yf, colmap(iSc,:), 'FaceAlpha',0.25, 'LineStyle','none', 'HandleVisibility','off');
    plot(ax(1), time, ens_mean_gl(iSc, :), 'Color',colmap(iSc, :), "LineWidth",2);

end

ax(1).YLim = [-10,200];
ax(1).XLim = [1750,2200];

ax(1).XLabel.String = 'time';
ax(1).YLabel.String = 'GL retreat (km)';
ax(1).XTick = 1800:100:2200;

% add obs
for i = 1:2
    if i == 1
        plot(ax(1), obs_time(i), observations_gl(i), 'ro', 'markersize', 10, 'markerfacecolor', 'r');
    else
        plot(ax(1), obs_time(i), observations_gl(i), 'ro', 'markersize', 10, 'markerfacecolor', 'r', 'HandleVisibility','off');

    end

    plot(ax(1), obs_time(i)*[1,1], observations_gl(i) + [-obs_noise_gl(i), obs_noise_gl(i)], 'r', 'linewidth', 2, 'HandleVisibility','off');
end

% add legend
l = legend(ax(1), ["All forcings", "No 1940s event", "No trend", "No trend nor 1940s event", "observations"], 'Location', 'northwest', "FontSize",fs, 'LineWidth',1);

% add box for figure 7
YLf7 = [-5,40];
XLf7 = [1750,2050];


plot(ax(1), XLf7, [YLf7(1), YLf7(1)],'k--', 'LineWidth',1, 'HandleVisibility','off');
plot(ax(1), XLf7, [YLf7(2), YLf7(2)],'k--', 'LineWidth',1, 'HandleVisibility','off');
plot(ax(1), [XLf7(1), XLf7(1)], YLf7,'k--', 'LineWidth',1, 'HandleVisibility','off');
plot(ax(1), [XLf7(2), XLf7(2)], YLf7,'k--', 'LineWidth',1, 'HandleVisibility','off');
%% Compute the attributable percentage stats
[~, idx1] = min(abs(time - 1930)); % fixed

np = length(time) - idx1; % number of time points


time_after_1930    = time(idx1:end);
total_retreat      = nan(1,length(time_after_1930));
total_retreat_ub   = nan(1,length(time_after_1930));
total_retreat_lb   = nan(1,length(time_after_1930));
noanth_retreat     = nan(1,length(time_after_1930));
noanth_retreat_ub  = nan(1,length(time_after_1930));
noanth_retreat_lb  = nan(1,length(time_after_1930));

for ip = 0:np
    idx2 = idx1 + ip; %index of the point in question
    total_retreat(ip+1)    =  ens_mean_gl(1,idx1) - ens_mean_gl(1,idx2); %total retreat in the all forcings
    total_retreat_ub(ip+1) =  ens_mean_gl(1,idx1) - (ens_mean_gl(1,idx2) + ens_std_gl(1,idx2));
    total_retreat_lb(ip+1) =  ens_mean_gl(1,idx1)  - (ens_mean_gl(1,idx2) - ens_std_gl(1,idx2));



    noanth_retreat(ip+1)    =  ens_mean_gl(3,idx1) - ens_mean_gl(3,idx2);
    noanth_retreat_lb(ip+1) =  ens_mean_gl(3,idx1) - (ens_mean_gl(3,idx2) - ens_std_gl(3,idx2));
    noanth_retreat_ub(ip+1) =  ens_mean_gl(3,idx1) - (ens_mean_gl(3,idx2) + ens_std_gl(3,idx2));



end
    

anthro_retreat_far = (total_retreat - noanth_retreat)./total_retreat;
anthro_retreat_far_ub = (total_retreat_lb - noanth_retreat_ub)./total_retreat_lb ;
anthro_retreat_far_lb = (total_retreat_ub - noanth_retreat_lb)./total_retreat_ub;

%% make second panel
xf = [time_after_1930; flip(time_after_1930)];
yf = [anthro_retreat_far_lb'; flip(anthro_retreat_far_ub')];

plot(ax(2), [1930, 2300], [0,0], 'k--', 'LineWidth',1.5)
fill(ax(2), xf, yf,  [0, 47, 167]/255, 'FaceAlpha', 0.2, "LineStyle",'none');
plot(ax(2), time_after_1930, anthro_retreat_far, 'color', [0, 47, 167]/255, 'LineWidth',2)
shg

ax(2).XLim = [1960, 2200];
ax(2).YLim = 0.01*[-50,100];
ax(2).XLabel.String = 'time';
ax(2).YLabel.String = 'fraction of attributable retreat';