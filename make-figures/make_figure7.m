% Make figure 7 of the manuscript: reconstructed grounding line positions
% and ice volume for the four cases: default, no trend, no 1940s ENSO, and
% no trend nor 1940s ENSO.

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.

%% Preliminaries
addpath("../functions")

realizations = 1:14;
samples      = 1:10;
outpath   = "../model-inputs-and-outputs/forward-runs/";
scenarios = ["default", "nobump", "notrend", "notrend_nobump"];

colmap  = linspecer(length(scenarios));

init_time = 1750;
fig =figure(1); clf; fig.Position(3:4) = [1150, 850];

for i = 1:4
    ax(i) = subplot(2,2,i);
    hold(ax(i), 'on')
    box(ax(i), 'on')
    grid on
    ax(i).FontSize = 14;
end

%% Pull data
ss_gl  = zeros(length(scenarios), length(realizations), length(samples), 276);
ss_grv = zeros(length(scenarios), length(realizations), length(samples), 276);
for iSc = 1:length(scenarios)
    for ir = 1:length(realizations)
        for is = 1:length(samples)
            fpath = strcat(outpath,scenarios(iSc), "/realization_", sprintf('%03d', realizations(ir)),"/sample_", sprintf('%03d',samples(is)), "/output_trajectory.mat");
            data = load(fpath);
            gl_retreat = data.gl_pos_cts;
            gl_retreat = (gl_retreat - gl_retreat(1))/1e3;% turn into retreat in km
            ss_gl(iSc, ir,is,:) = gl_retreat;
            ss_grv(iSc, ir, is, :) = data.grv;
       

        end
    end
end

%average over realizations (3rd entry) and members (4th entry) 
ens_mean_grv = squeeze(mean(mean(ss_grv,3),2));
ens_mean_gl  = squeeze(mean(mean(ss_gl,3),2));  

ens_std_gl   = squeeze(std(std(ss_gl, 0, 2), 0, 3)); 
ens_std_grv  = squeeze(std(std(ss_grv, 0, 2), 0, 3)); 

time = data.t;
time = time + init_time;

%% get the obs
observations = readmatrix("../observations/truth_actual.csv");
observations_gl = observations(1:2); %remove the gl constraint
observations_gl = (observations_gl - observations_gl(1))/1e3; % turn into retreat in km
observations_grv = observations(3); 

obs_noise     = readmatrix("../observations/noise_actual.csv");
obs_noise_gl  = obs_noise(1:2)/1e3;
obs_noise_grv = obs_noise(3);

obs_time     = readmatrix("../observations/truth_times.csv");
obs_time     = obs_time + init_time;


%% Make gl plot

for iSc = 1:length(scenarios)
    xf = [time; flip(time)];
    yf = [(ens_mean_gl(iSc, :) - ens_std_gl(iSc, :)),flip(ens_mean_gl(iSc, :) + ens_std_gl(iSc, :)) ]';
    fill(ax(1), xf, yf, colmap(iSc,:), 'FaceAlpha',0.25, 'LineStyle','none', 'HandleVisibility','off');
    plot(ax(1), time, ens_mean_gl(iSc, :), 'Color',colmap(iSc, :), "LineWidth",2);

end

ax(1).YLim = [-5,40];
ax(1).XLim = [1750,2050];

ax(1).XLabel.String = 'time';
ax(1).YLabel.String = 'GL retreat (km)';
ax(1).XTick = 1800:50:2000;

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
l = legend(ax(1), ["All forcings", "No 1940s event", "No anthro trend", "No anthro trend nor 1940s event", "observations"], 'Location', 'northwest', "FontSize",14, 'LineWidth',1);


% print out the attribution values for gl 
[~, idx1] = min(abs(time - 1930));
[~, idx2] = min(abs(time - 2015));

total_retreat  =  ens_mean_gl(1,idx1) - ens_mean_gl(1,idx2); %total retreat in the all forcings
noanth_retreat =  ens_mean_gl(3,idx1) - ens_mean_gl(3,idx2);
nobump_retreat =  ens_mean_gl(2,idx1) - ens_mean_gl(2,idx2);

total_retreat_ub =  ens_mean_gl(1,idx1) - (ens_mean_gl(1,idx2) + ens_std_gl(1,idx2));%maximum possible retreat (don't put sd in the 1930s position bc same for each)
noanth_retreat_ub =  ens_mean_gl(3,idx1) - (ens_mean_gl(3,idx2) + ens_std_gl(3,idx2));
nobump_retreat_ub =  ens_mean_gl(2,idx1) - (ens_mean_gl(2,idx2) + ens_std_gl(2,idx2));
total_retreat_lb =  ens_mean_gl(1,idx1)  - (ens_mean_gl(1,idx2) - ens_std_gl(1,idx2));
noanth_retreat_lb =  ens_mean_gl(3,idx1) - (ens_mean_gl(3,idx2) - ens_std_gl(3,idx2));
nobump_retreat_lb =  ens_mean_gl(2,idx1) - (ens_mean_gl(2,idx2) - ens_std_gl(2,idx2));

anthro_perc_retreat = (total_retreat - noanth_retreat)/total_retreat * 100;
bump_perc_retreat   = (total_retreat - nobump_retreat)/total_retreat * 100;

fprintf('anthropogenic forcing responsible for %.1f percent of gl retreat from 1930 to 2015 \n',anthro_perc_retreat )
fprintf('min anthropogenic forcing responsible for %.1f percent of gl retreat from 1930 to 2015 \n', (total_retreat_lb - noanth_retreat_ub)/total_retreat_lb * 100 )
fprintf('max anthropogenic forcing responsible for %.1f percent of gl retreat from 1930 to 2015 \n', (total_retreat_ub - noanth_retreat_lb)/total_retreat_ub * 100 )
fprintf('excess anthropogenic gl retreat is %.1f km \n', -(total_retreat - noanth_retreat))
fprintf('max anthropogenic gl retreat is %.1f km \n', -(total_retreat_lb - noanth_retreat_ub))
fprintf('min anthropogenic gl retreat is %.1f km \n', -(total_retreat_ub - noanth_retreat_lb))
fprintf('\n')

fprintf('1940s event responsible for %.1f percent of gl retreat from 1930 to 2015 \n',bump_perc_retreat )
fprintf('min 1940s event forcing responsible for %.1f percent of gl retreat from 1930 to 2015 \n', (total_retreat_lb - nobump_retreat_ub)/total_retreat_lb * 100 )
fprintf('max 1940s event forcing responsible for %.1f percent of gl retreat from 1930 to 2015 \n', (total_retreat_ub - nobump_retreat_lb)/total_retreat_ub * 100 )
fprintf('excess 1940s events gl retreat is %.1f km \n', -(total_retreat - nobump_retreat))
fprintf('max anthropogenic gl retreat is %.1f km \n', -(total_retreat_lb - nobump_retreat_ub))
fprintf('min anthropogenic gl retreat is %.1f km \n', -(total_retreat_ub - nobump_retreat_lb))
%% Make grv plot

for iSc = 1:length(scenarios)
    xf = [time; flip(time)];
    yf = [(ens_mean_grv(iSc, :) - ens_std_grv(iSc, :)),flip(ens_mean_grv(iSc, :) + ens_std_grv(iSc, :)) ]/1e14';
    fill(ax(2), xf, yf, colmap(iSc,:), 'FaceAlpha',0.25, 'LineStyle','none', 'HandleVisibility','off');
    plot(ax(2), time, ens_mean_grv(iSc, :)/1e14, 'Color',colmap(iSc, :), "LineWidth",2);

end

% add obs 
plot(ax(2), obs_time(3), observations_grv/1e14, 'ro', 'markersize', 10, 'markerfacecolor', 'r');
plot(ax(2), obs_time(3)*[1,1], (observations_grv + [-obs_noise_grv, obs_noise_grv])/1e14, 'r', 'linewidth', 2, 'HandleVisibility','off');


ax(2).YLim = [4.5,4.9];
ax(2).XLim = [1750,2050];

ax(2).XLabel.String = 'time';
ax(2).YLabel.String = 'Grounded volume (10^{14} m^3)';
ax(2).XTick = 1800:50:2000;
ax(2).YTick = 4.5:.1:5;

% attribution stuff
total_retreat_grv  =  ens_mean_grv(1,idx1) - ens_mean_grv(1,idx2); %total retreat in the all forcings
noanth_retreat_grv =  ens_mean_grv(3,idx1) - ens_mean_grv(3,idx2);
nobump_retreat_grv =  ens_mean_grv(2,idx1) - ens_mean_grv(2,idx2);

fprintf('excess anthropogenic grv loss is %.2f  x 10^12 m^3 \n', (total_retreat_grv - noanth_retreat_grv)/1e12)
fprintf('excess 1940s ENSO grv loss is %.2f x 10^12 m^3 \n', (total_retreat_grv - nobump_retreat_grv)/1e12)

%% Add distributions

[~,idx1] = min(abs(time - obs_time(1)));
[~,idx2] = min(abs(time - obs_time(2))); %get indices for 1930 and 2015

ss_gl_1930 = ss_gl(:,:,:, idx1);
ss_gl_2015 = ss_gl(:,:,:, idx2);
ss_gl_diff = ss_gl_2015 - ss_gl_1930;

sz = size(ss_gl_diff);

ss_gl_diff = reshape(ss_gl_diff, [sz(1), sz(2)*sz(3)]); %remove realization dependence

% for jj = 1:4
%     [fp,xfp]=kde(ss_gl_diff(jj,:)); %do KDE on the gl differences
%     plot(ax(3), xfp, fp, 'color', colmap(jj,:),"LineWidth",2);
% end


ss_grv_2015 = ss_grv(:,:,:, idx2);
ss_grv_2015 = reshape(ss_grv_2015, [sz(1), sz(2)*sz(3)]); %remove realization dependence

% for jj = 1:4
%     [fp,xfp]=kde(ss_grv_2015(jj,:)); %do KDE on the gl differences
%     plot(ax(4), xfp, fp, 'color', colmap(jj,:),"LineWidth",2);
% end


%add shaded band for observations
xf = [30,42,42,30]; %hard code lol
yf = [0,0,5,5];
fill(ax(3), xf, yf,'r', 'FaceAlpha', 0.15, 'LineStyle', 'none');
plot(ax(3),(observations(2) - observations(1))/1e3*[1,1], [0,5], 'r--', 'LineWidth',2);


for jj = 1:4
    Qs = std(ss_gl_diff(jj,:));
    Qc = mean(ss_gl_diff(jj,:));

    %plot the 5-95
    plot(ax(3), [Qc-Qs,Qc + Qs], [jj,jj], 'LineWidth',10, 'Color',colmap(jj,:));
    plot(ax(3), Qc,jj, 'ko', 'MarkerEdgeColor', 'k','MarkerSize', 8, "LineWidth",4)

end


ax(3).Position(4) = 0.08;
ax(3).XLabel.String = "GL change (km)";
grid(ax(3), 'off');
ax(3).YLim = [0,5];
ax(3).YTick = [];
ax(3).Position(2) = 0.43;



%% Repeat for the grounded volume

%add background obs
xf = [-obs_noise(3), obs_noise(3), obs_noise(3),- obs_noise(3)]/1e12; %hard code lol
yf = [0,0,5,5];
fill(ax(4), xf, yf,'r', 'FaceAlpha', 0.15, 'LineStyle', 'none');
plot(ax(4),[0,0], [0,5], 'r--', 'LineWidth',2);

ss_grv_2015 = observations(3) - ss_grv_2015;
for jj = 1:4
    Qs = std(ss_grv_2015(jj,:))/1e12;
    Qc = mean(ss_grv_2015(jj,:))/1e12;

    %plot the 5-95
    plot(ax(4), [Qc-Qs,Qc + Qs], [jj,jj], 'LineWidth',10, 'Color',colmap(jj,:));
    plot(ax(4), Qc,jj, 'ko', 'MarkerEdgeColor', 'k','MarkerSize', 8, "LineWidth",4)

end

ax(4).Position(4) = 0.08;
ax(4).XLabel.String = "2015 grounded volume relative to observed (10^{12} m^3)";
grid(ax(4), 'off');
ax(4).YLim = [0,5];
ax(4).YTick = [];
ax(4).Position(2) = 0.43;


