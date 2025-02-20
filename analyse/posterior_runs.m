% plot the evolution of the posterior mean runs

folder = "/data/hpcdata/users/aleey/projects/AttributionRealWorld/manual-EKI/ensembler/cases/manual-eki-posteriors-";

anthropogenic = 0:15;
counterfactual = 16:31; %indices of runs
counterfactual_nobump = 32:47; %indices of runs
counterfactual_nobumpnotrend = 48:63; %indices of runs

anthropogenic = [0,1,2,3,4,5,6,7,8,9,11,12,13,14,15]; %remove one??
anthropogenic = [0,1,2,5,6,7,8,9,11,12,13,14,15]; %remove one??
counterfactual = 16 + anthropogenic; %indices of runs
counterfactual_nobump = 32 + anthropogenic; %indices of runs
counterfactual_nobumpnotrend = 48 + anthropogenic; %indices of runs

gendata = 1; 
addpath('../functions')

anthro_colour = [0,0,1];
counter_colour = [0,1,0];
counter_nobump_colour = [1,0,0];



start_year = 1500; 
dx = 3e3;
dy = 3e3;


if gendata
trajectories = nan(4,length(anthropogenic),301);
gl           = nan(4,length(anthropogenic),301); 
x0 = 64; %index along which to measure gl

for i = 1:length(anthropogenic)
    fpath_anth = strcat(folder,num2str(anthropogenic(i)), "/run/outfile.nc");
    t = ncread(fpath_anth, "TIME");
    h = ncread(fpath_anth, "h");
    grfrac = ncread(fpath_anth, "grfrac");

    %loop over time points
    grv = nan(1,length(t));
    yy = ncread(fpath_anth, "y");
    gl_position = nan(1,length(t));

    for it = 1:length(t)
        grv(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
	    gl_position(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)), x0);
    end

    trajectories(1,i,:) = grv;
    gl(1,i,:) = gl_position; 


    % repeat for counter
    fpath_anth = strcat(folder,num2str(counterfactual(i)), "/run/outfile.nc");
    t = ncread(fpath_anth, "TIME");
    t = t + start_year;
    h = ncread(fpath_anth, "h");
    grfrac = ncread(fpath_anth, "grfrac");

    %loop over time points
    grv = nan(1,length(t));
    yy = ncread(fpath_anth, "y");
    gl_position = nan(1,length(t));


    for it = 1:length(t)
        grv(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
        gl_position(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)), x0);

    end

    trajectories(2,i,:) = grv;
    gl(2,i,:) = gl_position; 

    % repeat for no bump
    fpath_anth = strcat(folder,num2str(counterfactual_nobump(i)), "/run/outfile.nc");
    t = ncread(fpath_anth, "TIME");
    t = t + start_year;
    h = ncread(fpath_anth, "h");
    grfrac = ncread(fpath_anth, "grfrac");

    %loop over time points
    grv = nan(1,length(t));
    yy = ncread(fpath_anth, "y");
    gl_position = nan(1,length(t));


    for it = 1:length(t)
        grv(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
        gl_position(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)), x0);

    end

    trajectories(3,i,:) = grv;
    gl(3,i,:) = gl_position; 

    %repeat for no bump or trend
    fpath_anth = strcat(folder,num2str(counterfactual_nobumpnotrend(i)), "/run/outfile.nc");
    t = ncread(fpath_anth, "TIME");
    t = t + start_year;
    h = ncread(fpath_anth, "h");
    grfrac = ncread(fpath_anth, "grfrac");

    %loop over time points
    grv = nan(1,length(t));
    yy = ncread(fpath_anth, "y");
    gl_position = nan(1,length(t));


    for it = 1:length(t)
        grv(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
        gl_position(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)), x0);

    end

    trajectories(4,i,:) = grv;
    gl(4,i,:) = gl_position; 

end
end

%% make the plot
figure(1); clf; 
ax(1) = subplot(1,2,1); hold on;


% add the ensemble mean and std
anthro_mean = squeeze(mean(trajectories(1,:,:),2));
anthro_std  = squeeze(std(trajectories(1,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [anthro_mean - anthro_std; flip(anthro_mean + anthro_std)];
%fill(ax(1), xf, yf, 'b', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(1), t+250,anthro_mean, 'b', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(1), t + 250, squeeze(trajectories(1,i,:)), 'b',"HandleVisibility", "off");
    p.Color(4) = 0.2;
end

counter_mean = squeeze(mean(trajectories(2,:,:),2));
counter_std  = squeeze(std(trajectories(2,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_mean - counter_std; flip(counter_mean + counter_std)];
%fill(ax(1),xf, yf, 'g', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(1),t+250,counter_mean, 'g', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(1), t + 250, squeeze(trajectories(2,i,:)), 'g',"HandleVisibility", "off");
    p.Color(4) = 0.2;
end

counter_nobump_mean = squeeze(mean(trajectories(3,:,:),2));
counter_nobump_std  = squeeze(std(trajectories(3,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_nobump_mean - counter_nobump_std; flip(counter_nobump_mean + counter_nobump_std)];
%fill(ax(1),xf, yf, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(1),t+250,counter_nobump_mean, 'r', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(1), t + 250, squeeze(trajectories(3,i,:)), 'r', "HandleVisibility", "off");
    p.Color(4) = 0.2;
end

counter_nobumpnotrend_mean = squeeze(mean(trajectories(4,:,:),2));
counter_nobumpnotrend_std  = squeeze(std(trajectories(4,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_nobumpnotrend_mean - counter_nobumpnotrend_std; flip(counter_nobumpnotrend_mean + counter_nobumpnotrend_std)];
%fill(ax(1),xf, yf, 'c', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(1),t+250,counter_nobumpnotrend_mean, 'c', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(1), t + 250, squeeze(trajectories(4,i,:)), 'c', "HandleVisibility", "off");
    p.Color(4) = 0.2;
end

legend({"Anthropogenic", "Counterfactual no trend", "Counterfactual no bump", "Counterfactual no bump or trend"}, 'location', 'southwest');
%legend({"Anthropogenic", "Counterfactual no trend", "Counterfactual no bump"}, 'location', 'southwest');

% add the observations
obs = readmatrix("../observations/truth_actual_cts.csv");
obs_times = readmatrix("../observations/truth_times.csv") + 1750;
plot(ax(1),obs_times(3), obs(3), 'ko', 'markersize', 10, 'markerfacecolor', 'k',"HandleVisibility", "off");
plot(ax(1),obs_times(3)*[1,1], obs(3)*[0.99,1.01], 'color', 'k', 'linewidth', 1.5, "HandleVisibility", "off"); %fix this to be correct percentages
box on;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% repeat for gl position
ax(2) = subplot(1,2,2); hold on;


% add the ensemble mean and std
anthro_mean = squeeze(mean(gl(1,:,:),2));
anthro_std  = squeeze(std(gl(1,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [anthro_mean - anthro_std; flip(anthro_mean + anthro_std)];
%fill(ax(2),xf, yf, 'b', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(2),t+250,anthro_mean, 'b', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(2), t + 250, squeeze(gl(1,i,:)), 'b', "HandleVisibility", "off");
    p.Color(4) = 0.2;    
end

counter_mean = squeeze(mean(gl(2,:,:),2));
counter_std  = squeeze(std(gl(2,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_mean - counter_std; flip(counter_mean + counter_std)];
%fill(ax(2),xf, yf, 'g', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(2),t+250,counter_mean, 'g', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(2), t + 250, squeeze(gl(2,i,:)), 'r', "HandleVisibility", "off");
    p.Color(4) = 0.2;
end

counter_nobump_mean = squeeze(mean(gl(3,:,:),2));
counter_nobump_std  = squeeze(std(gl(3,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_nobump_mean - counter_nobump_std; flip(counter_nobump_mean + counter_nobump_std)];
%fill(ax(2),xf, yf, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(2),t+250,counter_nobump_mean, 'r', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(2), t + 250, squeeze(gl(3,i,:)), 'r', "HandleVisibility", "off");
    p.Color(4) = 0.2;
end

counter_nobumpnotrend_mean = squeeze(mean(gl(4,:,:),2));
counter_nobumpnotrend_std  = squeeze(std(gl(4,:,:),0,2));
xf = [t;flip(t)]+250; 
yf = [counter_nobumpnotrend_mean - counter_nobumpnotrend_std; flip(counter_nobumpnotrend_mean + counter_nobumpnotrend_std)];
%fill(ax(2),xf, yf, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2, "HandleVisibility", "off");
plot(ax(2),t+250,counter_nobumpnotrend_mean, 'c', 'linewidth', 1.5)
for i = 1:length(anthropogenic)
    p = plot(ax(2), t + 250, squeeze(gl(3,i,:)), 'c', "HandleVisibility", "off");
    p.Color(4) = 0.2;
end

legend({"Anthropogenic", "Counterfactual no trend", "Counterfactual no bump", "Counterfactual no bump no trend"}, 'location', 'northwest');

% add the observations
obs = readmatrix("../observations/truth_actual.csv");
%obs = readmatrix("../observations/truth_actual_cts.csv");
obs_times = readmatrix("../observations/truth_times.csv") + 1750;

plot(ax(2),obs_times(1), obs(1), 'ko', 'markersize', 10, 'markerfacecolor', 'k',"HandleVisibility", "off");
plot(ax(2),obs_times(2), obs(2), 'ko', 'markersize', 10, 'markerfacecolor', 'k',"HandleVisibility", "off");
obs_noise = readmatrix("../observations/noise_actual.csv");
plot(ax(2),obs_times(1)*[1,1],[ obs(1)- obs_noise(1),  obs(1) + obs_noise(1)], 'color', 'k', 'linewidth', 1.5, "HandleVisibility", "off"); 
plot(ax(2),obs_times(2)*[1,1],[ obs(2)- obs_noise(2),  obs(2) + obs_noise(2)], 'color', 'k', 'linewidth', 1.5, "HandleVisibility", "off"); 
box on;
grid on;
