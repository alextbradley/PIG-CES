% output the observations in csv format. Outputs two csv files:
% "truth_volume.csv": observations of total ice volume
% "truth_vaf.csv": observations of ice vaf
% "truth_times.csv": corresponding times at which they are measured
% "noise_volume.csv": the noise in the total volume observations. 
% "noise_vaf.csv": noise in the ice vaf observation

% 1930 grounding line position is that at the ridge crest, equal to the
% cold forcing steady state from PIG_306 (steady_ridge_state.mat).
%
% 2015 grounding line position is from the relaxation of the present day
% conditions from PIG_400

addpath('../functions/')
dx = 3e3; dy = 3e3;
truth_times = [1930,2015]; %times at which obs are made
truth_times = [180,265]; %times at which obs are made in model time

data_2015 = load("2015_state.mat");
volume_2015 = sum(sum(data_2015.h))*dx*dy;

data_1930 = load("steady_ridge_state.mat");
volume_1930 = sum(sum(data_1930.h))*dx*dy;

data = [volume_1930, volume_2015]/1e14;

noise = [0.1, 0.01];

writematrix(data, 'truth_volume.csv');
writematrix(noise, 'noise_volume.csv');
writematrix(truth_times, 'truth_times.csv');

%repeat for grounded volume
data_2015 = load("2015_state.mat");
volume_2015 = sum(sum(data_2015.h .* data_2015.grfrac))*dx*dy;

data_1930 = load("steady_ridge_state.mat");
volume_1930 = sum(sum(data_1930.h .* data_1930.grfrac))*dx*dy;

data = [volume_1930, volume_2015]/1e14;

noise = [0.1, 0.01];

writematrix(data, 'truth_grounded_volume.csv');
writematrix(noise, 'noise_grounded_volume.csv');

% repeat for vaf

haf_2015 = data_2015.h + (1028.0/918.0)*(data_2015.b); %height above floatation
vaf_2015 = sum(sum(haf_2015(haf_2015 > 0)))*dx*dy / 1e14; %vaf

haf_1930 = data_1930.h + (1028.0/918.0)*(data_1930.b); %height above floatation
vaf_1930 = sum(sum(haf_1930(haf_1930 > 0)))*dx*dy / 1e14; %vaf

data = [vaf_1930, vaf_2015];

noise = [0.1, 0.01];

writematrix(data, 'truth_vaf.csv');
writematrix(noise, 'noise_vaf.csv');

