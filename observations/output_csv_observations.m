% output the observations in csv format. Outputs two csv files:
% "truth.csv": = [0,0,0], array of observations
%                (1) of GL position in 1930 measured relative to the GL position in 1930
%                (2) of GL position in 2015 measured relative to the GL position in 2015
%                (3) of grounded volume in 2015, releative to mass in 2015
% "sd.csv": standard deviation of these observations
% "truth_actual.csv": raw values of the grounding line position and mass
% "truth_actual_cts.csv": raw values of the grounding line position and mass using the cts grounding line position
% "truth_times.csv": corresponding times at which they are measured
% "noise.csv": = [1,1,1] the noise in these observations, scaled by the sd

% 1930 grounding line position is that at the ridge crest, equal to the
% cold forcing steady state from PIG_306 (steady_ridge_state.mat).
%
% 2015 grounding line position is from the relaxation of the present day
% conditions from PIG_400

addpath('../functions/')
truth_times = [1930,2015,2015]; %times at which obs are made
truth_times = [180,265,265]; %times at which obs are made in model time
x0 = 64; %index along which to take grounding line

data_2015 = load("2015_state.mat");
grfrac_2015 = data_2015.grfrac;
yy          = data_2015.y;
yy          = yy(1,:);

dx = 3e3;
dy = 3e3;
volume_2015 = sum(sum(data_2015.h .* data_2015.grfrac))*dx*dy;
gl_pos_2015 = get_gl_pos(yy,grfrac_2015,x0);
gl_pos_2015_cts = get_gl_pos_cts(yy,grfrac_2015,x0);

data_1930 = load("steady_ridge_state.mat");
grfrac_1930 = data_1930.grfrac;
gl_pos_1930 = get_gl_pos(yy,grfrac_1930,x0);
gl_pos_1930_cts = get_gl_pos_cts(yy,grfrac_1930,x0);

data = [gl_pos_1930, gl_pos_2015, volume_2015];
data_cts = [gl_pos_1930_cts, gl_pos_2015_cts, volume_2015];
data_normalized = [0,0,0];

noise_actual = [4500,4500,1e13]; %stds in these obs
noise = ones(size(data));

writematrix(data, 'truth_actual.csv');
writematrix(data_cts, 'truth_actual_cts.csv');
writematrix(data_normalized, 'truth.csv')
writematrix(noise, 'noise.csv');
writematrix(noise_actual, 'noise.csv');
writematrix(truth_times, 'truth_times.csv');
