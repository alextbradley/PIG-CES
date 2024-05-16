% output the observations in csv format. Outputs two csv files:
% "truth.csv": observations of grounding line position
% "truth_times.csv": corresponding times at which they are measured
% "noise.csv": the noise in these observations. 

% 1930 grounding line position is that at the ridge crest, equal to the
% cold forcing steady state from PIG_306 (steady_ridge_state.mat).
%
% 2015 grounding line position is from the relaxation of the present day
% conditions from PIG_400

addpath('../functions/')
truth_times = [1930,2015]; %times at which obs are made
truth_times = [180,265]; %times at which obs are made in model time
x0 = 64; %index along which to take grounding line

data_2015 = load("2015_state.mat");
grfrac_2015 = data_2015.grfrac;
yy          = data_2015.y;
yy          = yy(1,:);

gl_pos_2015 = get_gl_pos(yy,grfrac_2015,x0);

data_1930 = load("steady_ridge_state.mat");
grfrac_1930 = data_1930.grfrac;
gl_pos_1930 = get_gl_pos(yy,grfrac_1930,x0);

data = [gl_pos_1930, gl_pos_2015];

noise = 4.5e3*ones(size(data));

writematrix(data, 'truth.csv');
writematrix(noise, 'noise.csv');
writematrix(truth_times, 'truth_times.csv');
