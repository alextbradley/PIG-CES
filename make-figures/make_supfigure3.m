% Make a plot showing the correlation between sliding pre-factor and melt
% rate prefactor.
clear
addpath("../functions");
realizations = 1:14;
iterations  = 1:5;
members     = 1:20;
start_year = 1750;

% Prep stuff
n_output = 3; %3 output variables
n_input  = 7; %7 input variables
n_meta   = 2; %number of meta data (iteration, member)

n_realizations = length(realizations);
n_iterations = length(iterations);
n_members    = length(members);
model_output = struct();
model_input  = struct();
meta_data    = struct(); %store iteration and member numbers
count = 1;
input_headers  = ["weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",...
    "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration"]; %note that ungrounded weertman c is redundant
output_headers =  ["dimensionless_gl_error_1930", "dimensionless_gl_error_2015", "dimensionless_grv_error_2015"];

for ir = 1:n_realizations
    realization = realizations(ir);
    for ii = 1:n_iterations
        iter = iterations(ii);
        padded_realization = sprintf("%03d", realization);
        padded_iteration  = sprintf("%03d", iter);

        % get model parameters
        params_path = strcat("../model-inputs-and-outputs/realization", padded_realization, "/iteration", padded_iteration, "/params.csv");
        params = readmatrix(params_path);

        for im = 1:n_members
            mem = members(im);
            padded_member  = sprintf("%03d", mem);
            output_path = strcat("../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/outputs.csv");
            trajectory_output_path = strcat("../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/output_trajectory.mat");
            outputs = readmatrix(output_path);
            trajectories = load(trajectory_output_path);

            %calibration quantities
            model_output(ir,ii, im).dimensionless_gl_1930_error  = outputs(1);
            model_output(ir,ii, im).dimensionless_gl_2015_error  = outputs(2);
            model_output(ir,ii, im).dimensionless_grv_2015_error = outputs(3);

            %time series bits
            if isfield(trajectories, 'time')
                model_output(ir,ii,im).time = trajectories.time + start_year;
            else
                model_output(ir,ii,im).time = trajectories.t + start_year;
            end

            model_output(ir,ii,im).grv = trajectories.grv;
            grv = trajectories.grv;
            idx = 266; %index of 2015
            model_output(ir,ii,im).grv2015 = grv(idx);
            model_output(ir,ii,im).gl_pos_discrete = trajectories.gl_pos_discrete;
            model_output(ir,ii,im).gl_pos_cts = trajectories.gl_pos_cts;
            gl_pos_cts =  trajectories.gl_pos_cts;

            model_output(ir,ii,im).gl_retreat = gl_pos_cts - gl_pos_cts(1);

            model_input(ir,ii,im).weertman_c_prefactor = params(im,1);
            model_input(ir,ii,im).ungrounded_weertmanC_prefactor = params(im,2);
            model_input(ir,ii,im).glen_a_ref_prefactor = params(im,3);
            model_input(ir,ii,im).melt_rate_prefactor_exponent = params(im,4);
            model_input(ir,ii,im).per_century_trend = params(im,5);
            model_input(ir,ii,im).bump_amplitude = params(im,6);
            model_input(ir,ii,im).bump_duration = params(im,7);

            meta_data(ir,ii,im).iteration = iter;
            meta_data(ir,ii,im).member = mem;
            meta_data(ir,ii,im).realization = realization;
        end %end loop over members
    end %end loop over iterations
end %end loop over realizations

%%
wc_all = [model_input(:).weertman_c_prefactor];
m_all = [model_input(:).melt_rate_prefactor_exponent];
grv2015_all = [model_output(:).grv2015];
obs       = readmatrix("../observations/truth_actual_cts.csv");
obs_grv2015 = obs(3);
relerr = (grv2015_all - obs_grv2015)/1e13;
relerr(relerr < -5) = -5;
relerr(relerr > 5) = 5; %normalise



figure(1); clf; scatter(m_all, wc_all, 30,relerr, 'filled', 'MarkerEdgeColor','k');
c = colorbar;
c.Label.String = '2015 grounded volume error (10^{12} m^3)';
c.FontSize = 14;
colormap(redblue(101));
box on;
ax = gca;
ax.FontSize = 12;
ax.XLabel.String = 'melt prefactor';
ax.YLabel.String = 'sliding prefactor';
clim([-4,4])
shg

% repeat with removed large errors and linear fit
m_all_red = m_all(abs(relerr)<0.1);
wc_all_red = wc_all(abs(relerr)<0.1);
relerr_red = relerr(abs(relerr)<0.1);

% do a linear fit 
figure(2); clf; scatter(m_all_red, wc_all_red, 30,relerr_red, 'filled', 'MarkerEdgeColor','k');
%c = colorbar;
%c.Label.String = '2015 grounded volume error (10^{12} m^3)';
%c.FontSize = 14;
colormap(redblue(101));
box on;
ax2 = gca;
ax2.FontSize = 12;
ax2.XLabel.String = 'melt prefactor';
ax2.YLabel.String = 'sliding prefactor';
clim([-4,4])
shg

ax2.XLim = ax.XLim;
ax2.YLim = ax.YLim;
