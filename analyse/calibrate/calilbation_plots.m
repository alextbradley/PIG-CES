% Master script for plots associated with calibration
%
% Plots:
% scatter_error_and_parameter_colour_by_iteration
%       make a scatter plot of the parameter values and normalized error
%
% plot_grv_and_ctsglpos_trajectories
%       plot time series of grounded volume nad grounding line position,
%       and color by iteration, using cts grounding line measurement
%
% plot_error_correlations
%       produce scatter plots of error in grounded volume against error in
%       both of the grounding line measurements
%

%specify plots
scatter_error_and_parameter_colour_by_iteration = 0;
plot_grv_and_ctsglpos_trajectories              = 0;
plot_error_correlations                         = 1; 


gendata = 0
%% Preliminaries
addpath('../../functions/');

%specify ensemble
realization = 30;
iterations  = 1:5;
members     = 1:20;

start_year = 1750;

numplot = 1;
fontsize = 14;
%% Get the data
if gendata

% Prep stuff
n_output = 3; %3 output variables
n_input  = 7; %7 input variables
n_meta   = 2; %number of meta data (iteration, member)

n_iterations = length(iterations);
n_members    = length(members);
model_output = struct();
model_input  = struct();
meta_data    = struct(); %store iteration and member numbers
count = 1;
input_headers  = ["weertman_c_prefactor", "ungrounded_weertmanC_prefactor" ,"glen_a_ref_prefactor",...
    "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration"];
output_headers =  ["dimensionless_gl_error_1930", "dimensionless_gl_error_2015", "dimensionless_grv_error_2015"];


for ii = 1:n_iterations
    iter = iterations(ii);
    padded_realization = sprintf("%03d", realization);
    padded_iteration  = sprintf("%03d", iter);

    % get model parameters
    params_path = strcat("../../model-inputs-and-outputs/realization", padded_realization, "/iteration", padded_iteration, "/params.csv");
    params = readmatrix(params_path);

    for im = 1:n_members
        mem = members(im);
        padded_member  = sprintf("%03d", mem);
        output_path = strcat("../../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/outputs.csv");
        trajectory_output_path = strcat("../../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/output_trajectory.mat");
        outputs = readmatrix(output_path);
        trajectories = load(trajectory_output_path);

        %calibration quantities
        model_output(ii, im).dimensionless_gl_1930_error  = outputs(1);
        model_output(ii, im).dimensionless_gl_2015_error  = outputs(2);
        model_output(ii, im).dimensionless_grv_2015_error = outputs(3);

        %time series bits
        model_output(ii,im).time = trajectories.time + start_year;
        model_output(ii,im).grv = trajectories.grv; %this will be changed to grv
        model_output(ii,im).gl_pos_discrete = trajectories.gl_pos_discrete;
        model_output(ii,im).gl_pos_cts = trajectories.gl_pos_cts;
        
        model_input(ii,im).weertman_c_prefactor = params(im,1);
        model_input(ii,im).ungrounded_weertmanC_prefactor = params(im,1);
        model_input(ii,im).glen_a_ref_prefactor = params(im,1);
        model_input(ii,im).melt_rate_prefactor_exponent = params(im,1);
        model_input(ii,im).per_century_trend = params(im,1);
        model_input(ii,im).bump_amplitude = params(im,1);
        model_input(ii,im).bump_duration = params(im,1);

        meta_data(ii,im).iteration = iter;
        meta_data(ii,im).member = mem;       
    end %end loop over members
end %end loop over iterations
end %end gendata flag

%% scatter_error_and_parameter_colour_by_iteration
if scatter_error_and_parameter_colour_by_iteration



end


%% plot_grv_and_glpos_trajectories
if plot_grv_and_ctsglpos_trajectories
   
    colmap = cmocean('ice',length(iterations)+2 );
    obscolor = [1,0.,0.];
    figure(numplot); clf;
    for i = 1:2
        ax(i) = subplot(1,2,i); 
        hold(ax(i), "on");
        box(ax(i), "on");
        ax(i).FontSize = fontsize;
    end

    for ii = 1:length(iterations)
        for im = 1:length(members)
            plot(ax(1), model_output(ii,im).time, smooth(model_output(ii,im).gl_pos_cts,1),'linewidth', 1.5, 'Color',colmap(ii,:) );

            subplot(1,2,2); hold on;
            plot(ax(2), model_output(ii,im).time, model_output(ii,im).grv,'linewidth', 1.5, 'Color',colmap(ii,:) );

        end
    end %end loop over iterations

    % add observations
    obs       = readmatrix("../../observations/truth_actual_cts.csv");
    obs_times = readmatrix("../../observations/truth_times.csv") + start_year;
    obs_err   = readmatrix('../../observations/noise_actual.csv');
   
   
    plot(ax(1),obs_times(1)*[1,1], obs(1)+obs_err(1)*[-1,1],'color', obscolor, 'linewidth', 1.5);
    plot(ax(1),obs_times(2)*[1,1], obs(2)+obs_err(2)*[-1,1], 'color', obscolor, 'linewidth', 1.5);
    plot(ax(1),obs_times(1:2), obs(1:2), 'ro', 'markersize', 10, 'markerfacecolor', obscolor);

    plot(ax(2),obs_times(3)*[1,1], obs(3)+obs_err(3)*[-1,1], 'color', obscolor, 'linewidth', 1.5);
    plot(ax(2),obs_times(3), obs(3), 'ro', 'markersize', 10, 'markerfacecolor', obscolor);



    %tidy stuff
    ax(1).XLabel.String = "time";
    ax(2).XLabel.String = "time";

    ax(1).YLabel.String = "grounding line position";
    ax(2).YLabel.String = "grounded volume";
    
    ax(1).YLim = [-3.1,-2.55]*1e5;
    ax(2).YLim = [4.3,5.1]*1e14;

    % add a colorbar
    c = colorbar(ax(2));
    c.Position(1) = 0.915;
    c.Position(3) = 0.01;
    set(c, 'Colormap',colmap(1:length(iterations), :))

    dtick = 1/length(iterations);
    c.Ticks = [dtick/2:dtick:(1-dtick/2)];
    c.TickLabels = 1:length(iterations);
    c.Label.String = 'iteration';
    c.FontSize = fontsize;
    
    numplot = numplot + 1;
    
end

%% plot_error_correlations
if plot_error_correlations

    figure(numplot); clf; hold on; box on;

    colmap = cmocean('ice',length(iterations)+2 );

    errs_1930 = [model_output(:).dimensionless_gl_1930_error];
    errs_2015 = [model_output(:).dimensionless_gl_2015_error];
    errs_grv  = [model_output(:).dimensionless_grv_2015_error];
    iters     = [meta_data(:).iteration];
    
    cmap = nan(length(errs_1930),3);
    for i = 1:length(errs_1930)
        cmap(i,:) = colmap(iters(i),:);
    end

    scatter(errs_grv, errs_1930+errs_2015,60, cmap, 'filled')

    ax(1) = gca;
    ax(1).YLim = [-20, 20];
    ax(1).XLim = [-50, 50];
    ax(1).FontSize = fontsize;
    ax(1).XLabel.String = 'dimensionless error in grounded volume in 2015';
    ax(1).YLabel.String = 'grounding line error 1930 + 2015';

    % fit linear model
    idx = (errs_grv > -50) & (errs_grv < 50);  %remove outliers
    mdl = fitlm(errs_grv,errs_1930+errs_2015);
    cfs = mdl.Coefficients.Estimate;
    xx = ax(1).XLim;
    plot(xx, cfs(1) + cfs(2)*xx, 'k--', 'linewidth', 1.5)
    
end
