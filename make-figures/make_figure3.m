% Make figure 3 of the manuscript, showing the convergence of the ensemble
% Kalman Iteration in (a) the grounding line position and (b) grounded
% volume. Also plot (g--i) the parameter values coloured by iteration and
% (j) average error as a function of iteration for each realization of
% forcing and (k) histograms of error in EKI vs LHC
%
% 20/02/2025, ATB (alex.bradley@kcl.ac.uk). MIT licence.

%% Preliminaries
gendata = 1;
addpath('../functions/');

%specify ensemble
realization = 3;
iterations  = 1:5;
members     = 1:20;

start_year = 1750;

%% Set up the figure
fig = figure(1); clf;
fig.Position(3:4) = [800, 900];

positions = [0.07, 0.06, 0.41, 0.24;
    0.56, 0.06, 0.41, 0.24;
    0.07, 0.36,  0.27, 0.13;
    0.385,  0.36,  0.27, 0.13;
    0.7,  0.36,  0.27, 0.13;
    0.07, 0.55, 0.27, 0.13;
    0.385, 0.55, 0.27, 0.13;
    0.7,  0.55, 0.27, 0.13;
    0.07, 0.74,  0.41, 0.25;
    0.56, 0.74,  0.41, 0.25];

for i = 1:10
    ax(i) = subplot('Position',positions(i,:));
    box(ax(i), 'on');
    hold(ax(i), 'on')
    ax(i).FontSize = 14;
    ax(i).FontName = 'arial';
end


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
        "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration"]; %note that ungrounded weertman c is redundant
    output_headers =  ["dimensionless_gl_error_1930", "dimensionless_gl_error_2015", "dimensionless_grv_error_2015"];


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
            model_output(ii, im).dimensionless_gl_1930_error  = outputs(1);
            model_output(ii, im).dimensionless_gl_2015_error  = outputs(2);
            model_output(ii, im).dimensionless_grv_2015_error = outputs(3);

            %time series bits
            if ((realization <= 30) && (realization >= 26))
                model_output(ii,im).time = trajectories.time + start_year;
            else
                model_output(ii,im).time = trajectories.t + start_year;
            end

            model_output(ii,im).grv = trajectories.grv;
            grv = trajectories.grv;
            [~,idx] = min(abs( trajectories.t + start_year - 2015));
            model_output(ii,im).grv2015 = grv(idx);
            model_output(ii,im).gl_pos_discrete = trajectories.gl_pos_discrete;
            model_output(ii,im).gl_pos_cts = trajectories.gl_pos_cts;
            gl_pos_cts =  trajectories.gl_pos_cts;

            model_output(ii,im).gl_retreat = gl_pos_cts - gl_pos_cts(1);

            model_input(ii,im).weertman_c_prefactor = params(im,1);
            model_input(ii,im).ungrounded_weertmanC_prefactor = params(im,2);
            model_input(ii,im).glen_a_ref_prefactor = params(im,3);
            model_input(ii,im).melt_rate_prefactor_exponent = params(im,4);
            model_input(ii,im).per_century_trend = params(im,5);
            model_input(ii,im).bump_amplitude = params(im,6);
            model_input(ii,im).bump_duration = params(im,7);

            meta_data(ii,im).iteration = iter;
            meta_data(ii,im).member = mem;
        end %end loop over members
    end %end loop over iterations

    % get the LOOCV data
    loocv_data = readmatrix(strcat('../mcmc_output/realization', padded_realization, '_LOOCV_grv2015.csv'));

    % get the mcmc data
    mcmc_data = readmatrix(strcat('../mcmc_output/mcmc_output_realization', padded_realization, '.csv'));

end %end gendata flag

%% Panels (a) and (b): convergence of the EKI
colmap = cmocean('ice',length(iterations)+2 );
obscolor = [1,0.,0.];

for ii = 1:length(iterations)
    for im = 1:length(members)
        plot(ax(9), model_output(ii,im).time, model_output(ii,im).grv/1e14,'linewidth', 1.5, 'Color',colmap(ii,:) );
        plot(ax(10), model_output(ii,im).time, model_output(ii,im).gl_retreat/1e3,'linewidth', 1.5, 'Color',colmap(ii,:) );

    end
end %end loop over iterations

% add observations
obs       = readmatrix("../observations/truth_actual_cts.csv");
obs_times = readmatrix("../observations/truth_times.csv") + start_year;
obs_err   = readmatrix('../observations/noise_actual.csv');


plot(ax(9),obs_times(3)*[1,1], (obs(3)+obs_err(3)*[-1,1])/1e14, 'color', obscolor, 'linewidth', 1.5);
plot(ax(9),obs_times(3), obs(3)/1e14, 'ro', 'markersize', 10, 'markerfacecolor', obscolor);


plot(ax(10),obs_times(1), (obs(1) - gl_pos_cts(1))/1e3, 'ro', 'markersize', 10, 'markerfacecolor', obscolor);
plot(ax(10),obs_times(1)*[1,1], (obs(1) - gl_pos_cts(1) + obs_err(1)*[-1,1])/1e3 , 'color', obscolor, 'linewidth', 1.5);
plot(ax(10),obs_times(2), (obs(2) - gl_pos_cts(1))/1e3, 'ro', 'markersize', 10, 'markerfacecolor', obscolor);
plot(ax(10),obs_times(2)*[1,1], (obs(2) - gl_pos_cts(1) + obs_err(2)*[-1,1])/1e3 , 'color', obscolor, 'linewidth', 1.5);



%tidy stuff
ax(9).XLabel.String = "year";
ax(9).YLabel.String = "grounded volume (10^{14} m^3)";
ax(10).XLabel.String = "year";
ax(10).YLabel.String = "grounding line retreat (km)";

ax(9).YLim = [4.3,5];
ax(10).YLim = [-4,50];

% add a colorbar
c = colorbar(ax(9));
c.Position(3) = 0.01;
c.Position(2) = 0.75;
c.Position(1) = 0.13;
c.Position(4) = 0.15;
set(c, 'Colormap',colmap(1:length(iterations), :))

dtick = 1/length(iterations);
c.Ticks = [dtick/2:dtick:(1-dtick/2)];
c.TickLabels = 1:length(iterations);
c.Label.String = 'iteration';
c.FontSize = 14;

%% Panels (c)--(h): convergence of parameters
sz = 30; %point size

% add the obs as background
xl = [0.75, 1.4;
    0.36, 2;
    -1,1.2;
    -500, 500;
    -200, 400;
    0,12];

for i = 1:6
    plot(ax(i+2), xl(i,:), [obs(3), obs(3)]/1e14,'k--',  'color', obscolor, 'linewidth', 1.5);
end

% add the scatter of points
for ii = 1:length(iterations)
    for im = 1:length(members)

        scatter(ax(3), model_input(ii,im).weertman_c_prefactor, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');
        scatter(ax(4), model_input(ii,im).glen_a_ref_prefactor, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');
        scatter(ax(5), model_input(ii,im).melt_rate_prefactor_exponent, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');
        scatter(ax(6), model_input(ii,im).per_century_trend, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');
        scatter(ax(7), model_input(ii,im).bump_amplitude, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');
        scatter(ax(8), model_input(ii,im).bump_duration, model_output(ii,im).grv2015/1e14, sz, colmap(ii,:), 'filled');

    end %end loop members
end %end loop over iterations


ax(3).XLabel.String = "sliding prefactor"';
ax(4).XLabel.String = "viscosity prefactor";
ax(5).XLabel.String = "melt prefactor";
ax(6).XLabel.String = "per century trend";
ax(7).XLabel.String = "1940s event amplitude";
ax(8).XLabel.String = "1940s event duration";

for i = 3:8
    ax(i).YLim = [4.3,5];
    ax(i).XLim = xl(i-2,:);
end

ax(3).YLabel.String = "grv (10^{14} m^3)";
ax(6).YLabel.String = "grv (10^{14} m^3)";

%% Panel (i): model-obs error as a function of iteration

allrel = 1:14;
allerrs = nan(length(allrel), length(iterations), length(members));
for ir = 1:length(allrel)
for ii = 1:n_iterations
    iter = iterations(ii);
    padded_realization = sprintf("%03d", allrel(ir));
    padded_iteration  = sprintf("%03d", iter);

    for im = 1:n_members
        mem = members(im);
        padded_member  = sprintf("%03d", mem);
        output_path = strcat("../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/outputs.csv");
        A = readmatrix(output_path);
        allerrs(ir,ii, im) = A(3)*1e12;
    end
end
end

mean_errs = mean(abs(allerrs), 3); %take the mean over the ensemble members

for ir = 1:length(allrel)
    p = plot(ax(1), iterations, mean_errs(ir, :)/1e14, 'k', 'linewidth', 1.5);
    p.Color(4) = 0.5; %set the plot alpha
end

ax(1).YLabel.String = 'MAE in grounded volume (10^{14} m^3)';
ax(1).XLabel.String = 'iteration';

%% Panel (j): histogram of EKI versus LHC

% get the LHC data
LHCpath = "../model-inputs-and-outputs/realization015_lhc/member";
LHCmem  = 0:99;

LHCerrs = nan(1,length(LHCmem));
for i = 1:length(LHCmem)
    fname = strcat(LHCpath, num2str(LHCmem(i)), '/outputs.csv' );
    A = readmatrix(fname); %outputs are in
    LHCerrs(i) = A(3)*1e12;

end

histogram(ax(2), abs(LHCerrs)/1e14, 'BinWidth', 0.02, 'FaceAlpha',0.5, 'FaceColor',  [0,47,167]/255 ) %lots of bins bc big errs 

% get the corresponding EKI data (note different realization potentially)
EKIpath = "../model-inputs-and-outputs/realization015_lhc/";
EKI_rel = 15;
EKI_errs = nan(1,length(LHCmem));
count =1 ;
for ii = 1:n_iterations
    iter = iterations(ii);
    padded_realization = sprintf("%03d", EKI_rel);
    padded_iteration  = sprintf("%03d", iter);

    for im = 1:n_members
        mem = members(im);
        padded_member  = sprintf("%03d", mem);
        output_path = strcat("../model-inputs-and-outputs/realization", padded_realization,"/iteration", padded_iteration, "/member", padded_member, "/outputs.csv");
        A = readmatrix(output_path);
        EKI_errs(count) = A(3)*1e12;
        count = count + 1;
    end
end

histogram(ax(2), abs(EKI_errs)/1e14, 'BinWidth', 0.02, 'FaceAlpha',0.5, 'FaceColor',  [167, 0,47]/255 )

legend(ax(2), {"LHC", "EKI"})
ax(2).XLim = [0, 0.5];
ax(2).YLim = [0,21];
ax(2).XLabel.String = 'absolute error in grounded volume (10^{14} m^3)';
ax(2).YLabel.String = 'count';