% Make combined plot showing the CES procedure.
% (a) Calibration of the ice loss
% (b) Emulation of the ice loss
% (c) Posterior distributions of model parameters
%

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.


%% Preliminaries
gendata = 1;
addpath('../functions/');

%specify ensemble
realization = 3;
iterations  = 1:5;
members     = 1:20;

start_year = 1750;

numplot = 1;
fontsize = 14;


% set up the figure
fig = figure(1); clf; 
fig.Position(3:4) = [860,860];

gap = 0.02;
h = (0.42-5*gap)/6;

positions = [0.1,0.56, 0.4, 0.4;
            0.58,0.56,0.4, 0.4;
            0.1, 0.06, 0.3, h;
            0.1, 0.06 + h + gap, 0.3, h;
            0.1, 0.06 + 2*h+2*gap, 0.3, h;
            0.1, 0.06 + 3*h+3*gap, 0.3, h;
            0.1, 0.06 + 4*h+4*gap, 0.3, h;
            0.1, 0.06 + 5*h+5*gap, 0.3, h;
            0.45, 0.06, 0.53, h;
            0.45, 0.06 + h+gap, 0.53, h;
            0.45, 0.06 + 2*h+2*gap, 0.53, h;
            0.45, 0.06 + 3*h+3*gap, 0.53, h;
            0.45, 0.06 + 4*h+4*gap, 0.53, h;
            0.45, 0.06 + 5*h+5*gap, 0.53, h];


for i = 1:14
    ax(i) = subplot('Position', positions(i,:));
    box(ax(i), 'on');
    hold(ax(i), 'on');
    ax(i).FontSize =fontsize;
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
        "melt_rate_prefactor_exponent","per_century_trend","bump_amplitude","bump_duration"];
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
            model_output(ii,im).gl_pos_discrete = trajectories.gl_pos_discrete;
            model_output(ii,im).gl_pos_cts = trajectories.gl_pos_cts;

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

%% Panel (a): calibration
colmap = cmocean('ice',length(iterations)+2 );
obscolor = [1,0.,0.];

for ii = 1:length(iterations)
    for im = 1:length(members)
        plot(ax(1), model_output(ii,im).time, model_output(ii,im).grv/1e14,'linewidth', 1.5, 'Color',colmap(ii,:) );

    end
end %end loop over iterations

% add observations
obs       = readmatrix("../observations/truth_actual_cts.csv");
obs_times = readmatrix("../observations/truth_times.csv") + start_year;
obs_err   = readmatrix('../observations/noise_actual.csv');


plot(ax(1),obs_times(3)*[1,1], (obs(3)+obs_err(3)*[-1,1])/1e14, 'color', obscolor, 'linewidth', 1.5);
plot(ax(1),obs_times(3), obs(3)/1e14, 'ro', 'markersize', 10, 'markerfacecolor', obscolor);



%tidy stuff
ax(1).XLabel.String = "year";
ax(1).YLabel.String = "grounded volume (10^{14} m^3)";

ax(1).YLim = [4.3,5];

% add a colorbar
c = colorbar(ax(1));
c.Position(3) = 0.01;
c.Position(2) = 0.57;
c.Position(1) = 0.15;
c.Position(4) = 0.15;
set(c, 'Colormap',colmap(1:length(iterations), :))

dtick = 1/length(iterations);
c.Ticks = [dtick/2:dtick:(1-dtick/2)];
c.TickLabels = 1:length(iterations);
c.Label.String = 'iteration';
c.FontSize = fontsize;

%% Panel (b): emulation
%we emulate the error in 1e12 m^3 -- check this has gone into the mcmc
%correctly...
inrange_color  = [42/255, 103/255, 131/255];
outrange_color = [228/255,128/255,111/255];

modelled = loocv_data(:,1);
emulated = loocv_data(:,2);
emulated_l95 = loocv_data(:,3);
emulated_u95 = loocv_data(:,4);

plot(ax(2), [4.3,5], [4.3,5] , 'k--', 'linewidth',1)
inrange_counter = 0;
for i =1:length(loocv_data)
    if (modelled(i) < emulated_l95(i)) || (modelled(i) > emulated_u95(i))
        %add points
        plot(ax(2), (modelled(i)*1e12 + obs(3))/1e14, (emulated(i)*1e12 + obs(3))/1e14, 'o','markerfacecolor', outrange_color, 'MarkerEdgeColor',outrange_color);

        % add error bars
        plot(ax(2), [(modelled(i)*1e12 + obs(3))/1e14,(modelled(i)*1e12 + obs(3))/1e14], [(emulated_l95(i)*1e12 + obs(3))/1e14,(emulated_u95(i)*1e12 + obs(3))/1e14], 'color', outrange_color, 'LineWidth',1.5)

    else
        %add points
        plot(ax(2), (modelled(i)*1e12+ obs(3))/1e14, (emulated(i)*1e12 + obs(3))/1e14, 'o','markerfacecolor', inrange_color, 'MarkerEdgeColor',inrange_color);
        
        %add error bars
        plot(ax(2), [(modelled(i)*1e12 + obs(3))/1e14,(modelled(i)*1e12 + obs(3))/1e14], [(emulated_l95(i)*1e12 + obs(3))/1e14,(emulated_u95(i)*1e12 + obs(3))/1e14], 'color', inrange_color, 'LineWidth',1.5)

        %add one to counter
        inrange_counter = inrange_counter + 1;
    end

end

% add the obs
plot(ax(2), obs(3)/1e14, obs(3)/1e14, 'ro', 'markersize', 10, 'MarkerFaceColor','r')
%add viertical and horizontal lines
yl = ax(2).YLim; xl = ax(2).XLim;
plot(ax(2),[min(yl), obs(3)/1e14],[obs(3)/1e14,obs(3)/1e14], 'r--', 'linewidth', 1.5);
plot(ax(2),[obs(3)/1e14,obs(3)/1e14], [min(xl), obs(3)/1e14],'r--', 'linewidth', 1.5);


% tidy stuff
ax(2).XLim = [4.3,5];
ax(2).YLim = [4.3,5];
ax(2).XLabel.String = 'modelled 2015 grounded volume (x 10^{14} m^3)';
ax(2).YLabel.String = 'emulated 2015 grounded volume (x 10^{14} m^3)';

% add text w/ coverage and R^2
mdl = fitlm(modelled, emulated);
txt = sprintf("coverage: %.1f%%, \nR^2: %.3f, \nRMSE (/1e14 m^3): %.3f", inrange_counter/length(loocv_data)*100, mdl.Rsquared.Ordinary,(mdl.RMSE)/1e2);
tt  = text(ax(2), 4.32, 4.91, txt, 'Fontsize',14);


%% Panels (c): Prior and posterior distributions
prior_mean = [1.0, 1.0, 0.0, 0.0,   200.0, 5.0];
prior_sd   = [0.3, 0.3, 1.2, 200.0, 100.0, 2.5]; %priors are normally distributed

prior_colour = [134, 194,234]/255;
naive_posterior_colour = [18,75,113]/255;
mcmc_posterior_colour = [253,127, 14]/255;

for i = 1:6
    if i ~=1 
       % ax(i+2).XTickLabels = [];
    end

    % work out the median, IQR (1.35 SD for Gaussian), and 95% interval (1.96 SD for gaussian)

    % plot the prior
    plot(ax(i+2),  [prior_mean(i) - 1.96*prior_sd(i),prior_mean(i) + 1.96*prior_sd(i)],[0.7, 0.7], 'linewidth', 2, 'color', prior_colour);
    plot(ax(i+2),  [prior_mean(i) - 1.35*prior_sd(i),prior_mean(i) + 1.35*prior_sd(i)],[0.7, 0.7], 'linewidth', 4, 'color', prior_colour);

    % plot the naive posterior
    if i == 1
        params = [model_input(5,:).weertman_c_prefactor];
    elseif i == 2
        params = [model_input(5,:).glen_a_ref_prefactor];
    elseif i == 3
        params = [model_input(5,:).melt_rate_prefactor_exponent];
    elseif i == 4
        params = [model_input(5,:).per_century_trend];
    elseif i == 5
        params = [model_input(5,:).bump_amplitude];
    elseif i == 6
        params = [model_input(5,:).bump_duration];
    end
    naive_posterior_sd = std(params);
    naive_posterior_mean = mean(params);
    plot(ax(i+2),  [naive_posterior_mean - 1.96*naive_posterior_sd,naive_posterior_mean + 1.96*naive_posterior_sd],[0.5, 0.5], 'linewidth', 2, 'color', naive_posterior_colour);
    plot(ax(i+2),  [naive_posterior_mean - 1.35*naive_posterior_sd,naive_posterior_mean + 1.35*naive_posterior_sd],[0.5, 0.5], 'linewidth', 4, 'color', naive_posterior_colour);

    % plot the mcmc
    mcmc_median = median(mcmc_data(:,i));
    mcmc_iqr    = quantile(mcmc_data(:,i), [0.25, 0.5, 0.75]); mcmc_iqr = mcmc_iqr([1,3]); %remove the median
    mcmc_95     = quantile(mcmc_data(:,i), [0.05, 0.95]);

    plot(ax(i+2),  mcmc_95,[0.3, 0.3], 'linewidth', 2, 'color', mcmc_posterior_colour);
    plot(ax(i+2),  mcmc_iqr,[0.3, 0.3], 'linewidth', 4, 'color', mcmc_posterior_colour);

    ax(i+2).YLim = [0,1];
    ax(i+2).YTick= [];
    ax(i+2).TickDir = 'in';
end

if realization == 14
    ax(3).XLim = [0.25, 1.75];
    ax(4).XLim = [0.25, 1.75];
    ax(5).XLim = [-2.5 2.5];
    ax(6).XLim = [-400,400];
    ax(8).XLim = [0,12];
    ax(7).XLim = [0,400];

end




% fix the limits


%% Panels (d): MCMC traceplots
xx = 1:length(mcmc_data);
for i = 1:6
    plot(ax(i+8), xx(2:20:end)/1e4, mcmc_data(2:20:end,i), 'k');

    if i ~=1 
        ax(i+8).XTickLabels = [];
    end
end

if realization == 14
ax(14).YLim = [8.1, 8.6];
ax(13).YLim = [110,240]; ax(13).YTick = 100:50:300;
ax(12).YLim = [100, 300];
ax(11).YLim = [0.5, 0.9];
ax(10).YLim = [0.8, 1.2];
ax(9).YLim = [0.65, 0.95];

end
ax(9).XLabel.String = 'MCMC entry number (/1e4)';




