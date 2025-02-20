function [ss,p] = plot_grounding_line_trajectory_function(realization, iteration, members, x0, col, obsflag, ax)
%make a plot of the grounding line trajectories for members. x0 is the line along which to measure the gl. col is the linecolor for the plots. obsflag is a flag to add the observations

ss = get_grounding_line_trajectories(realization, iteration, members, x0);

for im = 1:length(members)
	p(im) = plot(ax, ss(im).t + 1750, ss(im).gl_pos, 'color', col, 'linewidth', 1.5); 
end %end loop over members

%add the obs
if obsflag
obs = csvread('../../observations/truth_actual.csv');
obs_times = csvread('../../observations/truth_times.csv');
obs_times = obs_times+1750;
obs_err = csvread('../../observations/noise_actual.csv');
for i = 1:length(obs_times)
plot(obs_times(i)*[1,1], obs(i)+obs_err(i)*[-1,1], 'k', 'linewidth', 1.5);
end
plot(obs_times, obs, 'ko', 'markersize', 10, 'markerfacecolor', 'k');

end %end obsflag

ax.XLabel.String ='time';
ax.YLabel.String = 'gl position';
ax.YLim = [-3,-2.5]*1e5;
ax.XLim = [1750, 2050];


