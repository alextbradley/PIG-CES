function [ss,p] = plot_ice_volume_trajectory_function(realization, iteration, members, col, obsflag, ax)
%make a plot of grounding line trajectories for members

ss = get_ice_volume_trajectories(realization, iteration, members);

for im = 1:length(members)
	p(im) = plot(ax, ss(im).t + 1750, ss(im).ice_volume, 'color', col, 'linewidth', 1.5);
%	drawnow; pause
end %end loop over members

%add the obs
if obsflag
	obs = csvread('../observations/truth_actual.csv');
	obs = obs(3);
	obs_times = csvread('../observations/truth_times.csv');
	obs_times = obs_times(3);
	obs_err = csvread('../observations/noise_actual.csv');
	obs_err = obs_err(3);
	obs_times = obs_times+1750;
	
for i = 1
plot(obs_times(i)*[1,1], obs(i)+obs_err(i)*[-1,1], 'k', 'linewidth', 1.5);
end
plot(obs_times, obs, 'ko', 'markersize', 10, 'markerfacecolor', 'k');

end %end obsflag


ax.YLim = [4.3,5.0]*1e14;
ax.XLim = [1750, 2050];

ax.XLabel.String ='time';
ax.YLabel.String = 'grounded ice volume (/1e14) m^3';

