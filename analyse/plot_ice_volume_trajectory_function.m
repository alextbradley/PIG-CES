function [ss,p] = plot_ice_volume_trajectory_function(realization, iteration, members, col, obsflag, ax, mode)
%make a plot of grounding line trajectories for members

ss = get_ice_volume_trajectories(realization, iteration, members, mode);

for im = 1:length(members)
	p(im) = plot(ax, ss(im).t + 1750, ss(im).ice_volume/1e14, 'color', col, 'linewidth', 1.5);
%	drawnow; pause
end %end loop over members

%add the obs
if obsflag
	if strcmp(mode, 'volume')
		obs = csvread('../observations/truth_volume.csv');
		obs_times = csvread('../observations/truth_times.csv');
		obs_err = csvread('../observations/noise_volume.csv');
	elseif strcmp(mode, 'grounded_volume')
		obs = csvread('../observations/truth_grounded_volume.csv');
		obs_times = csvread('../observations/truth_times.csv');
		obs_err = csvread('../observations/noise_grounded_volume.csv');
	elseif strcmp(mode, 'vaf')
		obs = csvread('../observations/truth_vaf.csv');
		obs_times = csvread('../observations/truth_times.csv');
		obs_err = csvread('../observations/noise_vaf.csv');
	end



	obs_times = obs_times+1750;
for i = 1:length(obs_times)
plot(obs_times(i)*[1,1], obs(i)+obs_err(i)*[-1,1], 'k', 'linewidth', 1.5);
end
plot(obs_times, obs, 'ko', 'markersize', 10, 'markerfacecolor', 'k');

end %end obsflag


%ax.YLim = [7.5,8.7];
ax.XLim = [1750, 2050];

ax.XLabel.String ='time';
if strcmp(mode, 'volume')
ax.YLabel.String = 'ice volume (/1e14) m^3';
elseif strcmp(mode, 'grounded_volume');
ax.YLabel.String = 'grounded ice volume (/1e14) m^3';
elseif strcmp(mode, 'vaf')
ax.YLabel.String = 'vaf (/1e14) m^3';
end
