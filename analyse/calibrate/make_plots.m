% make a plot of (1) grounding line trajectory and (2) vaf thru time with observations

function make_plots(realization, max_iter, members)

subplot(1,2,1); hold on; box on;
plot_trajectories_by_iteration(realization, max_iter,members);

subplot(1,2,2); hold on; box on;
plot_ice_volume_by_iteration(realization, max_iter,members); %it would be nice to have like the ensemble mean and std on here too...

end
