function plot_ice_volume_by_iteration(realization, max_iter,members,mode)
% make a plot of the ice mass as a function of time for a given realization
% mode specifies the type of plot as follows:
% 'volume': ice grounded volume
% 'vaf': ice volume above floatation


realization = sprintf('%03d', realization);
iteration = strings(1,max_iter);
for i = 1:max_iter
	iteration(i) =  sprintf('%03d', i);
end

member = strings(1,max_iter);
for i = 1:length(members)
        member(i) =  sprintf('%03d', members(i));
end

hold on; ax = gca;
colmap = parula(max_iter+2);

for ii = 1:length(iteration)
	col = colmap(ii,:); %colour of the lines
        obsflag = 0;
        if ii == length(iteration)
                obsflag = 1;
        end
        [ss,p]= plot_ice_volume_trajectory_function(realization, iteration(ii), member, col, obsflag, ax, mode);

end

