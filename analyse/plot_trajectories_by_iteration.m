function plot_trajectories_by_iteration(realization, max_iter,members)

x0 = 64;
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
        [ss,p]= plot_grounding_line_trajectory_function(realization, iteration(ii), member, x0, col, obsflag, ax);


end
