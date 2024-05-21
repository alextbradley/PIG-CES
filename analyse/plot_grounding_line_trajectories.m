% Plot the trajectories for a set of simulations. Colour the trajectories by the iteration

x0 = 64;     %where to measure gl

%run spec
realization = ["001","002", "003","004", "005", "006", "007"];
realization = ["010"];
iteration   = ["001"];
member      = ["001","002","003","004","005","006","007","008","009","010"];

iter_colmap = parula(length(iteration)+2);
member_colmap = parula(length(member));

for ir = 1:length(realization)
fig = figure(ir); clf; hold on; box on;
ax = gca;
title(strcat("realization ", realization(ir)))
for ii = 1:length(iteration)
	col = iter_colmap(ii,:); %colour of the lines
	obsflag = 0; 
	if ii == length(iteration)
		obsflag = 1;
	end
	[ss,p]= plot_grounding_line_trajectory_function(realization(ir), iteration(ii), member, x0, col, obsflag, ax);
end
end

