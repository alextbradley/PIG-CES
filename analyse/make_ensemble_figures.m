% Make ensemble figures. This script produces the following
% figures:
%
% (1) a 5x3 subplots showing the trajectories of gl
% position, with each subplot corresponding to a different
% realization of forcing (1-15). This is saved as 
% ./figures/trajectories.pdf
%
% (2) for each realization of forcing, and for each 
% iteration, a plot of the 'slice' for each member (i.e.
% a 5 x 2 subplots). These are saved as 
% ./figures/realizationXXX_iterationXXX_slices.pdf
%

realization = ["001","002", "003","004", "005", "006", "007", "008", "009", "010", "011", "012", "013", "014", "015"]; %this shouldn't change, we'll figure them out as we go which iterations are available

iteration = ["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"];
colmap = parula(length(iteration)); %iterations

fig = figure(1);

is_worth_slicing = zeros(length(realization), length(iteration)); %boolean array to decide which slices to output 
for ir = 1:length(realization)
ax(ir) = subplot(3,5,ir); hold on; box on
title(strcat("realization ", realization(ir)))
ir

for ii = 1:length(iteration)
        col = colmap(ii,:); %colour of the lines
        obsflag = 0;
        if ii == length(iteration)
                obsflag = 1;
        end
        [ss,p]= plot_grounding_line_trajectory_function(realization(ir), iteration(ii), member, x0, col, obsflag, ax(ir));
	
	if any(~isnan([ss(:).t]))
		is_worth_slicing(ir, ii) = 1; 
	end	
	
end
end
fig = gcf;
fig.Position(3:4) = [1600, 800];
exportgraphics(gcf, './ensemble_figures/trajectories.pdf')

%%%%%%%%%%%%%%%%%% %individual slices %%%%%%%%%%%%%%%%%%%%%


 
