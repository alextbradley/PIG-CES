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

plot_trajectories = 1; %make a 5 x 4 plot of trajectories (one panel for each realization of forcing)
plot_slices 	  = 0; %make plots of slices for each realization and each iteration

realization = ["001","002", "003","004", "005", "006", "007", "008", "009", "010", "011", "012", "013", "014", "015", "016", "017", "018", "019", "020"]; %this shouldn't change, we'll figure them out as we go which iterations are available

iteration = ["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"];

member = ["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"];
colmap = parula(length(iteration)); %iterations
x0 = 64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_trajectories

fig = figure(1);

is_worth_slicing = zeros(length(realization), length(iteration)); %boolean array to decide which slices to output 
for ir = 1:length(realization)
ax(ir) = subplot(4,5,ir); hold on; box on
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

end %end flag

%%%%%%%%%%%%%%%%%% %individual slices %%%%%%%%%%%%%%%%%%%%%
if plot_slices
members = ["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"];
% specify where to take the slice
x0 = 64; y0 = 20;
x1 = 64; y1 = 60;

for ir = 1:length(realization)
for ii = 1:length(iteration)

[ir, ii]
if is_worth_slicing(ir,ii)
fig(ir) = figure(ir); clf;
fig(ir).Position(3:4) = [1540, 800];
figname = strcat("./ensemble_figures/", "slices_realization", realization(ir), "_iteration", iteration(ii), ".pdf");
for im = 1:length(members)
	axs(im) = subplot(2,5,im);
	fname = strcat("realization", realization(ir), "_iteration", iteration(ii),"_member", members(im), ".pdf");
	plot_slice_function(realization(ir), iteration(ii), members(im), [x0,x1], [y0,y1], axs(im));
	axs(im).Title.String = strcat("member", members(im));
	axs(im).Title.Interpreter = 'none';
end %end loop over members

s = suptitle(strcat("realization", realization(ir), ", iteration", iteration(ii)));
end %end is worth slicing flag

exportgraphics(fig(ir), figname);
close(fig(ir))
end %end loop over iterations
end %end loop over realizations 
end %end if plot_slices
