function plot_slice_function(realization, iteration, member, x, y, ax)
%plot a slice taken along co-ordinates x = [x1,x2], y = [y1,y2]. CURRENTLY ONLY COMPATIBLE WITH x1 = x2!! fpath is the path of the nc file

fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration,"-", member, "/run/outfile.nc");

x0 = x(1); x1 = x(2);
y0 = y(1); y1 = y(2);

t = ncread(fpath, "TIME");
s = ncread(fpath, "s");
h = ncread(fpath, "h");
b = ncread(fpath, "b"); b = squeeze(b(:,:,1));

tout_idx = 1:20:length(t); %time indices to output

hold on; box on;

title(strcat("realization",realization, ", iteration", iteration, ", member ", member))
plot(ax, b(x0,y0:y1),  'color', [165,42,42]/255, 'linewidth', 1.75);
colmap = parula(length(tout_idx));

for i = 1:length(tout_idx)
plot(ax, squeeze(s(x0,y0:y1,tout_idx(i))), 'color', colmap(i,:), 'linewidth', 1.75)
plot(ax, squeeze(s(x0,y0:y1,tout_idx(i))-h(x0,y0:y1,tout_idx(i))), 'color', colmap(i,:), 'linewidth', 1.75)
end

box on;
c = colorbar(ax);

c.Ticks = 0:1/(length(tout_idx)-1) : 1;
c.TickLabels = t(tout_idx);
c.Label.String = 'time (years)';
xlim([0, y1-y0])

