% Plot a slice through the ice

realization = "001";
iteration   = "001";
member      = "001";

% specify where to take the slice
x0 = 64; y0 = 20;
x1 = 64; y1 = 60;

fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration,"-", member, "/run/outfile.nc");

t = ncread(fpath, "TIME");
s = ncread(fpath, "s");
h = ncread(fpath, "h");
b = ncread(fpath, "b"); b = squeeze(b(:,:,1));


tout_idx = 1:20:length(t); %time indices to output

figure(1); clf; hold on; box on;
plot(b(x0,y0:y1),  'color', [165,42,42]/255, 'linewidth', 1.75);
colmap = parula(length(tout_idx));

for i = 1:length(tout_idx)
plot(squeeze(s(x0,y0:y1,tout_idx(i))), 'color', colmap(i,:), 'linewidth', 1.75)
plot(squeeze(s(x0,y0:y1,tout_idx(i))-h(x0,y0:y1,tout_idx(i))), 'color', colmap(i,:), 'linewidth', 1.75)
end

box on;
c = colorbar;

c.Ticks = 0:1/(length(tout_idx)-1) : 1;
c.TickLabels = t(tout_idx);
c.Label.String = 'time (years)';
xlim([0, y1-y0])
