function plot_slice_function_mat(prefix, x, y, ax)
%plot a slice taken along co-ordinates x = [x1,x2], y = [y1,y2] with inputs from mat files. Prefix is the folder prefix. CURRENTLY ONLY COMPATIBLE WITH x1 = x2!! 

jdir = dir(strcat(prefix, "/run/*.mat"));

x0 = x(1); x1 = x(2);
y0 = y(1); y1 = y(2);


t = nan(1,length(jdir));
s = nan(268,169, length(jdir));
h = nan(268,169, length(jdir));
b = nan(268,169, length(jdir));

prefix
for i = 1:length(jdir)
    fname = strcat(prefix, "/run/" ,jdir(i).name);
    data = load(fname);
    s(:,:,i) = data.s;
    h(:,:,i) = data.h;
    b(:,:,i) = data.b;
    t(i) = data.t;
    
end
b = squeeze(b(:,:,1));


tout_idx = 1:20:length(t); %time indices to output

hold on; box on;

%title(strcat("realization",realization, ", iteration", iteration, ", member ", member))
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

