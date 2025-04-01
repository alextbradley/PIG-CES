%make figure 2 of the manuscript, showing (a) the ambient temperature and
%(b) salinity profile and (c) individual realizations of forcing. 

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.

%% Set up figure

fig = figure(1); clf; 
fig.Position(3:4) = [900, 330];

positions = [0.07, 0.13, 0.15, 0.8;
             0.23, 0.13, 0.15, 0.8;
             0.46, 0.13, 0.48,0.8];

for i = 1:3
    ax(i) = subplot('Position', positions(i,:));
    box(ax(i), 'on');
    ax(i).FontSize = 14;
    grid(ax(i), 'on')
    hold(ax(i), 'on');

end

colS = [0,47,167]/255;
colT = [167, 0, 47]/255;

%% Make salinity and temperature profiles

zz = -1000:10:0;

T_lower = 1.2;
T_upper = -1;
S_lower = 34.6;
S_upper = 34.0;

d_low = -700;
d_hi  = -300;

plot(ax(1), two_layer_function(zz,T_lower,T_upper,d_low,d_hi),zz, "Color",colT, 'LineWidth',2);
plot(ax(2), two_layer_function(zz,S_lower,S_upper,d_low,d_hi),zz, "Color",colS, 'LineWidth',2);

%% Make the realizations of forcing plots
realizations = 1:14;


for ir = 1:length(realizations)
    fpath = strcat("../model-inputs-and-outputs/realization",  sprintf('%03d', realizations(ir)), "/realization.mat");
    rr = load(fpath);

    pc = rr.pycnocline_center;
    time = rr.time; 

    %restrict a bit for each of plotting
    pc = pc(1:10:end);
    time = time(1:10:end);
    time = time + 1750;

    if ir == 1
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5,'Marker','none');
    else        
        pp = plot(ax(3), time, pc, 'k', 'LineWidth',1.5, "HandleVisibility","off", 'Marker', 'none');

    end

    pp.Color(4) = 0.1; %set the alpha

    if ir == 1
        ens_mean = pc;
    else
        ens_mean = ens_mean + pc;
    end

end

ens_mean = ens_mean/length(realizations);
plot(ax(3), time, ens_mean, 'k', 'linewidth', 2);
ax(3).YLim = [-610, -390];
ax(3).YLabel.String = 'pycnocline depth (m)';


% add the trend and bump
colR = [0,47,230]/255; %right axis colour
yyaxis(ax(3), 'right')
bump =  50 * exp(-(time - 1945).^2 / 2 / 6^2);
trend = 0.5* (time - 1960).*(time > 1960);
plot(ax(3), time, bump, 'Color',colR, 'LineWidth',2)
plot(ax(3), time, trend,'k--',  'Color',colR, 'LineWidth',2)
ax(3).YLim = [-110, 110];
ax(3).YAxis(2).Color = colR;
ax(3).YLabel.String = 'pycnocline depth perturbation (m)';

legend({'Individual realizations of forcing', 'Ensemble mean', '1940s event', 'Anthropogenic trend'}, 'FontSize', 14, 'Location','SouthWest')
%% tidy stuff
for i = 1:3
    ax(i).FontSize = 14;
    grid(ax(i), 'on')
    ax(i).FontName = 'Arial';
end
ax(2).YTick = ax(1).YTick;
ax(2).YTickLabel = [];
ax(1).XLim = [-1.2, 1.3];
ax(2).XLim = [33.9, 34.7];
ax(1).XLabel.String = 'temp (C)';
ax(2).XLabel.String = 'salinity (PSU)';
ax(1).YLabel.String = 'depth (m)';
ax(3).XLabel.String = 'time';

%uncomment to export
%exportgraphics(fig, "figures/raw/figure2.pdf", 'ContentType','vector');

function Aa = two_layer_function(z,v_low,v_hi,d_low,d_hi)
    Aa = v_low .* (z < d_low) + v_hi .* (z > d_hi) + (v_low + (v_hi - v_low)/(d_hi - d_low) * (z - d_low)).* ((z >= d_low) & (z <= d_hi));

end