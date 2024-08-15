% make plots of the realizations of forcing

realizations = 21:40;

mean_pc = zeros(3001,1);
plot_individually = 0;
plot_ensemble_mean = 1;

if plot_individually
    for ir = 1:length(realizations)
        padded_realization = sprintf('%03d', realizations(ir));
        fname = strcat('../model-inputs-and-outputs/realization', padded_realization, "/realization.mat");
        data = load(fname);


        figure(ir); clf; hold on; box on;
        plot(data.time + 1750, data.pycnocline_center, 'k', 'linewidth', 1.5);
        plot(1750 + [min(data.time), max(data.time)], [-500,-500], 'k--', 'linewidth', 1.5)
        ax = gca;
        ax.FontSize = 14;
        ax.YLim = [-650,-350];
        ax.XLabel.String = 'year';
        ax.YLabel.String = 'pycnocline position (m)';
        title(strcat("realization ", padded_realization));

        fig = gcf; fig.Position(3:4) = [720,350];


    end
end

if plot_ensemble_mean
    mean_pc = zeros(3001,1);
    figure(1);clf; hold on; box on;
    for ir = 1:length(realizations)
        padded_realization = sprintf('%03d', realizations(ir));
        fname = strcat('../model-inputs-and-outputs/realization', padded_realization, "/realization.mat");
        data = load(fname);


        
        p = plot(data.time + 1750, data.pycnocline_center, 'linewidth', 1.5);
        p.Color(4) = 0.3;
        
        mean_pc = mean_pc + data.pycnocline_center;


    end
    mean_pc = mean_pc / ir;

    plot(1750 + [min(data.time), max(data.time)], [-500,-500], 'k--', 'linewidth', 1.5)
    p = plot(data.time + 1750, mean_pc,'k',  'linewidth', 1.5);



ax = gca;
ax.FontSize = 14;
ax.YLim = [-650,-350];
ax.XLabel.String = 'year';
ax.YLabel.String = 'pycnocline position (m)';
fig = gcf; fig.Position(3:4) = [720,350];
end