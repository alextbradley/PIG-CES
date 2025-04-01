% Make figure 4 of the manuscript, showing (a)--(c) the leave one out cross
% validation for each of the three emulated quantities and (d) a bar chart
% of the coverage over each emulator.

% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.


%% Preliminaries
%

realization = 1; 
fnames = ["grv2015", "gl1930", "gl2015"];
fpath = '../mcmc_output/';

fig = figure(1); clf;

% setup axis limits
axxlim = [4.3,5;
    -15, 45;
    -45, 15];
axylim = axxlim;

axylim = [4.2,5.2;
    -30, 60;
    -60, 30];


for i = 1:3
    ax(i) = subplot(1,3,i);
    hold(ax(i), 'on');
    box(ax(i), 'on');
    %grid on
    ax(i).XLim = axxlim(i,:);
    ax(i).YLim = axylim(i,:);
    %ax(i).YTick = ax(i).XTick;
    ax(i).FontSize = 14;
    grid(ax(i), 'on');
end
fig.Position(3:4) = [1300,330];


%setup plot colours
inrange_color = [42/255, 103/255, 131/255];
outrange_color = [228/255,128/255,111/255];


observations = readmatrix("../observations/truth_actual.csv");
observations = observations([3,1,2]); %permute for right order

%%
for i = 1:3

    %read the data in
    fname = strcat(fpath, "realization00", num2str(realization), "_LOOCV_", fnames(i), ".csv");
    data = readmatrix(fname);

    output_normalization_mean = mean(data(:,1));
    output_normalization_sd   = std(data(:,1));


    modelled = data(:,1);
    emulated_central = data(:,2);
    emulated_lower   = data(:,3);
    emulated_upper   = data(:,4);

    %fix the baselines
    if i == 1
        modelled = (modelled*1e12 + observations(1))/1e14;
        emulated_central = (emulated_central*1e12 + observations(1))/1e14; 
        emulated_lower   = (emulated_lower*1e12 + observations(1))/1e14;
        emulated_upper   = (emulated_upper*1e12 + observations(1))/1e14;
    else
        modelled = modelled * 3; % because these are in grid cells of size 3km (we're going to plot in km later)
        emulated_central = emulated_central * 3;
        emulated_lower   = emulated_lower * 3;
        emulated_upper   = emulated_upper * 3;
    end

    % add 1-1 line
    plot(ax(i), axxlim(i,:), axxlim(i,:), 'k--', 'LineWidth',1.5);

    coverage = 0;

    for j = 1:length(modelled)
        if (emulated_lower(j) < modelled(j)) && (emulated_upper(j) > modelled(j))
            coverage = coverage + 1;
            plot(ax(i), modelled(j), emulated_central(j),'o', 'markerfacecolor', inrange_color, 'MarkerEdgeColor',inrange_color);

            %add error bar
            plot(ax(i), modelled(j)*[1,1],[emulated_lower(j), emulated_upper(j)], 'color', inrange_color, 'LineWidth',1 );
        else
            plot(ax(i), modelled(j), emulated_central(j),'o', 'markerfacecolor', outrange_color, 'MarkerEdgeColor',outrange_color);

            %add error bar
            plot(ax(i), modelled(j)*[1,1],[emulated_lower(j), emulated_upper(j)], 'color', outrange_color, 'LineWidth',1 );

        end

    end
    coverages(i) = coverage/length(modelled);
    RMSEs(i)     = sqrt(sum((modelled - emulated_central).^2)/length(modelled));

    % add the observations
    ms = 15;
    obscolor = [207, 0, 47]/255; 
    if i == 1
        plot(ax(i), observations(i)/1e14,  observations(i)/1e14, 'ro', 'MarkerFaceColor', obscolor, 'markersize', ms, 'markeredgecolor', obscolor);

    elseif i == 2
        plot(ax(i), 0, 0 , 'ro', 'MarkerFaceColor',obscolor, 'markersize', ms,  'markeredgecolor', obscolor);
    else
        plot(ax(i), (observations(2) - observations(3))/1e3, (observations(2) - observations(3))/1e3 , 'ro', 'MarkerFaceColor',obscolor, 'markersize', ms,  'markeredgecolor', obscolor);

    end

end %end loop over 3 observation points



%% tidy stuff
ax(1).XLabel.String = 'modelled 2015 grounded volume (10^{14} m^3)';
ax(1).YLabel.String = 'emulated 2015 grounded volume (10^{14} m^3)';
ax(2).XLabel.String = 'modelled 1930 grounding line retreat (km)';
ax(2).YLabel.String = 'emulated 1930 grounding line retreat (km)';
ax(3).XLabel.String = 'modelled 2015 grounding line retreat (km)';
ax(3).YLabel.String = 'emulated 2015 grounding line retreat (km)';


%% print stuff
coverages
RMSEs
%exportgraphics(fig, "figures/raw/figure4abc.pdf", 'ContentType','vector');

%% Make panel (d): histogram of coverages for each realizations
realizations = 1:14;
coverages = nan(length(realizations), 3);
for ir = 1:length(realizations)
for i = 1:3

    %read the data in
    fname = strcat(fpath,  sprintf("realization%03d", realizations(ir)), "_LOOCV_", fnames(i), ".csv");
    data = readmatrix(fname);

    output_normalization_mean = mean(data(:,1));
    output_normalization_sd   = std(data(:,1));


    modelled = data(:,1);
    emulated_central = data(:,2);
    emulated_lower   = data(:,3);
    emulated_upper   = data(:,4);

    %fix the baselines
    if i == 1
        modelled = (modelled*1e12 + observations(1))/1e14;
        emulated_central = (emulated_central*1e12 + observations(1))/1e14; 
        emulated_lower   = (emulated_lower*1e12 + observations(1))/1e14;
        emulated_upper   = (emulated_upper*1e12 + observations(1))/1e14;
    else
        modelled = modelled * 3; % because these are in grid cells of size 3km (we're going to plot in km later)
        emulated_central = emulated_central * 3;
        emulated_lower   = emulated_lower * 3;
        emulated_upper   = emulated_upper * 3;
    end


    coverage = 0;

    for j = 1:length(modelled)
        if (emulated_lower(j) < modelled(j)) && (emulated_upper(j) > modelled(j))
            coverage = coverage + 1;
        end

    end
    coverages(ir,i) = coverage/length(modelled);
    RMSEs(i)     = sqrt(sum((modelled - emulated_central).^2)/length(modelled));



end %end loop over 3 observation points

end %end loop over realizations
fig2 = figure(2); clf; hold on
plot([0,15], [90, 90], 'k--', 'LineWidth',1., 'HandleVisibility','off')
plot([0,15], [80, 80], 'k--', 'LineWidth',1.,'HandleVisibility','off')
fig2.Position(3) = fig.Position(3);
b = bar(coverages*100 );
ax2 = gca;box on
ax2.FontSize = 14;
ax2.XLabel.String = 'realization';
ax2.YLabel.String = 'coverage (%)';
l=legend({'2015 grounded volume', '1930 grounding line retreat', '2015 grounding line retreat'});
%l.Orientation = 'horizontal';
l.Location = 'southeast';
ax2.XLim = [0.5, 14.5];
ax2.XTick = realizations;

%exportgraphics(fig2, "figures/raw/figure4d.pdf", 'ContentType','vector');
