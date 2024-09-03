% make a plot of the posterior climate

% load in the data
realization = "026";
fname = strcat("emulate/mcmc_output/mcmc_output_realization", realization, ".csv");
A = readmatrix(fname);

%
input_headers = ["weertman c prefactor" ,"glen a ref prefactor",...
    "melt rate prefactor exponent","per century trend","bump amplitude","bump duration"];

% make histogram and kde
lims = [0.5, 2;
    0.5, 2;
    -2,  2;
    -400, 400;
    -200, 600;
    0, 10];  %limits of the kde

np = 1000;
xx = zeros(np,6);

figure(1);clf;
ss = struct;

for i = 1:6
    ax(i) = subplot(2,3,i); hold(ax(i), 'on'); box(ax(i), 'on');
    histogram(ax(i), A(:,i), 50, 'Normalization','pdf', 'LineStyle','none');

    xx(:,i) = linspace(lims(i,1),lims(i,2), np);


    pts = xx(:,i);
    kde = fitdist(A(:,i),'kernel');
    y = pdf(kde,pts);
    ss(i).kde = kde;
   
    plot(ax(i), pts, y, 'k','LineWidth',1.5)

    ax(i).XLim = lims(i,:);

    ax(i).XLabel.String = input_headers(i);

    yl = ax(i).YLim;
    plot(ax(i), mean(kde)*[1,1],yl, 'k--', 'linewidth', 1.5);
    plot(ax(i), (mean(kde)-std(kde))*[1,1],yl, 'k--', 'linewidth', 1.5, 'color', [0.5,0.5,0.5]);
    plot(ax(i), (mean(kde)+std(kde))*[1,1],yl, 'k--', 'linewidth', 1.5, 'color', [0.5,0.5,0.5]);
    

end

%% Posterior climate
%random part (deterministic)
fname = strcat('../model-inputs-and-outputs/realization', realization, "/realization.mat");
data = load(fname);
t = data.time;
t = t + 1750;


%trend
kde_trend = ss(4).kde;
trend_upper   = mean(kde_trend)+std(kde_trend);
trend_central = mean(kde_trend);
trend_lower   = mean(kde_trend)-std(kde_trend);


%bump amplitude
kde_bumpamplitude = ss(5).kde;
bumpamplitude_upper   = mean(kde_bumpamplitude)+std(kde_bumpamplitude);
bumpamplitude_central = mean(kde_bumpamplitude);
bumpamplitude_lower   = mean(kde_bumpamplitude)-std(kde_bumpamplitude);

%bump duration
kde_bumpduration = ss(6).kde;
bumpduration_upper   = mean(kde_bumpduration)+std(kde_bumpduration);
bumpduration_central = mean(kde_bumpduration);
bumpduration_lower   = mean(kde_bumpduration)-std(kde_bumpduration);

figure(2); clf; hold on; box on;


%add the random part
p = plot(t, data.pycnocline_center, 'r', 'linewidth', 1.5);
p.Color = [1,0,0, 0.35];

% add the trend part
p = plot(t, -500+ (trend_central/100* (t-1960) .* (t > 1960)), 'g', 'linewidth', 1.5);
p.Color = [0,0,1, 0.35];

% add the bump part
p = plot(t, -500 + bumpamplitude_central*exp(-(t - 1945).^2 /2 /bumpduration_central^2), 'LineWidth',1.5);
p.Color = [0,1,0, 0.5];

rf_central =  get_forcing(t,data.pycnocline_center,trend_central,bumpduration_central,bumpamplitude_central);
rf_lower   =  get_forcing(t,data.pycnocline_center,trend_lower,bumpduration_lower,bumpamplitude_lower);
rf_upper   =  get_forcing(t,data.pycnocline_center,trend_upper,bumpduration_upper,bumpamplitude_upper);


xf = [t;flip(t)];
yf = [rf_lower ;flip(rf_upper)];

fill(xf, yf, 'k', 'FaceAlpha',0.2, 'LineStyle','none');
plot(t, rf_central, 'k', 'linewidth', 1.5);


plot([min(t), max(t)], [-500,-500], 'k--', 'linewidth', 1.5)
%plot(t, rf_lower, 'r', 'linewidth', 1.5);
%plot(t, rf_upper, 'g', 'linewidth', 1.5);



legend( {'random part', 'trend','posterior climate', 'bounds'}, 'location', 'northwest');






ax = gca;
ax.FontSize = 14;
ax.YLim = [-650,-250];
ax.XLabel.String = 'year';
ax.YLabel.String = 'pycnocline position (m)';

function rf = get_forcing(t,random_part,trend,bump_duration,bump_amplitude)

    rf = random_part; %seed to the random part

    %add trend part
    rf = rf + trend/100* (t-1960) .* (t > 1960);

    %add bump part
    bump_time = 1945; 
    rf = rf + bump_amplitude*exp(-(t - bump_time).^2 /2 /bump_duration^2);
end
