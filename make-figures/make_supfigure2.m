% Make figure showing the comparison between simulated velocities and
% observed for 2015
addpath("../functions/")

absolute = 1; %to set to absolute difference (1) or relative (0)
u = ncread("../model-inputs-and-outputs/forward-runs/combined_realizations.nc", 'u');
v =  ncread("../model-inputs-and-outputs/forward-runs/combined_realizations.nc", 'v');
grfrac = ncread("../model-inputs-and-outputs/forward-runs/combined_realizations.nc", 'grfrac');

umean = mean(u, 3);
vmean = mean(v, 3);
grmean = mean(grfrac, 3);
velocs = sqrt(umean.^2 + vmean.^2);
velocs(velocs ==0) = nan;

clf;
ax(1) = subplot(1,3,1);
h = imagesc(log10(velocs));
c(1) = colorbar;
set(h, 'AlphaData', ~isnan(velocs));
shg
hold on
grmean(isnan(velocs)) = nan;
contour(grmean, [0.5, 0.5], 'k', 'LineWidth',2)
clim([0,4]);
colormap(ax(1), cmocean('haline'))

% get the observations
u_obs = fopen("../driver_files/Inverse_3km_u_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin");
u_obs = fread(u_obs, 'real*8', 'b');
u_obs = reshape(u_obs, [269, 169]);
u_obs = (u_obs(2:end,:) +  u_obs(1:end-1,:))/2; %put on the approx h grid

v_obs = fopen("../driver_files/Inverse_3km_v_velocs_clip_noNan_BedmachineV3_FULL_stripe_fix.bin");
v_obs = fread(v_obs, 'real*8', 'b');
v_obs = reshape(v_obs, [268, 170]);
v_obs = (v_obs(:,2:end) +  v_obs(:,1:end-1))/2; %put on the approx h grid

%work out the gl position
h_obs = fopen("../driver_files/Inverse_3km_thickness_clip_noNan_BedmachineV3_FULL_stripe_fix.bin");
h_obs = fread(h_obs, 'real*8', 'b');
h_obs = reshape(h_obs, [268, 169]);
bed   = fopen("../driver_files/Inverse_3km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin");
bed = fread(bed, 'real*8', 'b');
bed = reshape(bed, [268, 169]);

gr_obs = ones(size(h_obs));
gr_obs(h_obs < -918/1028*bed) = 0;

velocs_obs = sqrt(u_obs.^2 + v_obs.^2);
ax(2) = subplot(1,3,2);
h = imagesc(log10(velocs_obs));
c(2) = colorbar;
set(h, 'AlphaData', ~isnan(velocs));
gr_obs(isnan(velocs)) = nan;
shg
hold on
contour(gr_obs, [0.5, 0.5], 'c', 'LineWidth',1.5)
clim([0,4]);
colormap(ax(2), cmocean('haline'))




% plot the difference
ax(3) = subplot(1,3,3); hold on;
cmap = redblue(101);

if absolute
h = imagesc(velocs - velocs_obs);
% or differences 
else

h = imagesc((velocs - velocs_obs)./velocs_obs * 100);
end
set(h, 'AlphaData', ~isnan(velocs));
contour(ax(3), grmean, [0.5, 0.5], 'k', 'LineWidth',2)
contour(ax(3), gr_obs, [0.5, 0.5], 'c', 'LineWidth',1.5)
set(ax(3), 'YDir', 'reverse');

colormap(gca, cmap);
c(3) = colorbar;
if absolute
clim(1000*[-1,1]);
else
    clim(100*[-1,1]);
end

for i = 1:3
    ax(i).FontSize = 12;
    c(i).FontSize = 12;

end

for i = 1:2
    c(i).Ticks = 0:4;
    c(i).TickLabels = {"10^0","10^1","10^2","10^3", "10^4" };
    c(i).Label.String = 'ice speed (m/yr)';
    c(i).FontSize = 12;
end

if absolute
c(3).Ticks = -1000:500:1000;
c(3).Label.String = 'ice speed difference (m/yr)';
box on;
else

  c(3).Ticks = -100:20:100;
c(3).Label.String = 'ice speed % difference'; %
box on;
end