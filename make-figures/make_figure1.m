% make figure 1 of the manuscript showing (a) the Pine Island domain with
% grounding line positions in present day and in equilibrium on top of the
% bed topography and (b) slices along the flowline
%
% 28/01/2025, ATB. alex.bradley@kcl.ac.uk. MIT license.
%
%% Preliminaries
addpath('../functions');
fig1 = figure(1);clf; 
for i = 1:3
    ax(i) = subplot(3,1,i);
    if i ~= 3
        axis equal
        ax(i).XLabel.String = 'xps (m)';
        ax(i).YLabel.String = 'yps (m)';
        
    end
    hold(ax(i), 'on');
    box(ax(i), 'on');
    ax(i).FontSize = 14;
end
fig1.Position(3:4) = [420, 1200];
klein_blue = [0, 47, 167 ]/255;
klein_blue = [0, 255, 255 ]/255;
klein_red =  [167, 0, 47]/255;

%% Set up the grid
x0 = -1792500.0;
y0 = -400500.0;
dx = 3000.0;
dy = 3000.0;
nx = 268;
ny = 169;
xx = x0 + (dx/2 : dx : nx*dx);
yy = y0 + (dy/2 : dy : ny*dy);


%% Get files and stuff

bathy_path = "../driver_files/Inverse_3km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin";
bathy_fid  = fopen(bathy_path);
bathy      = fread(bathy_fid, 'real*8', 'b');
bathy      = reshape(bathy,[268, 169] );


% add grounding line based on ice thickness
icethickness_2015_path = "../driver_files/Inverse_3km_thickness_clip_noNan_BedmachineV3_FULL_stripe_fix.bin";bathy_fid  = fopen(bathy_path);
icethickness_2015_fid  = fopen(icethickness_2015_path);
icethickness2015       = fread(icethickness_2015_fid, 'real*8', 'b');
icethickness2015       = reshape(icethickness2015,[268, 169] );

% work out the floating fraction 
height_above_floatation_2015 = icethickness2015 - (1028.0/918.0)*( - bathy);
grfrac_2015 = height_above_floatation_2015;
grfrac_2015(grfrac_2015 < 0) = 0;
grfrac_2015(grfrac_2015 > 0) = 1;



% repeat for steady grounding line 
icethickness_steady_path = "../driver_files/steadyThickness_3km_coldForcing_dTpt01.bin";bathy_fid  = fopen(bathy_path);
icethickness_steady_fid  = fopen(icethickness_steady_path);
icethickness_steady       = fread(icethickness_steady_fid, 'real*8', 'b');
icethickness_steady       = reshape(icethickness_steady,[268, 169] );

% work out the floating fraction 
height_above_floatation_steady = icethickness_steady - (1028.0/918.0)*( - bathy);
grfrac_steady = height_above_floatation_steady;
grfrac_steady(grfrac_steady < 0) = 0;
grfrac_steady(grfrac_steady > 0) = 1;

% get the ice mask
icemask_path = "../driver_files/Inverse_3km_h_mask_clip_BedmachineV3_FULL_stripe_fix.bin";bathy_fid  = fopen(bathy_path);
icemask_fid  = fopen(icemask_path);
icemask       = fread(icemask_fid, 'real*8', 'b');
icemask       = reshape(icemask,[268, 169] );


%% make first panel: whole domain
imagesc(ax(1),xx, yy, bathy');
contour(ax(1),xx,yy, grfrac_2015', [0.5, 0.5], 'Color',klein_blue, 'LineWidth',1.5 )
contour(ax(1),xx,yy, grfrac_steady', [0.5, 0.5], 'Color',klein_red, 'LineWidth',1.5 )
contour(ax(1),xx,yy, icemask', [0.5, 0.5], 'm--', 'LineWidth',1.5 )


c = colorbar(ax(1));
c.Position(1) = 0.83;
c.Position(3) = 0.02;
c.FontSize =14;
c.Label.String = 'elevation (m)';

clim(ax(1), [-1200, 1200]);
colormap(ax(1),cmocean('topo'))
ax(1).Position(3) = 0.7;
ax(1).Position(1) = 0.1;


%% make first panel: zoom in
contourf(ax(2),xx, yy,bathy',50, 'LineStyle','none', 'HandleVisibility','off' );
contour(ax(2),xx,yy, grfrac_2015', [0.5, 0.5], 'Color',klein_blue, 'LineWidth',2 )
contour(ax(2),xx,yy, grfrac_steady', [0.5, 0.5], 'Color',klein_red, 'LineWidth',2)
contour(ax(2),xx,yy, icemask', [0.5, 0.5], 'm--', 'LineWidth',1.5, 'HandleVisibility','off' );


c2 = colorbar(ax(2));
c2.Position(1) = 0.82;
c2.Position(3) = 0.02;
c2.FontSize =14;
c2.Label.String = 'elevation (m)';
clim(ax(2), [-1200, 0]);
cc = cmocean('topo', 101);
colormap(ax(2),cc(1:50, :));
ax(2).Position(3) = 0.7;
ax(2).Position(1) = 0.1;
ax(2).YTick = [-3,-2]*1e5;
ax(2).XLim = 1e6 *[-1.78, -1.5];
ax(2).YLim = 1e5 *[-3.8, -1.8];
l = legend(ax(2), {'2015 grounding line', 'pre-industrial grounding line'});

% add a box in panel a
xc = ax(2).XLim;
yc = ax(2).YLim;


plot(ax(1), xc, [yc(1), yc(1)],'w--', 'LineWidth',2, 'HandleVisibility','off');
plot(ax(1), xc, [yc(2), yc(2)],'w--', 'LineWidth',2, 'HandleVisibility','off');
plot(ax(1), [xc(1), xc(1)], yc,'w--', 'LineWidth',2, 'HandleVisibility','off');
plot(ax(1), [xc(2), xc(2)], yc,'w--', 'LineWidth',2, 'HandleVisibility','off');


%% get stuff along the line
xt = [64,64]; %index of the centre line
yt = [20,55];
dl = (0:dy:diff(yt)*dy)/1000;
dl = max(dl) - dl;

plot(ax(2), xx(xt), yy(yt), 'w', 'LineWidth',4, 'HandleVisibility','off')

bathy_line         = bathy(xt(1), yt(1):yt(2));
h2015_line         = icethickness2015(xt(1), yt(1):yt(2));
grfrac_2015_line   = grfrac_2015(xt(1), yt(1):yt(2));
hsteady_line       = icethickness_steady(xt(1), yt(1):yt(2));
grfrac_steady_line = grfrac_steady(xt(1), yt(1):yt(2));
mask_line          = icemask(xt(1), yt(1):yt(2));

% work out the surfaces
s2015_line = bathy_line + h2015_line;
s2015_line(grfrac_2015_line == 0) = (1-918/1028)*h2015_line(grfrac_2015_line == 0);
s2015_line(mask_line==0) = nan;
b2015_line = s2015_line - h2015_line;


ssteady_line = bathy_line + hsteady_line;
ssteady_line(grfrac_steady_line == 0) = (1-918/1028)*hsteady_line(grfrac_steady_line == 0);
ssteady_line(mask_line==0) = nan;
bsteady_line = ssteady_line - hsteady_line;
%% make panel 3
yl3 = [-1300, 800];

% add ocean fill background
xfb = [dl, flip(dl)];
yfb = [0*ones(size(dl)), flip(bathy_line)];
fill(ax(3), xfb, yfb, [133, 183, 249]/255, 'LineStyle','none')


%fill the bed
xfb = [dl, flip(dl)];
yfb = [-1300*ones(size(dl)), flip(bathy_line)];
fill(ax(3), xfb, yfb, [164, 126, 64]/255, 'LineStyle','none')
plot(ax(3), dl, bathy_line, 'k', 'LineWidth',2);


% fill present day geometry
xfb = [dl, flip(dl)];
yfb = [b2015_line, flip(s2015_line)];
idx = ~isnan(yfb);

fill(ax(3), xfb(idx), yfb(idx), [190,190,190]/255, 'LineStyle','none')

plot(ax(3), dl, s2015_line, 'Color', klein_blue, 'LineWidth',2);
plot(ax(3), dl, b2015_line, 'Color', klein_blue, 'LineWidth',2);

plot(ax(3), dl, ssteady_line, 'Color', klein_red, 'LineWidth',2);
plot(ax(3), dl, bsteady_line, 'Color', klein_red, 'LineWidth',2);

ax(3).YLim = yl3;
ax(3).XLim = [0,max(dl)];
ax(3).XLabel.String = 'distance along transect (km)';
ax(3).YLabel.String = 'elevation (m)';

% exportgraphics(gcf, "figures/raw/figure1.pdf", 'ContentType','vector');