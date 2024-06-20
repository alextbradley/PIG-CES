function ss = get_ice_volume_trajectories(realization, iteration, members)
%return the data for a single iteration

addpath('../functions')
ss = struct;
dx = 3e3;
dy = 3e3;
for im = 1:length(members)
im
fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration,"-", members(im), "/run/outfile.nc");
if exist(fpath)
        t = ncread(fpath, "TIME");
        grfrac = ncread(fpath, "grfrac");
        h = ncread(fpath, "h");
        b = ncread(fpath, "b");
        yy = ncread(fpath, 'y');
        ice_volume = nan(size(t));

        for it = 1:length(t)
                if mod(it,1) == 0
                        ice_volume(it) = sum(sum(h(:,:,it).*grfrac(:,:,it)))*dx*dy;
                else
                        ice_volume(it) = nan;
                end
        end %end loop over time points
        ice_volume = ice_volume(~isnan(ice_volume)); %remove any nans
        t = t(~isnan(ice_volume));

        ss(im).t = t;
        ss(im).ice_volume = ice_volume;

else
        ss(im).t = nan;
        ss(im).ice_volume = nan;

end
end %end loop over members
