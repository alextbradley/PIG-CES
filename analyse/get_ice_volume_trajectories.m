function ss = get_ice_volume_trajectories(realization, iteration, members,mode)
%return the data for a single iteration

addpath('../functions')
ss = struct;
dx = 3e3;
dy = 3e3;
for im = 1:length(members)
fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration,"-", members(im), "/run/outfile.nc");
if exist(fpath)
        t = ncread(fpath, "TIME");
        grfrac = ncread(fpath, "grfrac");
        h = ncread(fpath, "h");
        b = ncread(fpath, "b");
        yy = ncread(fpath, 'y');
        ice_volume = nan(size(t));

        for it = 1:length(t)
                if strcmp(mode, "grounded_volume")
                        ice_volume(it) = sum(sum(h(:,:,it).*grfrac(:,:,it)))*dx*dy;
		elseif strcmp(mode, "volume")
                        ice_volume(it) = sum(sum(h(:,:,it)))*dx*dy;
                elseif strcmp(mode, "vaf")
                        haf = h - (1028.0/918.0)*(- b); %height above floatation
                        ice_volume(it) = sum(sum(haf(haf > 0)))*dx*dy; %vaf
                else
                        error("mode must be of type volume or vaf")
                end

        end %end loop over time points

        ss(im).t = t;
        ss(im).ice_volume = ice_volume;

else
        ss(im).t = nan;
        ss(im).ice_volume = nan;

end
end %end loop over members
