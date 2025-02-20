function ss = get_grounding_line_trajectories(realization, iteration, members, x0)
%return the data for a single iteration

addpath('../../functions')
ss = struct;
for im = 1:length(members)
fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration,"-", members(im), "/run/outfile.nc");
if exist(fpath)
        t = ncread(fpath, "TIME");
        grfrac = ncread(fpath, "grfrac");
        yy = ncread(fpath, 'y');
        gl_pos = nan(size(t));
        for it = 1:length(t)
                gl_pos(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)),x0);
        end %end loop over time points

        ss(im).t = t;
        ss(im).gl_pos = gl_pos;

else
        ss(im).t = nan;
        ss(im).gl_pos = nan;

end
end %end loop over members
