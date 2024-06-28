%output ice sheet trajectories

addpath('../functions')

realization = ["026", "027", "028", "029", "030"];
iteration   = ["001", "002", "003", "004", "005"];
member      = ["001","002","003","004","005","006","007","008", "009", "010","011","012","013","014","015","016","017","018","019","020"];

dx = 3e3;
dy = 3e3;
x0 = 64;  %index along which to measure grounding line position

for ir = 1:length(realization)
for ii = 1:length(iteration)
for im = 1:length(member)

    %output file
    fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization(ir), "/EKI_EKI-", realization(ir),"-", iteration(ii),"-", member(im), "/run/outfile.nc");

    if exist(fpath)
        %read the info 
        t = ncread(fpath, "TIME");
        h    = ncread(fpath, "h");
        grfrac = ncread(fpath, "grfrac");

        %loop over time points
        vaf = nan(1,length(t));
        gl_pos_discrete = nan(1,length(t));
        gl_pos_cts = nan(1,length(t));
        yy = ncread(fpath, "y");

        for it = 1:length(t)
            vaf(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
            gl_pos_discrete(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)),x0);
            gl_pos_cts(it) = get_gl_pos_cts(yy,squeeze(grfrac(:,:,it)),x0);
        end

        outfolder = strcat("/data/hpcdata/users/aleey/projects/AttributionRealWorld/manual-EKI/model-inputs-and-outputs/realization", realization(ir), "/iteration" , iteration(ii), "/member", member(im));

        %check whether an observation exists
        filename = strcat(outfolder, "/output_trajectory.mat");
	    save(filename, "time", "vaf", "gl_pos_discrete", "gl_pos_cts");

    end %end file exists flag
end %end loop over members
end %end loop over iterations
end %end loop over realizations