%output trajectories for the forward runs

addpath('../functions')

realizations = ["011", "012", "013", "014", "015", "016", "017", "018", "019",...
                "020", "021", "022", "023", "024", "025", "026", "027", "028",...
                "029", "030", "031", "032", "033", "036", "038", "039", "040"];

samples      = ["001","002", "003","004", "005", "006", "007", "008", "009", "010"];
trend        = "continued_nobump";


dx = 3e3;
dy = 3e3;
x0 = 64;  %index along which to measure grounding line position

for ir = 1:length(realizations)
for is = 1:length(samples)


    %output file location on HPC
    outfolder = "/data/icesheet_output/aleey/wavi/ARCHER2_EKI/forward-runs/";   

    fpath   = strcat(outfolder,"trend_", trend, "/realization_" ,realizations(ir),"/sample_",samples(is), "/run/outfile.nc");

    if exist(fpath)
        %read the info 
        t = ncread(fpath, "TIME");
        h    = ncread(fpath, "h");
        grfrac = ncread(fpath, "grfrac");

        %loop over time points
        grv = nan(1,length(t));
        gl_pos_discrete = nan(1,length(t));
        gl_pos_cts = nan(1,length(t));
        yy = ncread(fpath, "y");

        for it = 1:length(t)
            grv(it) = sum(sum(squeeze(h(:,:,it).*grfrac(:,:,it)))) * dx *dy;
            gl_pos_discrete(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)),x0);
            gl_pos_cts(it) = get_gl_pos_cts(yy,squeeze(grfrac(:,:,it)),x0);
        end

        outfolder = strcat("/data/hpcdata/users/aleey/projects/AttributionRealWorld/manual-EKI/model-inputs-and-outputs/forward-runs/trend_", trend, "/realization_" ,realizations(ir),"/sample_",samples(is));

        %check whether an observation exists
        if ~isfolder(outfolder)
            disp(strcat("making ", outfolder));
            mkdir(outfolder);
        end
        filename = strcat(outfolder, "/output_trajectory.mat");
	    save(filename, "t", "grv", "gl_pos_discrete", "gl_pos_cts");

    end %end file exists flag
end %end loop over samples
end %end loop over realizations