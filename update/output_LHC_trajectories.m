% output the trajctories for the latin hypercube case

addpath('../functions')

realization = ["21"];

n_runs = 0:99; %indices of runs
dx = 3e3;
dy = 3e3;
x0 = 64;  %index along which to measure grounding line position

for ir = 1:length(realization)
for im = 1:length(n_runs)

    %output file
    fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization0",realization(ir), "_lhc/manual-eki-latinhypercube-realization", realization(ir), "-", num2str(n_runs(im)), "/run/outfile.nc");

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

        output = observe_ice_sheet(fpath)
        output_cts = observe_ice_sheet_cts(fpath)

        outfolder = strcat("/data/hpcdata/users/aleey/projects/AttributionRealWorld/manual-EKI/model-inputs-and-outputs/realization0", realization(ir), "_lhc/member", num2str(n_runs(im)));
        if ~exist(outfolder)
        mkdir(outfolder)
        end

        %check whether an observation exists
        filename = strcat(outfolder, "/output_trajectory.mat");
	    save(filename, "t", "grv", "gl_pos_discrete", "gl_pos_cts");

        filename = strcat(outfolder, "/outputs.csv");
        writematrix(output, filename);

        filename = strcat(outfolder, "/outputs_cts.csv");
        writematrix(output_cts, filename);


    end %end file exists flag
end %end loop over runs
end %end loop over members 
