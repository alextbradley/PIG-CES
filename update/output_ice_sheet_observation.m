% output the observation of the ice sheet
addpath('../functions')

realization = ["001"];
iteration   = ["001"];
member      = ["001"];

for ir = 1:length(realization)
for ii = 1:length(iteration)
for im = 1:length(member)

fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization(ir), "/EKI_EKI-", realization(ir),"-", iteration(ii),"-", member(im), "/run/outfile.nc");

%generate the observation
output = observe_ice_sheet(fpath);

%save in the appropriate place as a csv
if ~isnan(output) %if we don't find the nc file, we output a nan
outfolder = strcat("/data/hpcdata/users/aleey/projects/AttributionRealWorld/manual-EKI/realization", realization(ir), "/iteration" , iteration(ii), "/member", member(im));

if ~isdir(outfolder)
	mkdir(outfolder)
end

%check whether an observation exists
filename = strcat(outfolder, "/outputs.csv");
if exist(filename)
	fprintf(strcat("found an output at ", filename ,", so skipping \n"))
else
	%save the observation in csv format
	writematrix(output, filename);
end

end %end the nc file exists flag

end
end
end

