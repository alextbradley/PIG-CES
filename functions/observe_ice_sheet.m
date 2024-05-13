function observation = observe_ice_sheet(fpath)
%output the observation of ice sheet retreat from an nc file at fpath

t_obs = [180,270]; %time in the simulations of the observations
x0 = 64; %index along which to take the observation

%check that the file exists
if ~(exist(fpath))
fprintf(strcat("did not find a file at ", fpath, ", so I am skipping \n" ));
observation = nan;

else

observation = nan(1,length(t_obs)); %placeholder
grfrac = ncread(fpath, 'grfrac');
time   = ncread(fpath,  'TIME');
yy     = ncread(fpath, 'y');

for i = 1:length(t_obs)
	[~,idx] = min(abs(t_obs(i) - time)); %index of the closest time point
	observation(i) = get_gl_pos(yy,squeeze(grfrac(:,:,idx)), x0);

end


end
end
