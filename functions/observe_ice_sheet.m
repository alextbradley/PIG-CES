function observation = observe_ice_sheet(fpath)
%output the observation of ice sheet retreat from an nc file at fpath

t_obs = csvread('../observations/truth_times.csv');
truth_actual = csvread("../observations/truth_actual.csv")
noise_actual = csvread("../observations/noise_actual.csv")
x0 = 64; %index along which to take the observation

%check that the file exists
if ~(exist(fpath))
fprintf(strcat("did not find a file at ", fpath, ", so I am skipping \n" ));
observation = nan;

else

observation = nan(1,3); %two observations of grounding line position and one of ice volume
grfrac = ncread(fpath, 'grfrac');
time   = ncread(fpath,  'TIME');
yy     = ncread(fpath, 'y');
h      = ncread(fpath, 'h');

for i = 1:2
	[~,idx] = min(abs(t_obs(i) - time)); %index of the closest time point
	modelled_gl_position = get_gl_pos(yy,squeeze(grfrac(:,:,idx)), x0);
	observation(i) = (modelled_gl_position - truth_actual(i))/noise_actual(i);
end

% ice volume obs
[~,idx] = min(abs(t_obs(3) - time)); %index of the closest time point
dx = 3e3;
dy = 3e3;
modelled_vaf = sum(sum(squeeze(h(:,:,idx).*grfrac(:,:,idx)))) * dx *dy
observation(3) = (modelled_vaf - truth_actual(3))/noise_actual(3);

end
end
