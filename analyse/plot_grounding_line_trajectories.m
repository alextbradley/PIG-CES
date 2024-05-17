% Plot the trajectories for a set of simulations. Colour the trajectories by the iteration

addpath('../functions')
gendata = 0 %flag to generate data, set to 1 to pass thru gendata loop. Leave uncommented so that we can check each time!
x0 = 64;     %where to measure gl

%run spec
realization = ["001","002", "003","004", "005", "006", "007"];
iteration   = ["001"];
member      = ["001","002","003","004","005","006","007","008","009","010"];

iter_colmap = parula(length(iteration)+2);
member_colmap = parula(length(member));

if gendata
ss = struct; 
for ir = 1:length(realization)
for ii = 1:length(iteration)
for im = 1:length(member)


fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization(ir), "/EKI_EKI-", realization(ir),"-", iteration(ii),"-", member(im), "/run/outfile.nc");

if exist(fpath)
t = ncread(fpath, "TIME");
grfrac = ncread(fpath, "grfrac");
yy = ncread(fpath, 'y');
gl_pos = nan(size(t));
for it = 1:length(t)
gl_pos(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)),x0);
end

ss(ir,ii,im).t = t;
ss(ir,ii,im).gl_pos = gl_pos;
else
ss(ir,ii,im).t = nan;
ss(ir,ii,im).gl_pos = nan;
end %end file exists flag
end %end loop over members
end %end loop over iterations
end %end loop over realization
end %end gendata flag


% make the plot
for ir = 1:length(realization)
figure(ir);clf;hold on;
title(['realization ' realization(ir)])
for ii = 1:length(iteration)
for im = 1:length(member)

gl_pos = ss(ir,ii,im).gl_pos;
gl_retreat = gl_pos - gl_pos(1);
plot(ss(ir,ii,im).t + 1750, gl_pos, 'color', member_colmap(im,:), 'linewidth', 1.5)

end %end loop over members
end %end loop over iterations

% add the observations
obs = csvread('../observations/truth.csv');
%obs = obs -gl_pos(1); %obs have gl position, not retreat 
obs_times = csvread('../observations/truth_times.csv');
obs_times = obs_times+1750;
obs_err = csvread('../observations/noise.csv');
for i = 1:length(obs_times)
plot(obs_times(i)*[1,1], obs(i)+obs_err(i)*[-1,1], 'k', 'linewidth', 1.5);
end
plot(obs_times, obs, 'ko', 'markersize', 10, 'markerfacecolor', 'k');

xlabel('time');
ylabel('grounding line position')
ylim([-3,-2.5]*1e5)
box on
end %end loop over realizations
