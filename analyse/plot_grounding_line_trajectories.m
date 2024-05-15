% Plot the trajectories for a set of simulations. Colour the trajectories by the iteration

addpath('../functions')
gendata = 1; %flag to generate data, set to 1 to pass thru gendata loop
x0 = 64;     %where to measure gl

%run spec
realization = "002";
iteration   = "001";
member      = ["001","002","003","004","005","006","007","008","009","010"];

iter_colmap = parula(length(iteration)+2);

if gendata
ss = struct; 
for ii = 1:length(iteration)
for im = 1:length(member)


fpath = strcat("/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization",realization, "/EKI_EKI-", realization,"-", iteration(ii),"-", member(im), "/run/outfile.nc");

if exist(fpath)
t = ncread(fpath, "TIME");
grfrac = ncread(fpath, "grfrac");
yy = ncread(fpath, 'y');
gl_pos = nan(size(t));
for it = 1:length(t)
gl_pos(it) = get_gl_pos(yy,squeeze(grfrac(:,:,it)),x0);
end

ss(ii,im).t = t;
ss(ii,im).gl_pos = gl_pos;
else
ss(ii,im).t = nan;
ss(ii,im).gl_pos = nan;
end %end file exists flag
end %end loop over members
end %end loop over iterations
end %end gendata flag


% make the plot
figure(1);clf;hold on;
for ii = 1:length(iteration)
for im = 1:length(member)

gl_pos = ss(ii,im).gl_pos;
gl_retreat = gl_pos - gl_pos(1);
plot(ss(ii,im).t + 1750, gl_pos, 'color', iter_colmap(ii,:), 'linewidth', 1.5)

end %end loop over members
end %end loop over iterations

% add the observations
obs = csvread('../observations/truth.csv');
%obs = obs -gl_pos(1); %obs have gl position, not retreat 
obs_times = csvread('../observations/truth_times.csv');
obs_times = obs_times;
obs_err = csvread('../observations/truth_times.csv');
for i = 1:length(obs_times)
plot(obs_times(i)*[1,1], obs(i)+obs_err(i)*[-1,1], 'k', 'linewidth', 1.5);
end
plot(obs_times, obs, 'ko', 'markersize', 10, 'markerfacecolor', 'k');

xlabel('time');
ylabel('grounding line position')
box on
