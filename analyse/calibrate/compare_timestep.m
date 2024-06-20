% plot the grounding line position and ice volume for the simulations with different temporal resolutions
addpath('../functions')
labels = ["dT = 0.05"; "dT = 0.025"; "dT = 0.0125"];

prefix ="/data/icesheet_output/aleey/wavi/ARCHER2_EKI/realization001";
folders = ["EKI_EKI-001-001-004", "EKI_EKI-001-001-008";
        "EKI_EKI-001-001-004_pt5timestep", "EKI_EKI-001-001-008_pt5timestep";
        "EKI_EKI-001-001-004_pt25timestep", "EKI_EKI-001-001-008_pt25timestep"];

x0 = 64;
dx = 3e3;
dy = 3e3;
sz = size(folders);

gendata = 0;
if gendata
ss = struct;
for j = 1:sz(2)
for i = 1:sz(1)
    jdir = dir(fullfile(prefix, folders(i,j), "run", "*.mat"));

    time = nan(1,length(jdir));
    gl_pos = nan(1,length(jdir));
    ice_vol = nan(1,length(jdir));
    for k = 1:length(jdir)
        fpath = fullfile(prefix, folders(i,j), "run", jdir(k).name);       
        data = load(fpath);

        grfrac = data.grfrac;
        yy = data.y;
        yy = yy(1,:);
        gl_pos(k) = get_gl_pos(yy,grfrac,x0);
        time(k) = data.t;
        ice_vol(k) = sum(sum(data.h))*dx*dy;

    end
    ss(i,j).time = time;
    ss(i,j).gl_pos = gl_pos;
    ss(i,j).ice_vol = ice_vol;
end
end
end

%% make the plot
figure(1); clf;  

cols = [0,0,1; 1,0,0]; %colour by simulation
handlevis = ["on", "off"];
styles = ["-", "--", ":"]; %timestep

% plot of grounding line position
subplot(1,2,1); hold on; box on;
for i = 1:sz(1)
for j = 1:sz(2)
    plot(ss(i,j).time, ss(i,j).gl_pos, 'color', cols(j,:), 'linestyle', styles(i), 'HandleVisibility', handlevis(j));
end
end
xlabel('time');
ylabel('grounding line position');
legend(labels, 'location', 'northwest');

% plot of ice volume
subplot(1,2,2); hold on; box on

for i = 1:sz(1)
for j = 1:sz(2)
    plot(ss(i,j).time, ss(i,j).ice_vol, 'color', cols(j,:), 'linestyle', styles(i))
end
end

xlabel('time');
ylabel('ice volume (m^3)');


%% make plot of slices

figure(2); clf;
x0 = 64; y0 = 20;
x1 = 64; y1 = 60;

count = 1;
for j = 1:sz(2)
for i = 1:sz(1)
    ax(i,j) = subplot(2,3,count);
    plot_slice_function_mat(fullfile(prefix, folders(i,j)),  [x0,x1], [y0,y1], ax(i,j))
    count = count + 1;
end
end
