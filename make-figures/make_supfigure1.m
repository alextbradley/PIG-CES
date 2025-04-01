% Make supfigure 1 of the manuscript showing the grounding line position in
% each realization.
%
%
%% First make a plot for each individual realization, because some of them are a bit trash

realizations = ["011", "012", "013", "014", "015", "016", "017", "018", "019",...
    "020", "021", "022", "023", "024", "025", "026", "027", "028",...
    "029", "030", "031", "032", "033", "036", "038", "039", "040"];

realizations = ["011", "012", "014", "015", "016", "017", "019",...
    "020", "022", "023", "026",...
    "029", "030", "032", "036"];

samples      = ["001","002", "003","004", "005","006","007","008", "009", "010"];
samples      = ["001","002", "003","004", "005"];


%trends        = ["continued", "continued_notrend"];
trends       = ["continued", "continued_nobump", "continued_notrend","continued_notrendnobump"]

count = 1;
figure(1); clf; hold on;
colours = [0,0,1;
           1,0,0; 
           0,1,0;
           0.5, 0.5, 0.5]; %colours for the trends

outputs_all = zeros(276, length(samples)*length(realizations), length(trends));
for it = 1:length(trends)
    for ir = 1:length(realizations)
        outputs = zeros(276,length(samples));

        for is = 1:length(samples)

            %load the trajectory: with trend
            traj = load(strcat("../model-inputs-and-outputs/forward-runs/trend_",trends(it),"/realization_", realizations(ir),  "/sample_",samples(is) ,"/output_trajectory.mat"));

            % add to the array
            gl_pos = traj.gl_pos_cts;
            retreat = (gl_pos- gl_pos(1))/1e3; %grounding line retreat in km
            outputs(:, is) = retreat;
            t = traj.t + 1750;


        end

        ens_mean = mean(outputs,2);
        ens_spread = std(outputs,1,2);

        xf = [t; flip(t)];
        yf = [(ens_mean - ens_spread); flip(ens_mean + ens_spread)];

        ax(count) = subplot(6,5,count);
        hold(ax(count), 'on');
        box(ax(count), 'on');
        %plot(ax(count), traj.t, outputs)
        fill(ax(count), xf, yf, colours(it,:), 'FaceAlpha', 0.1, 'linestyle', 'none');
        plot(ax(count), t, ens_mean,'color',  colours(it,:),'linewidth', 2)
        %plot(ax(count), traj.t, ens_mean+ens_spread, 'k--','linewidth', 2)
        %plot(ax(count), traj.t, ens_mean-ens_spread, 'k--','linewidth', 2)



        ax(count).XLim = [1750, 2300];
        ax(count).YLim = [-5,60];
        title (ax(count), strcat('realization', realizations(ir)));


        %add the obs
        if it == length(trends)
            obs = readmatrix("../observations/truth_actual_cts.csv");
            obs = obs(1:2); %take only the gl positions
            obs = (obs - gl_pos(1))/1e3;
            obs_times = readmatrix("../observations/truth_times.csv") + 1750;
            obs_times = obs_times(1:2);
            obs_noise = readmatrix("../observations/noise_actual.csv");
            obs_noise = obs_noise(1:2);
            obs_noise = obs_noise/1e3;

            for i = 1:2
                plot(obs_times(i), obs(i), 'ko', 'markersize', 10, 'markerfacecolor', 'k');
                plot(obs_times(i)*[1,1], obs(i)+[-obs_noise(i), obs_noise(i)], 'k', 'linewidth', 1.5);

            end

        end

        %drawnow
        %pause



        outputs_all(:,(1 + length(samples)*(count-1)):(length(samples)*count), it ) =  outputs;

        count = count + 1;



    end %end loop over realizations
    count = count - length(realizations);
end %end loop over trends

%% make a plot of the average over every run
figure(2); clf; hold on
box on

for it = 1:length(trends)

    ens_mean =  mean(outputs_all(:,:,it), 2);
    ens_spread =  std(outputs_all(:,:,it),1,  2);

    plot(t,ens_mean, 'color', colours(it,:), 'linewidth', 2)


    xf = [t; flip(t)];
    yf = [(ens_mean - ens_spread); flip(ens_mean + ens_spread)];
   fill(xf, yf, colours(it,:), 'FaceAlpha', 0.1, 'linestyle', 'none');


    %add the obs
    if it == length(trends)
        obs = readmatrix("../observations/truth_actual_cts.csv");
        obs = obs(1:2); %take only the gl positions
        obs = (obs - gl_pos(1))/1e3;
        obs_times = readmatrix("../observations/truth_times.csv") + 1750;
        obs_times = obs_times(1:2);
        obs_noise = readmatrix("../observations/noise_actual.csv");
        obs_noise = obs_noise(1:2);
        obs_noise = obs_noise/1e3;

        for i = 1:2
            plot(obs_times(i), obs(i), 'ko', 'markersize', 10, 'markerfacecolor', 'k');
            plot(obs_times(i)*[1,1], obs(i)+[-obs_noise(i), obs_noise(i)], 'k', 'linewidth', 1.5);
        end

    end
end


