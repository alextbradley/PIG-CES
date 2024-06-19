count = 1;
for realization = 26:30
figure(count) ; clf; 
suptitle(strcat("realization ", num2str(realization)))
make_plots(realization, 3, 1:20)
fig = gcf;
fig.Position(3:4) = [1050,500];
drawnow
count = count + 1;
end
