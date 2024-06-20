function plot_slices_for_iteration(realization, iteration, members)
%plot slices through trajectory for a given realization and iteratio§n

figure(1); clf;

realization = sprintf('%03d', realization);
iteration   = sprintf('%03d', iteration);
member = strings(1,length(members));
for i = 1:length(members)
        member(i) =  sprintf('%03d', members(i));
end

x0 = 64; y0 = 20;
x1 = 64; y1 = 60;

for i = 1:length(members)
ax(i) = subplot(2,5,i);

plot_slice_function(realization, iteration, member(i), [x0,x1], [y0,y1], ax(i))

end
