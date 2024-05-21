% Plot a slice through the ice

realization = "003";
iteration   = "001";
member      = ["001","002","003","004","005","006","007","008","009","010"];
member = "001"

% specify where to take the slice
x0 = 64; y0 = 20;
x1 = 64; y1 = 60;

for im = 1:length(member)

figure(im); 
fig = plot_slice_function(realization, iteration, member(im), [x0,x1], [y0,y1])

end
