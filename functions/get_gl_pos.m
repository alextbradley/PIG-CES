function gl_pos = get_gl_pos(yy,grfrac,x0)
%return the grounding line position along line x = x0. yy are the y-co-ordinates, grfrac is the grounding fraction.
grfrac_slice = grfrac(x0,:);
idx = find(grfrac_slice == 1, 1, 'first'); %first fully grounded
float_at_gl_gridcell = grfrac_slice(idx);

if ~isempty(idx) %sometimes when fully deglaciates, we dont get a grounded point. In this case return a nan
gl_pos = yy(idx); %first grid point that is fully floating
else
%gl_pos = yy(end); %take the last point 
gl_pos = -100000; %set upper bound on gl retreat
end

%interpolate???
%gl_pos = yy(idx-1) * float_at_gl_gridcell + yy(idx+1)*(1-float_at_gl_gridcell);
%gl_pos = yy(idx) * float_at_gl_gridcell + yy(idx+1)*(1-float_at_gl_gridcell);


end

