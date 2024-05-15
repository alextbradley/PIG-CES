function gl_pos = get_gl_pos(yy,grfrac,x0)
%return the grounding line position along line x = x0. yy are the y-co-ordinates, grfrac is the grounding fraction.
grfrac_slice = grfrac(x0,:);
idx = find(grfrac_slice == 1, 1, 'first'); %first fully grounded
idx = find(grfrac_slice >0, 1, 'first'); %first not fully floating
float_at_gl_gridcell = grfrac_slice(idx);
gl_pos = yy(idx-1) * float_at_gl_gridcell + yy(idx+1)*(1-float_at_gl_gridcell);%interpolate
%gl_pos = yy(idx) * float_at_gl_gridcell + yy(idx+1)*(1-float_at_gl_gridcell);%interpolate
%gl_pos = yy(idx); %first grid point that is not fully floating

end
