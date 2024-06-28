function gl_pos = get_gl_pos_cts(yy,grfrac,x0)
%return the grounding line position along line x = x0. yy are the y-co-ordinates, grfrac is the grounding fraction.
% This is a continous version where we take interpolation across the last fully floating cell and the first partially grounded cell. 
grfrac_slice = grfrac(x0,:);
idx = find(grfrac_slice > 0, 1, 'first'); %first partially floating point

if isempty(idx); %sometimes when fully deglaciates, we dont get a grounded point. In this case return an upper bound
    gl_pos = -100000; %set upper bound on gl retreat

else
    float_at_gl_gridcell = grfrac_slice(idx); %first partially grounded fraction
    gl_pos = yy(idx-1)*(1-float_at_gl_gridcell) + yy(idx)*float_at_gl_gridcell; 
end


end

