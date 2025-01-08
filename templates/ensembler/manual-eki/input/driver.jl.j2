using WAVI, JLD2

function driver()


#
#Grid and boundary conditions
#
nx = 268
ny = 169
nσ = 12
x0 = -1792500.0
y0 = -400500.0
dx = 3000.0
dy = 3000.0

h_mask=Array{Float64}(undef,nx,ny);
read!("Inverse_3km_h_mask_clip_BedmachineV3_FULL_stripe_fix.bin",h_mask)
h_mask.=ntoh.(h_mask)

#
# boundary conditions
#
u_iszero=Array{Float64}(undef,nx+1,ny);
read!("Inverse_3km_uiszero_clip_BedmachineV3_FULL_stripe_fix.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("Inverse_3km_viszero_clip_BedmachineV3_FULL_stripe_fix.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("Inverse_3km_sigma_grid_BedmachineV3_FULL_stripe_fix.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

grid = Grid(nx = nx,
            ny = ny,
            nσ = nσ,
            x0 = x0,
            y0 = y0,
            dx = dx,
            dy = dy,
            h_mask = h_mask,
            u_iszero = u_iszero,
            v_iszero = v_iszero,
            σ = sigma_grid)

#
#input files
#
bed=Array{Float64}(undef,nx,ny);
read!("Inverse_3km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",bed)
bed.=ntoh.(bed)

h=Array{Float64}(undef,nx,ny);
read!("steadyThickness_3km_coldForcing_dTpt01.bin",h)
h.=ntoh.(h)

viscosity=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_3km_viscosity3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",viscosity)
viscosity.=ntoh.(viscosity)

temp=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_3km_3Dtemp_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",temp)
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_3km_damage3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",damage)
damage.=ntoh.(damage)

dhdt=Array{Float64}(undef,nx,ny);
read!("Inverse_3km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdt)
dhdt.=ntoh.(dhdt)

initial_conditions = InitialConditions(initial_thickness = h,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage)

#
# physical paramerter
# 
sec_per_year = 365.25*24*60^2
glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9 #standard value used in WAVI
glen_a_ref_prefactor = {{ run.glen_a_ref_prefactor }}
glen_a_ref_prefactor = max(glen_a_ref_prefactor, 1e-1) #avoid negative value
glen_a_ref = glen_a_ref_prefactor .* glen_a_ref

weertman_c=Array{Float64}(undef,nx,ny);
read!("Inverse_3km_WeertmanC_clip_adjusted_noNan_BedmachineV3_FULL_stripe_fix.bin",weertman_c)
weertman_c.=ntoh.(weertman_c)
weertman_c_prefactor = {{ run.weertman_c_prefactor }}
weertman_c_prefactor = max(weertman_c_prefactor, 1e-1) #avoid negative weertmanc
weertman_c = weertman_c .* weertman_c_prefactor

# adjust the weertman c in floating areas
ungrounded_weertmanC_prefactor = {{ run.ungrounded_weertmanC_prefactor }}
ungrounded_weertmanC_prefactor = max(ungrounded_weertmanC_prefactor, 1e-1) #avoid negative weertmanc
weertman_c[weertman_c .== 10000] .= 10000.0 * ungrounded_weertmanC_prefactor

accumulation_rate=Array{Float64}(undef,nx,ny);
read!("Inverse_3km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

params = Params(accumulation_rate = accumulation_rate,
				  weertman_c = weertman_c, 
				  glen_a_ref = glen_a_ref)

#
#solver parameters
#
maxiter_picard = 6
#parallel_spec = SharedMemorySpec(ngridsx=2,ngridsy=1,overlap=1,damping=0.0,niterations=1)
parallel_spec = BasicParallelSpec()
tol_picard = 1.0e-4;
solver_params = SolverParams(maxiter_picard = maxiter_picard, tol_picard = tol_picard)

#
# melt rate model
# 
bump_amplitude      = {{ run.bump_amplitude }}
bump_duration       = {{ run.bump_duration }}
bump_duration       = max(bump_duration, .1) #must be positive
melt_rate_prefactor_exponent = {{ run.melt_rate_prefactor_exponent }}
per_century_trend   = {{ run.per_century_trend }}
random_seed         = {{ run.random_seed }}

end_time = 300. #end in 2050, start is 1750
bump_time   = 195. #1945
trend_onset = 210. #1960
#pc_max = -399.0
#pc_min = -401.0 #warm conditions
#pc_max = -599.0
#pc_min = -601.0 #cold conditions
pc_max = -400.0
pc_min = -600.0 #variable
pw     = 400.0
rf_threshold =2.0 

idealized_anthro_melt_rate = IdealizedAnthroMeltRate(bump_amplitude = bump_amplitude,
bump_width = bump_duration/2,
bump_time = bump_time,
per_century_trend = per_century_trend,
trend_onset = trend_onset,
pc_max = pc_max,
pc_min = pc_min,
M = 5.0*(2.0^(melt_rate_prefactor_exponent)),
random_seed = random_seed,
rf_threshold = rf_threshold,
pw = pw, 
smooth_timescale = 1.0)


model = Model(grid = grid,
            bed_elevation = bed,
            params = params,
            solver_params = solver_params,
            initial_conditions= initial_conditions, 
            parallel_spec = parallel_spec, 
	    melt_rate = idealized_anthro_melt_rate)

#
#timestepping parameters (end time prescribed above)
#
niter0 = 0
dt = 0.05
chkpt_freq = 1000.0
pchkpt_freq = 100.0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                         dt = dt, 
                                         end_time = end_time, 
                                         chkpt_freq = chkpt_freq, 
                                         pchkpt_freq = pchkpt_freq)

#
#output parameters
#
outputs = (h = model.fields.gh.h,
           u = model.fields.gh.u,
           v = model.fields.gh.v,
           b = model.fields.gh.b,
           s = model.fields.gh.s,
           a = model.fields.gh.accumulation,
           grfrac = model.fields.gh.grounded_fraction,
           m = model.fields.gh.basal_melt)

output_freq = 1.0
output_params = OutputParams(outputs = outputs, 
                            output_freq = output_freq,
                            output_format = "mat",
                            zip_format = "nc",
                            output_start = true)

#
# assemble the simulation
#
simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)
                
#
#perform the simulation
#
run_simulation!(simulation)

return simulation
end

