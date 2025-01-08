using WAVI
function driver()

#
#Grid and boundary conditions
#
nx = 640
ny = 80
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 1000.0
dy = 1000.0
h_mask=trues(nx,ny)
u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true
grid = Grid(nx = nx,
            ny = ny,
            nσ = nσ,
            x0 = x0,
            y0 = y0,
            dx = dx,
            dy = dy,
            h_mask = h_mask,
            u_iszero = u_iszero,
            v_iszero = v_iszero)

#
#Bed
#
bed = WAVI.mismip_plus_bed #function definition

#
#solver parameters
#
maxiter_picard = 1
solver_params = SolverParams(maxiter_picard = maxiter_picard)

#
#Physical parameters
#
default_thickness = {{ run.thickness }}
accumulation_rate = {{ run.accumulation }}
params = Params(default_thickness = default_thickness,
                accumulation_rate = accumulation_rate)

#
#make the model
#
model = Model(grid = grid,
              bed_elevation = bed,
              params = params,
              solver_params = solver_params)

#
#timestepping parameters
niter0 = 0 #CHANGE ME FOR A PICKUP (!! must be niter0 !!)
step_thickness = false #default = true!
dt = 0.1
end_time = 1000.
chkpt_freq = 1.
pchkpt_freq = 1.
timestepping_params = TimesteppingParams(
                                        niter0 = niter0,
                                        dt = dt,
                                        end_time = end_time,
                                        chkpt_freq = chkpt_freq,
                                        pchkpt_freq = pchkpt_freq,)

#
#output parameters
#
outputs = (h = model.fields.gh.h,
           u = model.fields.gh.u,
           v = model.fields.gh.v)
output_freq = 50.
output_params = OutputParams(outputs = outputs,
                            output_freq = output_freq,
                            output_format = "mat",
                            dump_vel = true,
          zip_format = "nc")

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
