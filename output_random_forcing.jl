using WAVI, MAT, Printf

t = 0:0.1:300
random_seed = 11:20
smooth_timescale = 1.0
rf_threshold = 2.0
pc_max = -400
pc_min = -600

r = 0.8
for ir = 1:length(random_seed)
    rfa = WAVI.generate_random_forcing_anomaly.(t, r, random_seed[ir], smooth_timescale)
    pc = WAVI.get_random_pc_component.(t, r, random_seed[ir], rf_threshold, pc_max, pc_min, smooth_timescale)

    padded_realization = @sprintf("%03d", random_seed[ir]);

    fname =  "./model-inputs-and-outputs/realization" * padded_realization * "/realization.mat"
   # file = matopen(fname, "w")
   dd = Dict(
        "time" => collect(t),
        "random_forcing_anomaly" => rfa,
        "pycnocline_center" => pc)

    matwrite(fname,dd)
#        close(file)



end
