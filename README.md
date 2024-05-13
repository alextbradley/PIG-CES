- Priors.toml stores the information on the priors, which are the same for each realization of forcing.

- observations/truth.jld2 contains the observations, used to update the EKI:
     Γ = 1.0 * 3*10.0^3 * I #noisy observation of grounding line position, where I is a 2 x 2identity matrix
     y is a 2 x 1 noisy observation of grounding line position

- results are stored in /realizationXXX/iterationYYY/memberZZZ, with the eki stored in /realizationXXX/iterationYYY/eki.jld2