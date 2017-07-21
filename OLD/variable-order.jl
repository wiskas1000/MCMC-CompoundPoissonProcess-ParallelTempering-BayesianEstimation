# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # DESCRIPTION AND CALL ORDER OF VARIABLES # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
* `singlePlot`: Boolean variable that decides if the plots are combined and exported into a single file.

* `N` : Int64 -- Number of segments
* `J` : Int64 -- Number of jump types
* `I` : Int64 -- Number of nonzero segments
* `G` : Int64 -- Number of parallel runs
* `M` : Int64 -- Number of iterations in the algorithm
#
# labdaDelta
* `Delta` : Float64 -- The time between observations, fixed for all observations
* `real_labda` : Float64 --
* `T` :

* `real_psi` : Array{Float64,1} -- `J`-element array
* `real_mu` : Array{Float64,1} -- `J`-element array
* `real_tau` : Float64 --

* `running_labda` : `M`-element Array{Float64,1} containing the values of labda
* `running_psi` : `J`x`M` Array{Float64,`J`} containing the values of psi
* `running_mu` : `J`x`M` Array{Float64,`J`} containing the values of mu
* `running_tau` : `M`-element Array{Float64,1} containing the values of tau
* `acceptance` : `M`-element Array{Float64,1} containing the values of acceptance
* `figurePath` : String, path to the directory were the plots can be stored
#
#     labda::Float64
#     psi::Array{Float64, 1}
#     mu::Array{Float64, 1}
#     tau::Float64
#     ksi::Array{Float64, 1}
#
* `z` : Array{Float64,1} -- `I`-element array differences in a sample of a Compound Poisson Process
* `I` : Int64 -- Number of nonzero segments
# Z
* `i` : Int64 -- Segment. Warning: `i` must be between `1` and `I`. There is no check implemented for this.
#
* `alpha_0` : Float64
* `beta_0` : Float64
* `alpha_1` : Float64
* `beta_1` : Float64
* `ksi` : Array{Float64,1} -- `J`-element array
* `kappa` : Float64
#

* `intensity` : Array{Float64,1} -- `J`-element array giving the intensity for each element
* `psiDelta` : Array{Float64,1} -- `J`-element array with the jump-intensity

current_n
current_a

in_mu
in_tau
in_ksi
