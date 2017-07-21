# W. SEWLAL 1383337 # 2015-12-03
module ParametersCPPMH
using Distributions

export MHiterationParameters, MHiterationDistributions

type MHiterationParameters
    labda::Float64
    psi::Array{Float64, 1}
    mu::Array{Float64, 1}
    tau::Float64
    ksi::Array{Float64, 1}
end

type MHiterationDistributions
    dist_psi::Array{Distribution, 1}
    dist_mu::Distribution
    dist_tau::Distribution
end

end
# Array{Distributions.Distribution{F<:Distributions.VariateForm,S<:Distributions.ValueSupport},1} #psi
# Distributions.MvNormal{PDMats.PDMat,Array{Float64,1}} #mu
# Distributions.Gamma #tau