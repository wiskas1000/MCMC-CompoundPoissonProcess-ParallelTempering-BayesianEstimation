# W. SEWLAL 1383337 # 2015-12-03
module ParallelTemperingCPPMH

using fileCheck
using Distributions
using ParametersCPPMH
using JLD

include("Z.jl")
export getZ

include("MH-main-functions.jl")
export MetropolisHastings, MHproposeNiAi!, MHsaveResults, MHloadResults

include("PT-main-functions.jl")
export parallelTempering, parallelTemperingExchange, parallelExchange, parallelTemperingExchangePU, PTreport, PTsaveResults, PTloadResults

end
