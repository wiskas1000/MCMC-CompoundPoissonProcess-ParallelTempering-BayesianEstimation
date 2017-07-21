# W. SEWLAL 1383337 # 2015-12-03
include("ParallelTemperingCPP.jl")
using ParallelTemperingCPP

# @everywhere const id = myid()
# @everywhere const G = Int64(4)
 
# gamma_temperatures = zeros(Float64, G)
# gamma_temperatures = [1.0, 1 + 1e-4, 1 + 2e-4, 1 + 3e-4]
# I CAN JUST SEND THE GAMMA IN A FUNCTION CALL?

# @everywhere remotecall_fetch(getindex, 1, gamma_temperatures, 1)



# gamma = SharedArray(Int32, G, init = S -> S[Base.localindexes(S)] = myid())
# mygamma = SharedArray(Float64, G, pids=workers())
# @everywhere gamma = SharedArray(Int32, G)

@everywhere include("parallel-startup.jl")
@everywhere using PTstartup
@everywhere println(getMyGamma())
# @everywhere info(id)
# @everywhere info(G)
# @everywhere info(workerid)
# @everywhere info(mygamma)

@everywhere @time sleep(3)
# @everywhere info(gamma[2])
# @parallel for i=1:4
# gamma[2] = id
# end


# include("MH-main-functions.jl")
# export MetropolisHastings, MHproposeNiAi!, MHsaveResults, MHloadResults

# include("PT-main-functions.jl")
# export parallelTempering, PTreport, PTsaveResults, PTloadResults
# @everywhere a = SharedArray(Float64,10)
# @parallel for i=1:10
#   a[i] = i
# end
# info(a)
