 # Pkg.add("ParallelDataTransfer")
# @everywhere using ParallelTemperingCPP
@everywhere using ParallelDataTransfer
@everywhere using JLD2
# W. SEWLAL 1383337 # 2015-12-03
include("ParallelTemperingCPP.jl")

const gamma_temperatures = Array{Float64, 1}(4)
const G = length(gamma_temperatures)
gamma_temperatures[1:G] = [1.0, 1.0 + 1e-4, 1.0 + 2e-4, 1.0 + 3e-4]
@broadcast G = getfrom(1, :G)
@broadcast const gamma_temperatures = getfrom(1, :gamma_temperatures)

# @everywhere info(G)
# @everywhere info(gamma_temperatures)

pathTMPStream = "/home/wikash/Documents/Simulations/TMP/"

@everywhere include("parallel-startup.jl")
@everywhere using PTstartup

@broadcast setIndexMyGamma(idWorker)
@broadcast setMyPathTMPStream(getfrom(1, :pathTMPStream))
@broadcast println(gamma_temperatures[getIndexMyGamma()])
# @broadcast println(getIndexMyGamma())
# @broadcast println(getMyPathTMPStream())


# @everywhere @time sleep(3)
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


# mygamma::Array{Float64, 1}(4)

# gamma_temperatures[1] = 1.0
# gamma_temperatures[2] = 1 + 1e-4
# gamma_temperatures[3] = 1 + 2e-4
# gamma_temperatures[4] = 1 + 3e-4
# info(gamma_temperatures)

# gamma = SharedArray(Int32, G, init = S -> S[Base.localindexes(S)] = myid())
# mygamma = SharedArray(Float64, G, pids=workers())
# @everywhere gamma = SharedArray(Int32, G)

function fWriteText(stream::IOStream, param::Float64, M::Int64)
    for i = 1:M
        println(stream, i + param)
        println(stream, 3.14)
        # flush(stream)
    end
    return
end


function fWriteBinary(stream::IOStream, param::Float64, M::Int64)
    for i = 1:M
        write(stream, i + param)
        # flush(stream)
        write(stream, 3.14)
    end
    return
end

testfile1 = "testfile1.dat"
testfile2 = "testfile2.dat"
stream1 = open(testfile1, "a")
stream2 = open(testfile2, "a")
M = 9
param = 0.5

@time fWriteBinary(stream1, param, M)
# @time fWriteText(stream1, param, M)
flush(stream1)
# @time sleep(10)
# # @time fWriteText(stream1, param, M + 10)
# @time fWriteBinary(stream1, param, M)
# @time fWriteBinary(stream2, param, M)
close(stream1)
# close(stream2)

# readF = open(testfile1, "r")
# readF do 
# open(testfile1,"r") do f
    # asdf = read(f)
    # for i in asdf
    #     println(i)
    # end
    # seekstart(f)
    # for i in f
    #     println(typeof(i))
    # # for i in f
    # end
# end
   # readline(f)


# group LAMBDA
# iteration 1 ....