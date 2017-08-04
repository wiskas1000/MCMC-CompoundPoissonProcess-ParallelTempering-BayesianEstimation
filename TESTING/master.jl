# W. SEWLAL 1383337 # 2015-12-03
@everywhere using ParallelDataTransfer
@everywhere using JLD2

const gamma_temperatures = Array{Float64, 1}(4)
const G = length(gamma_temperatures)
gamma_temperatures[1:G] = [1.0, 1.0 + 1e-4, 1.0 + 2e-4, 1.0 + 3e-4]
@broadcast G = getfrom(1, :G)
@broadcast const gamma_temperatures = getfrom(1, :gamma_temperatures)

pathTMPStream = "/home/wikash/Documents/Simulations/TMP/"

@everywhere include("parallel-startup.jl")
@everywhere using PTstartup

@broadcast setIndexMyGamma(idWorker)
@broadcast setMyPathTMPStream(getfrom(1, :pathTMPStream))
@broadcast println(gamma_temperatures[getIndexMyGamma()])
# @broadcast println(getIndexMyGamma())
# @broadcast println(getMyPathTMPStream())



duration = 0.06 # Simulation duration (seconds)
function stop(a::Timer)
    bContinueRunning[1] = false
end

Timer(stop::Function, 60, 0)
thinning = 10 # Store every $(thinning)-th iteration of the simulation

bContinueRunning = Array{Bool, 1}(1)
bContinueRunning[1] = true
i = 0

Timer(stop::Function, 0.1, 0)
while(bContinueRunning[1])
    if ((i % thinning) == 0)
        println(i)
    end
    i += 1
end

@broadcast const bufferLength = 1000
@broadcast const N = 5000

@everywhere include("CPPParallelTempering.jl")
@everywhere using CPPParallelTempering

include("generate-data.jl")
z = generateCPP(100, 2, 1.0, 3.0, [1.2, 1.8], [-2.2, 3], 4.0)
println(z)


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
        write(stream, 3.14)
    end
    return
end

testfile1 = "testfile1.csv"
testfile2 = "testfile2.dat"
stream1 = open(testfile1, "w")
M = 9
param = 0.5

a = gamma_temperatures
writecsv(stream1, a')
flush(stream1)
close(stream1)
stream2 = open(testfile1, "a")
writecsv(stream2, a')
writecsv(stream2, [1.+a'; 1.+a'; 21.+a'; 31.+a'])

flush(stream2)
close(stream2)


# stream2 = open(testfile2, "a")
# @time fWriteBinary(stream1, param, M)
# @time fWriteText(stream1, param, M)
# @time sleep(10)
# # @time fWriteText(stream1, param, M + 10)
# @time fWriteBinary(stream1, param, M)
# @time fWriteBinary(stream2, param, M)
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