# W. SEWLAL 1383337 # 2015-12-03
using Distributions
using JLD2
"""
`generateCPP(N::Int64, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64)`

Returns a truncated array z of a Compound Poisson Process. The returned array z consists of `I` nonzero differences in a sample of a Compound Poisson Process. Assumed is that the jump size distribution is a weighted mixture of `J` Gaussians with weights `psi`/`labda`, means `mu` and variance 1/`tau` (equal for all jump types). The total number of jumps is drawn from a Poisson distribution with mean `labda` * `Delta`.

Keyword arguments
-----------------

* `N` : Int64 -- Number of segments
* `J` : Int64 -- Number of jump types
* `Delta` : Float64 -- The time between observations, fixed for all observations
"""
function generateCPPDataset(N::Int64, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64)::Array{Float64, 1}
    # Check size psi, mu
    length(psi) != J && throw(ArgumentError("generateCPP(): psi is not of length J"))
    length(mu) != J && throw(ArgumentError("generateCPP(): mu is not of length J"))
    # Check if labda = sum(psi)
    (abs(sum(psi) - labda) > 1e15) && throw(ArgumentError("generateCPP(): sum(psi) does not equal labda"))

    # Create z, zero values of z are on the end-side of the vector
    original_z = zeros(Float64, N);
    tmpai = zeros(Int64, J); tmpni::Int64 = 0; nonzeroCounter::Int64 = 0;

    for i in 1:N
        tmpni = 0
        for j in 1:J
            tmpai[j] = rand(Poisson(psi[j] * Delta))
            tmpni += tmpai[j]
        end
        if(tmpni > 0)
            nonzeroCounter += 1
            original_z[nonzeroCounter] = rand(Normal(dot(tmpai, mu), sqrt(tmpni / tau)))
        end
    end
    z = original_z[1:nonzeroCounter]
    return z
end
function generateCPPDataset(N::Int64, configuration::Int64)::Array{Float64, 1}
    J, Delta, labda, psi, mu, tau = getConfigParameters(configuration)
    return generateCPPDataset(N, J, Delta, labda, psi, mu, tau)
end

# CHANGE pwd() into basepath
function getConfigParameters(configuration::Int64)
    if configuration == 1
        J = 2;
        Delta = 1.0;
        psi = [0.3, 0.7];
        labda = 1.0;
        mu  = [-2.0, 1.5];
        tau = 9.0;
        return J, Delta, labda, psi, mu, tau
    elseif configuration == 2
        J = 2;
        Delta = 1.0;
        psi = [1.8, 4.2];
        labda = 6.0;
        mu  = [-2.0, 1.5];
        tau = 9.0;
        return J, Delta, labda, psi, mu, tau
    elseif configuration == 3
        J = 2;
        Delta = 1.0;
        psi = [0.8, 0.2];
        labda = 1.0;
        mu  = [2.0, -1.0];
        tau = 1.0;
        return J, Delta, labda, psi, mu, tau
    elseif configuration == 4
        J = 2;
        Delta = 1.0;
        psi = [4.8, 1.2];
        labda = 6.0;
        mu  = [2.0, -1.0];
        tau = 1.0;
        return J, Delta, labda, psi, mu, tau
    elseif configuration == 5
        J = 4;
        Delta = 1.0;
        psi = [0.3, 0.4, 0.2, 0.1];
        labda = 1.0;
        mu  = [-3.0, 0.5, 1.5, 5.0];
        tau = 9.0;
        return J, Delta, labda, psi, mu, tau
    elseif configuration == 6
        J = 4;
        Delta = 1.0;
        psi = [0.9, 1.2, 0.6, 0.3];
        labda = 3.0;
        mu  = [-3.0, 0.5, 1.5, 5.0];
        tau = 9.0;
        return J, Delta, labda, psi, mu, tau
    else
        errorMsg = string("Configuration ", configuration, " not found")
        throw(ArgumentError(errorMsg))
    end
end
    # basedir = joinpath(pwd(), "configurations")
    # fname = string("dataset-z-configuration-", configuration, ".jl")
    # include(joinpath(basedir, fname))


# function initializeZ(N::Int64, J::Int64, Delta::Float64, param::MHiterationParameters)
# param::MHiterationParameters

function saveCPPDataset(dir::AbstractString, N::Int64, z::Array{Float64, 1}, configuration::Int64)
    savepath = joinpath(dir, string("configuration", configuration))
    makeDirectory(savepath)
    fnameBase = getCPPDatasetFilenameBase(dir, N, configuration)
    fnameExtension = ".jld2"
    fnamepath = proposeFilenamePath(fnameBase, fnameExtension, savepath)
    save(fnamepath, "z", z)
    return
end

function saveCPPDataset(dir::AbstractString, N::Int64, z::Array{Float64, 1}, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64)
    savepath = joinpath(dir, "manual")
    makeDirectory(savepath)
    fnameBase = getCPPDatasetFilenameBase(dir, N, J, Delta, labda, psi, mu, tau)
    fnameExtension = ".jld2"
    fnamepath = proposeFilenamePath(fnameBase, fnameExtension, savepath)
    save(fnamepath, "z", z)
    return
end

function getCPPDatasetFilenameBase(dir::AbstractString, N::Int64, configuration::Int64)::AbstractString
    if (N < 1000)
        fnameBase = string("zCPP-c", configuration, "-N-", N)
    else
        fnameBase = string("zCPP-c", configuration, "-N-", round(Int, N / 1000), "k")
    end
    return fnameBase
end

function getCPPDatasetFilenameBase(dir::AbstractString, N::Int64, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64)::AbstractString
    if (N < 1000)
        fnameBase = string("zCPP-manual-N-", N)
    else
        fnameBase = string("zCPP-manual-N-", round(Int, N / 1000), "k")
    end
    fnameBase = string(fnameBase, "-J", J, "-", Delta, "-", labda, "-", psi, "-", mu,  "-", tau)
    return fnameBase
end


function loadCPPDataset(dir::AbstractString, N::Int64, configuration::Int64)::Array{Float64, 1}
    loadpath = joinpath(dir, string("configuration", configuration))
    fnameBase = getCPPDatasetFilenameBase(dir, N, configuration)
    fnameExtension = ".jld2"    
    loadedvars = loadFile(fnameBase, fnameExtension, loadpath)
    return loadedvars["z"]
end

function loadCPPDataset(dir::AbstractString, N::Int64, configuration::Int64, version::Int64)::Array{Float64, 1}
    loadpath = joinpath(dir, string("configuration", configuration))
    fnameBase = getCPPDatasetFilenameBase(dir, N, configuration)
    fnameExtension = ".jld2"    
    loadedvars = loadFile(fnameBase, fnameExtension, loadpath, version)
    return loadedvars["z"]
end
function loadCPPDataset(dir::AbstractString, N::Int64, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64)::Array{Float64, 1}
    loadpath = joinpath(dir, "manual")
    fnameBase = getCPPDatasetFilenameBase(dir, N, J, Delta, labda, psi, mu, tau)
    fnameExtension = ".jld2"
    loadedvars = loadFile(fnameBase, fnameExtension, loadpath)
    return loadedvars["z"]
end

function loadCPPDataset(dir::AbstractString, N::Int64, J::Int64, Delta::Float64, labda::Float64, psi::Array{Float64, 1}, mu::Array{Float64, 1}, tau::Float64, version::Int64)::Array{Float64, 1}
    loadpath = joinpath(dir, "manual")
    fnameBase = getCPPDatasetFilenameBase(dir, N, J, Delta, labda, psi, mu, tau)
    fnameExtension = ".jld2"
    loadedvars = loadFile(fnameBase, fnameExtension, loadpath, version)
    return loadedvars["z"]
end

######### DEPRECATED FUNCTIONS BELOW #########

# """
# `getZ(N::Int64, real_J::Int64, Delta::Float64, real_param::MHiterationParameters, savedZPath::AbstractString, loadZ::Bool)`

# Returns a truncated array `z` of a Compound Poisson Process. The returned array z consists of `I` nonzero differences in a sample of a Compound Poisson Process. If `loadZ` is set to true, the function will load `z` from a savefile. Else, it will call `initializeZ` to generate `z`.

# Keyword arguments
# -----------------

# * `N` : Int64 -- Number of segments
# * `real_J` : Int64 -- Number of jump types
# * `Delta` : Float64 -- The time between observations, fixed for all observations
# * `real_param` : MHiterationParameters -- True parameters
# """
# function getZ(N::Int64, real_J::Int64, Delta::Float64, real_param::MHiterationParameters, savedZPath::AbstractString, loadZ::Bool)
#     fnameZBase = string("Z-",real_param)
#     fnameZExtension = ".jld"

#     if(loadZ)
#         loadedvar = loadFile(fnameZBase, fnameZExtension, savedZPath)
#         return loadedvar["Z"]
#     else
#         newZ = initializeZ(N, real_J, Delta, real_param)
#         savefilePath = fnameNew(fnameZBase, fnameZExtension, savedZPath)
#         save(savefilePath, "Z", newZ, "param", real_param)
#         return newZ
#     end
# end


# """
# `initializeZ(N::Int64, J::Int64, Delta::Float64, real_param::MHiterationParameters)`

# Returns a truncated array z of a Compound Poisson Process. The returned array z consists of `I` nonzero differences in a sample of a Compound Poisson Process. Assumed is that the jump size distribution is a weighted mixture of `J` Gaussians with weights `real_param.psi`/`real_param.labda`, means `real_param.mu` and variance 1/`real_param.tau` (equal for all jump types). The total number of jumps is drawn from a Poisson distribution with mean `real_param.labda` * `Delta`.

# Keyword arguments
# -----------------

# * `N` : Int64 -- Number of segments
# * `J` : Int64 -- Number of jump types
# * `Delta` : Float64 -- The time between observations, fixed for all observations
# * `real_param` : MHiterationParameters
# """
# function initializeZ(N::Int64, J::Int64, Delta::Float64, real_param::MHiterationParameters)
#     # Check size psi, mu
#     length(real_param.psi) != J && throw(ArgumentError("Error: real_psi or real_mu is not of length J"))
#     length(real_param.mu) != J && throw(ArgumentError("Error: real_psi or real_mu is not of length J"))
#     # Check if labda = sum(psi)
#     (abs(sum(real_param.psi) - real_param.labda) > 1e15) && throw(ArgumentError("Error: sum(real_psi) does not equal real_labda"))

#     # Create z, zero values of z are on the end-side of the vector
#     original_z = zeros(Float64, N);
#     tmpai = zeros(Int64, J); tmpni::Int64 = 0; nonzeroCounter::Int64 = 0;

#     for i in 1:N
#         tmpni = 0
#         for j in 1:J
#             tmpai[j] = rand(Poisson(real_param.psi[j] * Delta))
#             tmpni += tmpai[j]
#         end
#         if(tmpni > 0)
#             nonzeroCounter += 1
#             original_z[nonzeroCounter] = rand(Normal(dot(tmpai, real_param.mu), sqrt(tmpni / real_param.tau)))
#         end
#     end
#     z = original_z[1:nonzeroCounter]
#     return z
# end

