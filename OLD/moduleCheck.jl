# W. SEWLAL 1383337
InfoMessages && info("Checking and importing modules")
usedModules = ["Distributions", "JLD", "Cairo", "Gadfly", "StatsBase"]

# Check if all packages are installed
missingPackages = ""
for moduleName in usedModules
    try
        Pkg.installed(moduleName)
    catch
        missingPackages = string(moduleName, ", ", missingPackages)
    end
end
!isempty(missingPackages) && error(string("Packages ", missingPackages, " not installed."))

# Load modules
using Distributions
using JLD
using Cairo
using Gadfly
using StatsBase

# Load custom modules
include(joinpath(pathProject, "ParametersCPPMH.jl"))
include(joinpath(pathProject, "ParallelTemperingCPPMH.jl"))
include(joinpath(pathProject, "PlottingCPPMH.jl"))

using ParametersCPPMH
using ParallelTemperingCPPMH
using PlottingCPPMH
InfoMessages && info("Modules loaded")