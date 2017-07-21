#= W. SEWLAL 1383337 # 2015-12-03
# MASTER FILE THAT RUNS THE PARALLEL TEMPERING ALGORITHM
=#


# TOGGLE SWITCHES FOR ALGORITHM
const loadZ = !true # Set to true to load a saved z. If false, z will be created and written to a file.
const saveResultsMH = !false # If false, results are loaded from a savefile.
const saveResultsPT = !false # If false, results are loaded from a savefile.

const savePlots = true # Store the plots in a pdf file.
const PTnoSwap = false # If true, swapping between tempered densities is not allowed
const InfoMessages = true # Show main info messages during the simulation
const verboseInfoMessages = !false # Show more info messages per iteration during the simulation. Used for debugging purposes.

const thinning = true # Save every 10th iteration in stead of every iteration

# ENVIRONMENT VARIABLES
pathProject = pwd() # Directory with all the .jl files for the CPP-simulations.
pathSavefiles = "/home/wikash/Documents/msc-thesis/Project/" # Directory where figures, z and results are saved.


# CHECK PATH FOR PROJECT FILES, CHECK IF SAVEFILES LOCATION IS SET.
!isfile(joinpath(pathProject, "filePathCheck.jl")) && throw(ArgumentError(string("filePathCheck.jl not found in ", pathProject , ".\nMake sure that pathProject is set correctly and that the file is present.\n")))
InfoMessages && info("Checking project files and setup of the saveFiles directory")
include(joinpath(pathProject, "filePathCheck.jl"))
using fileCheck; filePathCheck(pathSavefiles, pathProject, savePlots)

# IMPORT MODULES (USED AND CUSTOM MODULES)
include(joinpath(pathProject, "moduleCheck.jl"))
# cholfact!(A, [uplo::Symbol,] Val{false}) -> Cholesky

# LOAD/CREATE Z
# include(joinpath(pathProject, "configurations/zconfiguration1.jl"))
# include(joinpath(pathProject, "configurations/zconfiguration2.jl"))
include(joinpath(pathProject, "configurations/zconfiguration4.jl"))

const z = getZ(N, real_J, Delta, real_param, pathSavedZ, loadZ)
const I = length(z); const Z = N - I

# ALGORITHM SETUP
InfoMessages && verboseInfoMessages && info("Setting up algorithm parameters")
# SIMULATION SETUP -> see simulation files
# include(joinpath(pathProject, "configurations/simulationconfiguration1manual.jl"))
# include(joinpath(pathProject, "configurations/simulationconfiguration2manual.jl"))
# include(joinpath(pathProject, "configurations/simulationconfiguration2AA.jl"))
include(joinpath(pathProject, "configurations/simulationconfiguration4manual.jl"))
# include(joinpath(pathProject, "configurations/simulationconfiguration4AA.jl"))
# include(joinpath(pathProject, "configurations/simulationconfiguration4CA.jl"))

const G = length(gamma)

# Parallel tempering setup #iterations PT
# const M_PT = 500 
# const M_PT = 1000
# const M_PT = 1500
# const M_PT = 3000
# const M_PT = 5000
# const M_PT = 15000
# const M_PT = 25000
const M_PT = 50000
# const M_PT = 60000
# const M_PT = 100000
# const M_PT = 250000


# Metropolis-Hastings setup #iterations MH
# const M_MH = M_PT * G
# const M_MH = 1000
# const M_MH = 3000
# const M_MH = 5000
# const M_MH = 15000
# const M_MH = 20000
# const M_MH = 20000
# const M_MH = 60000
# const M_MH = 100000
const M_MH = 200000
# const M_MH = 1000000

# BURNIN SETUP
if M_PT < 25000
    const burninlabel = round(Int64, (M_PT / 5))
# elseif M_PT == 15000
#     const burninlabel = 5000
elseif M_PT < 50000
    const burninlabel = 15000
elseif M_PT < 100000
    const burninlabel = 20000
else
    const burninlabel = 20000
end



InfoMessages && verboseInfoMessages && info("Algorithm setup done")


# START MH
if saveResultsMH
    @time MHparamArray, MHacceptance = MetropolisHastings(z, N, J, I, M_MH, T, Delta, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa, thinning, InfoMessages, verboseInfoMessages)

    MHsaveResults(z, MHparamArray, MHacceptance, M_MH, real_param, real_J, J, pathSavedResults, InfoMessages, verboseInfoMessages)
else
    MHparamArray, MHacceptance = MHloadResults(M_MH, real_param, real_J, J, pathSavedResults, InfoMessages, verboseInfoMessages)
    InfoMessages && verboseInfoMessages && info("Metropolis-Hastings simulation results loaded")
end


# START PT
if saveResultsPT
    @time paramArray, distArray, acceptance, swapCounter = parallelTempering(z, N, J, I, G, M_PT, T, Delta, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa, PTnoSwap, bswap, gamma, UpdateModusMH, thinning, InfoMessages, verboseInfoMessages)
    @time PTsaveResults(z, paramArray, distArray, acceptance, swapCounter, gamma, M_PT, real_param, real_J, J, UpdateModusMH, PTnoSwap, pathSavedResults, InfoMessages, verboseInfoMessages)
else
    paramArray, distArray, acceptance, swapCounter = PTloadResults(gamma, M_PT, real_param, real_J, J, UpdateModusMH, PTnoSwap, pathSavedResults, InfoMessages, verboseInfoMessages)
    InfoMessages && verboseInfoMessages && info("Parallel Tempering simulation results loaded")
end

# writecsvMH(real_param, real_J, MHparamArray, MHacceptance, paramArray, acceptance, N, J, I, G, M_PT, M_MH, gamma, UpdateModusMH, PTnoSwap, pathFigures, InfoMessages, verboseInfoMessages)

# PLOTTING AND SAVING RESULTS
if savePlots
    plotFilePaths, ESSAC, normdiffMH, normdiffPT, lambdadiffRealMHPT = pdfPlotResults(real_param, real_J, MHparamArray, MHacceptance, paramArray, acceptance, N, J, I, G, M_PT, M_MH, gamma, burninlabel, UpdateModusMH, PTnoSwap, thinning, pathFigures, InfoMessages, verboseInfoMessages)

    (InfoMessages && !PTnoSwap) && PTreport(real_param, real_J, J, acceptance, swapCounter, gamma, UpdateModusMH, M_PT, ESSAC, normdiffMH, normdiffPT, lambdadiffRealMHPT, thinning, pathSavedLogs, plotFilePaths)
else
    (InfoMessages && !PTnoSwap) && PTreport(real_param, real_J, J, acceptance, swapCounter, gamma, UpdateModusMH, M_PT, ESSAC, normdiffMH, normdiffPT, lambdadiffRealMHPT, thinning, pathSavedLogs)
end





# # # const N = 10000; const real_J = 2
# # const N = 5000; const real_J = 4
# # const Delta = 1.0; const labdaDelta = 3.0; const T = (N * Delta)
# # # const real_psi = [1.8, 4.2]; const real_ksi = [0.0, 0]; const real_labda = labdaDelta / Delta;
# # # const real_mu  = [-2.0, 1.5]; const real_tau = 9.0;
# # # const real_psi = [0.3,0.4,0.2,0.1]; const real_ksi = [0.0, 0, 0, 0]; const real_labda = labdaDelta / Delta;
# # const real_psi = [0.9,1.2,0.6,0.3]; const real_ksi = [0.0, 0, 0, 0]; const real_labda = labdaDelta / Delta;
# # # const real_psi = [1.8,2.4,1.2,0.6]; const real_ksi = [0.0, 0, 0, 0]; const real_labda = labdaDelta / Delta;
# # const real_mu  = [-3.0,0.5,1.5,5.0]; const real_tau = 9.0;
# # const real_param = MHiterationParameters(real_labda, real_psi, real_mu, real_tau, real_ksi)
# # # CHECK SUM(PSI) == labda. This check in getZ() does not work correctly currently!



# # # # # const J = 4; # #jump-types
# # # # # const alpha_0 = 1.0; const beta_0 = 1.0;
# # # # # const alpha_1 = 1.0; const beta_1 = 1.0;
# # # # # # const alpha_1 = 1.0; const beta_1 = 1/std(z);
# # # # # const ksi = [0.0, 0, 0, 0]
# # # # # # const ksi = [0.0, 0]
# # # # # const kappa = 1.0
# # # # # 
# # # # # # const J = 2; # #jump-types
# # # # # # const UpdateModusMH = "PTA-TempDensity" # "PTA-TempDensity" # "PTB-TempPosterior"
# # # # # # const UpdateModusMH = "PTC-TempPrior"
# # # # # const UpdateModusMH = "PTD-TempTau"
# # # # # # const UpdateModusMH = "PTB-TempPosterior" # "PTA-TempDensity" # "PTB-TempPosterior"
# # # # # # const gamma = [1, 1 - (8e-4), 1 - (15e-4), 1 - (20e-4)]; const G = length(gamma) #temperatures
# # # # # # const gamma = [1.0, 1 - (4e-3), 1 - (8e-3), 1 - (12e-3)];
# # # # # # const gamma = [1.0, 1 - (5e-4), 1 - (10e-4), 1 - (15e-4)];
# # # # # const gamma = [1.0, 1.75, 2.5, 3];
# # # # # # const gamma = [1.0, 2.5, 5, 10];
# # # # # # const gamma = [1.0, 5, 0.6, 0.25];
# # # # # # const gamma = [1.0, 0.5, 0.25, 0.1];
# # # # # # const gamma = [1.0, 0.9, 0.5, 0.25];
# # # # # # const gamma = [1.0, 0.9995, 0.9990, 0.99985];
# # # # # # const gamma = [1.0, 0.6, 0.25, 0.1];
# # # # # # const gamma = [1.0, 0.8, 0.6, 0.25];
# # # # # # const gamma = [1.0, 0.99, 0.98, 0.97];
# # # # # # const gamma = [1.0, 0.95, 0.9, 0.8];
# # # # # # const gamma = [1.0, 1 - (1e-2), 1 - (2e-2), 1 - (3e-2)];
# # # # # # const gamma = [1.0, 1 + (5e-4), 1 + (10e-4), 1 + (15e-4)];
# # # # # # const gamma = [1.0, 1 + (1e-3), 1 + (2e-3), 1 + (3e-3)];
# # # # # # const gamma = [1.0, 1 + (8e-4), 1 + (16e-4), 1 + (24e-4)];
# # # # # const bswap = 0.8
