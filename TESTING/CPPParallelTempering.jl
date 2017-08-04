# W. SEWLAL 1383337 # 2015-12-03
module CPPParallelTempering
using Distributions
using FileIOExtensionsCPP

include("generate-dataset.jl")
export generateCPPDataset, saveCPPDataset

# include("parallel-startup.jl")
# include("MH-main-functions.jl")
# export MetropolisHastings, MHproposeNiAi!, MHsaveResults, MHloadResults

# include("PT-main-functions.jl")
# export parallelTempering, PTreport, PTsaveResults, PTloadResults

end
