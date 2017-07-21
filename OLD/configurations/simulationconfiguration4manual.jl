# #jump-types
# const J = 4;
# 
# # Prior setup
# const alpha_0 = 1.0; const beta_0 = 1.0;
# const alpha_1 = 1.0; const beta_1 = 1.0;
# const ksi = [0.0, 0.0, 0.0, 0.0]
# const kappa = 1.0;
# 
# # Parallel tempering setup; UpdateModusMH = "PTA-TempDensity" # "PTA-TempDensity" # "PTB-TempPosterior" # "PTC-TempPrior"
# const UpdateModusMH = "PTB-TempPosterior";
# # const gamma = [1.0, 1 - (5e-4), 1 - (10e-4), 1 - (15e-4)];
# const gamma = [1.0, 1 + (5e-4), 1 + (10e-4), 1 + (15e-4)];
# const bswap = 0.8;
# 
# # const alpha_1 = 1.0; const beta_1 = 1/std(z);

#jump-types
const J = 4;

# Prior setup
const alpha_0 = 1.0; const beta_0 = 1.0;
# const alpha_1 = 1.0; const beta_1 = 1.0;
const alpha_1 = 1.0; const beta_1 = std(z);
const ksi = [0.0, 0.0, 0.0, 0.0]
const kappa = 1.0;

# Parallel tempering setup; UpdateModusMH = "PTA-TempDensity" # "PTA-TempDensity" # "PTB-TempPosterior" # "PTC-TempPrior"
const UpdateModusMH = "PTA-TempDensity";
# const UpdateModusMH = "PTB-TempPosterior";
# const UpdateModusMH = "PTC-TempPrior";
# const UpdateModusMH = "PTD-TempTau";
# const gamma = [1.0, 1 - (5e-4), 1 - (10e-4), 1 - (15e-4)];
# const gamma = [1.0, 1 - (1e-4), 1 - (1e-3), 1 - (1e-2)];
# const gamma = [1.0, 1 - (5e-3), 1 - (10e-3), 1 - (15e-3)];
# const gamma = [1.0, 1 + (5e-3), 1 + (10e-3), 1 + (15e-3)];
# const gamma = [1.0, 0.6, 0.25, 0.1];
# const gamma = [1.0, 0.9, 0.75, 0.6, 0.5, 0.25, 0.15, 0.1];
# const gamma = [1.0, 1.5, 2.5, 5];
# const gamma = [1.0, 0.8, 0.6, 0.25];
# const gamma = [1.0, 1 + (5e-4), 1 + (10e-4), 1 + (15e-4)];
const gamma = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
const bswap = 0.95;

# const alpha_1 = 1.0; const beta_1 = 1/std(z);
