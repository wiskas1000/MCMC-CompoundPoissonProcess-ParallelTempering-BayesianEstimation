#jump-types
const J = 4;

# Prior setup
const alpha_0 = 1.0; const beta_0 = 1.0;
const alpha_1 = 1.0; const beta_1 = 1.0;
const ksi = [0.0, 0.0, 0.0, 0.0]
const kappa = 1.0;

# Parallel tempering setup; UpdateModusMH = "PTA-TempDensity" # "PTA-TempDensity" # "PTB-TempPosterior" # "PTC-TempPrior"
const UpdateModusMH = "PTA-TempDensity";
const gamma = [1.0, 1 - (7e-4), 1 - (14e-4), 1 - (21e-4)];
const bswap = 0.8;

# const alpha_1 = 1.0; const beta_1 = 1/std(z);
