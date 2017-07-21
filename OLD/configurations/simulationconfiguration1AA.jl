#jump-types
const J = 2;

# Prior setup
const alpha_0 = 1.0; const beta_0 = 1.0;
const alpha_1 = 1.0; const beta_1 = 1.0;
const ksi = [0.0, 0];
const kappa = 1.0;

# Parallel tempering setup; UpdateModusMH = "PTA-TempDensity" # "PTA-TempDensity" # "PTB-TempPosterior" # "PTC-TempPrior"
const UpdateModusMH = "PTA-TempDensity";
const gamma = [1.0, 1 - (5e-4), 1 - (10e-4), 1 - (15e-4)];
const bswap = 0.8;

# const alpha_1 = 1.0; const beta_1 = 1/std(z);
# const ksi = [0.0, 0, 0, 0]

# const gamma = [1, 1 - (8e-4), 1 - (15e-4), 1 - (20e-4)]; const G = length(gamma) #temperatures
# const gamma = [1.0, 1.75, 2.5, 3];
# const gamma = [1.0, 1 - (4e-3), 1 - (8e-3), 1 - (12e-3)];
# const gamma = [1.0, 1 - (5e-4), 1 - (10e-4), 1 - (15e-4)];
# const gamma = [1.0, 2.5, 5, 10];
# const gamma = [1.0, 5, 0.6, 0.25];
# const gamma = [1.0, 0.5, 0.25, 0.1];
# const gamma = [1.0, 0.9, 0.5, 0.25];
# const gamma = [1.0, 0.9995, 0.9990, 0.99985];
# const gamma = [1.0, 0.6, 0.25, 0.1];
# const gamma = [1.0, 0.8, 0.6, 0.25];
# const gamma = [1.0, 0.99, 0.98, 0.97];
# const gamma = [1.0, 0.95, 0.9, 0.8];
# const gamma = [1.0, 1 - (1e-2), 1 - (2e-2), 1 - (3e-2)];
# const gamma = [1.0, 1 + (5e-4), 1 + (10e-4), 1 + (15e-4)];
# const gamma = [1.0, 1 + (1e-3), 1 + (2e-3), 1 + (3e-3)];
# const gamma = [1.0, 1 + (8e-4), 1 + (16e-4), 1 + (24e-4)];
