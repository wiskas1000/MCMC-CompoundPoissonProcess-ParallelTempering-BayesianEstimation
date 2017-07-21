const N = 5000; 
const real_J = 4;
const Delta = 1.0; const labdaDelta = 1.0; const T = (N * Delta);

const real_psi = [0.3, 0.4, 0.2, 0.1]; const real_ksi = [0.0, 0.0, 0.0, 0.0]; const real_labda = labdaDelta / Delta;
const real_mu  = [-3.0, 0.5, 1.5, 5.0]; const real_tau = 9.0;
const real_param = MHiterationParameters(real_labda, real_psi, real_mu, real_tau, real_ksi);
