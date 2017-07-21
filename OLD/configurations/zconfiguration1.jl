const N = 5000; 
const real_J = 2;
const Delta = 1.0; const labdaDelta = 1.0; const T = (N * Delta);

const real_psi = [0.3, 0.7]; const real_ksi = [0.0, 0.0]; const real_labda = labdaDelta / Delta;
const real_mu  = [-2.0, 1.5]; const real_tau = 9.0;
const real_param = MHiterationParameters(real_labda, real_psi, real_mu, real_tau, real_ksi);
