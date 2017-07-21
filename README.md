# MCMC-CompoundPoissonProcess-BayesianEstimation
Multicore implementations of MCMC algorithms for Bayesian estimation of a compound Poisson process

This implementation is written for Julia v0.6. 

Fresh restart of the original implementations. Implemented are Parallel Tempering/Exchange algorithms with a Metropolis-Hastings update step for estimating the intensity and jump density of a Compound Poisson Process, given samples of the increments.

An older, single core, implementation exists for julia v0.4 and might be published.

Reference: thesis
