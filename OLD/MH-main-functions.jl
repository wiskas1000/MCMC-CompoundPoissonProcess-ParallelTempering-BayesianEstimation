# W. SEWLAL 1383337

"""
`MetropolisHastings(z, N, J, I, M, T, Delta, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa, infoMessages, infoMessagesVerbose)`

Decompounds a Compound Poisson Process, given by its increments `z`, using a Metropolis-Hastings MCMC scheme.

Keyword arguments
-----------------

* `z` : Array{Float64,1}
* `N` : Int64
* `J` : Int64 -- Number of jump types
* `I` : Int64
* `M` : Int64
* `T` : Float64
* `Delta` : Float64
* `alpha_0` : Float64
* `beta_0` : Float64
* `alpha_1` :Float64
* `beta_1` : Float64
* `ksi` : Array{Float64,1} -- `J`-element array
* `kappa` : Float64
* `infoMessages` : Bool
* `infoMessagesVerbose` : Bool
"""
function MetropolisHastings(z::Array{Float64,1}, N::Int64, J::Int64, I::Int64, M::Int64, T::Float64, Delta::Float64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64, thinning::Bool, infoMessages::Bool, infoMessagesVerbose::Bool)
    infoMessages && info("Simulation start: Single Metropolis-Hastings (gamma = 1)")
    # INITIALIZE VARIABLES AND INITIAL STATES OF A AND N
    if (thinning)
        MHparamArray::Array{MHiterationParameters, 1} = Array(MHiterationParameters, Int(round(Int, M / 10)))
        MHacceptance::Array{Float64, 1} = zeros(Float64, Int(round(Int, M / 10)))
        maxiter::Int32 = M - 10
    else
        MHparamArray = Array(MHiterationParameters, M)
        MHacceptance = zeros(Float64, M)
        maxiter = M - 1
    end
    
    current_a::Array{Int64, 2} = zeros(Int64, J, I)
    current_n::Array{Int64, 1} = zeros(Int64, I)    

    # Computation-efficient constants
    const zSq::Array{Float64, 1} = abs2(z)
    const beta_0_plus_T::Float64 = beta_0 + T
    const alpha_1_plus_half_I::Float64 = alpha_1 + (I / 2)
    const kappa_times_Eye_J::Array{Float64,2} = kappa * eye(J)

    # Variables for reduction of garbage collection
    proposal_a::Array{Int64, 1} = zeros(Int64, J); proposal_n::Array{Int64, 1} = zeros(Int64, 1);

    # START ALGORITHM
    MHinitializeParameters!(MHparamArray, J, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
    MHinitializeCurrentStates!(current_a, current_n, proposal_a, proposal_n, J, I, (MHparamArray[1].psi * Delta));
    current_iterationParameters::MHiterationParameters = MHparamArray[1]
    current_acceptance::Float64 = 0

    # UPDATING EACH STEP
    for m in 1:maxiter
        infoMessages && (m % 500 == 0 && info(string("Start iteration m = ", m, " of ", M, "\t(", round(m/M * 100, 1) ,"%)")))
        current_acceptance = 0
        # Update segments
        for i in 1:I
            # Setup proposal
            succes = MHsegmentUpdate!(proposal_n, proposal_a, current_n, current_a, J, z[i], i, Delta, current_iterationParameters)
            if (succes)
                current_acceptance += 1
            end
        end

        # Update parameters
        out_labda, out_psi = MHupdatePsiLabda(J, alpha_0, beta_0_plus_T, current_a)
        out_mu, out_tau, out_ksi = MHupdateTauKsiMu(kappa_times_Eye_J, z, zSq, J, I, alpha_1_plus_half_I, beta_1, kappa, current_n, current_a, current_iterationParameters.ksi)
        current_iterationParameters = MHiterationParameters(out_labda, out_psi, out_mu, out_tau, out_ksi)
        if (thinning)
            if (m%10 == 0)
               MHparamArray[Int(round(m/10)) + 1] = current_iterationParameters
               MHacceptance[Int(round(m/10)) + 1] = current_acceptance
            end
        else
            MHparamArray[m + 1] = current_iterationParameters
            MHacceptance[m + 1] = current_acceptance
        end
    end

    MHacceptance = MHacceptance / I
    return MHparamArray, MHacceptance
end

# """
# `MHinitializeVariables(J, M, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)`
#
# Initializes the main variables of the algorithm and assigns an initial value using prior distributions.
#
# """
function MHinitializeParameters!(paramArray::Array{MHiterationParameters, 1}, J::Int64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64)
    first_dist_psi = Array(Distribution, J)
    first_dist_psi[:] = Gamma(alpha_0, (1 / beta_0))
    first_dist_tau = Gamma(alpha_1, (1 / beta_1))

    first_psi = Array(Float64, J)
    for j in 1:J
        first_psi[j] = rand(first_dist_psi[j])
    end
    first_labda = sum(first_psi)
    first_tau = rand(first_dist_tau)

    first_dist_mu = MvNormal(ksi, (1 / (first_tau * kappa)) * eye(J))
    first_mu  = rand(first_dist_mu)

    paramArray[1] = MHiterationParameters(first_labda, first_psi, first_mu, first_tau, ksi)
    return
end


"""
`MHinitializeCurrentStates(J, I, psiDelta)`

Initializes `a` using information of the prior distribution. Returns an array with an initial value on the number of jumps per segment and a matrix with the number of jumps per jumptype for each nonzero segment (for each segment where `z` != 0).

Keyword arguments
-----------------

* `J` : Int64 -- Number of jump types
* `I` : Int64 -- Number of nonzero segments
* `psiDelta` : Array{Float64,1} -- `J`-element array with the jump-intensity
"""
function MHinitializeCurrentStates!(current_a::Array{Int64, 2}, current_n::Array{Int64, 1}, proposal_a::Array{Int64, 1}, proposal_n::Array{Int64, 1}, J::Int64, I::Int64, psiDelta::Array{Float64,1})
    for i in 1:I
        MHproposeNiAi!(proposal_n, proposal_a, J, psiDelta)
        current_n[i] = proposal_n[1]
        current_a[:, i] = proposal_a
    end
    return
end




"""
`MHproposeNiAi!(ni::Array{Int64, 1}, ai::Array{Int64, 1}, J::Int64, intensity::Array{Float64,1})`

Returns a proposal on the amount of jumps for a single segment. Returns the total number of jumps and an array with the amount of jumps for each type. The proposed number of jumps of type `j` is drawn from a Poisson distribution with mean `intensity[j]`, conditioned to have a positive outcome.

Keyword arguments
-----------------

* `J` : Int64 -- Number of jump types
* `intensity` : Array{Float64,1} -- `J`-element array giving the intensity for each element
"""
function MHproposeNiAi!(ni::Array{Int64, 1}, ai::Array{Int64, 1}, J::Int64, intensity::Array{Float64,1})
    ni[1] = 0
    while (ni[1] <= 0)
        for j in 1:J
            ai[j] = rand(Poisson(intensity[j]))
            ni[1] += ai[j]
        end
    end
    return
end




"""
`MHsegmentUpdate(J, z, i, psiDelta, current_n, current_a, in_mu, in_tau)`

Metropolis-Hastings update for segment `i`. This function MHwill modify current_a and current_n if the proposal is accepted.

Keyword arguments
-----------------

* `J` : Int64 -- Number of jump types
* `z` : Array{Float64,1} -- `I`-element array differences in a sample of a Compound Poisson Process
* `i` : Int64 -- Segment. Warning: `i` must be between `1` and `I`. There is no check implemented for this.
* `psiDelta` : Array{Float64,1} -- `J`-element array with the jump-intensity
* `current_n` :
* `current_a` :
* `in_mu` :
* `in_tau` :
"""
function MHsegmentUpdate!(proposal_n::Array{Int64, 1}, proposal_a::Array{Int64, 1}, current_n::Array{Int64, 1}, current_a::Array{Int64, 2}, J::Int64, in_z::Float64, i::Int64, Delta::Float64, in_param::MHiterationParameters)
    succes::Bool = false
    MHproposeNiAi!(proposal_n, proposal_a, J, in_param.psi * Delta)

    # Setup acceptance
    numer = pdf(Normal(dot(proposal_a     , in_param.mu), sqrt(proposal_n[1] / in_param.tau)), in_z)
    denom = pdf(Normal(dot(current_a[:, i], in_param.mu), sqrt(current_n[i]  / in_param.tau)), in_z)
    denom == 0 && error("in MHsegmentUpdate():\t Denominator is 0")

    A = numer / denom

    r1::Float64 = rand()

    # Update segment
    if(r1 < min(A, 1))
        current_a[:, i] = proposal_a;
        current_n[i] = proposal_n[1];
        succes = true;
    end
    return succes
end

function MHupdatePsiLabda(J::Int64, alpha_0::Float64, beta_0_plus_T::Float64, current_a::Array{Int64, 2})
    const s = sum(current_a, 2)
    psi = zeros(Float64, J)

    for j in 1:J
        psi[j] = rand(Gamma((alpha_0 + s[j]), (1 / beta_0_plus_T)))
    end
    return sum(psi), psi
end


function MHupdateTauKsiMu(kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, J, I, alpha_1_plus_half_I::Float64, beta_1, kappa, current_n, current_a, in_ksi)
    P::Array{Float64,2} = kappa_times_Eye_J
    q = kappa * in_ksi
    R = kappa * dot(in_ksi, in_ksi)
    for i in 1:I
        P += (current_a[:, i] * current_a[:, i]') / current_n[i]
        for j in 1:J
            q[j] += (current_a[j, i] / current_n[i]) * z[i]
        end
        R += zSq[i] / current_n[i]
    end
    Pinv = inv(Symmetric(P))

    tau::Float64 = rand(Gamma((alpha_1_plus_half_I), (1 / (beta_1 + sum(R - q' * Pinv * q) / 2))))
    ksi::Array{Float64,1} = Pinv * q
    mu::Array{Float64,1} = rand(MvNormal(Pinv * q, (1 / tau) .* Pinv))
    return mu, tau, ksi
end


function MHsaveResults(z::Array{Float64,1}, paramArray::Array{MHiterationParameters,1}, acceptance::Array{Float64,1}, M::Int64, real_param::MHiterationParameters, real_J::Int64, J::Int64, dir::AbstractString, infoMessages::Bool, infoMessagesVerbose::Bool)
    infoMessages && info("Saving results")
    dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J)); makeDirectory(dir)

    fnameResultsBase = string("MH_M=", round(Int, M / 1000), "k_", real_param)
    fnameResultsExtension = ".jld"
    savefilePath = fnameNew(fnameResultsBase, fnameResultsExtension, dir)
    save(savefilePath, "Z", z, "Parameters", paramArray, "Acceptance", acceptance, "M", M)
    infoMessages && infoMessagesVerbose && info("Saving complete")
    return
end


function MHloadResults(M::Int64, real_param::MHiterationParameters, real_J::Int64, J::Int64, dir::AbstractString, infoMessages::Bool, infoMessagesVerbose::Bool)
    infoMessages && info("No simulation performed (saveResultsMH is set to false). Loading Metropolis-Hastings simulation results")

    fnameResultsBase = string("MH_M=", round(Int, M / 1000), "k_", real_param)
    fnameResultsExtension = ".jld"
    dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J))
    loadedvar = loadFile(fnameResultsBase, fnameResultsExtension, dir)
#     z = loadedvar["Z"]
#     println(length(z))
#     getZ(10000, real_J, 1.0, real_param, "/home/wikash/Documents/msc-thesis/Project/SimulationVariables/", false)
    return loadedvar["Parameters"], loadedvar["Acceptance"]
end
