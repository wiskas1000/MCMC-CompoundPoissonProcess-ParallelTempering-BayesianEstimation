#= W. SEWLAL 1383337 # 2015-12-03
# Indices order: (J, I, G, M)
# N: #segments
# J: #jump-types
# I: #nonzero segments
# G: #parallel
# M: #runs
=#

function parallelTempering(z::Array{Float64,1}, N::Int64, J::Int64, I::Int64, G::Int64, M::Int64, T::Float64, Delta::Float64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64, noSwap::Bool, bswap::Float64, gamma::Array{Float64,1}, UpdateModusMH::AbstractString, thinning::Bool, InfoMessages::Bool, detailedInfoMessages::Bool)
    InfoMessages && info(string("Simulation start: Parallel tempering (gamma = ", gamma, ", mode = " , UpdateModusMH, ")"))

    # INITIALIZE VARIABLES AND INITIAL STATES OF A AND N
    if (thinning)
        PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, Int(round(Int, M / 10)))
        PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, Int(round(Int, M / 10)))
        PTacceptance::Array{Float64,2} = zeros(Float64, G, Int(round(Int, M / 10)))
        maxiter::Int32 = M - 10
    else
        PTparamArray = Array(MHiterationParameters, G, M)
        PTdistArray = Array(MHiterationDistributions, G, M)
        PTacceptance = zeros(Float64, G, M)
        maxiter = M - 1
    end
    
#     PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, M)
#     PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, M)
#     PTacceptance::Array{Float64,2} = zeros(Float64, G, M)
    PTcurrent_a::Array{Int64,3} = zeros(Int64, J, I, G);
    PTcurrent_n::Array{Int64, 2} = zeros(Int64, I, G);

    # Computation-efficient constants
    const mode1::Bool =  UpdateModusMH == "PTA-TempDensity"
    const mode2::Bool =  UpdateModusMH == "PTB-TempPosterior"
    const mode3::Bool =  UpdateModusMH == "PTC-TempPrior"
    const mode4::Bool =  UpdateModusMH == "PTD-TempTau"
    const zSq::Array{Float64,1} = abs2(z)
    if mode1
        const psi_rate::Array{Float64,1} = (gamma * (beta_0 + T)) .^ -1
        const tau_shape::Array{Float64,1} = gamma * (alpha_1 + (I / 2)) + ((1 - gamma) * (1 - (J / 2)))
    elseif mode2
        const psi_rate = (beta_0 + gamma * T) .^ -1
        const tau_shape = alpha_1 + gamma  * (I / 2)
    elseif mode3
        const psi_rate = (gamma * beta_0 + T) .^ -1
        const tau_shape = (gamma * alpha_1 + (1 - gamma) * (1 - (J/2)) + (I / 2))
    elseif mode4
        const psi_rate = zeros(Float64, G) + (beta_0 + T) .^ -1
        const tau_shape = (gamma * (alpha_1 - 1) + (I / 2) + 1)
    end
    const kappa_times_Eye_J::Array{Float64,2} = kappa * eye(J)
    const lengthA::Int64 = length(PTcurrent_a[:, :, 1])
    tmpInt::Array{Int64,1} = zeros(Int64, 4)
    tmpInt[1] = lengthA
    tmpFloats::Array{Float64,1} = zeros(Float64, 5)

    # Variables for reduction of garbage collection
    proposal_a::Array{Int64, 1} = zeros(Int64, J); proposal_n::Array{Int64, 1} = zeros(Int64, 1);

    # START ALGORITHM
    swapCounter = zeros(Int64, (G - 1), 2)
    PTinitializeParameters!(PTparamArray, PTdistArray, J, G, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
    
    for g in 1:G
        PTinitializeCurrentStates!(PTcurrent_a, PTcurrent_n, proposal_a, proposal_n, J, I, g, PTparamArray[g, 1].psi * Delta)
        correctOrder!(PTparamArray[g, 1].psi, PTparamArray[g, 1].mu, PTparamArray[g, 1].ksi, PTcurrent_a, g)
    end
    current_iterationParameters::Array{MHiterationParameters,1} = PTparamArray[:, 1]
    current_iterationDistributions::Array{MHiterationDistributions,1} = PTdistArray[:, 1]
    current_acceptance::Array{Float64, 1} = zeros(Float64, G)


    # Initial variables for the swap part of the algorithm
    if(noSwap)
        const Bswap = 1;
    else
        const Bswap = bswap;
    end

    r1 = 0; r2 = 0; swapIndexG = zeros(Int64, 2); succes = false;
    out_dist_psi = Array(Distribution, J); out_dist_mu = Distribution; out_dist_tau = Distribution;
    out_psi = zeros(Float64, J); out_labda = 0.0; out_ksi = zeros(Float64, J);
    out_mu = zeros(Float64, J); out_tau = 0.0;

    for m in 1:(M-1)
        InfoMessages && (m % 500 == 0 && info(string("Start iteration m = ", m, " of ", M, "\t(", round(m/M * 100, 1) ,"%)")))
        r1 = rand();
        if(r1 < Bswap)
            for g in 1:G
                # update segments
                current_acceptance[g] = 0
                for indexI in 1:I
                    mode1 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode2 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode3 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode4 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    if(succes)
                        current_acceptance[g] += 1
                    end
                end

                # calculate new parameters
                in_s = sum(PTcurrent_a[:, :, g], 2)
                if mode1
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode2
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda2(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu2(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode3
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda3(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu3(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                elseif mode4
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda4(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu4(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                end
                correctOrder!(out_psi, out_mu, out_ksi, PTcurrent_a, g)

                # update parameters
                current_iterationDistributions[g] = MHiterationDistributions(out_dist_psi, out_dist_mu, out_dist_tau)
                current_iterationParameters[g] = MHiterationParameters(out_labda, out_psi, out_mu, out_tau, out_ksi)
            end
        else
#             # Keep old values and distributions
#             PTdistArray[:, m + 1] = current_iterationDistributions[:]
#             PTparamArray[:, m + 1] = current_iterationParameters[:]

            # Propose a swap
            swapIndexG[1] = rand(1:(G-1));
            swapIndexG[2] = swapIndexG[1] + 1;
            swapCounter[swapIndexG[1], 2] += 1
            r2 = rand();

            # Calculate C
            mode1 && (C = calculateC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode2 && (C = calculateC2!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode3 && (C = calculateC3!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode4 && (C = calculateC4!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))

            # Perform swap
            if(r2 < C)
                swapCounter[swapIndexG[1], 1] += 1
                # swap values
                current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[1]]
                current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[1]]
                PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[1]]
                PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[1]]
            end
        end

        if (thinning)
            if (m%10 == 0)
                PTdistArray[:, Int(round(m/10)) + 1] = current_iterationDistributions
                PTparamArray[:, Int(round(m/10)) + 1] = current_iterationParameters
                PTacceptance[:, Int(round(m/10)) + 1] = current_acceptance
            end
        else
            PTdistArray[:, m + 1] = current_iterationDistributions
            PTparamArray[:, m + 1] = current_iterationParameters
            PTacceptance[:, m + 1] = current_acceptance
        end
    end
    PTacceptance = PTacceptance / I
    return PTparamArray, PTdistArray, PTacceptance, swapCounter
end

function parallelTemperingExchange(z::Array{Float64,1}, N::Int64, J::Int64, I::Int64, G::Int64, M::Int64, T::Float64, Delta::Float64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64, noSwap::Bool, bswap::Float64, gamma::Array{Float64,1}, UpdateModusMH::AbstractString, thinning::Bool, InfoMessages::Bool, detailedInfoMessages::Bool)
    InfoMessages && info(string("Simulation start: Parallel tempering (gamma = ", gamma, ", mode = " , UpdateModusMH, ")"))

    # INITIALIZE VARIABLES AND INITIAL STATES OF A AND N
    if (thinning)
        PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, Int(round(Int, M / 10)))
        PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, Int(round(Int, M / 10)))
        PTacceptance::Array{Float64,2} = zeros(Float64, G, Int(round(Int, M / 10)))
        maxiter::Int32 = M - 10
    else
        PTparamArray = Array(MHiterationParameters, G, M)
        PTdistArray = Array(MHiterationDistributions, G, M)
        PTacceptance = zeros(Float64, G, M)
        maxiter = M - 1
    end
    
#     PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, M)
#     PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, M)
#     PTacceptance::Array{Float64,2} = zeros(Float64, G, M)
    PTcurrent_a::Array{Int64,3} = zeros(Int64, J, I, G);
    PTcurrent_n::Array{Int64, 2} = zeros(Int64, I, G);

    # Computation-efficient constants
    const mode1::Bool =  UpdateModusMH == "PTA-TempDensity"
    const mode2::Bool =  UpdateModusMH == "PTB-TempPosterior"
    const mode3::Bool =  UpdateModusMH == "PTC-TempPrior"
    const mode4::Bool =  UpdateModusMH == "PTD-TempTau"
    const zSq::Array{Float64,1} = abs2(z)
    if mode1
        const psi_rate::Array{Float64,1} = (gamma * (beta_0 + T)) .^ -1
        const tau_shape::Array{Float64,1} = gamma * (alpha_1 + (I / 2)) + ((1 - gamma) * (1 - (J / 2)))
    elseif mode2
        const psi_rate = (beta_0 + gamma * T) .^ -1
        const tau_shape = alpha_1 + gamma  * (I / 2)
    elseif mode3
        const psi_rate = (gamma * beta_0 + T) .^ -1
        const tau_shape = (gamma * alpha_1 + (1 - gamma) * (1 - (J/2)) + (I / 2))
    elseif mode4
        const psi_rate = zeros(Float64, G) + (beta_0 + T) .^ -1
        const tau_shape = (gamma * (alpha_1 - 1) + (I / 2) + 1)
    end
    const kappa_times_Eye_J::Array{Float64,2} = kappa * eye(J)
    const lengthA::Int64 = length(PTcurrent_a[:, :, 1])
    tmpInt::Array{Int64,1} = zeros(Int64, 4)
    tmpInt[1] = lengthA
    tmpFloats::Array{Float64,1} = zeros(Float64, 5)

    # Variables for reduction of garbage collection
    proposal_a::Array{Int64, 1} = zeros(Int64, J); proposal_n::Array{Int64, 1} = zeros(Int64, 1);

    # START ALGORITHM
    swapCounter = zeros(Int64, (G - 1), 2)
    PTinitializeParameters!(PTparamArray, PTdistArray, J, G, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
    
    for g in 1:G
        PTinitializeCurrentStates!(PTcurrent_a, PTcurrent_n, proposal_a, proposal_n, J, I, g, PTparamArray[g, 1].psi * Delta)
        correctOrder!(PTparamArray[g, 1].psi, PTparamArray[g, 1].mu, PTparamArray[g, 1].ksi, PTcurrent_a, g)
    end
    current_iterationParameters::Array{MHiterationParameters,1} = PTparamArray[:, 1]
    current_iterationDistributions::Array{MHiterationDistributions,1} = PTdistArray[:, 1]
    current_acceptance::Array{Float64, 1} = zeros(Float64, G)


    # Initial variables for the swap part of the algorithm
    if(noSwap)
        const Bswap = 1.0;
    else
        const Bswap = bswap;
    end
    const BExchange = 0.05; r3 = 0;

    r1 = 0; r2 = 0; swapIndexG = zeros(Int64, 2); succes = false;
    out_dist_psi = Array(Distribution, J); out_dist_mu = Distribution; out_dist_tau = Distribution;
    out_psi = zeros(Float64, J); out_labda = 0.0; out_ksi = zeros(Float64, J);
    out_mu = zeros(Float64, J); out_tau = 0.0;

    for m in 1:(M-1)
        InfoMessages && (m % 500 == 0 && info(string("Start iteration m = ", m, " of ", M, "\t(", round(m/M * 100, 1) ,"%)")))
        r1 = rand();
        if(r1 < Bswap)
            for g in 1:G
                # update segments
                current_acceptance[g] = 0
                for indexI in 1:I
                    mode1 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode2 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode3 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode4 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    if(succes)
                        current_acceptance[g] += 1
                    end
                end

                # calculate new parameters
                in_s = sum(PTcurrent_a[:, :, g], 2)
                if mode1
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode2
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda2(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu2(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode3
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda3(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu3(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                elseif mode4
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda4(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu4(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                end
                correctOrder!(out_psi, out_mu, out_ksi, PTcurrent_a, g)

                # update parameters
                current_iterationDistributions[g] = MHiterationDistributions(out_dist_psi, out_dist_mu, out_dist_tau)
                current_iterationParameters[g] = MHiterationParameters(out_labda, out_psi, out_mu, out_tau, out_ksi)
            end
        else
            r2 = rand();
            r3 = rand();
            swapIndexG[1] = rand(1:(G-1));
            swapIndexG[2] = swapIndexG[1] + 1;
            swapCounter[swapIndexG[1], 2] += 1

            if(r3 < BExchange)
                # Propose an exchange
                # Calculate C
#                 mode1 && (C = calculateC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode1 && (C = calculateExchangeC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode2 && (C = calculateC2!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode3 && (C = calculateExchangeC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
#                 mode3 && (C = calculateC3!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))        
                mode4 && (C = calculateC4!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                
                # Perform exchange
#                 if(r2 < C)
                if(C > 1)
                    swapCounter[swapIndexG[1], 1] += 1
                    # swap values
                    current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[1]]
                    current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[1]]
                    PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[1]]
                    PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[1]]
#                     current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[2]]
#                     current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[2]]
#                     PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[2]]
#                     PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[2]]
                else
#                     current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[1]]
#                     current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[1]]
#                     PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[1]]
#                     PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[1]]
                    current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]]
                    current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]]
                    PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]]
                    PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]]
#                     current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[2]]
#                     current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[2]]
#                     PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[2]]
#                     PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[2]]
                end
            else
                # Propose a swap
                # Calculate C
                mode1 && (C = calculateC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode2 && (C = calculateC2!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode3 && (C = calculateC3!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
                mode4 && (C = calculateC4!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))

                # Perform swap
                if(r2 < C)
                    swapCounter[swapIndexG[1], 1] += 1
                    # swap values
                    current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[1]]
                    current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[1]]
                    PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[1]]
                    PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[1]]
                end
            end
        end

        if (thinning)
            if (m%10 == 0)
                PTdistArray[:, Int(round(m/10)) + 1] = current_iterationDistributions
                PTparamArray[:, Int(round(m/10)) + 1] = current_iterationParameters
                PTacceptance[:, Int(round(m/10)) + 1] = current_acceptance
            end
        else
            PTdistArray[:, m + 1] = current_iterationDistributions
            PTparamArray[:, m + 1] = current_iterationParameters
            PTacceptance[:, m + 1] = current_acceptance
        end
    end
    PTacceptance = PTacceptance / I
    return PTparamArray, PTdistArray, PTacceptance, swapCounter
end

function parallelExchange(z::Array{Float64,1}, N::Int64, J::Int64, I::Int64, G::Int64, M::Int64, T::Float64, Delta::Float64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64, noSwap::Bool, bswap::Float64, gamma::Array{Float64,1}, UpdateModusMH::AbstractString, thinning::Bool, InfoMessages::Bool, detailedInfoMessages::Bool)
    InfoMessages && info(string("Simulation start: Parallel tempering (gamma = ", gamma, ", mode = " , UpdateModusMH, ")"))

    # INITIALIZE VARIABLES AND INITIAL STATES OF A AND N
    if (thinning)
        PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, Int(round(Int, M / 10)))
        PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, Int(round(Int, M / 10)))
        PTacceptance::Array{Float64,2} = zeros(Float64, G, Int(round(Int, M / 10)))
        maxiter::Int32 = M - 10
    else
        PTparamArray = Array(MHiterationParameters, G, M)
        PTdistArray = Array(MHiterationDistributions, G, M)
        PTacceptance = zeros(Float64, G, M)
        maxiter = M - 1
    end
    
#     PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, M)
#     PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, M)
#     PTacceptance::Array{Float64,2} = zeros(Float64, G, M)
    PTcurrent_a::Array{Int64,3} = zeros(Int64, J, I, G);
    PTcurrent_n::Array{Int64, 2} = zeros(Int64, I, G);

    # Computation-efficient constants
    const mode1::Bool =  UpdateModusMH == "PTA-TempDensity"
    const mode2::Bool =  UpdateModusMH == "PTB-TempPosterior"
    const mode3::Bool =  UpdateModusMH == "PTC-TempPrior"
    const mode4::Bool =  UpdateModusMH == "PTD-TempTau"
    const zSq::Array{Float64,1} = abs2(z)
    if mode1
        const psi_rate::Array{Float64,1} = (gamma * (beta_0 + T)) .^ -1
        const tau_shape::Array{Float64,1} = gamma * (alpha_1 + (I / 2)) + ((1 - gamma) * (1 - (J / 2)))
    elseif mode2
        const psi_rate = (beta_0 + gamma * T) .^ -1
        const tau_shape = alpha_1 + gamma  * (I / 2)
    elseif mode3
        const psi_rate = (gamma * beta_0 + T) .^ -1
        const tau_shape = (gamma * alpha_1 + (1 - gamma) * (1 - (J/2)) + (I / 2))
    elseif mode4
        const psi_rate = zeros(Float64, G) + (beta_0 + T) .^ -1
        const tau_shape = (gamma * (alpha_1 - 1) + (I / 2) + 1)
    end
    const kappa_times_Eye_J::Array{Float64,2} = kappa * eye(J)
    const lengthA::Int64 = length(PTcurrent_a[:, :, 1])
    tmpInt::Array{Int64,1} = zeros(Int64, 4)
    tmpInt[1] = lengthA
    tmpFloats::Array{Float64,1} = zeros(Float64, 5)

    # Variables for reduction of garbage collection
    proposal_a::Array{Int64, 1} = zeros(Int64, J); proposal_n::Array{Int64, 1} = zeros(Int64, 1);

    # START ALGORITHM
    swapCounter = zeros(Int64, (G - 1), 2)
    PTinitializeParameters!(PTparamArray, PTdistArray, J, G, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
    
    for g in 1:G
        PTinitializeCurrentStates!(PTcurrent_a, PTcurrent_n, proposal_a, proposal_n, J, I, g, PTparamArray[g, 1].psi * Delta)
        correctOrder!(PTparamArray[g, 1].psi, PTparamArray[g, 1].mu, PTparamArray[g, 1].ksi, PTcurrent_a, g)
    end
    current_iterationParameters::Array{MHiterationParameters,1} = PTparamArray[:, 1]
    current_iterationDistributions::Array{MHiterationDistributions,1} = PTdistArray[:, 1]
    current_acceptance::Array{Float64, 1} = zeros(Float64, G)


    # Initial variables for the swap part of the algorithm
    if(noSwap)
        const Bswap = 1;
    else
        const Bswap = bswap;
    end

    r1 = 0; r2 = 0; swapIndexG = zeros(Int64, 2); succes = false;
    out_dist_psi = Array(Distribution, J); out_dist_mu = Distribution; out_dist_tau = Distribution;
    out_psi = zeros(Float64, J); out_labda = 0.0; out_ksi = zeros(Float64, J);
    out_mu = zeros(Float64, J); out_tau = 0.0;

    for m in 1:(M-1)
        InfoMessages && (m % 500 == 0 && info(string("Start iteration m = ", m, " of ", M, "\t(", round(m/M * 100, 1) ,"%)")))
        r1 = rand();
        if(r1 < Bswap)
            for g in 1:G
                # update segments
                current_acceptance[g] = 0
                for indexI in 1:I
                    mode1 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode2 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode3 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode4 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    if(succes)
                        current_acceptance[g] += 1
                    end
                end

                # calculate new parameters
                in_s = sum(PTcurrent_a[:, :, g], 2)
                if mode1
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode2
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda2(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu2(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode3
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda3(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu3(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                elseif mode4
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda4(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu4(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                end
                correctOrder!(out_psi, out_mu, out_ksi, PTcurrent_a, g)

                # update parameters
                current_iterationDistributions[g] = MHiterationDistributions(out_dist_psi, out_dist_mu, out_dist_tau)
                current_iterationParameters[g] = MHiterationParameters(out_labda, out_psi, out_mu, out_tau, out_ksi)
            end
        else
#             # Keep old values and distributions
#             PTdistArray[:, m + 1] = current_iterationDistributions[:]
#             PTparamArray[:, m + 1] = current_iterationParameters[:]

            # Propose a swap
            swapIndexG[1] = rand(1:(G-1));
            swapIndexG[2] = swapIndexG[1] + 1;
            swapCounter[swapIndexG[1], 2] += 1
            r2 = rand();

            # Calculate C
            mode1 && (C = calculateC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode2 && (C = calculateC2!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode3 && (C = calculateC3!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode4 && (C = calculateC4!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            
            # Perform swap
            if(r2 < C)
                swapCounter[swapIndexG[1], 1] += 1
                # swap values
                current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[2]]
                current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[2]]
                PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[2]]
                PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[2]]
            else
                current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[1]]
                current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[1]]
                PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[1]]
                PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[1]]
            end
        end

        if (thinning)
            if (m%10 == 0)
                PTdistArray[:, Int(round(m/10)) + 1] = current_iterationDistributions
                PTparamArray[:, Int(round(m/10)) + 1] = current_iterationParameters
                PTacceptance[:, Int(round(m/10)) + 1] = current_acceptance
            end
        else
            PTdistArray[:, m + 1] = current_iterationDistributions
            PTparamArray[:, m + 1] = current_iterationParameters
            PTacceptance[:, m + 1] = current_acceptance
        end
    end
    PTacceptance = PTacceptance / I
    return PTparamArray, PTdistArray, PTacceptance, swapCounter
end



function parallelTemperingExchangePU(z::Array{Float64,1}, N::Int64, J::Int64, I::Int64, G::Int64, M::Int64, T::Float64, Delta::Float64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64, noSwap::Bool, bswap::Float64, gamma::Array{Float64,1}, UpdateModusMH::AbstractString, thinning::Bool, InfoMessages::Bool, detailedInfoMessages::Bool)
    InfoMessages && info(string("Simulation start: Parallel tempering (gamma = ", gamma, ", mode = " , UpdateModusMH, ")"))

    # INITIALIZE VARIABLES AND INITIAL STATES OF A AND N
    if (thinning)
        PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, Int(round(Int, M / 10)))
        PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, Int(round(Int, M / 10)))
        PTacceptance::Array{Float64,2} = zeros(Float64, G, Int(round(Int, M / 10)))
        maxiter::Int32 = M - 10
    else
        PTparamArray = Array(MHiterationParameters, G, M)
        PTdistArray = Array(MHiterationDistributions, G, M)
        PTacceptance = zeros(Float64, G, M)
        maxiter = M - 1
    end
    
#     PTparamArray::Array{MHiterationParameters,2} = Array(MHiterationParameters, G, M)
#     PTdistArray::Array{MHiterationDistributions,2}  = Array(MHiterationDistributions, G, M)
#     PTacceptance::Array{Float64,2} = zeros(Float64, G, M)
    PTcurrent_a::Array{Int64,3} = zeros(Int64, J, I, G);
    PTcurrent_n::Array{Int64, 2} = zeros(Int64, I, G);

    # Computation-efficient constants
    const mode1::Bool =  UpdateModusMH == "PTA-TempDensity"
    const mode2::Bool =  UpdateModusMH == "PTB-TempPosterior"
    const mode3::Bool =  UpdateModusMH == "PTC-TempPrior"
    const mode4::Bool =  UpdateModusMH == "PTD-TempTau"
    const zSq::Array{Float64,1} = abs2(z)
    if mode1
        const psi_rate::Array{Float64,1} = (gamma * (beta_0 + T)) .^ -1
        const tau_shape::Array{Float64,1} = gamma * (alpha_1 + (I / 2)) + ((1 - gamma) * (1 - (J / 2)))
    elseif mode2
        const psi_rate = (beta_0 + gamma * T) .^ -1
        const tau_shape = alpha_1 + gamma  * (I / 2)
    elseif mode3
        const psi_rate = (gamma * beta_0 + T) .^ -1
        const tau_shape = (gamma * alpha_1 + (1 - gamma) * (1 - (J/2)) + (I / 2))
    elseif mode4
        const psi_rate = zeros(Float64, G) + (beta_0 + T) .^ -1
        const tau_shape = (gamma * (alpha_1 - 1) + (I / 2) + 1)
    end
    const kappa_times_Eye_J::Array{Float64,2} = kappa * eye(J)
    const lengthA::Int64 = length(PTcurrent_a[:, :, 1])
    tmpInt::Array{Int64,1} = zeros(Int64, 4)
    tmpInt[1] = lengthA
    tmpFloats::Array{Float64,1} = zeros(Float64, 5)

    # Variables for reduction of garbage collection
    proposal_a::Array{Int64, 1} = zeros(Int64, J); proposal_n::Array{Int64, 1} = zeros(Int64, 1);

    # START ALGORITHM
    swapCounter = zeros(Int64, (G - 1), 2)
    PTinitializeParameters!(PTparamArray, PTdistArray, J, G, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
    
    for g in 1:G
        PTinitializeCurrentStates!(PTcurrent_a, PTcurrent_n, proposal_a, proposal_n, J, I, g, PTparamArray[g, 1].psi * Delta)
        correctOrder!(PTparamArray[g, 1].psi, PTparamArray[g, 1].mu, PTparamArray[g, 1].ksi, PTcurrent_a, g)
    end
    current_iterationParameters::Array{MHiterationParameters,1} = PTparamArray[:, 1]
    current_iterationDistributions::Array{MHiterationDistributions,1} = PTdistArray[:, 1]
    current_acceptance::Array{Float64, 1} = zeros(Float64, G)


    # Initial variables for the swap part of the algorithm
    if(noSwap)
        const Bswap = 1;
    else
        const Bswap = bswap;
    end

    r1 = 0; r2 = 0; swapIndexG = zeros(Int64, 2); succes = false;
    out_dist_psi = Array(Distribution, J); out_dist_mu = Distribution; out_dist_tau = Distribution;
    out_psi = zeros(Float64, J); out_labda = 0.0; out_ksi = zeros(Float64, J);
    out_mu = zeros(Float64, J); out_tau = 0.0;

    for m in 1:(M-1)
        InfoMessages && (m % 500 == 0 && info(string("Start iteration m = ", m, " of ", M, "\t(", round(m/M * 100, 1) ,"%)")))
        r1 = rand();
        if(r1 < Bswap)
            for g in 1:G
                # update segments
                current_acceptance[g] = 0
                for indexI in 1:I
                    mode1 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode2 && (succes = singleSegmentUpdate!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode3 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    mode4 && (succes = singleSegmentUpdate3!(proposal_n, proposal_a, PTcurrent_n, PTcurrent_a, J, z[indexI], Delta, indexI, g, current_iterationParameters[g], gamma[g], UpdateModusMH))
                    if(succes)
                        current_acceptance[g] += 1
                    end
                end

                # calculate new parameters
                in_s = sum(PTcurrent_a[:, :, g], 2)
                if mode1
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode2
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda2(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu2(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)
                elseif mode3
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda3(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu3(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                elseif mode4
                    out_dist_psi, out_psi, out_labda = PTupdateDistPsiLabda4(J, in_s, gamma[g], alpha_0, psi_rate, g)
                    out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu = updateDistTauKsiMu4(kappa_times_Eye_J, z, zSq, J, I, g, current_iterationParameters[g].ksi, gamma[g], PTcurrent_n, PTcurrent_a, tau_shape, beta_1, kappa)                
                end
                correctOrder!(out_psi, out_mu, out_ksi, PTcurrent_a, g)

                # update parameters
                current_iterationDistributions[g] = MHiterationDistributions(out_dist_psi, out_dist_mu, out_dist_tau)
                current_iterationParameters[g] = MHiterationParameters(out_labda, out_psi, out_mu, out_tau, out_ksi)
            end
        else
#             # Keep old values and distributions
#             PTdistArray[:, m + 1] = current_iterationDistributions[:]
#             PTparamArray[:, m + 1] = current_iterationParameters[:]

            # Propose a swap
            swapIndexG[1] = rand(1:(G-1));
            swapIndexG[2] = swapIndexG[1] + 1;
            swapCounter[swapIndexG[1], 2] += 1
            r2 = rand();

            # Calculate C
            mode1 && (C = calculateExchangeC!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode2 && (C = calculateC2!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode3 && (C = calculateC3!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            mode4 && (C = calculateC4!(tmpFloats, tmpInt, z, T, J, I, G, swapIndexG, PTcurrent_a, PTcurrent_n, current_iterationParameters[:], current_iterationDistributions[:], gamma, m))
            
            # Perform swap
            if(C < 1)
                swapCounter[swapIndexG[1], 1] += 1
                # swap values
#                 current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[1]]
#                 current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[1]]
#                 PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[1]]
#                 PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[1]]
                current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[2]], current_iterationDistributions[swapIndexG[2]]
                current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[2]], current_iterationParameters[swapIndexG[2]]
                PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[2]], PTcurrent_a[:, :, swapIndexG[2]]
                PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[2]], PTcurrent_n[:, swapIndexG[2]]
            else
                current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[2]] = current_iterationDistributions[swapIndexG[1]], current_iterationDistributions[swapIndexG[1]]
                current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[2]] = current_iterationParameters[swapIndexG[1]], current_iterationParameters[swapIndexG[1]]
                PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[2]] = PTcurrent_a[:, :, swapIndexG[1]], PTcurrent_a[:, :, swapIndexG[1]]
                PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[2]] = PTcurrent_n[:, swapIndexG[1]], PTcurrent_n[:, swapIndexG[1]]
            end
        end

        if (thinning)
            if (m%10 == 0)
                PTdistArray[:, Int(round(m/10)) + 1] = current_iterationDistributions
                PTparamArray[:, Int(round(m/10)) + 1] = current_iterationParameters
                PTacceptance[:, Int(round(m/10)) + 1] = current_acceptance
            end
        else
            PTdistArray[:, m + 1] = current_iterationDistributions
            PTparamArray[:, m + 1] = current_iterationParameters
            PTacceptance[:, m + 1] = current_acceptance
        end
    end
    PTacceptance = PTacceptance / I
    return PTparamArray, PTdistArray, PTacceptance, swapCounter
end


function correctOrder!(psi, mu, ksi, current_a, g)
    v = sortperm(mu)
    mu = mu[v]
    psi = psi[v]
    ksi = ksi[v]
    current_a[:, :, g] = current_a[v, :, g]
    return
end



# function proposeAiNi!(ai::Array{Int, 1}, J::Int64, lambda::Array{Float64,1})
#     ni = 0
#     while (ni <= 0)
#         for j in 1:J
#             ai[j] = rand(Poisson(lambda[j]))
#             ni += ai[j]
#         end
#     end
#     return ni
# end


function PTinitializeParameters!(PTparamArray::Array{MHiterationParameters,2}, PTdistArray::Array{MHiterationDistributions, 2}, J::Int64, G::Int64, alpha_0::Float64, beta_0::Float64, alpha_1::Float64, beta_1::Float64, ksi::Array{Float64,1}, kappa::Float64)
    first_dist_psi = Array(Distribution, J)
    first_dist_psi[:] = Gamma(alpha_0, (1 / beta_0))
    first_dist_tau = Gamma(alpha_1, (1 / beta_1))    
    
    for g in 1:G
        first_psi = Array(Float64, J)    
        for j in 1:J
            first_psi[j] = rand(first_dist_psi[j])
        end
        first_labda = sum(first_psi)
        
        first_tau = rand(first_dist_tau)
        first_dist_mu = MvNormal(ksi, (1 / (first_tau * kappa)) * eye(J))
        first_mu  = rand(first_dist_mu)
        
        PTparamArray[g, 1] = MHiterationParameters(first_labda, first_psi, first_mu, first_tau, ksi)
        PTdistArray[g, 1]  = MHiterationDistributions(first_dist_psi, first_dist_mu, first_dist_tau)        
    end    
    return
end
function PTinitializeCurrentStates!(PTcurrent_a::Array{Int64, 3}, PTcurrent_n::Array{Int64,2}, proposal_a::Array{Int64, 1}, proposal_n::Array{Int64, 1}, J::Int64, I::Int64, g::Int64, psiDelta_g::Array{Float64,1})
    for i in 1:I
        MHproposeNiAi!(proposal_n, proposal_a, J, psiDelta_g)
        PTcurrent_n[i, g] = proposal_n[1]
        PTcurrent_a[:, i, g] = proposal_a
    end
    return
end


# function zMinMuSqDividedSigmaSq(J, z, ai, ni, mu)
#     for j in 1:J
#         ans = ai[j] * mu[j]
#     end
#     ans = (z - ans)^2 / ni
#     return ans
# end


function calculateA(J::Int64, in_z::Float64, Delta::Float64, original_a, original_n, proposal_a, proposal_n, in_param, in_gamma::Float64)
    numer = log(pdf(Normal(dot(proposal_a, in_param.mu), sqrt(proposal_n[1] / in_param.tau)), in_z))
    denom = log(pdf(Normal(dot(original_a, in_param.mu), sqrt(original_n / in_param.tau)), in_z))
    (original_n == -Inf ||  original_n == Inf) && throw(ArgumentError("calculateA(): Fault in generating a: original_n == +-Inf"))
#     (denom == -Inf || denom == Inf) && throw(ArgumentError("calculateA(): denom == +-Inf"))
    original_n == 0  && throw(ArgumentError("calculateA(): Fault in generating a: original_n == 0"))
    (denom == 0) && return 1
#     denom == 0 && throw(ArgumentError("calculateA(): denom == 0"))
    (numer - denom == NaN) && throw(ArgumentError("calculateA(): numer - denom == NaN"))

    D = in_gamma * (numer - denom)

    E = 0
    for j in 1:J
        if(proposal_a[j] != original_a[j])
            E += (proposal_a[j] - original_a[j]) * (log(in_param.psi[j]) + log(Delta))
            if(original_a[j] > proposal_a[j])
                E += sum(log((proposal_a[j] + 1):original_a[j]))
            else
                E -= sum(log((original_a[j] + 1):proposal_a[j]))
            end
        end
    end
    E *= (in_gamma - 1)
    return exp(D + E)
end


function calculateA3(J::Int64, in_z::Float64, Delta::Float64, original_a, original_n, proposal_a, proposal_n, in_param, in_gamma::Float64)
    numer = pdf(Normal(dot(proposal_a, in_param.mu), sqrt(proposal_n[1] / in_param.tau)), in_z)
    denom = pdf(Normal(dot(original_a, in_param.mu), sqrt(original_n / in_param.tau)), in_z)
    denom == 0 && throw(ArgumentError("calculateA3(): denom == 0"))
    D = min(numer / denom, 1)
    return D
end


# update a single element for a single segment for a single g
function singleSegmentUpdate!(proposal_n::Array{Int,1}, proposal_a::Array{Int,1}, PTcurrent_n::Array{Int,2}, PTcurrent_a::Array{Int,3}, J::Int64, in_z::Float64, Delta::Float64, indexI::Int64, indexG::Int64, in_param::MHiterationParameters, in_gamma::Float64, RWmodus::AbstractString)
    succes::Bool = false;
    MHproposeNiAi!(proposal_n, proposal_a, J, in_param.psi * Delta)

    # accept/reject proposal
    A = calculateA(J, in_z, Delta, PTcurrent_a[:, indexI, indexG], PTcurrent_n[indexI, indexG], proposal_a, proposal_n, in_param, in_gamma)

    r3::Float64 = rand();
    # update segments if accept
    if(r3 < A)
        PTcurrent_a[:, indexI, indexG] = proposal_a
        PTcurrent_n[indexI, indexG] = proposal_n[1]
        succes = true;
    end
    return succes
end


function singleSegmentUpdate3!(proposal_n::Array{Int,1}, proposal_a::Array{Int,1}, PTcurrent_n::Array{Int,2}, PTcurrent_a::Array{Int,3}, J::Int64, in_z::Float64, Delta::Float64, indexI::Int64, indexG::Int64, in_param::MHiterationParameters, in_gamma::Float64, RWmodus::AbstractString)
    succes::Bool = false;
    MHproposeNiAi!(proposal_n, proposal_a, J, in_param.psi * Delta)

    # accept/reject proposal
    A = calculateA3(J, in_z, Delta, PTcurrent_a[:, indexI, indexG], PTcurrent_n[indexI, indexG], proposal_a, proposal_n, in_param, in_gamma)

    r3::Float64 = rand();
    # update segments if accept
    if(r3 < A)
        PTcurrent_a[:, indexI, indexG] = proposal_a
        PTcurrent_n[indexI, indexG] = proposal_n[1]
        succes = true;
    end
    return succes
end



function PTupdateDistPsiLabda(J::Int64, in_s, in_gamma::Float64, alpha_0::Float64, psi_rate::Array{Float64,1}, indexG::Int64)
    out_dist_psi = Array(Distribution, J); out_psi = zeros(J);
    for j in 1:J
        out_dist_psi[j] = Gamma((1 + in_gamma * (alpha_0 + in_s[j] - 1)), psi_rate[indexG])
        out_psi[j] = rand(out_dist_psi[j])
    end
    return out_dist_psi, out_psi, sum(out_psi)
end

function updateDistTauKsiMu(kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, J::Int64, I::Int64, g::Int64, in_ksi, in_gamma::Float64, PTcurrent_n, PTcurrent_a, tau_shape::Array{Float64,1}, beta_1::Float64, kappa::Float64)
    P::Array{Float64,2} = kappa_times_Eye_J
    q = kappa * in_ksi
    R = kappa * dot(in_ksi, in_ksi)
    for i in 1:I
        P += (PTcurrent_a[:, i, g] * PTcurrent_a[:, i, g]') / PTcurrent_n[i, g]
        for j in 1:J
            q[j] += (PTcurrent_a[j, i, g] / PTcurrent_n[i, g]) * z[i]
        end
        R += (zSq[i]) / PTcurrent_n[i, g]
    end
    Pinv = inv(P)

    out_ksi = Pinv * q
    out_dist_tau = Gamma(tau_shape[g], (1 / (in_gamma * (beta_1 + sum(R - q' * Pinv * q) / 2))))
    out_tau = rand(out_dist_tau)
    out_dist_mu = MvNormal(Pinv * q, Pinv / (in_gamma * out_tau))
    out_mu = rand(out_dist_mu)

    return out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu
#     info(out_dist_tau)
end


function PTupdateDistPsiLabda2(J::Int64, in_s, in_gamma::Float64, alpha_0::Float64, psi_rate::Array{Float64,1}, indexG::Int64)
    out_dist_psi = Array(Distribution, J); out_psi = zeros(J);
    for j in 1:J
        out_dist_psi[j] = Gamma(alpha_0 + in_gamma * in_s[j], psi_rate[indexG])
        out_psi[j] = rand(out_dist_psi[j])
    end
    return out_dist_psi, out_psi, sum(out_psi)
end
function updateDistTauKsiMu2(kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, J::Int64, I::Int64, g::Int64, in_ksi, in_gamma::Float64, PTcurrent_n, PTcurrent_a, tau_shape::Array{Float64,1}, beta_1::Float64, kappa::Float64)
    P::Array{Float64,2} = kappa_times_Eye_J
    q = kappa * in_ksi
    R = kappa * dot(in_ksi, in_ksi)
    for i in 1:I
        P += in_gamma * (PTcurrent_a[:, i, g] * PTcurrent_a[:, i, g]') / PTcurrent_n[i, g]
        for j in 1:J
            q[j] += in_gamma * (PTcurrent_a[j, i, g] / PTcurrent_n[i, g]) * z[i]
        end
        R += in_gamma * zSq[i] / PTcurrent_n[i, g]
    end
    Pinv = inv(P)

    out_ksi = Pinv * q
    out_dist_tau = Gamma(tau_shape[g], 1 / (beta_1 + sum(R - q' * Pinv * q) / 2))
    out_tau = rand(out_dist_tau)
    out_dist_mu = MvNormal(Pinv * q, Pinv / out_tau)
    out_mu = rand(out_dist_mu)

    return out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu
end



function PTupdateDistPsiLabda3(J::Int64, in_s, in_gamma::Float64, alpha_0::Float64, psi_rate::Array{Float64,1}, indexG::Int64)
    out_dist_psi = Array(Distribution, J); out_psi = zeros(J);
    for j in 1:J
        out_dist_psi[j] = Gamma(in_gamma * (alpha_0 - 1 ) + in_s[j] + 1, psi_rate[indexG])
        out_psi[j] = rand(out_dist_psi[j])
    end
    return out_dist_psi, out_psi, sum(out_psi)
end
function updateDistTauKsiMu3(kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, J::Int64, I::Int64, g::Int64, in_ksi, in_gamma::Float64, PTcurrent_n, PTcurrent_a, tau_shape::Array{Float64,1}, beta_1::Float64, kappa::Float64)
    P::Array{Float64,2} = in_gamma * kappa_times_Eye_J
    q = in_gamma * kappa * in_ksi
    R = in_gamma * kappa * dot(in_ksi, in_ksi)
    for i in 1:I
        P += (PTcurrent_a[:, i, g] * PTcurrent_a[:, i, g]') / PTcurrent_n[i, g]
        for j in 1:J
            q[j] += (PTcurrent_a[j, i, g] / PTcurrent_n[i, g]) * z[i]
        end
        R += zSq[i] / PTcurrent_n[i, g]
    end
    Pinv = inv(P)

    out_ksi = Pinv * q
    out_dist_tau = Gamma(tau_shape[g], 1 / ((in_gamma * beta_1) + sum(R - q' * Pinv * q) / 2))
    out_tau = rand(out_dist_tau)
    out_dist_mu = MvNormal(Pinv * q, Pinv / out_tau)
    out_mu = rand(out_dist_mu)

    return out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu
end


function PTupdateDistPsiLabda4(J::Int64, in_s, in_gamma::Float64, alpha_0::Float64, psi_rate::Array{Float64,1}, indexG::Int64)
    out_dist_psi = Array(Distribution, J); out_psi = zeros(J);
    for j in 1:J
        out_dist_psi[j] = Gamma(alpha_0 + in_s[j], psi_rate[indexG])
        out_psi[j] = rand(out_dist_psi[j])
    end
    return out_dist_psi, out_psi, sum(out_psi)
end
function updateDistTauKsiMu4(kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, J::Int64, I::Int64, g::Int64, in_ksi, in_gamma::Float64, PTcurrent_n, PTcurrent_a, tau_shape::Array{Float64,1}, beta_1::Float64, kappa::Float64)
    P::Array{Float64,2} = in_gamma * kappa_times_Eye_J
    q = kappa * in_ksi
    R = kappa * dot(in_ksi, in_ksi)
    for i in 1:I
        P += (PTcurrent_a[:, i, g] * PTcurrent_a[:, i, g]') / PTcurrent_n[i, g]
        for j in 1:J
            q[j] += (PTcurrent_a[j, i, g] / PTcurrent_n[i, g]) * z[i]
        end
        R += zSq[i] / PTcurrent_n[i, g]
    end
    Pinv = inv(P)

    out_ksi = Pinv * q
    out_dist_tau = Gamma(tau_shape[g], 1 / ((in_gamma * beta_1) + sum(R - q' * Pinv * q) / 2))
    out_tau = rand(out_dist_tau)
    out_dist_mu = MvNormal(Pinv * q, Pinv / out_tau)
    out_mu = rand(out_dist_mu)

    return out_ksi, out_dist_tau, out_tau, out_dist_mu, out_mu
end


# tmpInt[1] : nrElements
# tmpInt[2] : tmpValue
# tmpInt[3] : tmpGreater
function sumLogFaculty!(tmpInt::Array{Int64,1}, input)
    tmpInt[1] = length(input)
    tmpReturn::Float64 = 0
    tmpInt[3] = tmpInt[1] - sum(input .== 0) - sum(input .== 1) # #elements > 1.
    tmpInt[2] = 1

    while tmpInt[3] > 0
        tmpInt[2] += 1
        tmpReturn += tmpInt[3] * log(tmpInt[2])
        tmpInt[3] -= sum(input .== tmpInt[2])
    end
    return tmpReturn
end

# # # # tmpInt[1] : nrElements
# # # # tmpInt[2] : tmpValue
# # # # tmpInt[3] : tmpCounter
# # # # tmpInt[4] : tmpHits
function sumLog!(tmpInt::Array{Int64,1}, input)
    tmpInt[1] = length(input)
    tmpReturn::Float64 = 0
    tmpInt[3] = sum(input .== 0) + sum(input .== 1) # #elements > 1.
    tmpInt[2] = 1

    while tmpInt[3] < tmpInt[1]
        tmpInt[2] += 1
        tmpInt[4] = sum(input .== tmpInt[2])
        tmpReturn += tmpInt[4] * log(tmpInt[2])
        tmpInt[3] += tmpInt[4]
    end
    return tmpReturn
end

function sumLogPDFNormal!(tmpInt::Array{Int64,1}, z::Array{Float64,1}, I::Int64, tau::Float64, mu, indexG::Int64, current_n, current_a)
    sum::Float64 = 0
    for i in 1:I
        sum += (z[i] - dot(current_a[:, i, indexG], mu))^2 / current_n[i, indexG]
    end
    sum *= tau
    sum += I * (log(2) + log(pi) - log(tau))
    sum += sumLog!(tmpInt, current_n[:, indexG])
    sum *= -0.5
end

# # # # # OPTIMIZE : YOU CAN LOSE A LOT OF MEMORY OVERHEAD HERE
#                 for i in 1:I
#                     tmpTerm1 += log(pdf(Normal(dot(PTcurrent_a[:, i, swapGamma], PTparamArray[swapGamma].mu), sqrt(PTcurrent_n[i, swapGamma] / PTparamArray[swapGamma].tau)), z[i]))
#                 end
        #         println("\n")
        #         println(tmpTerm1)
#     tmpInt = zeros(Int64, 3)
#     tmpInt[1] = lengthA
# , lengthA::Int64
# tmpFloats::Array{Float64,1} = zeros(Float64, 5)
# # for x in 1:5
# #     tmpFloats[x] = 0
# # end
# #
# #
#         tmpTerm1 = 0; tmpTerm2 = 0; tmpTerm3 = 0; tmpPrior = 0;
#         tmpTerm1 = sumLogPDFNormal(z, I, PTparamArray[swapGamma].tau, PTparamArray[swapGamma].mu, swapGamma, PTcurrent_n, PTcurrent_a)
#
#         s = reshape(sum(PTcurrent_a, 2), J, G)
#         for j in 1:J
#             tmpTerm2 += s[j, swapGamma] * log(PTparamArray[swapGamma].psi[j]) # * Delta!
#         end
#
#         tmpTerm3 = sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
#
#         for j in 1:J
#             tmpPrior += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
#         end
#         tmpPrior += log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
#         tmpPrior += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
#         logCFrac[swapIndex] = tmpTerm1 + tmpTerm2 - tmpTerm3 - (PTparamArray[swapGamma].labda * T) + tmpPrior

function calculateExchangeC!(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
    logCFrac = zeros(Float64, 2);

    for swapIndex in 1:2
        swapGamma = swapIndexG[swapIndex]
        tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
        tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
        for j in 1:J
            tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
        end
        tmpFloats[2] = -1 * PTparamArray[swapGamma].labda * T
        tmpFloats[3] = sumLogPDFNormal!(tmpInt, z, I, PTparamArray[swapGamma].tau, PTparamArray[swapGamma].mu, swapGamma, PTcurrent_n, PTcurrent_a)
        tmpFloats[4] = 0
        for j in 1:J
            tmpFloats[4] += sum(PTcurrent_a[j, :, swapGamma]) * log(PTparamArray[swapGamma].psi[j]) # * Delta!
        end
        tmpFloats[5] = -1 * sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
        logCFrac[swapIndex] = sum(tmpFloats)
    end
    return exp(logCFrac[1] - logCFrac[2])
end



function calculateC!(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
    logCFrac = zeros(Float64, 2);

    for swapIndex in 1:2
        swapGamma = swapIndexG[swapIndex]
        tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
        tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
        for j in 1:J
            tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
        end
        tmpFloats[2] = -1 * PTparamArray[swapGamma].labda * T
        tmpFloats[3] = sumLogPDFNormal!(tmpInt, z, I, PTparamArray[swapGamma].tau, PTparamArray[swapGamma].mu, swapGamma, PTcurrent_n, PTcurrent_a)
        tmpFloats[4] = 0
        for j in 1:J
            tmpFloats[4] += sum(PTcurrent_a[j, :, swapGamma]) * log(PTparamArray[swapGamma].psi[j]) # * Delta!
        end
        tmpFloats[5] = -1 * sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
        logCFrac[swapIndex] = sum(tmpFloats)
    end
    return exp((gamma[swapIndexG[2]] - gamma[swapIndexG[1]]) * (logCFrac[1] - logCFrac[2]))
end

# function calculateC2(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
#     logCFrac = zeros(Float64, 2);
#     for swapIndex in 1:2
#         swapGamma = swapIndexG[swapIndex]
# #         tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
# #         tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
# #         for j in 1:J
# #             tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
# #         end
#         s = reshape(sum(PTcurrent_a, 2), J, G)
#         for j in 1:J
#             tmpTerm2 += s[j, swapGamma] * log(PTparamArray[swapGamma].psi[j])
#         end
#         for i in 1:I
#             tmpTerm1 += log(pdf(Normal(dot(PTcurrent_a[:, i, swapGamma], PTparamArray[swapGamma].mu), sqrt(PTcurrent_n[i, swapGamma] / PTparamArray[swapGamma].tau)), z[i]))
#         end
#         #         nrElements = length(PTcurrent_a[:, :, swapGamma])
#         #         tmpCounter = 0; tmpHits = 0; tmpValue = 0;
#         #         while(tmpCounter < nrElements)
#         #             tmpHits = sum(PTcurrent_a[:, :, swapGamma] .== tmpValue)
#         #             tmpCounter += tmpHits
#         #             tmpValue += 1
#         #             tmpTerm3 += tmpHits * sum(log(1:tmpValue))
#         #         end
#         tmpTerm3 = sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
#         #
#         for j in 1:J
#             tmpPrior[swapIndex] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
#         end
#         tmpPrior[swapIndex] += log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
#         tmpPrior[swapIndex] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
#         logCFrac[swapIndex] = tmpTerm1 + tmpTerm2 - tmpTerm3 - (PTparamArray[swapGamma].labda * T)
#     end
#     logC = (tmpPrior[1] - tmpPrior[2]) + (gamma[swapIndexG[2]] - gamma[swapIndexG[1]]) * (logCFrac[1] - logCFrac[2])
#     return exp(logC)
# end
function calculateC2!(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
    logCFrac = zeros(Float64, 2);

    for swapIndex in 1:2
        swapGamma = swapIndexG[swapIndex]
#         tmpFloats[1] = 0
#         tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
#         tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
#         for j in 1:J
#             tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
#         end
        tmpFloats[2] = -1 * PTparamArray[swapGamma].labda * T
        tmpFloats[3] = sumLogPDFNormal!(tmpInt, z, I, PTparamArray[swapGamma].tau, PTparamArray[swapGamma].mu, swapGamma, PTcurrent_n, PTcurrent_a)
        tmpFloats[4] = 0
        for j in 1:J
            tmpFloats[4] += sum(PTcurrent_a[j, :, swapGamma]) * log(PTparamArray[swapGamma].psi[j]) # * Delta!
        end
        tmpFloats[5] = -1 * sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
        logCFrac[swapIndex] = (tmpFloats[2] + tmpFloats[3] + tmpFloats[4] + tmpFloats[5])
    end
    return exp((gamma[swapIndexG[2]] - gamma[swapIndexG[1]]) * (logCFrac[1] - logCFrac[2]))
end


function calculateC3!(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
    logCFrac = zeros(Float64, 2);

    for swapIndex in 1:2
        swapGamma = swapIndexG[swapIndex]
        tmpFloats[1] = 0
        tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
        tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
        for j in 1:J
            tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
        end
        logCFrac[swapIndex] = tmpFloats[1]
    end
    return exp((gamma[swapIndexG[2]] - gamma[swapIndexG[1]]) * (logCFrac[1] - logCFrac[2]))
end

function calculateC4!(tmpFloats::Array{Float64,1}, tmpInt::Array{Int64,1}, z::Array{Float64,1}, T::Float64, J::Int64, I::Int64, G::Int64, swapIndexG, PTcurrent_a, PTcurrent_n, PTparamArray, PTdistArray, gamma, m::Int64)
    logCFrac = zeros(Float64, 2);

    for swapIndex in 1:2
        swapGamma = swapIndexG[swapIndex]
#         tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
#         tmpFloats[1] = 0
#         tmpFloats[1] = log(pdf(PTdistArray[swapGamma].dist_mu, PTparamArray[swapGamma].mu))
#         tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
#         for j in 1:J
#             tmpFloats[1] += log(pdf(PTdistArray[swapGamma].dist_psi[j], PTparamArray[swapGamma].psi[j]))
#         end
#         tmpFloats[2] = -1 * PTparamArray[swapGamma].labda * T
#         tmpFloats[3] = sumLogPDFNormal!(tmpInt, z, I, PTparamArray[swapGamma].tau, PTparamArray[swapGamma].mu, swapGamma, PTcurrent_n, PTcurrent_a)
#         tmpFloats[4] = 0
#         for j in 1:J
#             tmpFloats[4] += sum(PTcurrent_a[j, :, swapGamma]) * log(PTparamArray[swapGamma].psi[j]) # * Delta!
#         end
#         tmpFloats[5] = -1 * sumLogFaculty!(tmpInt, PTcurrent_a[:, :, swapGamma])
        logCFrac[swapIndex] = log(pdf(PTdistArray[swapGamma].dist_tau, PTparamArray[swapGamma].tau))
    end
    return exp((gamma[swapIndexG[2]] - gamma[swapIndexG[1]]) * (logCFrac[1] - logCFrac[2]))
end


function PTreport(real_param, real_J, J, acceptance, swapCounter, gamma, UpdateModusMH, M, ESSAC, normdiffMH, normdiffPT, lambdadiffRealMHPT, thinning, dir, plotFilePaths="")
    if (sum(swapCounter[:, 2] .== 0) == 0)
        info("")
        info(string("Simulation results for gamma = ", gamma, " , m = ", M))
        info(string("Total number of swaps: " , sum(swapCounter[:, 1]), " (", round(100 * sum(swapCounter[:, 1]) / sum(swapCounter[:, 2]), 2), "%)"))

        swapRatio = round(100 * swapCounter[:, 1] ./ swapCounter[:, 2], 2)
        for g in 1:length(swapRatio)
            infoStr = string("Swaps between tempered density ", g, " and ", g + 1, ": ")
            infoStr = string(infoStr, swapCounter[g, 1], " of ", swapCounter[g, 2], "\t (", swapRatio[g], "%)")
            info(infoStr)
        end
    else
        info("No swaps commenced")
    end

    (plotFilePaths == "") && return

    dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J)); makeDirectory(dir)
    fnameLogBase = string("PT", "_", UpdateModusMH, "_gamma-", gamma, "_M=", round(Int, M / 1000), "k_", real_param)
    fnameLogExtension = ".log"
    logFilePath = fnameNew(fnameLogBase, fnameLogExtension, dir)

    outFile = open(logFilePath, "w")
    try
        write(outFile, string("Simulation results \t", Dates.today(), "\n"))
        write(outFile, string("\nParallel tempering modus = ", UpdateModusMH))
        write(outFile, string("\ngamma = ", gamma))
        write(outFile, string("\nPT_M = ", M))

        write(outFile, string("\n\nPlotting files\n"))
        for fName in plotFilePaths
            write(outFile, string(fName, "\n"))
        end

        write(outFile, string("\nSwapping results"))
        if (sum(swapCounter[:, 2] .== 0) == 0)
            write(outFile, "\n")
            write(outFile, string("Total number of swaps: " , sum(swapCounter[:, 1]), " (", round(100 * sum(swapCounter[:, 1]) / sum(swapCounter[:, 2]), 2), "%)"))
            for g in 1:length(swapRatio)
                infoStr = string("\nSwaps between tempered density ", g, " and ", g + 1, ": ")
                infoStr = string(infoStr, swapCounter[g, 1], " of ", swapCounter[g, 2], "\t (", swapRatio[g], "%)")
                write(outFile, infoStr)
            end
        end

        write(outFile, string("\n\nLaTeX Code\n"))
        write(outFile, string("\n\nExchange ratio\n"))
        # SWAP RATIO
#         Tempered densities 
        nCols = 3 + (J - 1) - 1
        nCols2 = nCols + 1
        stringLaTeX = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{r\*{$nCols}{c}}\n\\toprule\n\\itshape Method & \\itshape Total swaps & \\itshape Exchange ratio (total) "
        for g in 1:length(swapRatio)
            stringLaTeX = string(stringLaTeX, " & \\itshape Density ", g, " \\& ", g + 1)
        end
        stringLaTeX = string(stringLaTeX, "\\\\ \n \\midrule\n")
        stringLaTeX = string(stringLaTeX, UpdateModusMH)
        stringLaTeX = string(stringLaTeX, " & \$", sum(swapCounter[:, 1]), "\$ & \$", round(100 * sum(swapCounter[:, 1]) / sum(swapCounter[:, 2]), 2), "\$\\%")
#         stringLaTeX = string(stringLaTeX, " & \$", sum(swapCounter[:, 1]), "\$ & \$", sum(swapCounter[:, 2]), "\$ & \$", round(100 * sum(swapCounter[:, 1]) / sum(swapCounter[:, 2]), 2), "\$\\%")
        for g in 1:length(swapRatio)
            stringLaTeX = string(stringLaTeX, " & \$", swapRatio[g], "\$\\%")
        end
        stringLaTeX = string(stringLaTeX, "\\\\ \n \\cmidrule{2-$nCols2}\n")
        stringLaTeX = string(stringLaTeX, "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
        write(outFile, stringLaTeX)
        
        # EFFECTIVE SAMPLE SIZE
        if J < 3
            write(outFile, string("\n\nEffective Sample Size\n"))
            nCols = 3 + 2 * J - 1
            nCols2 = nCols + 1
            stringLaTeX2 = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{r\*{$nCols}{c}}\n\\toprule\n"
            stringLaTeX2 = string(stringLaTeX2, "\\itshape Parameter & \$\\tau\$")
            for j in 1:J
                stringLaTeX2 = string(stringLaTeX2, " & \$\\mu_", j, "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, " & \$\\lambda\$")
            for j in 1:J
                stringLaTeX2 = string(stringLaTeX2, " & \$\\psi_", j, "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\midrule\n")
            stringLaTeX2 = string(stringLaTeX2, UpdateModusMH)

            for j in 1:(1 + J + 1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][1], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n")
            stringLaTeX2 = string(stringLaTeX2, "Gibbs Sampler (solo)")
            for j in 1:(1 + J + 1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][2], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n")
            stringLaTeX2 = string(stringLaTeX2, "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
            write(outFile, stringLaTeX2)
        else
            write(outFile, string("\n\nEffective Sample Size\n"))
            nCols = 1 + 1 + J - 1
            nCols2 = nCols + 1
            stringLaTeX2 = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{r\*{$nCols}{c}}\n\\toprule\n"
            stringLaTeX2 = string(stringLaTeX2, "\\itshape Parameter & \$\\tau\$")
            for j in 1:J
                stringLaTeX2 = string(stringLaTeX2, " & \$\\mu_", j, "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\midrule\n")
            stringLaTeX2 = string(stringLaTeX2, UpdateModusMH)
            for j in 1:(1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][1], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n")
            stringLaTeX2 = string(stringLaTeX2, "Gibbs Sampler (solo)")
            for j in 1:(1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][2], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n\\toprule\n")
            stringLaTeX2 = string(stringLaTeX2, "\\itshape Parameter")
            stringLaTeX2 = string(stringLaTeX2, " & \$\\lambda\$")
            for j in 1:J
                stringLaTeX2 = string(stringLaTeX2, " & \$\\psi_", j, "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\midrule\n")
            stringLaTeX2 = string(stringLaTeX2, UpdateModusMH)
            for j in (1 + J + 1):(1 + J + 1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][1], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n")
            stringLaTeX2 = string(stringLaTeX2, "Gibbs Sampler (solo)")
            for j in (1 + J + 1):(1 + J + 1 + J)
                stringLaTeX2 = string(stringLaTeX2, " & \$", round(ESSAC[j][2], 2) , "\$ ")
            end
            stringLaTeX2 = string(stringLaTeX2, "\\\\ \n \\cmidrule{2-$nCols2}\n")
            stringLaTeX2 = string(stringLaTeX2, "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
            write(outFile, stringLaTeX2)
        end
        
        # NORM DIFFERENCE F
#         info(normdiffMH)
#         info(normdiffPT)
        write(outFile, string("\n\nNorm difference\n"))
        nCols = 2 - 1
        nCols2 = nCols + 1
        stringLaTeX3 = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{r\*{$nCols}{c}}\n\\toprule\n"        
        stringLaTeX3 = string(stringLaTeX3, "\\itshape Method & \\itshape Norm difference")
        stringLaTeX3 = string(stringLaTeX3, "\\\\ \n \\midrule\n")
        stringLaTeX3 = string(stringLaTeX3, UpdateModusMH, " & \$", round(normdiffPT, 4), "\$")
        stringLaTeX3 = string(stringLaTeX3, "\\\\ \n \\cmidrule{2-$nCols2}\n")
        stringLaTeX3 = string(stringLaTeX3, "Gibbs Sampler (solo) & \$", round(normdiffMH, 4), "\$")
        stringLaTeX3 = string(stringLaTeX3, "\\\\ \n \\cmidrule{2-$nCols2}\n")
        stringLaTeX3 = string(stringLaTeX3, "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
        write(outFile, stringLaTeX3)
        
        # DIFFERENCE ESTIMATES OF LAMBDA
#         lambdadiffRealMHPT
        write(outFile, string("\n\Lambda estimates\n"))
        nCols = 3 - 1
        nCols2 = nCols + 1
        stringLaTeX4 = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{r\*{$nCols}{c}}\n\\toprule\n"        
        stringLaTeX4 = string(stringLaTeX4, "\\itshape Method & \\itshape Estimate  & \\itshape Absolute difference ")
        stringLaTeX4 = string(stringLaTeX4, "\\\\ \n \\midrule\n")
        stringLaTeX4 = string(stringLaTeX4, "Gibbs Sampler (solo) & \$", round(lambdadiffRealMHPT[2], 4), "\$", " & \$", round(lambdadiffRealMHPT[3], 4), "\$")
        stringLaTeX4 = string(stringLaTeX4, "\\\\ \n \\cmidrule{2-$nCols2}\n")
        stringLaTeX4 = string(stringLaTeX4, UpdateModusMH, " & \$", round(lambdadiffRealMHPT[4], 4), "\$", " & \$", round(lambdadiffRealMHPT[5], 4), "\$")
        stringLaTeX4 = string(stringLaTeX4, "\\\\ \n \\cmidrule{2-$nCols2}\n")
        stringLaTeX4 = string(stringLaTeX4, "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
        write(outFile, stringLaTeX4)
        
        finally
        close(outFile)
    end

    return
end

function PTsaveResults(z::Array{Float64,1}, PTparamArray::Array{MHiterationParameters,2}, PTdistArray::Array{MHiterationDistributions,2}, acceptance::Array{Float64,2}, swapCounter, gamma::Array{Float64,1}, M::Int64, real_param::MHiterationParameters, real_J::Int64, J::Int64, UpdateModusMH::AbstractString, PTnoSwap::Bool, dir::AbstractString, infoMessages::Bool, infoMessagesVerbose::Bool)
    infoMessages && info("Saving results")
    dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J)); makeDirectory(dir)

    fnameResultsBase = string("PT", "_", UpdateModusMH, "_gamma-", gamma, "_M=", round(Int, M / 1000), "k_", real_param)
    PTnoSwap && (fnameResultsBase = string(fnameResultsBase, "_noSwap"))
    fnameResultsExtension = ".jld"
    savefilePath = fnameNew(fnameResultsBase, fnameResultsExtension, dir)

    save(savefilePath, "Z", z, "Parameters", PTparamArray, "Acceptance", acceptance, "Distributions", PTdistArray, "Swaps", swapCounter, "Gamma", gamma, "M", M)
    infoMessages && infoMessagesVerbose && info("Saving complete")
    return
end

function PTloadResults(gamma::Array{Float64,1}, M::Int64, real_param::MHiterationParameters, real_J::Int64, J::Int64, UpdateModusMH::AbstractString, PTnoSwap::Bool, dir::AbstractString, infoMessages::Bool, infoMessagesVerbose::Bool)
    infoMessages && info("No simulation performed (saveResultsPT is set to false). Loading Parallel Tempering simulation results")

    fnameResultsBase = string("PT", "_", UpdateModusMH, "_gamma-", gamma, "_M=", round(Int, M / 1000), "k_", real_param)
    PTnoSwap && (fnameResultsBase = string(fnameResultsBase, "_noSwap"))
    fnameResultsExtension = ".jld"
    dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J))
    loadedvar = loadFile(fnameResultsBase, fnameResultsExtension, dir)
    return loadedvar["Parameters"], loadedvar["Distributions"], loadedvar["Acceptance"], loadedvar["Swaps"]
end
