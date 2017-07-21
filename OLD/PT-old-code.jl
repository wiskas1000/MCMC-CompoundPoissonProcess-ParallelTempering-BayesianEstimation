#     if mean(MH_mu[1, :]) < mean(MH_mu[2, :])
#         MH_mu[1, :], MH_mu[2, :] = MH_mu[2, :], MH_mu[1, :]
#         MH_psi[1, :], MH_psi[2, :] = MH_psi[2, :], MH_psi[1, :]
#     end


# # # # # # # # # # # # # # #  PT
#     PTparamArray::Array{Array{MHiterationParameters, 1}, 1} = Array(Array{MHiterationParameters, 1}, G)
#     PTdistArray::Array{Array{MHiterationDistributions, 1}, 1} = Array(Array{MHiterationDistributions, 1}, G)
#     PTacceptance::Array{Array{Float64,1},1} = Array(Array{Float64,1}, G)
#     PTcurrent_a::Array{Array{Int64, 2},1} = Array(Array{Int64, 2}, G)
#     PTcurrent_n::Array{Array{Int64, 1},1} = Array(Array{Int64, 1}, G)
#     for g in 1:G
#         PTparamArray[g] = Array(MHiterationParameters, M)
#         PTdistArray[g]  = Array(MHiterationDistributions, M)
#         PTacceptance[g] = zeros(Float64, M)
#         PTcurrent_a[g] = zeros(Int64, J, I)
#         PTcurrent_n[g] = zeros(Int64, I)
#     end
#
#     PTinitializeParameters!(PTparamArray, PTdistArray, J, G, alpha_0, beta_0, alpha_1, beta_1, ksi, kappa)
#     PTinitializeCurrentStates!(PTcurrent_a, PTcurrent_n, proposal_a, proposal_n, J, I, G, PTparamArray[1][1].psi * Delta)
#
#
# PTMHStep!(PTparamArray[g], PTdistArray[g], PTacceptance[g], PTcurrent_n[g], PTcurrent_a[g], gamma[g], proposal_n, proposal_a, m, I, J, z, zSq, Delta, kappa_times_Eye_J, alpha_0, beta_0_plus_T, alpha_1_plus_half_I, beta_1, kappa, UpdateModusMH)
#




# # # # # # # # # # # mh
#         updateq!(q, current_a[:, i], current_n[i], z[i], J)

# function updateq!(q::Array{Float64,1}, in_a::Array{Int64, 1}, in_n::Int64, in_z::Float64, J::Int64)
#     for j in 1:J
#         q[j] += (in_z / in_n) * in_a[j]
#     end
# end


#         for j in 1:J
#             for k in j:J
#                 P[j, k] += current_a[j, i] * current_a[k, i] / current_n[i]
#                 (j != k) && (P[k, j] = P[j, k])
#             end
#         end
#         updateP!(P, current_a[:, i], current_n[i], J)
# function updateP!(P::Array{Float64,2}, in_a::Array{Int64,1}, in_n::Int64, J::Int64)
#     for j in 1:J
#         for k in j:J
#             P[j, k] += in_a[j] * in_a[k] / in_n
#             (j != k) && (P[k, j] = P[j, k])
#         end
#     end
#     return
# end


#     out_labda::Array{Float64, 1} = zeros(Float64, 1); out_psi::Array{Float64, 1} = zeros(Float64, J); out_mu::Array{Float64, 1} = zeros(Float64, J); out_tau::Array{Float64, 1} = zeros(Float64, 1); out_ksi::Array{Float64, 1} = zeros(Float64, J);
#         MHupdatePsiLabda!(out_labda, out_psi, J, alpha_0, beta_0_plus_T, current_a)
#         MHupdateTauKsiMu!(out_mu, out_tau, out_ksi, kappa_times_Eye_J, z, zSq, I, alpha_1_plus_half_I, beta_1, kappa, current_n, current_a, MHparamArray[m].ksi)
#         MHparamArray[m + 1] = MHiterationParameters(out_labda[1], out_psi[:], out_mu[:], out_tau[1], out_ksi[:])

















# """
# `MHupdatePsiLabda(J, T, alpha_0, beta_0, current_a)`
#
# Draw psi and labda based on the latest state of a.
#
# Keyword arguments
# -----------------
#
# * `J` : Int64 -- Number of jump types
# * `T` :
# * `alpha_0` : Float64
# * `beta_0` : Float64
# * `current_a` :
# """
# function MHupdatePsiLabda!(labda::Array{Float64,1}, psi::Array{Float64,1}, J::Int64, alpha_0::Float64, beta_0_plus_T::Float64, current_a::Array{Int64, 2})
#     const s = sum(current_a, 2)
#
#     for j in 1:J
#         psi[j] = rand(Gamma((alpha_0 + s[j]), (1 / beta_0_plus_T)))
#     end
#     labda[1] = sum(psi)
#     return
# end
# """
# `MHupdateTauKsiMu(J, z, I, alpha_1, beta_1, kappa, current_n, current_a, in_ksi)`
#
# Draw tau, ksi and mu based on the latest state of a.
#
# Keyword arguments
# -----------------
#
# * `J` : Int64 -- Number of jump types
# * `z` : Array{Float64,1} -- `I`-element array differences in a sample of a Compound Poisson Process
# * `I` : Int64 -- Number of nonzero segments
# * `beta_1` : Float64
# * `kappa` : Float64
# * `current_n` :
# * `current_a` :
# * `in_ksi` :
# """
# function MHupdateTauKsiMu!(mu::Array{Float64,1}, tau::Array{Float64,1}, ksi::Array{Float64,1}, kappa_times_Eye_J::Array{Float64,2}, z::Array{Float64,1}, zSq::Array{Float64,1}, I, alpha_1_plus_half_I::Float64, beta_1, kappa, current_n, current_a, in_ksi)
#     P::Array{Float64,2} = kappa_times_Eye_J
#     for i in 1:I
#         P += (current_a[:, i] * current_a[:, i]') / current_n[i]
#     end
#     q = kappa * in_ksi + current_a ./ current_n * z
#     R = kappa * dot(in_ksi, in_ksi) + (1 ./ current_n) * zSq
#     Pinv = inv(P)
#
#     tau[1] = rand(Gamma((alpha_1_plus_half_I), (1 / (beta_1 + sum(R - q' * Pinv * q) / 2))))
#     ksi[:] = Pinv * q
#     mu[:] = rand(MvNormal(Pinv * q, Pinv / tau[1]))
#     return
# end


# running_psi   = zeros(J, G, M);
# running_mu    = zeros(J, G, M);
# running_ksi   = zeros(J, G, M);
# running_tau   = zeros(G, M);
# running_labda = zeros(G, M);
# dist_psi      = Array(Distribution, J, G, M);
# dist_mu       = Array(Distribution, G, M);
# dist_tau      = Array(Distribution, G, M);
# acceptance    = zeros(G, M);
#
# running_psi[:, 1, 1] = rand(Gamma(alpha_0, (1 / beta_0)), J);
# running_tau[1, 1] = rand(Gamma(alpha_1, (1 / beta_1)));
# running_mu[:, 1, 1] = rand(MvNormal(ksi, (1 / (running_tau[1, 1] * kappa)) * eye(J)));
# running_labda[1, 1] = sum(running_psi[:, 1, 1]);
#
# for g in 2:G
#     running_psi[:, g, 1] = running_psi[:, 1, 1]
#     running_mu[:, g, 1]  = running_mu[:, 1, 1]
#     running_tau[g, 1]    = running_tau[1, 1]
#     running_labda[g, 1]  = running_labda[1, 1]
# end
# for g in 1:G
#     dist_tau[g, 1] = Gamma(alpha_1, (1 / beta_1))
#     dist_mu[g, 1] = MvNormal(ksi, (1 / (running_tau[1, 1] * kappa)) * eye(J))
#     for j in 1:J
#         dist_psi[j, g, 1] = Gamma(alpha_0, (1 / beta_0))
#     end
# end


# Create initial states
# current_a = zeros(Int64, J, N, G); ai = zeros(Int64, J);
# for i in 1:N
#     tmpni = 0
#     while (tmpni <= 0)
#         tmpni = 0
#         for j in 1:J
#             ai[j] = rand(Poisson(((running_psi[j, 1, 1] * Delta) / running_labda[1, 1])))
#             tmpni += ai[j]
#         end
#     end
#     for g in 1:G
#         current_a[:, i, g] = ai
#     end
# end
#
# current_n = zeros(N, G);
# for g in 1:G
#     current_n[:, g] = sum(current_a[:, :, g], 1);
# end


# # pdfPlotAllPT(false, running_tau, running_labda, running_psi, running_mu, acceptance, N, figurePath)
#
#
# running_psi
# running_mu
# running_tau
# running_labda
# #

# indexI = rand(1:N);
#
# if(lambertW)
#         out_labda = out_labda + lambertw(-1 * exp(-1 * out_labda) * out_labda)
#         out_psi = out_psi / out_labda
#     end
#
#
#         original_a = current_a[:, indexI, indexG]
#         original_n = current_n[indexI, indexG]
#         D1 = 0.5 * log(original_n) - log(proposal_n)
#         D2 = 0.5 * in_tau^2 * (((z[indexI] - dot(original_a, in_mu)) ^ 2) / (original_n ^ 2)) - (((z[indexI] - dot(proposal_a, in_mu)) ^ 2) / (proposal_n ^ 2))
#         D = in_gamma * (D1 + D2)
#         E = 0
#         for j in 1:J
#             if(proposal_a[j] != original_a[j])
#                 E += (proposal_a[j] - original_a[j]) * (log(in_psi[j]) + log(Delta))
#                 if(original_a[j] > proposal_a[j])
#                     E += sum(log((proposal_a[j] + 1):original_a[j]))
#                 else
#                     E -= sum(log((original_a[j] + 1):proposal_a[j]))
#                 end
#             end
#         end
#         E = E * (in_gamma - 1)
#         A2 = exp(D + E)
#         toc()
#
#         tic()
#         numer = pdf(Normal(dot(proposal_a, in_mu), sqrt(proposal_n / in_tau)), z[indexI])
#         denom = pdf(Normal(dot(current_a[:, indexI, indexG], in_mu), sqrt(current_n[indexI, indexG] / in_tau)), z[indexI])
#         numer2 = 1
#         denom2 = 1
#         for j in 1:J
#             numer2 *= ((in_psi[j] * Delta) ^ (proposal_a[j])) / factorial(proposal_a[j])
#             denom2 *= ((in_psi[j] * Delta) ^ (current_a[j, indexI, indexG])) / factorial(current_a[j, indexI, indexG])
#         end
#         if((denom != 0) && (denom2 != 0))
#             A = ((numer / denom) ^ in_gamma) * ((numer2 / denom2) ^ (in_gamma - 1))
#         end
#         toc()
#         testa1 = exp(D)
#         testa2 = ((numer / denom) ^ in_gamma)
#         testb1 = exp(E)
#         testb2 = ((numer2 / denom2) ^ (in_gamma - 1))
#         println("test D results:", testa1 - testa2)
#         println("test E results:", testb1 - testb2)
# #         println(testb1)
# #         println(testb2)  testb1 = testb2 !!! TRUE
#
# #         println(A)
# #         println(A2)
#         # numer can be denom!
#
#
#
# # update a single element for a single segment for a single g
# function singleSegmentUpdate!(succes::Bool, current_n::Array{Int64,2}, current_a::Array{Int64,3}, J::Int64, I::Int64, G::Int64, z, Delta::Float64, indexI::Int64, indexG::Int64, in_param::MHiterationParameters, in_gamma::Float64, RWmodus::AbstractString, infoMessages::Bool)
#     succes = false;
#     if(RWmodus == "RW 1 segment MH")
#         infoMessages && info("Update segment: proposal RW 1 segment MH")
#         proposal_n, proposal_a = proposeAiNi(J, in_param.psi * Delta)
#
#         # accept/reject proposal
#         A = calculateA(J, z[indexI], Delta, current_a[:, indexI, indexG], current_n[indexI, indexG], proposal_a, proposal_n, in_param, in_gamma)
#         denom = 1;
#     elseif(RWmodus == "RW 1 segment MH-orig")
#         infoMessages && info("Update segment: proposal RW 1 segment MH-orig")
#         proposal_n, proposal_a = proposeAiNi(J, in_param.psi * Delta)
#
#         # accept/reject proposal
#         numer = pdf(Normal(dot(proposal_a, in_param.mu), sqrt(proposal_n / in_param.tau)), z[indexI])
#         denom = pdf(Normal(dot(current_a[:, indexI, indexG], in_param.mu), sqrt(current_n[indexI, indexG] / in_param.tau)), z[indexI])
#         (denom != 0) && (A = (numer / denom) ^ in_gamma)
#     elseif(RWmodus == "RW 1 element")
#         indexJ = rand(1:J); updatePlus = rand(Bool);
#
#         proposal_a = current_a[:, indexI, indexG]
#         proposal_n = current_n[indexI, indexG]
#         if(updatePlus)
#             proposal_a[indexJ] += 1
#             proposal_n += 1
#         elseif((current_a[indexJ, indexI, indexG] > 0) && (sum(current_a[:, indexI, indexG]) > 1))
#                 proposal_a[indexJ] -= 1
#                 proposal_n -= 1
#         end
#
#         # accept/reject proposal
#         numer = pdf(Normal(dot(proposal_a, in_param.mu), sqrt(proposal_n / in_param.tau)), z[indexI])
#         denom = pdf(Normal(dot(current_a[:, indexI, indexG], in_param.mu), sqrt(current_n[indexI, indexG] / in_param.tau)), z[indexI])
#         for j in 1:J
#             numer *= ((in_param.psi[j] * Delta) ^ (proposal_a[j])) / factorial(proposal_a[j])
#             denom *= ((in_param.psi[j] * Delta) ^ (current_a[j, indexI, indexG])) / factorial(current_a[j, indexI, indexG])
#         end
#         (denom != 0) && (A = (numer / denom) ^ in_gamma)
#     end
#
#     denom == 0 && throw(ArgumentError("denom == 0"))
#     r3 = rand();
#     # update segments if accept
#     if(r3 < A)
#         infoMessages && info("Proposal accept: segment update")
#         succes = true;
#         current_a[:, indexI, indexG] = proposal_a
#         current_n[indexI, indexG] = proposal_n
#     end
#     return succes
# end
