module PTstartup
using Distributions

const idWorker = myid() - 1
const myGamma = zeros(Float64, 1)
const idMyGamma = zeros(Int64, 1)
# myPathTMPStream::AbstractString

export idWorker, idMyGamma, myGamma
export getMyGamma, setMyGamma, getIndexMyGamma, setIndexMyGamma, setMyPathTMPStream, getMyPathTMPStream

function getIndexMyGamma()::Int64
    return idMyGamma[1]
end

function setIndexMyGamma(newIndex::Int64)
    idMyGamma[1] = newIndex
    return
end

function getMyGamma()::Float64
    return myGamma[1]
end

function setMyGamma(newGamma::Float64)
    myGamma[1] = newGamma
    return
end

function setMyPathTMPStream(pathTMPStream::AbstractString)
    global myPathTMPStream = joinpath(pathTMPStream, string(idWorker))
    # myPathTMPStreams = Array{AbstractString, 1}(G)
end

function getMyPathTMPStream()::AbstractString
    return myPathTMPStream
end

end


# const G = 3
# 
# const J = 4
# const current_tau = zeros(Float64, J)
# const current_psi = zeros(Float64, J)
# const current_mu  = zeros(Float64, J)
# 
# dist_psi = Array(Distribution, J)
# dist_tau = Array(Distribution, J)
# dist_mu = Array(Distribution, J)
# 
# const current_labda = zeros(Float64, 1)
# const current_acceptance = zeros(Float64, 1)
# # current_labda = Float64(0.0)
# # current_acceptance = Float64(0.0)
# info(current_mu)
# info(current_acceptance)
# 
# gamma_temperatures = SharedArray(Float64, G)
# proc_workers = workers()
# info(gamma_temperatures)
# info(proc_workers)
# 
# current_mu[3] = 2.0
# info(current_mu)




    # running_tau = zeros(M);
    # running_psi = zeros(J, M);
    # running_mu = zeros(J, M);
    # running_labda = zeros(M);
    # running_ksi = zeros(J, M);
    # acceptance = zeros(M);


# first_dist_psi = Array(Distribution, J)
# first_dist_psi[:] = Gamma(alpha_0, (1 / beta_0))
# first_dist_tau = Gamma(alpha_1, (1 / beta_1))

# first_psi = Array(Float64, J)
# for j in 1:J
#     first_psi[j] = rand(first_dist_psi[j])
# end
# first_labda = sum(first_psi)
# first_tau = rand(first_dist_tau)

# first_dist_mu = MvNormal(ksi, (1 / (first_tau * kappa)) * eye(J))
# first_mu  = rand(first_dist_mu)

