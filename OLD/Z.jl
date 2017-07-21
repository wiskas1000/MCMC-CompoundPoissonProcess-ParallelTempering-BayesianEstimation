# W. SEWLAL 1383337 # 2015-12-03
"""
`getZ(N::Int64, real_J::Int64, Delta::Float64, real_param::MHiterationParameters, savedZPath::AbstractString, loadZ::Bool)`

Returns a truncated array `z` of a Compound Poisson Process. The returned array z consists of `I` nonzero differences in a sample of a Compound Poisson Process. If `loadZ` is set to true, the function will load `z` from a savefile. Else, it will call `initializeZ` to generate `z`.

Keyword arguments
-----------------

* `N` : Int64 -- Number of segments
* `real_J` : Int64 -- Number of jump types
* `Delta` : Float64 -- The time between observations, fixed for all observations
* `real_param` : MHiterationParameters -- True parameters
"""
function getZ(N::Int64, real_J::Int64, Delta::Float64, real_param::MHiterationParameters, savedZPath::AbstractString, loadZ::Bool)
    fnameZBase = string("Z-",real_param)
    fnameZExtension = ".jld"

    if(loadZ)
        loadedvar = loadFile(fnameZBase, fnameZExtension, savedZPath)
        return loadedvar["Z"]
    else
        newZ = initializeZ(N, real_J, Delta, real_param)
        savefilePath = fnameNew(fnameZBase, fnameZExtension, savedZPath)
        save(savefilePath, "Z", newZ, "param", real_param)
        return newZ
    end
end


"""
`initializeZ(N::Int64, J::Int64, Delta::Float64, real_param::MHiterationParameters)`

Returns a truncated array z of a Compound Poisson Process. The returned array z consists of `I` nonzero differences in a sample of a Compound Poisson Process. Assumed is that the jump size distribution is a weighted mixture of `J` Gaussians with weights `real_param.psi`/`real_param.labda`, means `real_param.mu` and variance 1/`real_param.tau` (equal for all jump types). The total number of jumps is drawn from a Poisson distribution with mean `real_param.labda` * `Delta`.

Keyword arguments
-----------------

* `N` : Int64 -- Number of segments
* `J` : Int64 -- Number of jump types
* `Delta` : Float64 -- The time between observations, fixed for all observations
* `real_param` : MHiterationParameters
"""
function initializeZ(N::Int64, J::Int64, Delta::Float64, real_param::MHiterationParameters)
    # Check size psi, mu
    length(real_param.psi) != J && throw(ArgumentError("Error: real_psi or real_mu is not of length J"))
    length(real_param.mu) != J && throw(ArgumentError("Error: real_psi or real_mu is not of length J"))
    # Check if labda = sum(psi)
    (abs(sum(real_param.psi) - real_param.labda) > 1e15) && throw(ArgumentError("Error: sum(real_psi) does not equal real_labda"))

    # Create z, zero values of z are on the end-side of the vector
    original_z = zeros(Float64, N);
    tmpai = zeros(Int64, J); tmpni::Int64 = 0; nonzeroCounter::Int64 = 0;

    for i in 1:N
        tmpni = 0
        for j in 1:J
            tmpai[j] = rand(Poisson(real_param.psi[j] * Delta))
            tmpni += tmpai[j]
        end
        if(tmpni > 0)
            nonzeroCounter += 1
            original_z[nonzeroCounter] = rand(Normal(dot(tmpai, real_param.mu), sqrt(tmpni / real_param.tau)))
        end
    end
    z = original_z[1:nonzeroCounter]
    return z
end
