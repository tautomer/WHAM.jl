module WHAM
using DelimitedFiles
using Statistics: mean, varm
using NumericalIntegration: integrate, cumul_integrate, SimpsonEven

struct WHAMParamaters1D
    β::Float64
    nWindow::Int64
    nBin::Int64
    lowerBound::Float64
    k::Float64
    binWidth::Float64
    windowCenter::Vector{Float64}
    binCenter::Vector{Float64}
end

mutable struct WHAMArrays1D
    pBiased::Matrix{Float64}
    vBiased::Matrix{Float64}
    nPoints::Vector{Float64}
end

struct UIArrays1D
    mean::Vector{Float64}
    var::Vector{Float64}
end
function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T, nBin::Integer, binWidth::T) where T <: AbstractFloat

    if nw != length(windowCenter)
        error("Number of windows and size of window centers mismatch.")
    end
    lowerBound = bound[1]
    tmp = lowerBound - 0.5*binWidth
    binCenter = [i*binWidth + tmp for i in 1:nBin]
    whamParam = WHAMParamaters1D(1.0/temp, nw, nBin, lowerBound, k, binWidth,
        windowCenter, binCenter)
    whamArray = WHAMArrays1D(Matrix{Float64}(undef, nBin, nw),
        Matrix{Float64}(undef, nBin, nw), Vector{Int32}(undef, nw))
    uiArray = UIArrays1D(Vector{Float64}(undef, nw), Vector{Float64}(undef, nw))
    return whamParam, whamArray, uiArray
end

function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T) where T <: AbstractFloat

    println(
    """Warning: no number of bins or bin width provided in the arguments.
    Will use default value of nBin = 100. It's strongly advised to provide
    one of those values.""")
    nBin = 100
    sort!(bound)
    binWidth = (bound[2]-bound[1]) / nBin
    return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
end

function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T; nBin::Integer=nothing) where T <: AbstractFloat
    sort!(bound)
    binWidth = (bound[2]-bound[1]) / nBin
    return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
end

# function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
#     k::T; binWidth::T=nothing) where T <: AbstractFloat
#     sort!(bound)
#     tmp = (bound[2]-bound[1]) / binWidth
#     if ! isinteger(tmp)
#         println(
#         """The number of bins calculated based the bin size and boundaries
#         isn't integer. Setup will use the next integer value based on the
#         computed number of bins.
#         """)
#         nBin = floor(tmp) + 1
#     end
#     binWidth = (bound[2]-bound[1]) / nBin
#     return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
# end
# 
# function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
#     k::T; nBin::Integer=nothing, binWidth::T=nothing) where T <: AbstractFloat
#     sort!(bound)
#     if abs(nBin*binWidth - (bound[2]-bound[1])) > eps()
#         println(
#         """Both number of bins and bin width are provided in the arguments.
#         The total size of n_bin × w_bin isn't consistent with the boundaries.
#         Setup will be based on the number of bins to avoid possible conflicts.
#         """)
#     end
#     binWidth = (bound[2]-bound[1]) / nBin
#     return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
# end

"""
    function biasedDistibution(traj::Vector)
    
compute the biased distibutions and corresponding baised potentials
"""
function biasedDistibution(traj::Vector{T}, iw::Integer,
    param::WHAMParamaters1D, array::WHAMArrays1D) where T <: AbstractFloat

    nb = param.nBin
    lb = param.lowerBound
    bw = param.binWidth
    fc = param.k
    hist = zeros(nb)
    nCollected = 0.0
    println("Processing window number $iw")
    for rc in traj
        dx = (rc-lb) / bw
        if dx < 0 || dx >= nb
            continue
        end
        index = floor(Int16, dx+1.0)
        hist[index] += 1.0
        nCollected += 1.0
    end
    array.pBiased[:, iw] = hist ./ nCollected
    # file = string("hist_", iw, ".txt")
    # open(file, "w") do io
    #     writedlm(io, hist)
    # end
    array.nPoints[iw] = nCollected
    @. array.vBiased[:, iw] = exp(-param.β*fc*(param.binCenter-param.windowCenter[iw])^2)
    return array
end 
 
"""
    function unbias(param::WHAMParamaters1D, array::WHAMArrays1D, tol::Float64=10e-9)

The procedure and notations mostly follow [this Nwchem document](http://www.nwschem-sw.org/images/Nwchem-new-pmf.pdf)
"""
function unbias(param::WHAMParamaters1D, array::WHAMArrays1D,
    tol::Float64=1e-12, maxCycle::Integer=20000000)
    
    interval = maxCycle / 100
    nb = param.nBin + 0
    nw = param.nWindow + 0
    beta = param.β
    println("""Starting unbias procedure.
    Total number of windows: $nw
    Total number of bins: $nb
    Self-consistent criteria: $tol
    Max number of iterations: $maxCycle""")

    pUnbiased = Vector{Float64}(undef, nb)
    fi = ones(nw)
    # to pre compile the function
    unbiasedProbability!(pUnbiased, fi, array)
    computeError(pUnbiased, fi, array)
    eps = tol

    iter = 0
    while eps >= tol
        unbiasedProbability!(pUnbiased, fi, array)

        eps, fi = computeError(pUnbiased, fi, array)
        iter += 1
        if iter % interval == 0
            println("Current number of iterations: $iter ",
            "Current total error: $eps")
        end
    end
    if iter > maxCycle
        error("Failed to converge after $maxCycle cycles")
    else
        println("Converged in $iter cycles")
    end

    pmf = Vector{Float64}(undef, nb)
    @. pmf = - 1.0 / beta * log(pUnbiased)
    pmin = minimum(pmf)
    pmf .-= pmin
    return param.binCenter, pmf
end

function unbiasedProbability!(pUnbiased::Vector{T}, fi::Vector{T},
    array::WHAMArrays1D) where T <: AbstractFloat

    pUnbiased[:] = (array.pBiased * array.nPoints) ./ (array.vBiased *
        (array.nPoints ./ fi))
    index = findall(x->x==0.0, pUnbiased)
    pUnbiased[index] .= 1e-15
    return pUnbiased
end

function computeError(pUnbiased::Vector{T}, fi::Vector{T},
    array::WHAMArrays1D) where T <: AbstractFloat

    fi_old = deepcopy(fi)
    fi[:] = transpose(array.vBiased) * pUnbiased
    eps = sum( (1.0 .- fi./fi_old).^2 )
    return eps, fi
end

function windowStats(x::AbstractVector{T}) where T <: AbstractFloat
    meanval = mean(x)
    sqdev = varm(x, meanval, corrected=false)
    return meanval, sqdev
end

using TimerOutputs

"""
    function integration(array::UIArrays1D, param::WHAMParamaters1D)

Umbrella integration scheme. Reference https://aip.scitation.org/doi/10.1063/1.2052648
"""
function integration(param::WHAMParamaters1D, array::UIArrays1D)
    println("""Performing numerical integration.
    Total number of windows: $(param.nWindow)
    Total number of bins: $(param.nBin)""")
    c1 = 1 / sqrt(2π)
    temp = 1.0 / param.β
    k = 2.0param.k
    function pbiased(σ2, ξbar, x)
        p = c1 / sqrt(σ2) * exp(-0.5*(x-ξbar)^2/σ2) 
        return p
    end
    function unbiasDeriv(σ2, ξbar, ξᵢ, x)
        ∇ = temp * (x-ξbar) / σ2 - k * (x-ξᵢ)
        return ∇
    end
    reset_timer!()
    @timeit "get p" p = pbiased.(array.var, array.mean,
        reshape(param.binCenter, 1, param.nBin))
    @timeit "get der" der = unbiasDeriv.(array.var, array.mean, param.windowCenter,
        reshape(param.binCenter, 1, param.nBin))
    @timeit "normlize p" p ./= sum(p, dims=1)
    @timeit "get deriv" ∂A∂ξ = vec(sum(p .* der, dims=1))
    @timeit "integration" pmf = cumul_integrate(param.binCenter, ∂A∂ξ)
    print_timer()
    pmin = minimum(pmf)
    pmf .-= pmin
    return param.binCenter, pmf
end

function TSTRate(bin::AbstractArray{T}, pmf::AbstractArray{T}, β::T, μ::T
    ) where T <: AbstractFloat

    globalMin = findmin(pmf)[2]
    middle = floor(Int64, length(pmf)/2)
    if globalMin >= middle
        localMin = findmin(pmf[1:middle])[2]
        ts = findmax(pmf[localMin:globalMin])[2] + localMin - 1
        reactant = pmf[ts:end]
        bin = bin[ts:end]
    else
        localMin = findmin(pmf[middle:end])[2] + middle - 1
        ts = findmax(pmf[globalMin:localMin])[2] + globalMin - 1
        reactant = pmf[ts:-1:1]
        bin = bin[ts:-1:1]
    end
    reactant = @. exp(-β*reactant) 
    prefactor = 1.0 / (2π*β*μ)
    qr = integrate(bin, reactant, SimpsonEven())
    return prefactor * reactant[1] / abs(qr)
end

end # module