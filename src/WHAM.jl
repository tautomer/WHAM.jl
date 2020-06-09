module WHAM
using StaticArrays

struct WHAMParamaters1D
    temp::Float64
    nWindow::Int16
    nBin::Int16
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

function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T, nBin::Integer, binWidth::T) where T <: AbstractFloat

    if nw != length(windowCenter)
        error("Number of windows and size of window centers mismatch.")
    end
    lowerBound = bound[1]
    tmp = lowerBound - 0.5*binWidth
    binCenter = [i*binWidth + tmp for i in 1:nBin]
    whamParam = WHAMParamaters1D(temp, nw, nBin, lowerBound, k, binWidth,
        windowCenter, binCenter)
    whamArray = WHAMArrays1D(Matrix{Float64}(undef, nBin, nw),
        Matrix{Float64}(undef, nBin, nw), Vector{Int32}(undef, nw))
    return whamParam, whamArray
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

function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T; binWidth::T=nothing) where T <: AbstractFloat
    sort!(bound)
    tmp = (bound[2]-bound[1]) / binWidth
    if ! isinteger(tmp)
        println(
        """The number of bins calculated based the bin size and boundaries
        isn't integer. Setup will use the next integer value based on the
        computed number of bins.
        """)
        nBin = floor(tmp) + 1
    end
    binWidth = (bound[2]-bound[1]) / nBin
    return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
end

function setup(temp::T, nw::Integer, bound::Vector{T}, windowCenter::Vector{T},
    k::T; nBin::Integer=nothing, binWidth::T=nothing) where T <: AbstractFloat
    sort!(bound)
    if abs(nBin*binWidth - (bound[2]-bound[1])) > eps()
        println(
        """Both number of bins and bin width are provided in the arguments.
        The total size of n_bin × w_bin isn't consistent with the boundaries.
        Setup will be based on the number of bins to avoid possible conflicts.
        """)
    end
    binWidth = (bound[2]-bound[1]) / nBin
    return setup(temp, nw, bound, windowCenter, k, nBin, binWidth)
end

"""
    function biasedDistibution(traj::Vector)
    
compute the biased distibutions and corresponding baised potentials
"""
const x = 1
function biasedDistibution(traj::Vector{T}, iw::Integer,
    param::WHAMParamaters1D, array::WHAMArrays1D) where T <: AbstractFloat

    nb = param.nBin
    lb = param.lowerBound
    bw = param.binWidth
    fc = param.k
    hist = zeros(nb)
    nCollected = 0.0
    println("Processing window number $iw")
    for rc in 1:traj
        dx = (rc-lb) / bw
        if dx < 0 || dx >= nb
            continue
        end
        index = Int16(dx + 1)
        hist[index] += 1.0
        nCollected += 1.0
    end
    array.pBiased[:, iw] = hist ./ nCollected
    array.nPoints[iw] = nCollected
    @. array.vBiased[iw] = exp(-param.temp*fc*(param.binCenter-param.windowCenter[iw])^2)
    return array
end

"""
    function unbias(param::WHAMParamaters1D, array::WHAMArrays1D, tol::Float64=10e-9)

The procedure and notations mostly follow [this Nwchem document](http://www.nwchem-sw.org/images/Nwchem-new-pmf.pdf)
"""
function unbias(param::WHAMParamaters1D, array::WHAMArrays1D,
    tol::Float64=10e-9, maxCycle::Integer=10000000)
    
    nb = param.nBin
    nw = param.nWindow
    beta = param.temp
    println("""Starting unbias procedure.
    Total number of windows: $nw
    Total number of bins: $nb
    Self-consistent criteria: $tol
    Max number of iterations: $maxCycle""")

    pUnbiased = MVector{n, Float64}(undef)
    fi = @MVector fill(1.0, nw)
    pBiased = @MMatrix array.pBiased
    vBiased = @MMatrix array.vBiased
    nPoints = @MVector array.nPoints
    eps = tol

    iter = 0
    while eps >= tol && iter <= maxCycle
        for i in 1:nb
            numer = 0.0; denom = 0.0
            denom += nPoints * vBiased[i, :] / fi 
            numer += nPoints * pBiased[i, :]
            numer /= denom
            if numer == 0.0
                numer = 1e-15
            end
            pUnbiased[i] = numer
        end

        eps = 0.0
        for i in 1:nw
            fi_old = fi
            fi_new = sum(vBiased[:, i] .* pUnbiased)
            fi[i] = fi_new
            eps += (1.0 - fi_new/fi_old)^2
        end
        iter += 1
        if iter % 100000 == 0
            println("""Current number of iterations: $iter
            Current total error: $eps""")
        end
    end

    pmf = - 1.0 / beta * log(pUnbiased)
    pmin = minimum(pmf)
    pmf .-= pmin
    return pmf
end

end # module