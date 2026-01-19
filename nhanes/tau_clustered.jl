using StatsBase

# Extend this method
import Statistics: var


mutable struct KendallTau

    # The focus is the tau correlation between X and Y
    X::Vector{Float64}
    Y::Vector{Float64}

    # Cluster labels
    I::Vector{Int64}

    # Sampling weights
    wt::Vector{Float64}

    # Ranges of values for each cluster
    ix::Vector{UnitRange{Int64}}

    # The estimated tau values for all distinct cluster pairs
    tauij::Matrix{Float64}

    # The average of the between-cluster tau estimates
    taubar::Float64

    # The estimated marginal CDF for X
    Fx::Vector{Float64}

    # The estimated  marginal CDF for Y
    Fy::Vector{Float64}

    # The estimated joint CDF for X and Y, under independence
    Fxy::Vector{Float64}

    # Quantities defined in section 3.1 of Hunsberger et al.
    h1::Vector{Float64}
    h2::Matrix{Float64}
end

"""
    KendallTau(X, Y, I, wt)

Estimate the weighted Kendalls-Tau correlation betwen X and Y, using
weights wt, based on a cluster sample defined by the labels in I.

Hunsberger, S., Long, L., Reese, S., Hong, G., Myles, I., Zerbe, C.,
Chetchotisakd, P. & Shih, J. (2022). Rank correlation inferences for
clustered data with small sample size. Statistica Neerlandica
"""
function KendallTau(X, Y, I, wt)

    # Sort the group label vector
    if !issorted(I)
        ii = sortperm(I)
        X = X[ii]
        Y = Y[ii]
        I = I[ii]
        wt = wt[ii]
    end

    # Find the position ranges for each group
    u = unique(I)
    ix = []
    for v in u
        push!(ix, searchsorted(I, v))
    end

    n = length(X) # number of observations
    m = length(u) # number of groups

    KT = KendallTau(X, Y, I, wt, ix, zeros(m, m), 0.0, zeros(n), zeros(n), zeros(n), zeros(m), zeros(m, m))
    return KT
end

"""
Estimate Kendall's tau using pairs of observations in two independent
clusters.
"""
function get_tauij(X1, X2, Y1, Y2, wt1, wt2)

    @assert length(X1) == length(Y1) == length(wt1)
    @assert length(X2) == length(Y2) == length(wt2)

    concordant = 0
    ntotal = 0
    for i in eachindex(X1)
        for j in eachindex(X2)
            f1 = X1[i] > X2[j]
            f2 = Y1[i] > Y2[j]
            m = wt1[i] * wt2[j]
            concordant += (f1 == f2) ? m : -m
            ntotal += m
        end
    end

    return concordant / ntotal
end

"""
Compute all between-cluster Kendall's tau estimates, and their mean.
"""
function fill_tauij(KT::KendallTau)

    (; X, Y, ix, wt, tauij) = KT

    m = length(ix)
    taubar = 0.0
    n = div(m * (m  - 1), 2)
    for i in 1:m
        ii = ix[i]
        for j in 1:i-1
            jj = ix[j]
            tauij[i, j] = get_tauij(X[ii], X[jj], Y[ii], Y[jj], wt[ii], wt[jj])
            taubar += tauij[i, j]
        end
    end
    taubar /= n

    KT.taubar = taubar
end

"""
Estimate the CDF using a weighted sample.
"""
function wcdf!(x, w, F)

    @assert length(x) == length(w) == length(F)

    ii = sortperm(x)
    xs = x[ii]
    ws = w[ii]
    cw = cumsum(ws)
    tot = sum(ws)

    for i in eachindex(x)
        j = searchsorted(xs, x[i])
        F[i] = mean(cw[j]) / tot
    end
end

"""
Calculate all CDF estimates.
"""
function fill_F(KT::KendallTau)

    (; X, Y, Fx, Fy, Fxy, wt) = KT

    wcdf!(X, wt, Fx)
    wcdf!(Y, wt, Fy)

    # Estimate the joint CDF under independence
    Fxy .= Fx .* Fy
end

function fill_h1(KT::KendallTau)

    (; Fx, Fy, Fxy, h1, ix, taubar) = KT

    for i in eachindex(ix)
        for j in ix[i]
            h1[i] += 1 - 2*Fx[j] - 2*Fy[j] + 4*Fxy[j]
        end
        h1[i] /= length(ix[i])
        h1[i] -= taubar
    end
end

function fill_h2(KT::KendallTau)

    (; Fx, Fy, Fxy, h1, h2, ix, tauij, taubar) = KT

    for i in eachindex(ix)
        for j in 1:i-1
            h2[i, j] = tauij[i, j] - taubar - h1[i] - h1[j]
        end
    end
end

"""
Carry out all calculations needed for estimation and inference.
"""
function fit!(KT::KendallTau)
    fill_tauij(KT)
    fill_F(KT)
    fill_h1(KT)
    fill_h2(KT)
end

"""
Returns the estimated sampling variance of the tau statistic.  If type
is 'adjusted' account for the clustered sample, if type is
'unadjusted' treat the observations as independent.
"""
function var(KT::KendallTau; type=:adjusted)

    (; X, Y, ix, h1, h2) = KT

    n = length(X) # Sample size
    m = length(ix) # number of groups

    if type == :unadjusted
        return 2*(2n+5) / (9*n*(n-1))
    elseif type == :adjusted
        va = 4 * sum(abs2, h1) / m^2
        v0 = 0.0
        m0 = div(m * (m - 1), 2)

        # Sum of squares of the lower triangle of h2.
        for i in 1:m
            for j in 1:i-1
                v0 += h2[i, j]^2
            end
        end
        v0 /= (m0*n)
        va += v0
    else
        error("!!")
    end
end

function test1()

    r = 0.3   # Pearson correlation between X and Y
    n = 1000  # overall sample size
    g = 10    # number of groups
    icc = 0.5 # group ICC
    wt = 1 .+ rand(n) # sampling weights

    I = sample(1:g, n; replace=true)
    u = randn(g)
    E = sqrt(icc)*u[I] + sqrt(1 - icc)*randn(n)
    X = randn(n)
    Y = r*X + sqrt(1-r^2)*E
    KT = KendallTau(X, Y, I, wt)
    fit!(KT)
    return KT
end

#KT = test1()
