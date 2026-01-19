using StatsBase
using LinearAlgebra
using Random
using Distributions
using Serialization

include("configure.jl")
include("prep_nhanes.jl")
include("tau_clustered.jl")

"""
Test the null hypothesis of indepdendence based on the data (X, Y),
which are a weighted, clustered sample.

The approach uses nproj random projections, calculating the tau
correlation and p-value for the data following each projection, then
combining the p-values using the Pareto combination test.
"""
function pareto_independence(X, Y, strat, wt; nproj=100)

    @assert size(X, 1) == size(Y, 1) == length(wt)

    n, p = size(X)
    q = size(Y, 2)

    pvals = zeros(nproj)
    for i in 1:nproj

        # Random projections
        XX = X * randn(p)
        YY = Y * randn(q)

        # Estimate the tau correlation and its standard error
        kt = KendallTau(XX, YY, strat, wt)
        fit!(kt)
        r0 = kt.taubar
        va = var(kt, type=:adjusted)
        se = sqrt(va)

        # Calculate a p-value using a normal approximation for the Z-score
        pvals[i] = 2*cdf(Normal(), -abs(r0)/se)
    end

    return pvals
end

function corrproj(X, Y, strat, wt; nrep=1000)
    f = 1 # Proportion of total sample to retain
    n = size(X, 1) # Total sample size
    m = n # Current sample size
    rslt = [] # All results

    # Repeatedly reduce the sample size by 20%
    for k in 1:10
        ap = zeros(nrep) # Pareto aggregated p-values
        bonf = zeros(nrep) # Bonferroni aggregated p-values

        # Repeat to get the variation over random subsamples and random projections
        for j in 1:nrep
            XX, YY, stratx, wtx = X, Y, strat, wt
            if f < 1
                # Obtain a random subsample of the desired size
                m = Int(ceil(f * n))
                ii = sample(1:n, m)
                XX, YY, stratx, wtx = X[ii, :], Y[ii, :], stratx[ii], wtx[ii]
            end

            # Test independence
            pv = pareto_independence(XX, YY, stratx, wtx)

            # Harmonic mean p-value
            ap[j] = 1 / mean(1 ./ (1e-5 .+ pv))

            # Bonferroni p-values
            bonf[j] = clamp(minimum(pv) * length(pv), 0, 1)
        end

        # Save all results for this sample size
        push!(rslt, [m, median(ap), quantile(ap, 0.1), quantile(ap, 0.9), median(bonf)])
        f *= 0.8
    end

    return rslt
end

"""
Assess independence betwen the variables in v1 and v2 at a sequence of
sample sizes, then store the results in the file named 'fn'.ser.
Results are stratified by sex.
"""
function main(v1, v2, fn)

    # Repeat the analysis for each sex
    # RIAGENDR 1=male, 2=female
    rr = []
    for riagendr in [1, 2]
        dat0 = filter(r->r[:RIAGENDR] == riagendr, dat)
        X1 = dat0[:, v1] |> Matrix
        X2 = dat0[:, v2] |> Matrix

        # Create a stratum indicator by combining the survey stratum and PSU
        strat = 2*dat0[:, :SDMVSTRA] + dat0[:, :SDMVPSU]

        # The weight variable
        wt = dat0[:, :WTINT2YR]
        r1 = corrproj(X1, X2, strat, wt)
        push!(rr, r1)
    end

    rx = (male=rr[1], female=rr[2])

    serialize("$(fn).ser", rx)
end

# Test independence between dental and lab variables
main(denvars, labvars, "den_lab")

# Test independence betwen body measurement and lab variables
main(bmxvars, labvars, "bmx_lab")

# Test independence between body composition and lab variables
main(dexavars, labvars, "dexa_lab")
