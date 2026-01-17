using DataFrames
using CSV
using LinearAlgebra
using ReadStatTables
using Statistics

"""
Before running the script, download data files from here:

https://wwwn.cdc.gov/nchs/nhanes

Set the path variable 'pa' to point to the location of the data files.
"""

yr = split(basename(pa), "-")[1]
sfx = Dict("2013"=>"H", "2015"=>"I", "2017"=>"J")[yr]
wave = basename(pa)

function zscore(x)
    return (x .- mean(x)) ./ std(x)
end

"""
Use rank-revealing QR to identify a subset of columns that are not strongly collinear.
"""
function drop_redundant(dm; tol=1e-5)
    dm = dm[completecases(dm), :]
    dd = select(dm, Not(:SEQN))
    dz = Matrix(dd)
    qq = qr(dz, ColumnNorm())
    r = abs.(diag(qq.R))
    ii = qq.jpvt
    i1 = ii[r .> tol]
    i0 = ii[r .<= tol]
    na = names(dd)
    println("Dropping: ", na[i0])
    dd = dd[:, na[i1]]
    dd[:, :SEQN] = dm[:, :SEQN]
    dd = dd[completecases(dd), :]
    return dd
end

# Read data files
demog = readstat(joinpath(pa, "DEMO_$(sfx).XPT")) |> DataFrame
bio = readstat(joinpath(pa, "BIOPRO_$(sfx).XPT")) |> DataFrame
bpx = readstat(joinpath(pa, "BPX_$(sfx).XPT")) |> DataFrame
bmx = readstat(joinpath(pa, "BMX_$(sfx).XPT")) |> DataFrame
dxx = readstat(joinpath(pa, "DXX_$(sfx).XPT")) |> DataFrame
dxx = readstat(joinpath(pa, "DXX_$(sfx).XPT")) |> DataFrame
dden = readstat(joinpath(pa, "OHXDEN_$(sfx).XPT")) |> DataFrame

# Select a few demographic variables
demog = select(demog, [:SEQN, :RIAGENDR, :RIDAGEYR, :SDMVSTRA, :SDMVPSU, :WTINT2YR])

# Select a few body dimension variables
bmx = select(bmx, [:SEQN, :BMXWT, :BMXHT, :BMXARMC, :BMXWAIST])

# Select a few DEXA variables
dxx = select(dxx, [:SEQN, :DXXTRFAT, :DXDTOLE, :DXDTOFAT, :DXDTOBMC])

# The response is systolic blood pressure.
bpx = select(bpx, [:SEQN, :BPXSY1, :BPXSY2, :BPXDI1, :BPXDI2])

# Generate dental variables
ctc = [x for x in names(dden) if endswith(x, "CTC")]
dden = select(dden, vcat("SEQN", ctc))
xden = Matrix(select(dden, ctc))
dden[:, :dentZ] = sum(xden .== "Z"; dims=2)[:]
dden[:, :dentE] = sum(xden .== "E"; dims=2)[:]
dden[:, :dentP] = sum(xden .== "P"; dims=2)[:]
dden[:, :dentF] = sum(xden .== "F"; dims=2)[:]
dden = select(dden, [:SEQN, :dentZ, :dentE, :dentP, :dentF])

# Merge
dat = leftjoin(demog, bio, on=:SEQN)
dat = leftjoin(dat, bpx, on=:SEQN)
dat = leftjoin(dat, bmx, on=:SEQN)
dat = leftjoin(dat, dxx, on=:SEQN)
dat = leftjoin(dat, dden, on=:SEQN)
dat = dat[completecases(dat), :]
dat = filter(r->30<=r[:RIDAGEYR]<=50, dat)
dat = drop_redundant(dat)

# Z-score most variables
for v in names(dat)
    if !(v in ["SEQN", "SDMVSTRA", "SDMVPSU", "RIAGENDR", "WTINT2YR"])
        dat[!, v] = zscore(dat[:, v])
    end
end

# Blocks of variables for correlation testing
labvars = [x for x in names(dat) if startswith(x, "LB")]
bmxvars = ["BMXWT", "BMXHT", "BMXARMC", "BMXWAIST"]
demogvars = ["RIAGENDR", "RIDAGEYR"]
depvars = ["deppc1", "deppc2", "deppc3"]
dexavars = ["DXXTRFAT", "DXDTOLE", "DXDTOFAT", "DXDTOBMC"]
denvars = ["dentZ", "dentE", "dentP", "dentF"]

# This variable has a strange distribution
labvars = labvars[labvars .!= "LBDSATLC"]

# Generate lists of variable names for each block
labvars = intersect(labvars, names(dat))
bmxvars = intersect(bmxvars, names(dat))
demogvars = intersect(demogvars, names(dat))
depvars = intersect(depvars, names(dat))
dexavars = intersect(dexavars, names(dat))
denvars = intersect(denvars, names(dat))

"""
Identify the variable in 'va' to drop so as to maximize the minimum singular value.
"""
function svd_drop(va, dat)
    sl = []
    for v in va
        # Drop one variable
        vad = va[va .!= v]
        xx = Matrix(dat[:, vad])
        _,ss,_ = svd(xx)
        push!(sl, minimum(ss))
    end

    _, ii = findmax(sl)
    va = va[va .!= va[ii]]
    return va
end

# Drop nearly redundant lab variables
labvars2 = svd_drop(labvars, dat)
labvars3 = svd_drop(labvars2, dat)
labvars4 = svd_drop(labvars3, dat)
labvars5 = svd_drop(labvars4, dat)
labvars = labvars5
