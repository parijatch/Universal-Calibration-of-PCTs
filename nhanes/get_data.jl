using Printf

include("configure.jl")

fn = ["DEMO", "BIOPRO", "BPX", "BMX", "DXX", "OHXDEN"]

urlbase = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/YYYY/DataFiles/"

for f in fn
    url = @sprintf("%s%s_%s.xpt", urlbase, f, sfx[year])
    url = replace(url, "YYYY"=>string(year))

    println("Downloading ", url)
    download(url, joinpath(pa, basename(url)))
end
