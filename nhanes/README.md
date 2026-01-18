# NHANES

This directory contains all code needed to reproduce the NHANES data analysis, illustrating the use
of the Pareto Combination Test to assess the independence between multivariate health phenotypes
using NHANES data.

Set the year and data directory in the file **configure.jl**.  Then run **get_data.jl** to download the files,
and run **corr_proj.jl** to run the full analysis (this will take approximately one hour).  Finally, run the
**tables.jl** script to generate the latex output.

The raw data files and documentation are available [here](https://wwwn.cdc.gov/nchs/nhanes).

