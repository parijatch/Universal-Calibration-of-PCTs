"""
Generate latex tables of the results.
"""

using Serialization
using Printf

den_lab = deserialize("den_lab.ser")
dexa_lab = deserialize("dexa_lab.ser")
bmx_lab = deserialize("bmx_lab.ser")

function make_table(rr, rn)
    for j in eachindex(rr[1])
        row = rr[1][j]
        print("$(rn) & ")
        print(@sprintf("%4d & %5.2f & %5.2f & %5.2f & %5.2f &&", row...))
        row = rr[2][j]
        print(@sprintf("%4d & %5.2f & %5.2f & %5.2f & %5.2f \\\\", row...))
        if j == length(rr[1])
            println("\\hline")
        else
            println("")
        end
    end
end

println(raw"\begin{tabular}{lrrrrrrrrrrr}")
println(raw"& \multicolumn{5}{c}{Female} && \multicolumn{5}{c}{Male}\\\cline{2-6}\cline{8-12}")
println(raw"& $n\;$ & $q_{50}$ & $q_{10}$ & $q_{90}$ & Bonf && $n\;$ & $q_{50}$ & $q_{10}$ & $q_{90}$ & Bonf\\\\")

make_table(den_lab, "den/lab")
make_table(bmx_lab, "bmx/lab")
make_table(dexa_lab, "dexa/lab")

println("\\end{tabular}")
