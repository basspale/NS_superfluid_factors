push!(LOAD_PATH, "./")

using MurcaNoneqNumericalScan
using DelimitedFiles

function FM(xi::Float64)
    return 1.0 + (22020.0*xi^2)/(11513.0*pi^2) + (5670.0*xi^4)/(11513.0*pi^4) + (420.0*xi^6)/(11513.0*pi^6) + (9.0*xi^8)/(11513.0*pi^8)
end

function HM(xi::Float64)
    #=
    Fernandez and Reisenegger (2005), apj 625, 291-306.
    =#
    return (14680.0*xi)/(11513.0*pi^2) + (7560.0*xi^3)/(11513.0*pi^4) + (840*xi^5)/(11513.0*pi^6) + (24.0*xi^7)/(11513.0*pi^8)
end


"""
Neutron branch
"""

# define vN/xi
length1 = 20
length2 = 20
vn_xis = range(0, stop=1/3.0, length=length1)
vn_xis = vcat(vn_xis, range(vn_xis[end]+vn_xis[2]-vn_xis[1], stop=10., length=length2))
vp_xis = range(0, stop=1.0, length=length1)
vp_xis = vcat(vp_xis, range(vp_xis[end]+vp_xis[2]-vp_xis[1], stop=10., length=length2))

# Store vN/xi
open("../output_data/murca_n_vn_over_xi.dat", "w") do f
    writedlm(f, vn_xis, " ")
end

open("../output_data/murca_n_vp_over_xi.dat", "w") do f
    writedlm(f, vp_xis, " ")
end

#logxis = 0:0.1:2
logxis = 1.1:0.1:2.0
n = 15
for logxi in logxis
    @show logxi
    xi = exp10(logxi)

    # Rate
    res_tmp = map(vn_xi->map(vp_xi->Irate_n_SFnp(vn_xi * xi, vp_xi * xi, xi, n) / HM(xi), vp_xis), vn_xis)
    res_tmp = hcat(res_tmp...) # M[i,j], i=neutron gap, j=proton gap
    open("../output_data/Rrate_murca_n_SFnp_nonzeroT_logxi_$(logxi).dat", "w") do f
        writedlm(f, res_tmp, " ")
    end

    # Emissivity
    res_tmp = map(vn_xi->map(vp_xi->Iemis_n_SFnp(vn_xi * xi, vp_xi * xi, xi, n) / FM(xi), vp_xis), vn_xis)
    res_tmp = hcat(res_tmp...) # M[i,j], i=neutron gap, j=proton gap
    open("../output_data/Remis_murca_n_SFnp_nonzeroT_logxi_$(logxi).dat", "w") do f
        writedlm(f, res_tmp, " ")
    end
end
