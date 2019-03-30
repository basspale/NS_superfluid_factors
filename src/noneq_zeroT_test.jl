push!(LOAD_PATH, "./")

using MurcaNoneqZeroTemp
using DelimitedFiles

#=
Fernandez and Reisenegger (2005), apj 625, 291-306.
=#

function FM(xi::Float64)
    return 1.0 + (22020.0*xi^2)/(11513.0*pi^2) + (5670.0*xi^4)/(11513.0*pi^4) + (420.0*xi^6)/(11513.0*pi^6) + (9.0*xi^8)/(11513.0*pi^8)
end

function HM(xi::Float64)
    return (14680.0*xi)/(11513.0*pi^2) + (7560.0*xi^3)/(11513.0*pi^4) + (840*xi^5)/(11513.0*pi^6) + (24.0*xi^7)/(11513.0*pi^8)
end

xi = 10000.
rn = 0.001
rp = 0.004020
@show In_rate_SFnp_trp(rn*xi, rp*xi, xi, 10, 10) / HM(xi)
@show Ip_rate_SFnp_trp(rn*xi, rp*xi, xi, 20) / HM(xi)
@show Ip_emis_SFnp_trp(rn*xi, rp*xi, xi, 20) / FM(xi)
# Crosscheck
xi = 10000.0
dpns = [[100., 100.],
        [200., 200.],
        [10.,900.],
        [300.,10.]]

for dpn in dpns .* 10
    dp = dpn[1]
    dn = dpn[2]
    Rpe = Ip_emis_SFnp_trp(dn, dp, xi, 40) / FM(xi)
    Rpr = Ip_rate_SFnp_trp(dn, dp, xi, 40) / HM(xi)
    @show (dn/xi, dp/xi, Rpr, Rpe)
end

# xi = 10000.
# dn = 0.011090909090909092
# dp = 0.00402020202020202

# for n = [5,10,15,20,40]
#     Rrate_p = Ip_rate_SFnp_trp(dn, dp, xi, n) / HM(xi)
#     @show n, Rrate_p
# end

