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

xi = 1000.0
l = 100
dp = range(1.0,stop=336.0,length=l)
dn = range(1.0,stop=1000.0,length=l)
xsize, = size(dn)
ysize, = size(dp)

Rs = zeros(ysize,xsize)

for iy in 1:ysize
    for ix in 1:xsize
        Rs[iy,ix] = Ip_emis_SFnp_trp(dn[ix], dp[iy], xi, 10) / FM(xi)
    end
    println(iy, ", ")
end

open("Remis_murca_p_xarray.dat", "w") do f
    println(f, "#range of x (neutron gap over xi)")
    for ix in 1:xsize
        print(f, dn[ix]/xi, " ")
    end
end

open("Remis_murca_p_yarray.dat", "w") do f
    println(f, "#range of y (proton gap over xi)")
    for iy in 1:ysize
        print(f, dp[iy]/xi, " ")
    end
end

open("Remis_murca_p.dat", "w") do f
    println(f, "#RnEmis [x:ngap, y:pgap]")
    for ix in 1:xsize
        for iy in 1:ysize
            print(f, Rs[iy,ix], " ")
        end
        print(f, "\n")
    end
end

