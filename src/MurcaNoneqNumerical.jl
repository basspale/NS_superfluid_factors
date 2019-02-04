module MurcaNoneqNumerical

"""
Superfluid reduction factor for modified Urca process.
Numerical integration is done with Gauss-Laguerre quadrature and trapezotal rule.

Reduction factor is named I{n,p}_SF{n,p}{n,p}, meaning that
- In (Ip): reduction factor for neutron (proton) branch.
- SFn (SFp): only neutron (proton) superfluidity is turned on 
- SFnp: both proton and neutron are superfluid.
"""

export In_SFp, Ip_SFp, In_SFn, Ip_SFn, In_SFnp, Ip_SFnp

using FastGaussQuadrature
#using Dierckx

f(x::Float64) = 1.0/(exp(x)+1.0)

"""
Both proton and neutron superfluidity
"""

# function integrand_In_SFnp(tnu::Float64, tn::Float64, tp::Float64, tn1::Float64, tn2::Float64,
#                            yn::Float64, yp::Float64, yn1::Float64, yn2::Float64, xi::Float64)
#     prefactor = 60480.0/11513.0/pi^8
#     integrand_tmp_emis = 0.0
#     integrand_tmp_rate = 0.0

#     Snu = xi + yn + yp + yn1 + yn2
#     Sn = xi + sqrt(2*yn) + yp + yn1 + yn2
#     Sp = xi + yn + sqrt(2*yp) +yn1 + yn2
#     Sn1 = xi + yn + yp + sqrt(2*yn1) + yn2
#     Sn1 = xi + yn + yp + yn1 + sqrt(2*yn2)
#     xnu = Snu*tnu
#     xn = Sn*tn
#     xp = Sp*tp
#     xn1 = Sn1*tn1
#     xn2 = Sn2*tn2
#     J = Snu*Sn*Sp*Sn1*Sn2

#     zp = sqrt(xp^2 + yp^2)
#     zn = sqrt(xn^2 + yn^2)
#     zn1 = sqrt(xn1^2 + yn1^2)
#     zn2 = sqrt(xn2^2 + yn2^2)
#     for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
#         integrand_tmp_emis += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zn1)*f(j4*zn2) * (f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2-xi)+f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2+xi)) * J
#         integrand_tmp_rate += xnu^2*f(j1*zn)*f(j2*zp)*f(j3*zn1)*f(j4*zn2) * (f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2-xi)-f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2+xi)) * J
#     end
    
#     return [prefactor * integrand_tmp_rate, prefactor * integrand_tmp_emis]
# end

# function integrand_Ip_SFnp(tnu::Float64, tn::Float64, tp::Float64, tp1::Float64, tp2::Float64, yn::Float64, yp::Float64, xi::Float64)
#     prefactor = 60480.0/11513.0/pi^8
#     integrand_tmp_emis = 0.0
#     integrand_tmp_rate = 0.0
    
#     Snu = xi + yn + yp*3
#     Sn = xi + sqrt(2*yn) + yp*3
#     Sp = xi + yn + sqrt(2*yp) + 2*yp
#     Sp1 = xi + yn + sqrt(2*yp) + 2*yp
#     Sp2 = xi + yn + sqrt(2*yp) + 2*yp
#     xnu = Snu*tnu
#     xn = Sn*tn
#     xp = Sp*tp
#     xp1 = Sp1*tp1
#     xp2 = Sp2*tp2
#     J = Snu*Sn*Sp*Sp1*Sp2
    
#     zn = sqrt(xn^2 + yn^2)
#     zp = sqrt(xp^2 + yp^2)
#     zp1 = sqrt(xp1^2 + yp^2)
#     zp2 = sqrt(xp2^2 + yp^2)
#     for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
#         integrand_tmp_emis += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * (f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2-xi)+f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2+xi)) * J
#         integrand_tmp_rate += xnu^2*f(j1*zn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * (f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2-xi)-f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2+xi)) * J
#     end

#     return [prefactor * integrand_tmp_rate, prefactor * integrand_tmp_emis]
# end

function integrand_In_SFnp(xnu::Float64, xn::Float64, xp::Float64, xn1::Float64, xn2::Float64,
        yn::Float64, yp::Float64, yn1::Float64, yn2::Float64, xi::Float64)
    prefactor = 60480.0/11513.0/pi^8
    integrand_tmp_emis = 0.0
    integrand_tmp_rate = 0.0
    zp = sqrt(xp^2 + yp^2)
    zn = sqrt(xn^2 + yn^2)
    zn1 = sqrt(xn1^2 + yn1^2)
    zn2 = sqrt(xn2^2 + yn2^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp_emis += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zn1)*f(j4*zn2) * (f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2-xi)+f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2+xi))
        integrand_tmp_rate += xnu^2*f(j1*zn)*f(j2*zp)*f(j3*zn1)*f(j4*zn2) * (f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2-xi)-f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2+xi))
    end
    
    return [prefactor * integrand_tmp_rate, prefactor * integrand_tmp_emis]
end

function integrand_Ip_SFnp(xnu::Float64, xn::Float64, xp::Float64, xp1::Float64, xp2::Float64, yn::Float64, yp::Float64, xi::Float64)
    prefactor = 60480.0/11513.0/pi^8
    integrand_tmp_emis = 0.0
    integrand_tmp_rate = 0.0
    zn = sqrt(xn^2 + yn^2)
    zp = sqrt(xp^2 + yp^2)
    zp1 = sqrt(xp1^2 + yp^2)
    zp2 = sqrt(xp2^2 + yp^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp_emis += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * (f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2-xi)+f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2+xi))
        integrand_tmp_rate += xnu^2*f(j1*zn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * (f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2-xi)-f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2+xi))
    end

    return [prefactor * integrand_tmp_rate, prefactor * integrand_tmp_emis]
end

function In_SFnp(vn::Float64, vp::Float64, xi::Float64, n::Int64, h1::Float64, h2::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h1:h1:1.0-h1
    phi1s = 0.0-h2:h2:2*pi-h2
    
    yp=vp
    
    integ_tmp = zeros(2)
    
    for ctn=ctns, phi1=phi1s
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2,xi) .*h1*h2/2/pi
        end
    end
    
    println(integ_tmp)
    
    for ctn=ctns, phi1=[0, 2*pi]
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (0.5 * wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2,xi) .*h1*h2/2/pi
        end
    end
    
    for ctn=[0.0,1.0], phi1=phi1s
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (0.5 * wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2,xi) .*h1*h2/2/pi
            #integ_tmp .+= wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi/2
        end
    end
    
    for ctn=[0.0,1.0], phi1=[0, 2*pi]
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (0.25 * wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2,xi) .*h1*h2/2/pi
            #integ_tmp .+= wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi/4
        end
    end
    
    return integ_tmp
end

function Ip_SFnp(vn::Float64, vp::Float64, xi::Float64, n::Int64, h::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h:h:1.0-h
    yp = vp
    
    integ_tmp = zeros(2)
    
    for ctn=ctns
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_Ip_SFnp(xnu, x1, x2, x3, x4, yn, yp, xi) .* h
        end
    end
    
    println(integ_tmp)
    
    for ctn=[0.0,1.0]
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp .+= (0.5 * wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4)) .* integrand_Ip_SFnp(xnu, x1, x2, x3, x4, yn, yp, xi) .* h
            #integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFnp(xnu, x1, x2, x3, x4, yn, yp) *h/2
        end
    end
    
    return integ_tmp
end

end
