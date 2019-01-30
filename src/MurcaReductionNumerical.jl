module MurcaReductionNumerical

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
Only proton superfluidity
"""

function integrand_In_SFp(xnu::Float64, xn::Float64, xp::Float64, xn1::Float64, xn2::Float64, y::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zp = sqrt(xp^2 + y^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*xn)*f(j2*zp)*f(j3*xn1)*f(j4*xn2) * f(xnu-j1*xn-j2*zp-j3*xn1-j4*xn2)
    end
    
    return prefactor * integrand_tmp
end

function integrand_Ip_SFp(xnu::Float64, xn::Float64, xp::Float64, xp1::Float64, xp2::Float64, y::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zp = sqrt(xp^2 + y^2)
    zp1 = sqrt(xp1^2 + y^2)
    zp2 = sqrt(xp2^2 + y^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*xn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * f(xnu-j1*xn-j2*zp-j3*zp1-j4*zp2)
    end
    
    return prefactor * integrand_tmp
end

function In_SFp(y::Float64, n::Int64)
    nodes, weights = gausslaguerre(n)
    
    integ_tmp = 0.0
    for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
        wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
        integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFp(xnu, x1, x2, x3, x4, y)
    end
    
    return integ_tmp
end

function Ip_SFp(y::Float64, n::Int64)
    nodes, weights = gausslaguerre(n)
    
    integ_tmp = 0.0
    for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
        wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
        integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFp(xnu, x1, x2, x3, x4, y)
    end
    
    return integ_tmp
end


"""
Only neutron superfluidity
"""

function integrand_In_SFn(xnu::Float64, xn::Float64, xp::Float64, xn1::Float64, xn2::Float64,
        yn::Float64, yn1::Float64, yn2::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zn = sqrt(xn^2 + yn^2)
    zn1 = sqrt(xn1^2 + yn1^2)
    zn2 = sqrt(xn2^2 + yn2^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*zn)*f(j2*xp)*f(j3*zn1)*f(j4*zn2) * f(xnu-j1*zn-j2*xp-j3*zn1-j4*zn2)
    end
    
    return prefactor * integrand_tmp
end

function integrand_Ip_SFn(xnu::Float64, xn::Float64, xp::Float64, xp1::Float64, xp2::Float64, yn::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zn = sqrt(xn^2 + yn^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*zn)*f(j2*xp)*f(j3*xp1)*f(j4*xp2) * f(xnu-j1*zn-j2*xp-j3*xp1-j4*xp2)
    end
    
    return prefactor * integrand_tmp
end

function In_SFn(vn::Float64, n::Int64, h1::Float64, h2::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h1:h1:1.0-h1
    phi1s = 0.0-h2:h2:2*pi-h2
    
    integ_tmp = 0.0
    
    for ctn=ctns, phi1=phi1s
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFn(xnu, x1, x2, x3, x4, yn, yn1, yn2) *h1*h2/2/pi
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFn(xnu, x1, x2, x3, x4, yn, yn1, yn2) *h1*h2/2/pi/2
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFn(xnu, x1, x2, x3, x4, yn, yn1, yn2) *h1*h2/2/pi/2
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFn(xnu, x1, x2, x3, x4, yn, yn1, yn2) *h1*h2/2/pi/4
        end
    end
    
    return integ_tmp
end

function Ip_SFn(vn::Float64, n::Int64, h::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h:h:1.0-h
    
    integ_tmp = 0.0
    
    for ctn=ctns
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFn(xnu, x1, x2, x3, x4, yn) *h
        end
    end
    
    println(integ_tmp)
    
    for ctn=[0.0,1.0]
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFn(xnu, x1, x2, x3, x4, yn) *h/2
        end
    end
    
    return integ_tmp
end


"""
Both proton and neutron superfluidity
"""

function integrand_In_SFnp(xnu::Float64, xn::Float64, xp::Float64, xn1::Float64, xn2::Float64,
        yn::Float64, yp::Float64, yn1::Float64, yn2::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zp = sqrt(xp^2 + yp^2)
    zn = sqrt(xn^2 + yn^2)
    zn1 = sqrt(xn1^2 + yn1^2)
    zn2 = sqrt(xn2^2 + yn2^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zn1)*f(j4*zn2) * f(xnu-j1*zn-j2*zp-j3*zn1-j4*zn2)
    end
    
    return prefactor * integrand_tmp
end

function integrand_Ip_SFnp(xnu::Float64, xn::Float64, xp::Float64, xp1::Float64, xp2::Float64, yn::Float64, yp::Float64)
    prefactor = 120960.0/11513.0/pi^8
    integrand_tmp = 0.0
    zn = sqrt(xn^2 + yn^2)
    zp = sqrt(xp^2 + yp^2)
    zp1 = sqrt(xp1^2 + yp^2)
    zp2 = sqrt(xp2^2 + yp^2)
    for j1=[-1,1], j2=[-1,1], j3=[-1,1], j4=[-1,1]
        integrand_tmp += xnu^3*f(j1*zn)*f(j2*zp)*f(j3*zp1)*f(j4*zp2) * f(xnu-j1*zn-j2*zp-j3*zp1-j4*zp2)
    end
    
    return prefactor * integrand_tmp
end

function In_SFnp(vn::Float64, vp::Float64, n::Int64, h1::Float64, h2::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h1:h1:1.0-h1
    phi1s = 0.0-h2:h2:2*pi-h2
    
    yp=vp
    
    integ_tmp = 0.0
    
    for ctn=ctns, phi1=phi1s
        ctn1 = -0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        ctn2 = +0.5*sqrt(3) * sqrt(1-ctn^2)* cos(phi1) -0.5*ctn
        yn = vn*sqrt(1+3*ctn^2)
        yn1 = vn*sqrt(1+3*ctn1^2)
        yn2 = vn*sqrt(1+3*ctn2^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
        
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi/2
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi/2
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
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_In_SFnp(xnu, x1, x2, x3, x4, yn, yp, yn1, yn2) *h1*h2/2/pi/4
        end
    end
    
    return integ_tmp
end

function Ip_SFnp(vn::Float64, vp::Float64, n::Int64, h::Float64)
    nodes, weights = gausslaguerre(n)
    
    ctns = 0.0-h:h:1.0-h
    yp = vp
    
    integ_tmp = 0.0
    
    for ctn=ctns
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFnp(xnu, x1, x2, x3, x4, yn, yp) *h
        end
    end
    
    println(integ_tmp)
    
    for ctn=[0.0,1.0]
        yn = vn*sqrt(1+3*ctn^2)
        for inu=1:n, i1=1:n, i2=1:n, i3=1:n, i4=1:n
            xnu, x1, x2, x3, x4 = nodes[[inu, i1, i2, i3, i4]]
            wnu, w1, w2, w3, w4 = weights[[inu, i1, i2, i3, i4]]
            integ_tmp += wnu*w1*w2*w3*w4 * exp(xnu+x1+x2+x3+x4) * integrand_Ip_SFnp(xnu, x1, x2, x3, x4, yn, yp) *h/2
        end
    end
    
    return integ_tmp
end

end
