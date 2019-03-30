module MurcaNoneqZeroTemp

export In_rate_SFnp_trp, Ip_rate_SFnp_trp, In_emis_SFnp_trp, Ip_emis_SFnp_trp

function Kr(u::Float64, w::Float64)::Float64
    # w=Delta_p/eta
    eps = 1.e-10
    return -3*w^2/2 * (w^2/4 + (1-u)^2) * log((1-u+sqrt((1-u)^2-w^2+eps))/w) + sqrt((1-u)^2-w^2+eps) * (13*w^2*(1-u)/8 + (1-u)^3/4)
end

function Ke(u::Float64, w::Float64)::Float64
    # w=Delta_p/eta
    eps = 1.e-10
    return -2*w^2 * (3*w^2/4.0*(1-u) + (1-u)^3) * log((1-u+sqrt((1-u)^2-w^2+eps))/w) + sqrt((1-u)^2-w^2+eps)/5.0 * ((1-u)^4 + w^2/6.0*(83*(1-u)^2 + 16*w^2))
end

function integrandMr(tn::Float64, ti::Float64, tf::Float64, dn::Float64, dni::Float64, dnf::Float64, dp::Float64, xsi::Float64)::Float64
    if xsi <= dp + dn + dni + dnf
        return 0.0
    else
    eps = 1.e-10
    rp = dp/xsi
    rn = dn/xsi
    rni = dni/xsi
    rnf = dnf/xsi
    xn = tanh(pi/2*sinh(tn))
    xi = tanh(pi/2*sinh(ti))
    xf = tanh(pi/2*sinh(tf))
    an = (1-rp-rni-rnf-rn)/2
    bn = (1-rp-rni-rnf+rn)/2
    un = an*xn + bn
    ai = (1-rnf-rp-un-rni)/2
    bi = (1-rnf-rp-un+rni)/2
    ui = ai*xi + bi
    af = (1-rp-un-ui-rnf)/2
    bf = (1-rp-un-ui+rnf)/2
    uf = af*xf + bf
    J = an*ai*af * (pi/2*cosh(tn))/cosh(pi/2*sinh(tn))^2 * (pi/2*cosh(ti))/cosh(pi/2*sinh(ti))^2 * (pi/2*cosh(tf))/cosh(pi/2*sinh(tf))^2
    return un*ui*uf/sqrt(un^2 - rn^2+eps)/sqrt(ui^2 - rni^2+eps)/sqrt(uf^2 - rnf^2+eps) * J * Kr(un+ui+uf, rp)
    end
end

function integrandMe(tn::Float64, ti::Float64, tf::Float64, dn::Float64, dni::Float64, dnf::Float64, dp::Float64, xsi::Float64)::Float64
    if xsi <= dp + dn + dni + dnf
        return 0.0
    else
    eps = 1.e-10
    rp = dp/xsi
    rn = dn/xsi
    rni = dni/xsi
    rnf = dnf/xsi
    xn = tanh(pi/2*sinh(tn))
    xi = tanh(pi/2*sinh(ti))
    xf = tanh(pi/2*sinh(tf))
    an = (1-rp-rni-rnf-rn)/2
    bn = (1-rp-rni-rnf+rn)/2
    un = an*xn + bn
    ai = (1-rnf-rp-un-rni)/2
    bi = (1-rnf-rp-un+rni)/2
    ui = ai*xi + bi
    af = (1-rp-un-ui-rnf)/2
    bf = (1-rp-un-ui+rnf)/2
    uf = af*xf + bf
    J = an*ai*af * (pi/2*cosh(tn))/cosh(pi/2*sinh(tn))^2 * (pi/2*cosh(ti))/cosh(pi/2*sinh(ti))^2 * (pi/2*cosh(tf))/cosh(pi/2*sinh(tf))^2
    return un*ui*uf/sqrt(un^2 - rn^2+eps)/sqrt(ui^2 - rni^2+eps)/sqrt(uf^2 - rnf^2+eps) * J * Ke(un+ui+uf, rp)
    end
end


function Itilde_rate(dn::Float64, dni::Float64, dnf::Float64, dp::Float64, xsi::Float64, xf=4.0, xi=-4.0, n=40.0)::Float64
    h = (xf - xi)/n
    res = 0.0
    prefactor = xsi^7/3 *60480/(11513*pi^8)
    for i1 in 1:n
        t1 = xi + h*i1
        for i2 in 1:n
            t2 = xi + h*i2
            for i3 in 1:n
                t3 = xi + h*i3
                res += integrandMr(t1, t2, t3, dn, dni, dnf, dp, xsi)*h^3
            end
        end
    end
    return res * prefactor
end

function Itilde_emis(dn::Float64, dni::Float64, dnf::Float64, dp::Float64, xsi::Float64, xf=4.0, xi=-4.0, n=40.0)::Float64
    h = (xf - xi)/n
    res = 0.0
    prefactor = xsi^8/4.0 *60480/(11513*pi^8)
    for i1 in 1:n
        t1 = xi + h*i1
        for i2 in 1:n
            t2 = xi + h*i2
            for i3 in 1:n
                t3 = xi + h*i3
                res += integrandMe(t1, t2, t3, dn, dni, dnf, dp, xsi)*h^3
            end
        end
    end
    return res * prefactor
end

function In_rate_SFnp_trp(dn0::Float64,dp::Float64, xi::Float64, nct=10.0, nphi=10.0)::Float64
    dct = 1.0/nct
    dphi = 2*pi/nphi
    
    ctns = dct:dct:1.0-dct
    phin1s = dphi:dphi:2*pi-dphi
    
    res = 0.0
    
    for cosn in [0.0, 1.0]
        sinn = sqrt(1-cosn^2)
        for phin1 in [0.0, 2*pi]
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.25 * Itilde_rate(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi)
        end
    end
            
    for cosn in [0.0, 1.0]
        sinn = sqrt(1-cosn^2)
        for phin1 in phin1s
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.5 * Itilde_rate(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi)
        end
    end
        
    for phin1 in [0.0, 2*pi]
        for cosn in ctns
            sinn = sqrt(1-cosn^2)
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.5 * Itilde_rate(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi)
        end
    end
    
    for cosn in ctns
        dn = dn0 * sqrt(1.0+3.0*cosn^2)
        sinn = sqrt(1-cosn^2)
        for phin1 in phin1s
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += Itilde_rate(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi) 
        end
        #println(res)
    end
    
    return res
end

function Ip_rate_SFnp_trp(dn0::Float64,dp::Float64, xi::Float64, nct=10.0)::Float64
    dct = 1.0/nct
    ctns = dct:dct:1.0-dct

    res = 0.0
    for cosn in [0.0, 1.0]
        dn = dn0 * sqrt(1.0 + 3.0*cosn^2)
        res += 0.5 * Itilde_rate(dn, dp, dp, dp, xi) * dct
    end

    for cosn in ctns
        dn = dn0 * sqrt(1.0 + 3.0*cosn^2)
        res += Itilde_rate(dn, dp, dp, dp, xi) * dct
    end

    return res
end

function In_emis_SFnp_trp(dn0::Float64,dp::Float64, xi::Float64, nct=10.0, nphi=10.0)::Float64
    #=
    Reduction factor; I/FM(xi)
    =#
    dct = 1.0/nct
    dphi = 2*pi/nphi
    
    ctns = dct:dct:1.0-dct
    phin1s = dphi:dphi:2*pi-dphi
    
    res = 0.0
    
    for cosn in [0.0, 1.0]
        sinn = sqrt(1-cosn^2)
        for phin1 in [0.0, 2*pi]
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.25 * Itilde_emis(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi) 
        end
    end
            
    for cosn in [0.0, 1.0]
        sinn = sqrt(1-cosn^2)
        for phin1 in phin1s
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.5 * Itilde_emis(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi) 
        end
    end
        
    for phin1 in [0.0, 2*pi]
        for cosn in ctns
            sinn = sqrt(1-cosn^2)
            dn = dn0 * sqrt(1.0+3.0*cosn^2)
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += 0.5 * Itilde_emis(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi) 
        end
    end
    
    for cosn in ctns
        dn = dn0 * sqrt(1.0+3.0*cosn^2)
        sinn = sqrt(1-cosn^2)
        for phin1 in phin1s
            cosn1 = -sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            cosn2 = sqrt(3)/2.0*sinn*cos(phin1) - 0.5*cosn
            dn1 = dn0 * sqrt(1.0+3.0*cosn1^2)
            dn2 = dn0 * sqrt(1.0+3.0*cosn2^2)
            res += Itilde_emis(dn, dn1, dn2, dp, xi) * dct * dphi/(2*pi) 
        end
        #println(res)
    end
    
    return res
end

function Ip_emis_SFnp_trp(dn0::Float64,dp::Float64, xi::Float64, nct=10.0)::Float64
    dct = 1.0/nct
    ctns = dct:dct:1.0-dct

    res = 0.0
    for cosn in [0.0, 1.0]
        dn = dn0 * sqrt(1.0 + 3.0*cosn^2)
        res += 0.5 * Itilde_emis(dn, dp, dp, dp, xi) * dct
    end

    for cosn in ctns
        dn = dn0 * sqrt(1.0 + 3.0*cosn^2)
        res += Itilde_emis(dn, dp, dp, dp, xi) * dct
    end

    return res
end
    
end
