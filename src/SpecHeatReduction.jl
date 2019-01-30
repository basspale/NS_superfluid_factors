module SpecHeatReduction

"""
Superfluid reduction factors for nucleon specfic heat

A: 1S0
B: 3P2(mJ=0)
C: 3P2(mJ=2)

Refs:
Levenfish and Yakovlev, Astronomy Reports, Volume 38, Issue 2, March 1994, pp.247-251
Yakovlev and Levenfish, Phys. Usp., 42, 737-778, 1999
"""

using FastGaussQuadrature

export vA, vB, vC, RA, RB, RC, RA_fit, RB_fit, RC_fit

# Fitting fromulas for gap as a function of temperature; t = T/T_c
# vp{} is its derivative for t
vA(t) = sqrt(1-t)*(1.456 - 0.157/sqrt(t) + 1.764/t)
vpA(t) = -0.5/sqrt(1-t)*(1.456 - 0.157/sqrt(t) + 1.764/t) + sqrt(1-t)*(0.5*0.157/t^(1.5) - 1.764/t^2)

vB(t) = sqrt(1-t)*(0.7893 + 1.188/t)
vpB(t) = -0.5/sqrt(1-t)*(0.7893 + 1.188/t) - sqrt(1-t)*1.188/t^2

vC(t) = sqrt(1-t^4)/t*(2.03 - 0.4903*t^4 + 0.1727*t^8)
vpC(t) = -2*t^3/sqrt(1-t^4)/t*(2.03 - 0.4903*t^4 + 0.1727*t^8) + sqrt(1-t^4)*(-2.03/t^2 - 3*0.4903*t^2 + 7*0.1727*t^6)


# Nuerical integration by Gauss-Laguerre quadrature
function RA(t, n)
    factor = -6/pi^2
    
    nodes, weights = gausslaguerre(n)
    
    y = vA(t)
    I = 0.0
    
    for i=1:n
        x = nodes[i]
        w = weights[i]
        
        z = sqrt(x^2 + y^2)
        I += w * exp(z+x)/(exp(z) + 1)^2 * (-z^2 + vA(t)^2 + t*vA(t)*vpA(t))
    end
    #println(factor)
    return I * factor
end

function RB(t, n, h)
    factor = -6/pi^2
    nodes, weights = gausslaguerre(n)
    
    v = vB(t)
    I = 0.0
    
    for ct=0.0:h:1.0, i=1:n
        x = nodes[i]
        w = weights[i]
        F = 1+3*ct^2
        
        z = sqrt(x^2 + v^2*F)
        I += w * exp(z+x)/(exp(z) + 1)^2 *  (-z^2 + vB(t)^2*F + t*vB(t)*vpB(t)*F)
    end
    return I * factor * h
end

function RC(t, n, h)
    factor = -6/pi^2
    nodes, weights = gausslaguerre(n)
    
    v = vC(t)
    I = 0.0
    
    for ct=0.0:h:1.0, i=1:n
        x = nodes[i]
        w = weights[i]
        F = 1-ct^2
        
        z = sqrt(x^2 + v^2*F)
        I += w * exp(z+x)/(exp(z) + 1)^2 * (-z^2 + vC(t)^2*F + t*vC(t)*vpC(t)*F)
    end
    return I * factor * h
end

# Fitting formula for reduction factors
RA_fit(v) = (0.4186 + sqrt(1.007^2+0.5010^2*v^2))^(2.5) * exp(1.456 - sqrt(1.456^2 + v^2))
RB_fit(v) = (0.6893 + sqrt(0.790^2+0.2824^2*v^2))^2 * exp(1.934 - sqrt(1.934^2+v^2))
RC_fit(v) = (2.188 - (9.537e-5)^2*v^2 + 0.1491^4*v^4)/(1.0 + 0.2846^2*v^2 + 0.01335^4*v^4 + 0.1815^6*v^6)

end
