# Crosschecking fitting formulas of various superfluid reduction factors of nucleons.

## What is reduction factor?
The thermal property of neutron star (NS) is greatly influenced by the nucleon (neutron and/or proton) superfluidity, which makes energy gap above the Fermi surface.

The main effects are
- Reduction of specific heat

- Reduction of emissivity of modified Urca procss (and of direct Urca process for vary massive NSs)

These are calculated and summarized as simple fitting functions in the literature.

The purpose is this repositry is to crosscheck these formulas by the direct numerical integration.

## Method

In performing the energy integration, the Gauss-Laguerre quadrature method is particularly useful.
See [wikipedia](https://en.wikipedia.org/wiki/Gaussian_quadrature) for the detail.

The calculation is done by Julia. In the code, I use [FastGaussQuadrature](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) for energy integration.

The neutron gap is also dependent on the angle around the quantization axis.
In integrating over the angle, I use simple trapezoidal rule.

## Codes

-```src```: Julia codes for numerical integration or computation of fitting formulas.

-```notebook```: Jupyter notebooks to compare the numerical integration to the fitting formulas.

-```output_data```: store the results of numerical integration that takes a long run time.

-```input_data```: data of coefficients for fitting formulas.

## Refs.
- [Yakovlev et al., 2001](http://www.sciencedirect.com/science/article/pii/S0370157300001319)

- [Yakovlev et al., 1999](http://mr.crossref.org/iPage?doi=10.1070%2FPU1999v042n08ABEH000556)

- [Levenfish et al., 1994](http://ads.nao.ac.jp/abs/1994AstL...20...43L)
