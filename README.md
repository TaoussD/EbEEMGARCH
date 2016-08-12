# EbEEMGARCH

R package to estimate MGARCH(1,1) models equation by equation

This package provides various tools to simulate and estimate MGARCH(1,1) models. After, I will include tools to estimate CVaR of financial series who follow MGARCH models.


This package is based on these papers :
- C. Francq. & J.M. Zakoian, *Estimating multivariate GARCH and Stochastic Correlation models equation by equation*
- C. Francq. & J.M. Zakoian, *Joint inference on market and estimation risks in dynamic portfolios* 

All these papers are available at [http://perso.univ-lille3.fr/~cfrancq](http://perso.univ-lille3.fr/~cfrancq)

## Installation

The package can be installed from the sources available on the repo via this command in a R console ([devtools](https://github.com/hadley/devtools), [Rcpp](https://github.com/RcppCore/Rcpp/) and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) are required)
```R
library(Rcpp)
library(devtools)
install_github("TaoussD/EbEEMGARCH")
``` 

devtools and Rcpp can be easily installed in a R console (available on CRAN)
```R
install.packages("devtools")
install.packages("Rcpp")
```
If you still have troubles, please update R and install the associated Rtools.

## Methods

- estimCCC.EbEE : Estimation of the parameters of a MGARCH(1,1) CCC-diagonal or semi-diagonal model equation by equation
- estimDCC.EbEE : Estimation of the parameters of a Engle or Aielli MGARCH(1,1) DCC equation by equation
- residuals_DCC : Compute residuals, Ht and Rt for Engle & Aielli DCC models (usual notations)
- GarchCCC.sim : Simulation of MGARCH(1,1) CCC-diagonal or semi-diagonal 
- GarchDCC.sim : Simulation of Aielli or Engle MGARCH(1,1) DCC semi-diagonal
- MSD.CCC.EbEE : Compute mean and variance of the estimator through Monte Carlo methods for CCC models
- MSD.DCC.EbEE : Compute mean and variance of the estimator through Monte Carlo methods for DCC models
- VaR.Spherical : Compute CVaR of portfolios where yields follow a MGARCH(1,1) DCC and residuals have a spherical distribution
- VaR.FHS : Compute CVaR of portfolios where yields follow a MGARCH(1,1) DCC using FHS method
- estim.Markowitz : Compute the assets' weights and the associated CVaR using Markowitz methodology
- vech0 : vech0 operator
- inv.vech0 : Inverse of vech0 operator
- Sqrt : Square root of a symetric semi-definite positive matrix
- InvSqrt : Generalized inverse of a square root of a symetric semi-definite positive matrix

Homepage of the documentation available in R with

```R
?EbEEMGARCH
```

*More methods are under development*

## Authors

D. Taouss & C. Francq
