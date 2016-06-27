# EbEEMGARCH
*Package under development*

R package to estimate MGARCH(1,1) model equation by equation

This package provides various tools to simulate and estimate MGARCH(1,1) models. After, I will include tools to estimate VaR of financial series who follow MGARCH models.


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


## Methods

- estim.EBEE : Estimation of the parameters of a MGARCH(1,1) CCC-diagonal or semi-diagonal model equation by equation
- Mgarch.sim : Simulation of MGARCH(1,1) CCC-diagonal or semi-diagonal 
- vech0 : vech0 operator
- Sqrt : Square root of a symetric semi-definite positive matrix

Homepage of the documentation available in R with

```R
?EbEEMGARCH
```

*More methods are under development*

## Authors

D. Taouss & C. Francq
