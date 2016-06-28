\name{estimCCC.EbEE}
\alias{estimCCC.EbEE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of multivariate GARCH CCC equation by equation
%%  ~~function to do ... ~~
}
\description{Estimation of a MGARCH(1,1) CCC-diagonal or semi-diagonal equation by equation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimCCC.EbEE(Omega, Alpha, Beta, eps, r = 10, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Omega}{Initialisation parameter for Omega (vector)
%%     ~~Describe \code{x} here~~
}
  \item{Alpha}{Initialisation parameter for Alpha (vector), must be a matrix if you are using the "sdiagonal" model
%%     ~~Describe \code{x} here~~
}
  \item{Beta}{Initialisation parameter for Beta (vector)
%%     ~~Describe \code{x} here~~
}
  \item{eps}{m*n matrix of the serie to be estimated
%%     ~~Describe \code{x} here~~
}
\item{r}{.
%%     ~~Describe \code{x} here~~
}
\item{model}{model="diagonal" if MGARCH(1,1) diagonal
\cr model="sdiagonal" if MGARCH(1,1) is semi-diagonal i.e. Alpha is plain and Beta diagonal
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Omega}{Estimation of Omega}
\item{Alpha}{Estimation of Alpha}
\item{Beta}{Estimation of Beta}
\item{R}{Estimation of R, the correlation matrix}
}
\references{
%% ~put references to the literature/web site here ~
C. Francq and J.M. Zakoian, Estimating multivariate GARCH and Stochastic Correlation models equation by equation, October 2014
}
\author{
%%  ~~who you are~~
D. Taouss, C. Francq
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EbEEMGARCH}}{  Homepage of the documentation}}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
####
# Sampling a diagonal model
####

m <- 3 #Number of assets in the portfolio
Omega0 <- rep(0.01, m)
Alpha0 <- rep(0.05, m)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))

Epsi <- GarchCCC.sim(2500, Omega0, Alpha0, Beta0,"diagonal", R0, "normal")


####
# Estimation of the parameters
####

Omegainit <- rep(0.1, m)
Alphainit <- rep(0.5, m)
Betainit <- rep(0.7, m)

estimCCC.EbEE(Omegainit,Alphainit,Betainit,Epsi,r=10,model="diagonal")

####
# If we allow a semi-diagonal model
####

omegainit <- rep(0.2, m)
alphainit <- matrix(rep(0.2 / m, m ^ 2), nrow = m)
betainit <- rep(0.70, m)

estimCCC.EbEE(omegainit,alphainit,betainit,Epsi,r=10,model="sdiagonal")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~MGARCH }
\keyword{ ~CCC model }% __ONLY ONE__ keyword per line
\keyword{ ~DCC model }% __ONLY ONE__ keyword per line

