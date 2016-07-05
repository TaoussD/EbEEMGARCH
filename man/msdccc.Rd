\name{MSD.CCC.EbEE}
\alias{MSD.CCC.EbEE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute mean and variance of estimators through Monte-Carlo methods
%%  ~~function to do ... ~~
}
\description{Generate a sample of the estimator to compute mean and standard error
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSD.CCC.EbEE(theta0, init, nobs, iter, type, noise, nu = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta0}{List of the real parameters}
  \item{init}{List of initialisation parameters}
  \item{nobs}{Number of observations in the sample}
  \item{iter}{Number of iterations}
  \item{type}{type="diagonal" for estimation as an MGARCH(1,1) CCC-diagonal
\cr type="sdiagonal" for estimation as a MGARCH(1,1) CCC semi-diagonal}
  \item{noise}{"normal" or "student"}
  \item{nu}{Degrees of freedom of the t-distribution, leave blank if normal-noise}
%%     ~~Describe \code{x} here~~
}

\details{Check example to see how to create the lists in argument
%%  ~~ If necessary, more details than the description above ~~
}
\value{With usual notations of MGARCH(1,1) CCC models\cr
\item{Omega.mean}{Mean of Omega}
\item{Omega.sd}{Standard deviation of Omega}
\item{A.mean}{Mean of A}
\item{A.sd}{Standard deviation of A}
\item{B.mean}{Mean of B}
\item{B.sd}{Standard deviation of B}
\item{R.mean}{Mean of R, correlation matrix}
\item{R.sd}{Standard deviation of R}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{C. Francq & J.M. Zakoian, Estimating multivariate GARCH and Stochastic Correlation models equation by equation, October 2014
%% ~put references to the literature/web site here ~
}
\author{D. Taouss & C. Francq
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EbEEMGARCH}}{  Homepage of the documentation}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m<-3
Omega0 <- rep(0.01, m)
Alpha0 <- rep(0.05, m)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))
theta0<-list(Omega=Omega0,A=Alpha0,B=Beta0,R=R0)


Omegainit <- rep(0.1, m)
Alphainit <- rep(0.5, m)
Betainit <- rep(0.7, m)
init<-list(Omega=Omegainit,A=Alphainit,B=Betainit)

MSD.CCC.EbEE(theta0,init,2000,10,"diagonal","normal")
}
