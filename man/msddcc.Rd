\name{MSD.DCC.EbEE}
\alias{MSD.DCC.EbEE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute mean and variance of estimators through Monte-Carlo methods
%%  ~~function to do ... ~~
}
\description{Generate a sample of the estimator to compute mean and standard error
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSD.DCC.EbEE(theta0, init, nobs, iter, type, noise, nu=Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta0}{List of the real parameters}
  \item{init}{List of initialisation parameters}
  \item{nobs}{Number of observations in the sample}
  \item{iter}{Number of iterations}
  \item{type}{type="Engel" for estimation as an Engle-DCC
\cr type="Aielli" for estimation as an Aielli-DCC}
  \item{noise}{"normal" or "student"}
  \item{nu}{Degrees of freedom of the t-distribution, leave blank if normal-noise}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{With usual notations of GARCH(1,1) DCC models\cr
\item{Omega.mean}{Mean of Omega}
\item{Omega.sd}{Standard deviation of Omega}
\item{A.mean}{Mean of A}
\item{A.sd}{Standard deviation of A}
\item{B.mean}{Mean of B}
\item{B.sd}{Standard deviation of B}
\item{S.mean}{Mean of S, correlation matrix}
\item{S.sd}{Standard deviation of S}
\item{alpha.mean}{Mean of alpha}
\item{alpha.sd}{Standard deviation of alpha}
\item{beta.mean}{Mean of beta}
\item{beta.sd}{Standard deviation of beta}
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
\note{Can take a lot of time if iter is big
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EbEEMGARCH}}{  Homepage of the documentation}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Sampling some data
m<-2
Omega <- c(1, 1);
A <- matrix(rep(0.025, m ^ 2), ncol = m)
B <- c(0.8, 0.8);
S <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
alpha <- 0.05;
beta <- 0.99 - alpha
n <- 2500
nu <-7
eps <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = nu, model="Aielli",noise = "student")

m<-2
omegainit <- rep(0.02, m)
Ainit <- matrix(rep(0.03, m ^ 2), ncol = m)
Binit <- rep(0.7, m)
Sinit<-diag(rep(1,m))
alphainit <- 0.05
betainit <- 0.90 - alphainit


init<-list(A=Ainit,B=Binit,S=Sinit,alpha=alphainit,beta=betainit,Omega=omegainit)
theta0<-list(Omega=Omega,A=A,B=B,alpha=alpha,beta=beta,S=S)
MSD.DCC.EbEE(theta0,init,2000,5,type="Aielli",noise="student",nu=7)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

