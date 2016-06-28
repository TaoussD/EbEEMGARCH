\name{GarchDCC.sim}
\alias{GarchDCC.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of a MGARCH(1,1) DCC semi-diagonal
%%  ~~function to do ... ~~
}
\description{Simulation of a Engle or Aielli MGARCH(1,1) DCC semi-diagonal with student or normal noise
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GarchDCC.sim(n, omega, Alpha, beta, aalpha, bbeta, S, nu = Inf, valinit = 500, Aielli = TRUE, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observation}
  \item{omega}{Vector Omega}
  \item{Alpha}{Vector of the diagonal of Alpha}
  \item{beta}{Vector of the diagonal of Beta}
  \item{aalpha}{Scalar alpha in Aielli's notation}
  \item{bbeta}{Scalar beta in Aielli's notation}
  \item{S}{Variance of the noise (matrix)}
  \item{nu}{Degrees of freedom of the t-distribution, leave blank if normal-noise}
  \item{valinit}{Burn-in}
  \item{Aielli}{Aielli or Engle DCC-sample}
  \item{noise}{"normal" or "student"}

%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{eps}{Simulations}
\item{cor}{Correlation Matrix}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
###
#Simulation of a Aielli DCC semi-diagonal with student noise
###

omega <- c(0.01, 0.01);
Alpha <- matrix(c(0.03, 0.01, 0.01, 0.03), nrow = 2)
beta <- c(0.8, 0.8);
S <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
aalpha <- 0.01;
bbeta <- 0.99 - aalpha
n <- 800
nu <-7

eps <- GarchDCC.sim(n, omega, Alpha, beta, aalpha, bbeta, S, nu = nu, noise = "student")

###
#Simulation of a Engle DCC semi-diagonal with normal noise
###

eps <- GarchDCC.sim(n, omega, Alpha, beta, aalpha, bbeta, S, noise = "normal",Aielli=FALSE)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
