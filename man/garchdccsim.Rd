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
GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = Inf, valinit = 500, model, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{With usual notations of MGARCH(1,1) DCCC
  \item{n}{Number of observation}
  \item{Omega}{Vector Omega (constant)}
  \item{A}{Matrix A}
  \item{B}{Vector of the diagonal of B}
  \item{alpha}{Scalar alpha in Aielli's notation}
  \item{beta}{Scalar beta in Aielli's notation}
  \item{S}{Variance of the noise (matrix)}
  \item{nu}{Degrees of freedom of the t-distribution, leave blank if normal-noise}
  \item{valinit}{Burn-in}
  \item{model}{type="Engle" for estimation as an Engle-DCC
\cr type="Aielli" for estimation as an Aielli-DCC
}
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

Omega <- c(0.01, 0.01);
A <- matrix(c(0.03, 0.01, 0.01, 0.03), nrow = 2)
B <- c(0.8, 0.8);
S <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
alpha <- 0.05;
beta <- 0.99 - alpha
n <- 2500
nu <-14

eps <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = nu, noise = "student",model="Aielli")

###
#Simulation of a Engle DCC semi-diagonal with normal noise
###

eps <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = Inf, noise = "normal",model="Engle")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
