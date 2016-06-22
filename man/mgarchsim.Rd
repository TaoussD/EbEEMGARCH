\name{Mgarch.sim}
\alias{Mgarch.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of diagonal Multivariate GARCH
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
Fonction to simulate data from full diagonal MGARCH(1,1) model
\usage{
Mgarch.sim(n, omega, alpha, beta, R, nu = Inf, valinit = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observation}
  \item{omega}{Vector Omega}
  \item{alpha}{Vector of the diagonal of Alpha}
  \item{beta}{Vector of the diagonal of Beta}
  \item{R}{Correlation matrix}
  \item{nu}{??}
  \item{valinit}{Burn-in}
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{D. Taouss & C. Francq
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Omega0 <- rep(0.01, m)
Alpha0 <- rep(0.05, m)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))

Epsi <- Mgarch.sim(n = 2500, omega = Omega0, alpha = Alpha0, beta = Beta0, R.mat = R0, nu = Inf, valinit = 500)
Epsi <- as.data.frame(Epsi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ MGARCH }% __ONLY ONE__ keyword per line
