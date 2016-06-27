\name{Mgarch.sim}
\alias{Mgarch.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of a MGARCH(1,1) CCC-diagonal or semi-diagonal
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
Fonction to simulate data from full diagonal MGARCH(1,1) model
\usage{
Mgarch.sim(n, omega, alpha, beta, R, noise, nu = Inf, valinit = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observation}
  \item{omega}{Vector Omega}
  \item{alpha}{Vector of the diagonal of Alpha}
  \item{beta}{Vector of the diagonal of Beta}
  \item{model}{model="diagonal" if MGARCH(1,1) diagonal
  \cr model="sdiagonal" if MGARCH(1,1) is semi-diagonal i.e. Alpha is plain and Beta diagonal}
  \item{R}{Variance of the noise (matrix)}
  \item{noise}{"normal" or "student"}
  \item{nu}{Degrees of freedom of the t-distribution, leave blank if normal-noise}
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
dataframe of the observations
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
\code{\link{EbEEMGARCH}}{  Homepage of the documentation}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

####
# Sampling a diagonal model with normal noise
####

m <- 3 #Number of assets in the portfolio
Omega0 <- rep(0.01, m)
Alpha0 <- rep(0.05, m)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))

Epsi <- Mgarch.sim(2500, Omega0, Alpha0, Beta0,"diagonal", R0, "normal")

####
# Sampling a semi-d1iagonal model with student noise
####

m <- 3 #Number of assets in the portfolio
Omega0 <- rep(0.01, m)
Alpha0 <- matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))

Epsi <- Mgarch.sim(2500, Omega0, Alpha0, Beta0,"sdiagonal", R0, "student",7)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ MGARCH }% __ONLY ONE__ keyword per line
