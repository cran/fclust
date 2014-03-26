\name{FKM.noise} 
\alias{FKM.noise}
\title{Fuzzy \emph{k}-means with noise cluster}
\description{
 Performs the fuzzy \emph{k}-means clustering algorithm with noise cluster.
} 
\usage{ 
 FKM.noise (X, k, m, delta, RS, stand, startU, conv, maxit)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{Number of clusters (default: 2)}
 \item{m}{Parameter of fuzziness (default: 2)}
 \item{delta}{Noise distance (default: average squared Euclidean distance between objects and prototypes from \code{FKM} using the same values of \code{k} and \code{m})}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational starting point for the membership degree matrix \code{U} (default: no rational start)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+6)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.noise})}
 \item{clus}{Matrix containing the indices of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indices of the medoid objects (\code{NULL} for \code{FKM.noise})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{cput}{Vector containing the computational times (user times) for the \code{RS} starts}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.noise})}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.noise})} 
 \item{delta}{Noise distance}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{call}{Matched call}
} 
\note{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational starting point.
} 
\references{
 Dave' R.N., 1991. Characterization and detection of noise in clustering. Pattern Recognition Letters, 12, 657-664.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{Fclust}}, \code{\link{Fclust.index}}}
\examples{
data(butterfly)
clust=FKM.noise(butterfly,RS=5,delta=3)
}
\keyword{multivariate}