\name{FKM.med.noise} 
\alias{FKM.med.noise}
\title{Fuzzy \emph{k}-medoids with noise cluster}
\description{
 Performs the fuzzy \emph{k}-medoids clustering algorithm with noise cluster.
 \cr Differently from fuzzy \emph{k}-means where the cluster prototypes (centroids) are artificial objects computed as weighted means, in the fuzzy \emph{k}-medoids the cluster prototypes (medoids) are a subset of the observed objects.
 \cr The noise cluster is an additional cluster (with respect to the \emph{k} standard clusters) such that objects recognized to be outliers are assigned to it with high membership degrees.
} 
\usage{ 
 FKM.med.noise (X, k, m, delta, RS, stand, startU, conv, maxit)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{Number of clusters (default: 2)}
 \item{m}{Parameter of fuzziness (default: 1.5)}
 \item{delta}{Noise distance (default: average squared Euclidean distance between objects and prototypes from \code{FKM.med} using the same values of \code{k} and \code{m})}
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
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.med.noise})}
 \item{clus}{Matrix containing the indices of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indices of the medoid objects}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{cput}{Vector containing the computational times (user times) for the \code{RS} starts}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.med.noise})}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.med.noise})}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.med.noise})} 
 \item{delta}{Noise distance}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{call}{Matched call}
} 
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational starting point.
 \cr As for \code{FKM.med}, in \code{FKM.med.noise} the parameter of fuzziness is usually lower than the one used in \code{FKM}.
} 
\references{
 Dave' R.N., 1991. Characterization and detection of noise in clustering. Pattern Recognition Letters, 12, 657-664.
 \cr Krishnapuram R., Joshi A., Nasraoui O., Yi L., 2001. Low-complexity fuzzy relational clustering algorithms for web mining. IEEE Transactions on Fuzzy Systems, 9, 595-607.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{FKM.med}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{butterfly}}}
\examples{
## butterfly data
data(butterfly)
## fuzzy k-medoids with noise cluster
clust=FKM.med.noise(butterfly,RS=5,delta=3)
}
\keyword{multivariate}


