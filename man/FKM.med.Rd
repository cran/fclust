\name{FKM.med} 
\alias{FKM.med}
\title{Fuzzy \emph{k}-medoids}
\description{
 Performs the fuzzy \emph{k}-medoids clustering algorithm.
} 
\usage{ 
 FKM.med (X, k, m, RS, stand, startU, conv, maxit)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{Number of clusters (default: 2)}
 \item{m}{Parameter of fuzziness (default: 1.5)}
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
 \item{F}{Array containing the covariance matrices of all the clusters (\code{NULL} for \code{FKM.med})}
 \item{clus}{Matrix containing the indices of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indices of the medoid objects}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{cput}{Vector containing the computational times (user times) for the \code{RS} starts}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.med})}
 \item{vp}{Volume parameter (\code{NULL} for \code{FKM.med})} 
 \item{delta}{Noise distance (\code{NULL} for \code{FKM.med})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{call}{Matched call}
} 
\note{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational starting point.
 \cr In \code{FKM.med} the parameter of fuzziness is usually lower than the one used in \code{FKM}.
} 
\references{
 Krishnapuram R., Joshi A., Nasraoui O., Yi L., 2001. Low-complexity fuzzy relational clustering algorithms for web mining. IEEE Transactions on Fuzzy Systems, 9, 595-607.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{Fclust}}, \code{\link{Fclust.index}}}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
\dontrun{
## It may take more than a few seconds
clust=FKM.med(Mc[,1:(ncol(Mc)-1)],k=6,m=1.1,RS=10,stand=1)
}
}
\keyword{multivariate}