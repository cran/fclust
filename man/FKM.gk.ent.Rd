\name{FKM.gk.ent} 
\alias{FKM.gk.ent}
\title{Gustafson and Kessel - like fuzzy \emph{k}-means with entropy regularization}
\description{
 Performs the Gustafson and Kessel - like fuzzy \emph{k}-means clustering algorithm with entropy regularization.
} 
\usage{ 
 FKM.gk.ent (X, k, ent, vp, RS, stand, startU, conv, maxit)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{Number of clusters (default: 2)}
 \item{ent}{Degree of fuzzy entropy (default: 1)}
 \item{vp}{Volume parameter (default: rep(1,k))}
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
 \item{F}{Array containing the covariance matrices of all the clusters}
 \item{clus}{Matrix containing the indices of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indices of the medoid objects (\code{NULL} for \code{FKM.gk.ent})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{cput}{Vector containing the computational times (user times) for the \code{RS} starts}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness (\code{NULL} for \code{FKM.gk.ent})}
 \item{ent}{Degree of fuzzy entropy}
 \item{vp}{Volume parameter} 
 \item{delta}{Noise distance (\code{NULL} for \code{FKM.gk.ent})}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{call}{Matched call}
} 
\note{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the first element of \code{value}, \code{cput} and \code{iter} refer to the rational starting point.
 \cr If a cluster covariance matrix becomes singular, the algorithm stops and the element of \code{value} is {NaN}.
 \cr The default value for \code{ent} is in general not reasonable if \code{FKM.gk.ent} is run using raw data. 
 \cr The update of the membership degrees requires the computation of exponential functions. In some cases, this may produce \code{NaN} values and the algorithm stops. Such a problem is usually solved by running \code{FKM.gk.ent} using standardized data (\code{stand=1}).
} 
\references{
 Ferraro M.B., Giordani P., 2013. A new fuzzy clustering algorithm with entropy regularization. Proceedings of the meeting on Classification and Data Analysis (CLADAG).
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{Fclust}}, \code{\link{Fclust.index}}}
\examples{
data(unemployment)
\dontrun{
## It may take more than a few seconds
clust=FKM.gk.ent(unemployment,k=3,ent=0.2,RS=10,stand=1)
}
}
\keyword{multivariate}