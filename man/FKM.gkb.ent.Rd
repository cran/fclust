\name{FKM.gkb.ent}
\alias{FKM.gkb.ent}
\title{Gustafson, Kessel and Babuska - like fuzzy k-means with entropy regularization}
\description{
 Performs the Gustafson, Kessel and Babuska - like fuzzy \emph{k}-means clustering algorithm with entropy regularization.
 \cr Differently from fuzzy \emph{k}-means, it is able to discover non-spherical clusters.
 \cr The Babuska et al. variant improves the computation of the fuzzy covariance matrices in the standard Gustafson and Kessel clustering algorithm.
 \cr The entropy regularization allows us to avoid using the artificial fuzziness parameter \emph{m}. This is replaced by the degree of fuzzy entropy \emph{ent}, related to the concept of temperature in statistical physics.
 An interesting property of the fuzzy \emph{k}-means with entropy regularization is that the prototypes are obtained as weighted means with weights equal to the membership degrees (rather than to the membership degrees at
 the power of \emph{m} as is for the fuzzy k-means).
}
\usage{
 FKM.gkb.ent (X, k, ent, vp, gam, mcn, RS, stand, startU, index, alpha, conv, maxit, seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{ent}{Degree of fuzzy entropy (default: 1)}
 \item{vp}{Volume parameter (default: rep(1,k))}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices (default: 0)}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices (default: 1e+15)}
 \item{RS}{Number of (random) starts (default: 1)}
 \item{stand}{Standardization: if \code{stand=1}, the clustering algorithm is run using standardized data (default: no standardization)}
 \item{startU}{Rational start for the membership degree matrix \code{U} (default: no rational start)}
 \item{index}{Cluster validity index to select the number of clusters: \code{"PC"} (partition coefficient), \code{"PE"} (partition entropy), \code{"MPC"} (modified partition coefficient), \code{"SIL"} (silhouette), \code{"SIL.F"} (fuzzy silhouette), \code{"XB"} (Xie and Beni) (default: "SIL.F")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
 \item{conv}{Convergence criterion (default: 1e-9)}
 \item{maxit}{Maximum number of iterations (default: 1e+2)}
 \item{seed}{Seed value for random number generation (default: NULL)}
}
\value{
 Object of class \code{fclust}, which is a list with the following components:
 \item{U}{Membership degree matrix}
 \item{H}{Prototype matrix}
 \item{F}{Array containing the covariance matrices of all the clusters}
 \item{clus}{Matrix containing the indexes of the clusters where the objects are assigned (column 1) and the associated membership degrees (column 2)}
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.gkb.ent})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{A integer value or vector indicating the number of clusters. (default: 2:6)}
 \item{m}{Parameter of fuzziness (\code{NULL} for \code{FKM.gkb.ent})}
 \item{ent}{Degree of fuzzy entropy}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.gkb.ent})}
 \item{vp}{Volume parameter (default: \code{rep(1,max(k)}). If \code{k} is a vector, for each group the first \code{k} element of \code{vp}are considered.}
 \item{delta}{Noise distance (\code{NULL} for \code{FKM.gkb.ent})}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.gkb.ent})}
 \item{call}{Matched call}
}
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the argument \code{RS} is ignored (the algorithm is run using the rational start) and therefore \code{value}, \code{cput} and \code{iter} refer to such a rational start.
 \cr If a cluster covariance matrix becomes singular, the algorithm stops and the element of \code{value} is \code{NaN}.
 \cr The default value for \code{ent} is in general not reasonable if \code{FKM.gk.ent} is run using raw data.
 \cr The update of the membership degrees requires the computation of exponential functions. In some cases, this may produce \code{NaN} values and the algorithm stops. Such a problem is usually solved by running \code{FKM.gk.ent} using standardized data (\code{stand=1}).
}
\references{
 Babuska R., van der Veen P.J., Kaymak U., 2002.  Improved covariance estimation for Gustafson-Kessel clustering. Proceedings of the IEEE International Conference on Fuzzy Systems (FUZZ-IEEE), 1081-1085.
 \cr Ferraro M.B., Giordani P., 2013. A new fuzzy clustering algorithm with entropy regularization. Proceedings of the meeting on Classification and Data Analysis (CLADAG).
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM.gk.ent}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{unemployment}}}
\examples{
\dontrun{
## unemployment data
data(unemployment)
## Gustafson, Kessel and Babuska-like fuzzy k-means with entropy regularization,
##fixing the number of clusters
clust=FKM.gkb.ent(unemployment,k=3,ent=0.2,RS=10,stand=1)
## Gustafson, Kessel and Babuska-like fuzzy k-means with entropy regularization,
##selecting the number of clusters
clust=FKM.gkb.ent(unemployment,k=2:6,ent=0.2,RS=10,stand=1)}
}
\keyword{multivariate}
