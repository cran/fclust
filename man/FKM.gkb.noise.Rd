\name{FKM.gkb.noise}
\alias{FKM.gkb.noise}
\title{Gustafson, Kessel and Babuska - like fuzzy k-means with noise cluster}
\description{
 Performs the Gustafson, Kessel and Babuska - like fuzzy \emph{k}-means clustering algorithm with noise cluster.
 \cr Differently from fuzzy \emph{k}-means, it is able to discover non-spherical clusters.
 \cr The Babuska et al. variant improves the computation of the fuzzy covariance matrices in the standard Gustafson and Kessel clustering algorithm.
 \cr The noise cluster is an additional cluster (with respect to the \emph{k} standard clusters) such that objects recognized to be outliers are assigned to it with high membership degrees.
}
\usage{
 FKM.gkb.noise (X,k,m,vp,delta,gam,mcn,RS,stand,startU,index,alpha,conv,maxit,seed)
}
\arguments{
 \item{X}{Matrix or data.frame}
 \item{k}{An integer value or vector specifying the number of clusters for which the \code{index} is to be calculated (default: 2:6)}
 \item{m}{Parameter of fuzziness (default: 2)}
 \item{vp}{Volume parameter (default: rep(1,k))}
 \item{delta}{Noise distance (default: average Euclidean distance between objects and prototypes from \code{FKM.gk} using the same values of \code{k} and \code{m})}
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
 \item{medoid}{Vector containing the indexes of the medoid objects (\code{NULL} for \code{FKM.gkb.noise})}
 \item{value}{Vector containing the loss function values for the \code{RS} starts}
 \item{criterion}{Vector containing the values of the cluster validity index}
 \item{iter}{Vector containing the numbers of iterations for the \code{RS} starts}
 \item{k}{Number of clusters}
 \item{m}{Parameter of fuzziness}
 \item{ent}{Degree of fuzzy entropy (\code{NULL} for \code{FKM.gkb.noise})}
 \item{b}{Parameter of the polynomial fuzzifier (\code{NULL} for \code{FKM.gkb.noise})}
 \item{vp}{Volume parameter (default: \code{rep(1,max(k)}). If \code{k} is a vector, for each group the first \code{k} element of \code{vp}are considered.}
 \item{delta}{Noise distance}
 \item{gam}{Weighting parameter for the fuzzy covariance matrices}
 \item{mcn}{Maximum condition number for the fuzzy covariance matrices}
 \item{stand}{Standardization (Yes if \code{stand=1}, No if \code{stand=0})}
 \item{Xca}{Data used in the clustering algorithm (standardized data if \code{stand=1})}
 \item{X}{Raw data}
 \item{D}{Dissimilarity matrix (\code{NULL} for \code{FKM.gkb.noise})}
 \item{call}{Matched call}
}
\details{
 If \code{startU} is given, the argument \code{k} is ignored (the number of clusters is \code{ncol(startU)}).
 \cr If \code{startU} is given, the argument \code{RS} is ignored (the algorithm is run using the rational start) and therefore \code{value}, \code{cput} and \code{iter} refer to such a rational start.
 \cr If a cluster covariance matrix becomes singular, then the algorithm stops and the element of \code{value} is \code{NaN}.
}
\references{
 Babuska R., van der Veen P.J., Kaymak U., 2002.  Improved covariance estimation for Gustafson-Kessel clustering. Proceedings of the IEEE International Conference on Fuzzy Systems (FUZZ-IEEE), 1081-1085.
 \cr Dave' R.N., 1991. Characterization and detection of noise in clustering. Pattern Recognition Letters, 12, 657-664.
 \cr Gustafson E.E., Kessel W.C., 1978. Fuzzy clustering with a fuzzy covariance matrix. Proceedings of the IEEE Conference on Decision and Control, pp. 761-766.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{FKM.gk.noise}}, \code{\link{Fclust}}, \code{\link{Fclust.index}}, \code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{unemployment}}}
\examples{
\dontrun{
## unemployment data
data(unemployment)
## Gustafson, Kessel and Babuska-like fuzzy k-means with noise cluster,
##fixing the number of clusters
clust=FKM.gkb.noise(unemployment,k=3,delta=20,RS=10)
## Gustafson, Kessel and Babuska-like fuzzy k-means with noise cluster,
##selecting the number of clusters
clust=FKM.gkb.noise(unemployment,k=2:6,delta=20,RS=10)}
}
\keyword{multivariate}
