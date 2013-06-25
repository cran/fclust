\name{Fclust.index} 
\alias{Fclust.index}
\title{Clustering indices}
\description{
 Performs some clustering indices for choosing the optimal number of cluster \emph{k}.
} 
\usage{ 
 Fclust.index (fclust.obj, index, alpha)
}
\arguments{
 \item{fclust.obj}{Object of class \code{fclust}}
 \item{index}{Clustering indices: "PC" (partition coefficient), "PE" (partition entropy), "MPC" (modified partition coefficient), "SIL" (silhouette), "SIL.F" (fuzzy silhouette), "XB" (Xie and Beni index), "ALL" for all the indices (default: "ALL")}
 \item{alpha}{Weighting coefficient for the fuzzy silhouette index \code{SIL.F} (default: 1)}
}
\value{
 \item{out.index}{Vector containing the index values}
} 
\note{
 \code{index} is not case-sensitive.
} 
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
out.index=Fclust.index(clust)
}
\keyword{multivariate}