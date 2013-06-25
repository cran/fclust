\name{SIL} 
\alias{SIL}
\title{Silhouette index}
\description{
 Produces the silhouette index. The optimal number of cluster \emph{k} is achieved when the index value is maximized.
} 
\usage{ 
 SIL (Xca, U)
}
\arguments{
 \item{Xca}{Matrix or data.frame}
 \item{U}{Membership degree matrix}
}
\value{
 \item{sil.obj}{Vector containing the silhouette indices for all the objects}
 \item{sil}{Value of the silhouette index (mean of \code{sil.obj})}
}
\note{
 \code{Xca} should contain the same dataset used in the clustering algorithm, i.e., if the clustering algorithm is run using standardized data, then \code{SIL} should be computed using the same standardized data.
}  
\references{
 Kaufman L., Rousseeuw P.J., 1990. Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York. 
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
sil=SIL(clust$Xca,clust$U)
}
\keyword{multivariate}