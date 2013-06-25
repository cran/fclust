\name{PC} 
\alias{PC}
\title{Partition coefficient}
\description{
 Produces the partition coefficient index. The optimal number of cluster \emph{k} is achieved when the index value is maximized.
} 
\usage{ 
 PC (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{pc}{Value of the partition coefficient index}
} 
\references{
 Bezdek J.C., 1974. Cluster validity with fuzzy sets. Journal of Cybernetics, 3, 58-73.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{PE}}, \code{\link{MPC}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
pc=PC(clust$U)
}
\keyword{multivariate}