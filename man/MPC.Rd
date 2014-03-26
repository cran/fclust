\name{MPC} 
\alias{MPC}
\title{Modified partition coefficient}
\description{
 Produces the modified partition coefficient index. The optimal number of cluster \emph{k} is achieved when the index value is maximized.
} 
\usage{ 
 MPC (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{mpc}{Value of the modified partition coefficient index}
} 
\references{
 Dave' R.N., 1996. Validating fuzzy partitions obtained through \emph{c}-shells clustering. Pattern Recognition Letters, 17, 613-623.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{PC}}, \code{\link{PE}}, \code{\link{SIL}}, \code{\link{SIL.F}}, \code{\link{XB}}, \code{\link{Fclust}}}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
mpc=MPC(clust$U)
}
\keyword{multivariate}