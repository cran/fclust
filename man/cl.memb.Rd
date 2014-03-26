\name{cl.memb} 
\alias{cl.memb}
\title{Cluster membership}
\description{
 Produces a summary of the membership degree information.
} 
\usage{ 
 cl.memb (U)
}
\arguments{
 \item{U}{Membership degree matrix}
}
\value{
 \item{info.U}{Matrix containing the indices of the clusters where the objects are assigned (row 1) and the associated membership degrees (row 2)}
} 
\note{
 An object is assigned to a cluster according to the maximal membership degree.
} 
\author{Paolo Giordani, Maria Brigida Ferraro}
\seealso{\code{\link{cl.memb.H}}, \code{\link{cl.memb.t}}}
\examples{
n=20
k=3
U=matrix(runif(n*k,0,1), nrow=n, ncol=k)
U=U/apply(U,1,sum)
info.U=cl.memb(U)
}
\keyword{multivariate}