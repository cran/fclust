\name{Fclust} 
\alias{Fclust}
\title{Interactive fuzzy clustering analysis}
\description{
 Performs a fuzzy clustering analysis using the algorithms available in the package.
} 
\usage{ 
 Fclust (X)
}
\arguments{
 \item{X}{Matrix or data.frame}
}
\value{
 \item{clust}{Object of class \code{fclust}}
} 
\author{Paolo Giordani, Maria Brigida Ferraro}\seealso{\code{\link{FKM}}, \code{\link{FKM.ent}}, \code{\link{FKM.gk}}, \code{\link{FKM.gk.ent}}, \code{\link{FKM.noise}}, \code{\link{FKM.med}}, \code{\link{Fclust.index}}}
\examples{
\dontrun{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
clust=Fclust(Mc)
}
}
\keyword{multivariate}