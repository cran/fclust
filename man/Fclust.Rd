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
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{print.fclust}}, \code{\link{summary.fclust}}, \code{\link{plot.fclust}}, \code{\link{FKM}}, \code{\link{FKM.ent}}, \code{\link{FKM.gk}}, \code{\link{FKM.gk.ent}}, \code{\link{FKM.gkb}}, \code{\link{FKM.gkb.ent}}, \code{\link{FKM.med}}, \code{\link{FKM.pf}}, \code{\link{FKM.noise}}, \code{\link{FKM.ent.noise}}, \code{\link{FKM.gk.noise}}, \code{\link{FKM.gkb.ent.noise}},  \code{\link{FKM.gkb.noise}}, \code{\link{FKM.gk.ent.noise}},\code{\link{FKM.med.noise}}, \code{\link{FKM.pf.noise}}, \code{\link{Fclust.index}}, \code{\link{Mc}}}
\examples{
\dontrun{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
## interactive fuzzy clustering
## (excluded the factor column Type (last column))
##clust=Fclust(Mc[,1:(ncol(Mc)-1)])}
}
\keyword{multivariate}
