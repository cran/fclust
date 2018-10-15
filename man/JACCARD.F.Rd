\name{JACCARD.F}
\alias{JACCARD.F}
\title{
Fuzzy Jaccard index}
\description{
Produces the fuzzy version of the Jaccard index between a hard (reference) partition and a fuzzy partition.
}
\usage{
JACCARD.F(partHard, partFuzzy, t_norm)
}
\arguments{
\item{partHard}{(Reference) partition matrix or data.frame}
\item{partFuzzy}{Fuzzy membership degree matrix or data.frame}
\item{t_norm}{Type of the triangular norm: "minimum" (minimum triangular norm), "triangular product" (product norm) (default: "minimum")}
}
\value{
\code{jaccard.f}{Value of the fuzzy Jaccard index}
}
\details{
\code{partHard} and \code{partFuzzy} must have the same number of rows, but, not necessarily, the same number of columns (i.e., the same number of clusters in the two partitions).
}
\references{
Campello, R.J., 2007. A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment. Pattern Recognition Letters, 28, 833-841.
\cr Jaccard, P., 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles, 37, 547-579.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{
\code{\link{ARI.F}}, \code{\link{RI.F}}, \code{\link{Fclust.compare}}
}
\examples{
\dontrun{
require("mclust")
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
## fuzzy k-means
## (excluded the factor column Type (last column))
clust=FKM(Mc[,1:(ncol(Mc)-1)],k=6,m=1.5,stand=1)
## fuzzy Jaccard index
classTrue=unmap(classification=Mc$Type)
classEst=unmap(classification=clust$clus[,1])
jaccard.f=JACCARD.F(partHard=classTrue,partFuzzy=classEst)
}
}
\keyword{multivariate}
