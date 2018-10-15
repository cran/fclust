\name{Fclust.compare}
\alias{Fclust.compare}
\title{
Similarity between partitions}
\description{
Performs some measures of similarity between a hard (reference) partition and a fuzzy partition.}
\usage{
Fclust.compare(partHard, partFuzzy, index, tnorm)
}
\arguments{
\item{partHard}{(Reference) partition matrix or data.frame}
\item{partFuzzy}{Fuzzy membership degree matrix or data.frame}
\item{index}{Measures of similarity: "ARI.F" (fuzzy version of the adjuster Rand index), "RI.F" (fuzzy version of the Rand index), "JACCARD.F" (fuzzy version of the Jaccard index), "ALL" for all the indices (default: "ALL")}
\item{tnorm}{Type of the triangular norm: "minimum" (minimum triangular norm), "triangular product" (product norm) (default: "minimum")}
}
\value{
 \code{out.index}{Vector containing the similarity measures}
}
\details{
\code{partHard} and \code{partFuzzy} must have the same number of rows, but, not necessarily, the same number of columns (i.e., the same number of clusters in the two partitions).
\cr \code{index} is not case-sensitive.
\cr All the measures of similarity share the same properties of their non-fuzzy counterpart.
}
\references{
Campello, R.J., 2007. A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment. Pattern Recognition Letters, 28, 833-841.
\cr Hubert, L., Arabie, P., 1985. Comparing partitions. Journal of Classification, 2, 193-218.
\cr Jaccard, P., 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles, 37, 547-579.
\cr Rand, W.M., 1971. Objective criteria for the evaluation of clustering methods. Journal of the American Statistical Association, 66, 846-850.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{
\code{\link{RI.F}}, \code{\link{ARI.F}}, \code{\link{JACCARD.F}}
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
## all measures of similarity
classTrue=unmap(classification=Mc$Type)
classEst=unmap(classification=clust$clus[,1])
all.indexes=Fclust.compare(partHard=classTrue,partFuzzy=classEst)
## fuzzy adjusted Rand index
Fari.index=Fclust.compare(partHard=classTrue,partFuzzy=classEst,index="ARI.F")
}
}
\keyword{multivariate}
