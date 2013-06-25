\name{butterfly}
\alias{butterfly}
\docType{data}
\title{Butterfly data}
\description{
  Synthetic dataset with 2 clusters and some outliers.
}
\usage{data(butterfly)}
\format{
 A matrix with 17 rows and 2 columns.
} 
\details{
FKM with noise cluster should be used for finding outliers. In fact, differently from FKM, the membership degrees of the outliers are low for all the clusters.}
\author{Paolo Giordani, Maria Brigida Ferraro}\examples{
data(butterfly)
plot(butterfly,type='n')
text(butterfly[,1],butterfly[,2],labels=rownames(butterfly),cex=0.7,lwd=2)
round(FKM(butterfly)$U,2)
round(FKM.noise(butterfly,delta=3)$U,2)
}
\keyword{data}
\keyword{multivariate}