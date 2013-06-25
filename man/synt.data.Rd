\name{synt.data}
\alias{synt.data}
\docType{data}
\title{Synthetic data}
\description{
  Synthetic dataset with 2 non-spherical clusters.
}
\usage{data(synt.data)}
\format{
 A matrix with 302 rows and 2 columns.
} 
\details{
Although two clusters are clearly visible, FKM fails to discover them. 
The Gustafson - Kessel variant of FKM should be used for finding the known-in-advance clusters.}
\author{Paolo Giordani, Maria Brigida Ferraro}
\examples{
data(synt.data)
plot(synt.data)
cor(synt.data)
syntFKM=FKM(synt.data)
syntFKM.gk=FKM.gk(synt.data)
par(mfcol = c(2,1))
plot(synt.data,col=syntFKM$clus[,1])
plot(synt.data,col=syntFKM.gk$clus[,1])
}
\keyword{data}
\keyword{multivariate}