\name{Mc}
\alias{Mc}
\docType{data}
\title{McDonald's data}
\description{
 We provide a Nutrition analysis of McDonald's menu items.
}
\usage{data(Mc)}
\format{
 A data.frame with 81 rows and 16 columns.
} 
\details{
Data are from McDonald's USA Nutrition Facts for Popular Menu Items. A subset of menu items is reported. Beverages are excluded.
In case of duplications, regular size or medium size information is reported. The variable Type is a factor the levels of which specify the kind of the menu items.
Although some menu items could be well described by more than one level, only one level of the variable Type specifies each menu item. 
Percent Daily Values (\%DV) are based on a 2,000 calorie diet. Some menu items are registered trademarks.
}
\author{Paolo Giordani, Maria Brigida Ferraro}
\examples{
data(Mc)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
Mc=Mc[,-1]
p=(ncol(Mc)-1)
## FKM
clust.FKM=FKM(Mc[,1:p],k=6,m=1.5,stand=1)
Mc[,ncol(Mc)+1]=factor(clust.FKM$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM")
levels(Mc$Cluster.FKM)=paste("Clus FKM",1:clust.FKM$k,sep=" ")
table(Mc$Type,Mc$Cluster.FKM)
clust.FKM$Hraw=Hraw(clust.FKM$X,clust.FKM$H)
clust.FKM$Hraw
\dontrun{
## It may take more than a few seconds
## FKM.ent
clust.FKM.ent=FKM.ent(Mc[,1:p],k=6,ent=3,RS=10,stand=1)
Mc[,ncol(Mc)+1]=factor(clust.FKM.ent$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM.ent")
levels(Mc$Cluster.FKM.ent)=paste("Clus FKM.ent",1:clust.FKM.ent$k,sep=" ")
table(Mc$Type,Mc$Cluster.FKM.ent)
clust.FKM.ent$Hraw=Hraw(clust.FKM.ent$X,clust.FKM.ent$H)
clust.FKM.ent$Hraw
## FKM.med
clust.FKM.med=FKM.med(Mc[,1:p],k=6,m=1.1,RS=10,stand=1)
Mc[,ncol(Mc)+1]=factor(clust.FKM.med$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM.med")
levels(Mc$Cluster.FKM.med)=paste("Clus FKM.med",1:clust.FKM.med$k,sep=" ")
table(Mc$Type,Mc$Cluster.FKM.med)
clust.FKM.med$Hraw=Hraw(clust.FKM.med$X,clust.FKM.med$H)
clust.FKM.med$Hraw
## or, equivalently,
Mc[clust.FKM.med$medoid,1:p]
}
}
\keyword{data}
\keyword{multivariate}
