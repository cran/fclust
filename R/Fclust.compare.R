Fclust.compare <- function(partHard, partFuzzy, index, tnorm = c("minimum","product"))
{

  if (missing(partHard))
    stop("The hard partitions partHard must be given")
  if (missing(partFuzzy))
    stop("The fuzzy partitions partFuzzy must be given")

  if (is.null(partHard))
    stop("The hard partitions partHard is empty")

  if (is.null(partFuzzy))
    stop("The fuzzy partitions partFuzzy is empty")

  partHard = as.matrix(partHard)
  partFuzzy = as.matrix(partFuzzy)

  if((dim(partHard)[1] != dim(partFuzzy)[1]))
    stop("partHard and partFuzzy must have the same number of observations")

  lg=length(c("ARI.F","RI.F","JACCARD.F","ALL"))

  if (missing(index))
  {
    indexN=lg
  }
  else
  {
    if (is.null(index))
    {
      indexN=lg
    }
    else
    {
      indexN=match(toupper(index),c("ARI.F","RI.F","JACCARD.F","ALL"))
    }
  }
  if (any(is.na(indexN)))
  {
    indexN=lg
    cat("(At least one) no match is found: all the indices will be computed ",fill=TRUE)
  }
  index=indexN
  out.index=rep(0,lg-1)

  tnorm <- match.arg(tnorm, choices = eval(formals(Fclust.compare)$tnorm))

  partHard  = as.matrix(partHard)
  partFuzzy = as.matrix(partFuzzy)


  out = partition_comp(HardClust = partHard,Fuzzy = partFuzzy, t_norm = tnorm)

  if (any(index==1) || any(index==lg))
    out.index[1]=out$adjRand.F
  if (any(index==2) || any(index==lg))
    out.index[2]=out$Rand.F
  if (any(index==3) || any(index==lg))
    out.index[3]=out$Jaccard.F

  names(out.index)=c("ARI.F","RI.F","JACCARD.F")
  if (max(index)<lg)
    out.index=out.index[index]
  return(out.index)

}
