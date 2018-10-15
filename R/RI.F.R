RI.F <- function(partHard,partFuzzy, t_norm = c("minimum","product"))
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

  if(any(dim(partHard) != dim(partFuzzy)))
    stop("partHard and partFuzzy must be matrix with the same dimension")

  t_norm <- match.arg(t_norm, choices = eval(formals(RI.F)$t_norm))
  out = partition_comp(HardClust = partHard,Fuzzy = partFuzzy, t_norm = t_norm)
  return(out$Rand.F)

}
