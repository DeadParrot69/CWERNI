#' Generate neutral community
#'
#' @param theta a number
#' @param Ivec  a vector
#' @param Jvec  a vector
#'
#' @return a vector of species abundances
#' @export
#'
#' @examples
#' generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000))

generate_spat_abund = function(theta,Ivec,Jvec)
{
  numsam = length(Jvec)
  J = sum(Jvec)
  locspecnum = matrix(0,nrow = numsam,ncol = J)
  mcspecnum = rep(0,J)
  abund = matrix(0,nrow = numsam,ncol = J)
  k = 0
  n = 0
  for(i in 1:numsam)
  {
    for(j in 1:Jvec[i])
    {
      bnd = Ivec[i]/(Ivec[i] + j - 1)
      if(stats::runif(1) > bnd)
      {
        locspecnum[i,j] = locspecnum[i,sample(1:j,1)]
      } else
      {
        k = k + 1
        if(stats::runif(1) <= theta/(theta + k - 1))
        {
          n = n + 1
          mcspecnum[k] = n
        } else
        {
          mcspecnum[k] = mcspecnum[sample(1:k,1)]
        }
        locspecnum[i,j] = mcspecnum[k]
      }
      abund[i,locspecnum[i,j]] = abund[i,locspecnum[i,j]] + 1
    }
  }
  zeros = 1
  numspec = J
  while(zeros == 1)
  {
    sumsites = sum(abund[1:numsam,numspec])
    if(sumsites == 0)
    {
      numspec = numspec - 1
    } else
    {
      zeros = 0
    }
  }
  abund = abund[1:numsam,1:numspec]
  return(abund)
}
