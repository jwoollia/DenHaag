#' @title Constructs the inverse of the numerator relationship matrix.
#' @description Does what it says!
#' @param pedf A potentially reduced data frame for the population but must contain "id", "sire" and "dam".
#' @return The inverse of the numerator relationship matrix.
#' @export
#
a_inv <- function(pedf)
# pedf needs to contain "id","sire","dam","f"
{
nn <- nrow(pedf)
ainv <- as.double(matrix(nrow=nn, ncol=nn))
ainv <- diag(nn)
ahlp <- matrix(c(0.5,0.5,-1,0.5,0.5,-1,-1,-1,2), nrow=3, ncol=3)
#
for (i in 1:nn)
  {if (pedf$sire[i]>0)
    {ainv[i,i] = 0
    ihlp <- c(pedf$sire[i],pedf$dam[i],i)
    delta = 1 - (pedf$f[pedf$sire[i]]+pedf$f[pedf$dam[i]])/2
    ainv[ihlp,ihlp] <- ainv[ihlp,ihlp]+ahlp/delta
    }
  }
return(ainv)
}
