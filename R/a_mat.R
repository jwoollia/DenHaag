#' @title Supplies the numerator relationship matrix for a pedigree.
#' @description Does what it says!
#' @param pedf A potentially reduced data frame for the population but must contain "id", "sire" and "dam".
#' @return The numerator relationship matrix for the population.
#' @export
#
a_mat <- function(pedf)
# pedf needs to contain "id","sire","dam"
{
nn <- nrow(pedf)
nrm <- as.double(matrix(nrow=nn,ncol=nn))
nrm <-diag(nn)
for (i in 1:nn)
  {if (!(pedf$sire[i]==0) & !(pedf$dam[i]==0))
    {nrm[i,i]=1.0+nrm[pedf$sire[i],pedf$dam[i]]/2
    }
  for (j in 1:i-1)
    {if (!(pedf$sire[i]==0) & !(pedf$dam[i]==0))
      {nrm[i,j]=nrm[j,pedf$sire[i]]/2 + nrm[j,pedf$dam[i]]/2
      nrm[j,i]=nrm[i,j]
      }
    else if (!(pedf$dam[i]==0))
      {nrm[i,j]=nrm[j,pedf$dam[i]]
      nrm[j,i]=nrm[i,j]
      }
    else if (!(pedf$sire[i]==0))
      {nrm[i,j]=nrm[j,pedf$sire[i]]
      nrm[j,i]=nrm[i,j]
      }
    }
  }
return(nrm)
}
