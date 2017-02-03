#' @title Recommends contributions for achieving the target group coancestry.
#' @description Carries out the standard optimum contribution algorithm, including constraints on totals contributions
#'    of males and females, and non-negative contributions. The approximate solution is obtained by the algorithm of
#'    Meuwissen. It prints out the recommendations after scaling the contributions to give the required total offspring.
#'    The function checks whether the number of parents to be considered exceeds the population, whether there are both
#'    males and female parents available, and whether the group coancestry is achievable. If these checks fail, error
#'    messages are printed.
#' @param ddf The data frame of 9 columns holding the population data.
#' @param tp The number of parents to be considered for contributions, assumed to be the most recent.
#' @param tc The number of offspring to be created.
#' @param gcoa The target group coancestry, which needs to be supplied by the user.
#' @return A logical indicating the success of the algorithm or otherwise.
#' @export
#
oc_sel <- function(ddf,tp,tc,gcoa)
{
abar <- 2.0*gcoa
tn <- nrow(ddf)
if(tn<tp)
  {print("Number of parents requested exceeds population!")
  return(FALSE)
  }
if( !(sum(ddf$sex==1)>0))
  {print("No males among the parents!")
  return(FALSE)
  }
if( !(sum(ddf$sex==2)>0))
  {print("No females among the parents!")
  return(FALSE)
  }
print(paste0("Recommendations for ",tc, " offspring from the most recent ",tp," parents"))
# cand is list of candidates
# q is matrix of sexes of candidates
# ebv is ebv of candidates
# nrm contains their numerator relationships
cand <- ddf[ddf$id>tn-tp,"id"]
q <- matrix(nrow=tp,ncol=2)
q[,1] <- as.double((ddf$sex[cand] < 1.5))
q[,2] <- as.double((ddf$sex[cand] > 1.5))
ebv <- ddf[cand,"ebv"]
nrm <- a_mat(ddf[c(1,6,7)])
nrm <- nrm[cand,cand]
ind <- c(1:tp)
#
nn <- 0
tt <- tp
while (!(nn==tt))
  {nn <- tt
  ai <- solve(nrm[ind,ind])
  qq <- q[ind,]
  ebvv <- ebv[ind]
  qaqi <- solve(t(qq)%*%ai%*%qq)
  ter1 <- t(ebvv)%*% (ai-ai%*%qq%*%qaqi%*%t(qq)%*%ai) %*%ebvv
  ter2 <- 4*abar-sum(sum(qaqi))
  if(ter2 <= 0)
    {
    print("Constraint cannot be achieved: increase group coancestry")
    return(FALSE)
    }
  lgm0 <- sqrt(ter1/ter2)
  lgm1 <- qaqi%*% (t(qq)%*%ai%*%ebvv-c(1,1)*lgm0)
  cc <- ai%*% (ebvv-qq%*%lgm1)
  cc <- cc /matrix(2.0*lgm0,nrow=nrow(cc),ncol=ncol(cc))
  # check for negative c elements
  tt <- 0
  for (ii in 1:nn)
    {if (cc[ii]>=0)
      {tt <- tt+1
      ind[tt] <- ind[ii]
      }
    }
  ind <- ind[1:tt]
  }
ccdf <- data.frame("id"=cand[ind],"sex"=ddf$sex[cand[ind]],"ebv"=ebv[ind],"c"=cc,"noff"= 2.0*tc*cc)
print(ccdf,row.names=FALSE)
return(TRUE)
}
