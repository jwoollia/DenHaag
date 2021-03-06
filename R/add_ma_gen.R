#' @title Adds a generation of offspring after selecting mates by maximum avoidance (i.e. minimum coancestry)
#' @description The maximum avoidance is found by a call to amat and using it to evaluate 1000 random swaps. The average
#'    inbreeding of the offspring and the group coancestry are calculated and printed, together with the deviation from
#'    random mating. In creating the offspring, the function makes calls to the blup function for calculating the EBV, and
#'    the a_mat function for calculating the inbreeding coefficient from the leading diagonal.
#' @param ddf A data frame with the 9 columns as generated by new_pop
#' @param hh Heritability of the trait in the data frame
#' @return An data frame of the same form as the input but extended by a number of offspring deteremined by noff.
#'    If noff has unequal numbers of offspring from male and female parents, or gives numbers of offspring are not even,
#'    or the heritability is out of bounds then errors are printed and the input data frame is returned unchanged.
#' @export
#' @importFrom stats rnorm
#' @importFrom stats runif
#
add_ma_gen <- function(ddf,hh)
{
#
# errors in setting noff
#
if(!(sum((ddf$sex==1)*ddf$noff)==sum((ddf$sex==2)*ddf$noff)))
  {print("Numbers of male and female parents unequal in matings!")
  return(ddf)
  }
toff <- as.integer(sum(ddf$noff)/2)
if(toff<=0)
  {print("Offspring numbers are <= zero!")
  return(ddf)
  }
if(toff%%2==1)
  {print(paste0("Offspring numbers are expected to be even not odd! ",toff,"!"))
  return(ddf)
  }
if( (hh >= 1)|(hh <= 0) )
  {print(paste0("Heritability out of bounds! ",hh,"!"))
  return(ddf)
  }
print(paste0("Creating ",toff," offspring with maximum avoidance"))
#
# group coancestry
#
cc <- ddf$noff/(2.0*toff)
nrm <- a_mat(ddf[,c(1,6,7)])
gc <- 0.5*(t(cc)%*%nrm%*%cc)
#
# setting up a random sire list from noff
#
isir <- numeric(0)
while (max(ddf$noff[ddf$sex==1])>0)
  {isir <- append(isir,ddf$id[ddf$noff>0 & ddf$sex==1],after=length(isir))
  ddf$noff[ddf$noff>0 & ddf$sex==1] <- ddf$noff[ddf$noff>0 & ddf$sex==1]-1
  }
isir <- sample(isir)
#
# setting up a random dam list from noff
#
idam <- numeric(0)
while (max(ddf$noff[ddf$sex==2])>0)
  {idam <- append(idam,ddf$id[ddf$noff>0 & ddf$sex==2],after=length(idam))
  ddf$noff[ddf$noff>0 & ddf$sex==2] <- ddf$noff[ddf$noff>0 & ddf$sex==2]-1
  }
idam <- sample(idam)
#
# sampling for maximum avoidance
#
rnu <- stats::runif(2000,min=0,max=1)
rnu <- floor(rnu*toff)+1
rnu <- matrix(rnu,nrow=1000,ncol=2)
for (i in 1:1000)
  {js1 <- isir[rnu[i,1]]
  js2 <- isir[rnu[i,2]]
  jd1 <- idam[rnu[i,1]]
  jd2 <- idam[rnu[i,2]]
  # current mating
  v1 <- nrm[js1,jd1]+nrm[js2,jd2]
  # alternative
  v2 <- nrm[js1,jd2]+nrm[js2,jd1]
  if(v2<v1)
    # then swap
    {idam[rnu[i,1]] <- jd2
    idam[rnu[i,2]] <- jd1
    }
  }
avf <- 0
for (i in 1:length(isir))
  {avf <- avf + nrm[isir[i],idam[i]]
  }
avf <- avf/(2.0*toff)
alph <- 1 - (1-avf)/(1-gc)
print(paste0("Group coancestry ",round(gc,digits=4)," and average offspring F ",
             round(avf,digits=4)," with alpha ",round(alph,digits=4)))
#
# creating the offspring
#
dum <- data.frame(id=integer(toff),
                  sex=integer(toff),
                  noff=integer(toff),
                  ptype=double(toff),
                  ebv=double(toff),
                  sire=integer(toff),
                  dam=integer(toff),
                  tbv=double(toff),
                  f=double(toff))
e <- stats::rnorm(n=toff,mean=0,sd=1)
ms <- stats::rnorm(n=toff,mean=0,sd=1)
ms <- sqrt(hh/2)*ms*sqrt(c(rep(1, times=toff))-(ddf$f[isir]+ddf$f[idam])/2)
dum$id <- c(1:toff) + nrow(ddf)
dum$sex <- c(rep(c(1,2), each=toff/2))
dum$noff <- c(rep(0, times=toff))
dum$sire <- isir
dum$dam <- idam
dum$tbv <- (ddf$tbv[isir]+ddf$tbv[idam])/2 + ms
dum$ptype <- dum$tbv+sqrt(1-hh)*e
dum$ebv <- numeric(toff)
dum$f <- numeric(toff)
ddf <- rbind(ddf,dum)
#
# ebvs and fs
#
ddf$ebv <- blup(ddf,hh)
ddf$f <- diag(a_mat(ddf[,c(1,6,7)]))-c(rep(1,times=nrow(ddf)))
print(ddf, row.names=FALSE)
return(ddf)
}
