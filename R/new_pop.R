#' @title Sets up a base population for selection
#' @description The first step in the simulation.
#' @param nn Size of base, an even integer
#' @param hh Heritability of the simulated trait
#' @return A data frame of 9 labelled columns of integers and doubles. Errors are printed if the size
#'    of the base is not an even integer, and the heritability is out of bounds.
#' @export
#' @importFrom stats rnorm
#' @importFrom stats runif
#
new_pop <- function(nn,hh)
{
#
# errors in setting nn
#
if(nn<=0)
  {print(paste0("Population are extinct!",nn," <= 0!"))
  return()
  }
if(nn%%2==1)
  {print(paste0("Census size is expected to be even not odd! ",nn,"!"))
  return()
  }
if( (hh >= 1)|(hh <= 0) )
  {print(paste0("Heritability out of bounds! ",hh,"!"))
  return()
  }
print(paste0("Creating ",nn," individuals for base generation"))
burn_in <- stats::runif(n=1000)
#
xe <- stats::rnorm(n=nn,mean=0,sd=1)
xbv <- stats::rnorm(n=nn,mean=0,sd=1)
xbv <- sqrt(hh)*xbv
xp <- xbv+sqrt(1-hh)*xe
xebv <- hh*xp
jsir <- integer(nn)
jdam <- integer(nn)
xf <- numeric(nn)
jno <- numeric(nn)
jid <- c(1:nn)
jsex <- c(rep(c(1,2),each=nn/2))
ddf <- data.frame("id"=as.integer(jid),
                  "sex"=as.integer(jsex),
                  "noff"=as.integer(jno),
                  "ptype"=as.double(xp),
                  "ebv"=as.double(xebv),
                  "sire"=as.integer(jsir),
                  "dam"=as.integer(jdam),
                  "tbv"=as.double(xbv),
                  "f"=as.double(xf))
print(ddf,row.names=FALSE)
return(ddf)
}
