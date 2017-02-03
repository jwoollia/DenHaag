#' @title Calculates BLUP EBVs
#' @description Calculates BLUP EBVs after constructing the inverse of the numerator relationship matrix.
#' @param ddf The standard data frame of 9 columns containing population data.
#' @param hh The heritability of the trait.
#' @return A vector of EBVs
#' @export
#
blup <- function(ddf,hh)
{
nn <- nrow(ddf)
ainv <- a_inv(ddf[,c("id","sire","dam","f")])
alfa <- (1-hh)/hh
mme <- as.double(diag(nn)) + alfa*ainv
ddf$ebv <- solve(mme,ddf$ptype)
return(ddf$ebv)
}
