############################################################################
## Returns a monthly climatology from monthly data covering multiple years
## Beni Stocker, 16.06.2013

mean.bymonth <- function(vec) {

  ## Calculate mean annual cycle by months
  nmonths <- 12
  vecm <- rep(NA,nmonths)
  ismonth <- rep(FALSE,nmonths)
  for (m in seq(nmonths)){
    ismonth[]  <- FALSE
    ismonth[m] <- TRUE
    mask <- rep(ismonth,length(vec)/nmonths)
    nadd <- length(vec)-length(mask)
    if (nadd!=0) {
      mask <- c(mask,mask[1:nadd])
    }
    vecm[m] <- mean( vec[mask], na.rm=TRUE )
  }
  return(vecm)
  
}
