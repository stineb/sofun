daily2monthly <- function( dvals ){
  ## /////////////////////////////////////
  ## get monthly values
  ## dvals is vector of daily values
  ## -------------------------------------

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  ndayyear <- sum(ndaymonth)
  nmonth <- length(ndaymonth)

  ## can only treat even years, i.e. 365 days in a year
  if (length(dvals)%%ndayyear!=0) {print("length problem");stop}
  nyears <- length(dvals)/ndayyear

  mvals <- rep( 0.0, nmonth*nyears )
  day <- 0
  for (yr in 1:nyears){
    for (moy in 1:nmonth){
      for (dm in 1:ndaymonth[moy]){
        day <- day+1 
        mvals[moy] <- mvals[moy] + dvals[day]/ndaymonth[moy]
      }
    }
  }
  return(mvals)
}

monthly2daily <- function( mvals, method, mvals_pvy=NA, mvals_nxy=NA ) {
  ##/////////////////////////////////////////////////////////////////////////
  ## Returns daily values based on monthly values, using a defined method.
  ##-------------------------------------------------------------------------
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350,381) # day of year of middle-month-day
  ndayyear <- sum(ndaymonth)

  dvals <- rep(NA,ndayyear)

  if (method=="interpol"){
    ##--------------------------------------------------------------------
    ## LINEAR INTERPOLATION
    ## of monthly to quasi-daily values.
    ## If optional argument 'mvals_pvy' is provided, take December-value
    ## of previous year to interpolate to first 15 days of January,
    ## otherwise, use the same year's December value to get first 15 days.
    ##--------------------------------------------------------------------
    doy <- 0
    for (moy in 1:nmonth) {
      for (dm in 1:ndaymonth[moy]){
        doy <- doy+1
        if (doy > middaymonth[moy]) {
          dvals[doy] <- mvals[moy] + (doy-middaymonth[moy])/ndaymonth[moy] * (mvals[moy+1]-mvals[moy])
        } else if (doy < middaymonth[moy]) {
          if (moy==1) {
            if (is.na(mvals_pvy)) {
              dvals[doy] <- mvals_pvy[nmonth] + (middaymonth[nmonth]+doy)/ndaymonth[nmonth] * (mvals[moy]-mvals_pvy[nmonth])
            } else {
              dvals[doy] <- mvals[nmonth] + (middaymonth[nmonth]+doy)/ndaymonth[nmonth] * (mvals[moy]-mvals[nmonth])
            }
          } else {
            dvals[doy] <- mvals[moy-1] + (doy-middaymonth[moy-1])/ndaymonth[moy-1] * (mvals[moy]-mvals[moy-1])
          }
        } else {
          dvals[doy] <- mvals[moy]
        }          
      }
    }

  } else if (method=="polynom") {
    ##--------------------------------------------------------------------
    ## In addition to tempdaily daily values are calculated using a polynom of second
    ## order through the middpoints between months. Additionally, average of daily 
    ## values is identical to the monthly input data. That's crucial for modelling
    ## soil heat diffusion and freezing/thawing processes. 
    ##--------------------------------------------------------------------
    
    ## Starting conditons of december in previous year
    startt <- -30.5               ## midpoint between Nov-Dec of previous year
    endt <- 0.5                   ## midpoint between Dec-Jan of this year
    dt <- ndaymonth[nmonth]       ## number of Dec days
    if (is.na(mvals_pvy)){
      lastmonthtemp <- mvals[nmonth]
    } else{
      lastmonthtemp <- mvals_pvy[nmonth]
    }
    doy <- 0                      ## initialisation of this years days
    
    for (month in 1:nmonth) {
      dtold <- dt
      startt <- endt
      endt <- endt + dt
      if (month<nmonth) {
        dtnew <- ndaymonth[month+1]
        nextmonthtemp <- mvals[month+1]
      } else {
        dtnew <- ndaymonth[1]
        if (is.na(mvals_nxy)){
          nextmonthtemp <- mvals[1]          
        } else {
          nextmonthtemp <- mvals_nxy[1]
        }
      }

      starttemp <- (mvals[month]*dt+lastmonthtemp*dtold)/(dt+dtold)
      endtemp <- (nextmonthtemp*dtnew+mvals[month]*dt)/(dtnew+dt)
      deltatemp <- endtemp-starttemp
      
      ## Calculate vars for a,b,c coefficients in polynom y <- ax^2 +bx + c
      d2t <- endt**2.0 - startt**2.0
      d3t <- endt**3.0 - startt**3.0

      ## Take a sheet of paper and try solve the polynom, well here is the outcome
      polya <- (mvals[month]*dt - deltatemp*d2t/dt/2.0 - starttemp*dt + deltatemp*startt) / (d3t/3.0 - d2t**2.0/dt/2.0 - dt*startt**2.0 + startt*d2t)
      polyb <- deltatemp/dt - polya*(startt+endt)
      polyc <- starttemp - polya*startt**2.0 - polyb*startt

      ## Calculate daily values with the polynom function
      for (dm in 1:ndaymonth[month]){
        doy <- doy + 1
        dvals[doy] <- polya * doy**2.0 + polyb * doy + polyc        
      }
      lastmonthtemp <- mvals[month]
    }

    ## Calculate monthly means after interpolation - not absolutely identical to input
    doy<-0
    mtempint <- rep(NA,nmonth)
    for (m in 1:nmonth) {
      mtempint[m] <- 0.0
      for (dm in 1:ndaymonth[m]){
        doy <- doy + 1
        mtempint[m] <- mtempint[m]+dvals[doy]/ndaymonth[m]        
      }
    }

  } else if (method=="step") {
    ##--------------------------------------------------------------------
    ## STEP-WISE FILLING
    ## Use monthly values for all days in repsective month.
    ##--------------------------------------------------------------------
    doy <- 0
    for (moy in 1:nmonth) {
      for (dm in 1:ndaymonth[moy]){
        doy <- doy+1
        dvals[doy] <- mvals[moy]
      }
    }
  
  }
  return(dvals)
}

## read daily values
dtemp <- read.table( '/alphadata01/bstocker/sofun/trunk/input/dtemp_CH-Oe1_2000_365.txt' )[,1]

## get monthly values as a mean over days in respective month
mtemp <- daily2monthly( dtemp )

## expand to daily from monthly
dtemp_interpol <- monthly2daily( mtemp, "interpol" )
dtemp_polynom <- monthly2daily( mtemp, "polynom" )

## polynom method returns daily values with which monthly means should be identical to original monthly values
mtemp_polynom <- daily2monthly(dtemp_polynom)


plot( 1:length(dtemp), dtemp, type="l" )
lines( middaymonth, mtemp, col="red" )
# lines( 1:length(dtemp), dtemp_interpol, col="blue" )
# lines( 1:length(dtemp), dtemp_polynom, col="green" )
lines( middaymonth, mtemp_polynom, col="magenta")
