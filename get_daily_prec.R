get_daily_prec <- function( mval_prec, mval_wet, set_seed=FALSE ){
    #--------------------------------------------------------------------
    # Distributes monthly total precipitation to days, given number of 
    # monthly wet days. Adopted from LPX.
    #--------------------------------------------------------------------
    ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    ndayyear <- sum(ndaymonth)
    nmonth <- length(ndaymonth)

    c1 <- 1.0
    c2 <- 1.2

    if (set_seed) {set.seed(0)}
    prdaily_random <- array( NA, dim=c(ndayyear,2))
    for (doy in 1:ndayyear){
      prdaily_random[doy,] <- runif(2)
    }

    dval_prec <- rep(NA,ndayyear)
    doy <- 0
    prob <- 0.0
    prob_rain <- rep(NA,nmonth)
    mprecave <- rep(NA,nmonth)
    mprecip <- rep(NA,nmonth)
    for (moy in 1:nmonth){
      prob_rain[moy] <- 0.0
      mprecave[moy] <- 0.0
      mprecip[moy] <- 0.0      
    }
    daysum <- 0

    set.seed( prdaily_random[1,1] * 1e7 )
   
    for (moy in 1:nmonth){
      if ( mval_wet[moy]<=1.0 ) {mval_wet[moy] <- 1.0}
      prob_rain[moy] <- mval_wet[moy] / ndaymonth[moy]
      mprecave[moy] <- mval_prec[moy] / mval_wet[moy]
      dry <- TRUE
      iloop <- 0


      while( dry ){
        iloop <- iloop + 1
        nwet <- 0
        for (dm in 1:ndaymonth[moy]){
          doy <- doy + 1
   
          # Transitional probabilities (Geng et al. 1986)
          if (doy>1) {
            if (dval_prec[doy-1] < 0.1) {
              prob <- 0.75 * prob_rain[moy]
            } else { 
              prob <- 0.25 + (0.75 * prob_rain[moy])
            }
          }        
          # Determine we randomly and use Krysanova / Cramer estimates of 
          # parameter values (c1,c2) for an exponential distribution
          if (iloop==1) { 
            vv <- prdaily_random[doy,1]
          } else {
            # xxx problem: rand() generates a random number that leads to floating point exception
            vv <- runif(1)
          }
      
      
          if (vv>prob) {        
            dval_prec[doy] <- 0.0
          } else {        
            nwet <- nwet + 1        
            v1 <- prdaily_random[doy,2]        
            dval_prec[doy] <- ((-log(v1))^c2) * mprecave[moy] * c1         
            if (dval_prec[doy] < 0.1) dval_prec[doy] <- 0.0        
          }
      
          mprecip[moy] <- mprecip[moy] + dval_prec[doy]
        }
    
        # If it never rained this month and mprec[moy]>0 and mval_wet[moy]>0, do
        # again
        dry <- (nwet==0 && iloop<50 && mval_prec[moy]>0.1)
        if (iloop>50) {
          print('Daily.F, prdaily: Warning stopped after 50 tries in cell')
        }

        # Reset counter to start of month          
        if (dry) {
          doy <- doy - ndaymonth[moy]
        }

      } #while
  
      # normalise generated precipitation by monthly CRU values
      if ( moy > 1 ) {daysum <- daysum + ndaymonth[moy-1]}
      if ( mprecip[moy] < 1.0 ) {mprecip[moy] <- 1.0}
      for (dm in 1:ndaymonth[moy]){
        doy <- daysum + dm
        dval_prec[doy] <- dval_prec[doy] * (mval_prec[moy] / mprecip[moy])
        if ( dval_prec[doy] < 0.1 ) {dval_prec[doy] <- 0.0}
        # dval_prec[doy] <- mval_prec[moy] / ndaymonth[moy]  #no generator
      }
         
      # Alternative: equal distribution of rain for fixed number of wet days
      # prob <- prob_rain[moy] + prob
      # if (prob.ge.1.0) then   
      #   dval_prec[doy] <- mprec[moy]
      #   prob <- prob-1.0
      # } else {
      #   dval_prec[doy] <- 0.0
      #   prob <- prob
      # }
                      
    } 
    
    return( dval_prec )
}

mprec_val <- 50.0
mwetd_val <- 1.0
mprec <- c( mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val, mprec_val)
mwetd <- c( mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val, mwetd_val)
dprec <- get_daily_prec( mprec, mwetd )

# plot( 1:ndayyear, dprec, type="l" )
# print( paste( "required annual precip:", sum(mprec), "generated annual precip", sum(dprec) ) )
# print( paste( "required annual wetdays:", sum(mwetd), "generated annual wetdays", sum(dprec>0.0) ) )

