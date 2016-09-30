calc_vpd <- function( vap, temp, tmin=NA, tmax=NA ){
  ##-----------------------------------------------------------------------
  ## Output:   mean monthly vapor pressure deficit, Pa (vpd)
  ## Features: Returns mean monthly vapor pressure deficit
  ## Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure 
  ##           Calculation Methods, in Evaporation and Evapotranspiration: 
  ##           Measurements and Estimations, Springer, London.
  ##             vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
  ##             where:
  ##                 tc = average daily air temperature, deg C
  ##                 vap  = actual vapor pressure, kPa
  ##-----------------------------------------------------------------------
  ## arguments
  ## tc    ## mean monthly temperature, deg C
  ## vap   ## mean monthly vapor pressure, hPa (because CRU data is in hPa)
  ## tmin  ## (optional) mean monthly min daily air temp, deg C 
  ## tmax  ## (optional) mean monthly max daily air temp, deg C 
  ##
  ## function return variable
  ## vpd   ##  mean monthly vapor pressure deficit, Pa  
  ##-----------------------------------------------------------------------

  if ( !is.na(tmin) && !is.na(tmax) ) {
    my_tc <- 0.5 * (tmin + tmax)
  } else {
    my_tc <- tc      
  }

  ## calculate VPD in units of kPa
  vpd <- ( 0.611 * exp( (17.27 * my_tc)/(my_tc + 237.3) ) - 0.10 * vap )    

  ## this empirical equation may lead to negative values for VPD (happens very rarely). assume positive...
  vpd <- max( 0.0, vpd )

  ## convert to Pa
  vpd <- vpd * 1.0e3

  return( vpd )

}