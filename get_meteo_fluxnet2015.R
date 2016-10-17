get_meteo_fluxnet2015 <- function( path ){
  ##--------------------------------------------------------------------
  ## Function returns a dataframe containing all the data of flux-derived
  ## GPP for the station implicitly given by path (argument).
  ## Specific for FLUXNET 2015 data
  ##--------------------------------------------------------------------
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  ndayyear <- sum(ndaymonth)
  nmonth <- length(ndaymonth)

  ## from flux to energy conversion, umol/J (Meek et al., 1984), same as used in SPLASH (see Eq.50 in spash_doc.pdf)
  kfFEC <- 2.04

  ## get daily meteo data
  meteo <- read.csv( path, na.strings="-9999" )  

  ## add three columns: year, month, day in month
  meteo$year <- as.numeric(substr(meteo$TIMESTAMP,start=1,stop=4))
  meteo$moy  <- as.numeric(substr(meteo$TIMESTAMP,start=5,stop=6))
  meteo$dom  <- as.numeric(substr(meteo$TIMESTAMP,start=7,stop=8))

  meteo$year_dec <- meteo$year + (meteo$moy-1)/nmonth + (meteo$dom-1)/sum(ndaymonth)

  ## rename variables (columns)
  if (!is.null(meteo$NETRAD)){
    meteo <- rename( meteo, c( "TA_F"="temp", "VPD_F"="vpd", "P_F"="prec", "NETRAD"="nrad" ) ) #, "LW_IN_F"="lwin", "SW_IN_F"="swin" ) )    
    meteo$nrad <- meteo$nrad * 60 * 60 * 24  # given in W m-2 (avg.), required in J m-2 (daily total)
  } else {
    meteo <- rename( meteo, c( "TA_F"="temp", "VPD_F"="vpd", "P_F"="prec" ) ) #, "LW_IN_F"="lwin", "SW_IN_F"="swin" ) )   
    meteo$nrad <- rep( NA, dim(meteo)[1] ) 
  }

  ## Take only net PPFD (=in - out)
  if (!is.null(meteo$PPFD_IN) && !is.null(meteo$PPFD_OUT)){
    meteo$ppfd <- meteo$PPFD_IN - meteo$PPFD_OUT    
  } else {
    meteo$ppfd <- rep( NA, dim(meteo)[1] )
  }

  ## convert units
  meteo$vpd  <- meteo$vpd  * 1e2  # given in hPa, required in Pa
  meteo$ppfd <- meteo$ppfd * 1.0e-6 * kfFEC * 60 * 60 * 24  # given in W m-2, required in mol m-2 d-1 

  # meteo <- select( meteo, year, moy, dom, year_dec, temp, prec, vpd, ppfd, nrad ) 
  meteo <- subset( meteo, select=c( year, moy, dom, year_dec, temp, prec, vpd, ppfd, nrad ) )

  return( meteo )

}
