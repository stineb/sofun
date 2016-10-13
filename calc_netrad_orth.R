##----------------------------------------------------------------------
## This calculation of net radiation as a function of irradiance is used 
## in the SWBM by Orth et al. (2013)
## No reference is given therein.
##
## Arguments:
## - irrad: irradiance, W m-2
##
## Function return variable:
## - netrad: net radiation, W m-2 
##----------------------------------------------------------------------
calc_netrad_orth <- function( irrad ){

  # Scaling solar radiation to approximate net radiation (xxx ref? xxx)
  netrad <- irrad * 0.65 - 35.0

  return( netrad )

}