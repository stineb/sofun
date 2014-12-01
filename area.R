############################################################################
## Returns grid cell area in m2 on a spherical Earth
## MUST BE CONSISTENT WITH CALCULATION AND EARTH RADIUS AS DEFINED IN
## outannual.F, modelpara.inc !!!
## Beni Stocker, 25.03.2011.
## test
## -------------------------------------------------------------------------

## dx and dy are in units of degrees

area <- function(lat,dx,dy){
  r_earth <- 6370000

  area <- 4*r_earth^2*0.5*dx*pi/180*
    cos(abs(lat)*pi/180)*
    sin(0.5*dy*pi/180)
  return(area)
}

