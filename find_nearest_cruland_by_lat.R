find_nearest_cruland_by_lat <- function( lon, lat ){

  library(ncdf)
  source("/home/bstocker/.Rprofile")

  nc <- open.ncdf( paste( myhome, "data/cru/cru_ts_3.23/cru_ts3.23.190101.tmp.dat.nc", sep=""), readunlim=FALSE )
  crufield <- get.var.ncdf( nc, varid="TMP" )
  lon_vec <- get.var.ncdf( nc, varid="LON" )
  lat_vec <- get.var.ncdf( nc, varid="LAT" )
  crufield[crufield==-9999] <- NA
  close.ncdf(nc)

  ilon <- which.min( abs( lon_vec - lon ) )
  ilat <- which.min( abs( lat_vec - lat ) )

  if (!is.na(crufield[ilon,ilat])) {print("WRONG: THIS SHOULD BE NA!!!")}
  for (n in seq(2*length(lon_vec))){
    ilon_look <- (-1)^(n+1)*round((n+0.1)/2)+ilon
    if (ilon_look > length(lon_vec)) {ilon_look <- ilon_look %% length(lon_vec)} ## Wrap search around globe in latitudinal direction
    if (ilon_look < 1)               {ilon_look <- ilon_look + length(lon_vec) }
    print(paste("ilon_look",ilon_look))
    if (!is.na(crufield[ilon_look,ilat])) {
      break
    }
  }
  # if (!is.na(crufield[ilon_look,ilat])) {print("SUCCESSFULLY FOUND DATA")}
  return( lon_vec[ ilon_look ] )
  
}
