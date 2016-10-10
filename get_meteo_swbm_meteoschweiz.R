get_meteo_swbm_meteoschweiz <- function( path ){

  meteo <- read.table( path, skip=10, header=TRUE )

  # Reads the following columns of which the last three are represented by their codes
  # STA JAHR MO TG HH MM   211   224   237

  # Codes:
  # 211  Lufttemperatur 2 m über Boden; Tagesmittel  [°C]
  # 224  Globalstrahlung; Tagesmittel  [W/m²]
  # 237  Niederschlag; Kalendertag  [mm]

  colnames(meteo) <- c( "station_no", "year", "moy", "dom", "hr", "min", "temp", "rad", "prec" )

  return( meteo )

} 

