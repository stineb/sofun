add_year_dec <- function( df, year, moy, doy=NA ){

  if (is.na(doy)){
    doy <- rep( NA, dim(df)[1] )
    middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350)
    for (imo in 1:12){
      doy[ which(moy == imo) ] <- middaymonth[ imo ]
    }
  }

  year_dec <- year + (doy-1)/365
  df$year_dec <- year_dec

  return( df )

}