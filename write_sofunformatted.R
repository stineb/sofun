write_sofunformatted <- function( filnam, data ){
  ## write output with standard formatting
  if ( is.vector( data ) ){
    len <- length(data)
    formatted <- vector( "character", len )
    for (i in 1:len){
      formatted[i] <- sprintf("%16.6f", data[i] )
    }
    writeLines( formatted, filnam )
  } else if ( is.data.frame( data ) && length( dim( data ) )==2 ){
    len <- dim( data )[1]
    formatted <- vector( "character", len )
    for (i in 1:len){
      formatted[i] <- sprintf("%16.6f    %f", data[i,1], data[i,2] )
    }
    writeLines( formatted, filnam )    
  }
}

