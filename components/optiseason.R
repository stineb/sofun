sin_season <- function(day){
  sea <- 1/2* (sin( (2* pi* day/365) + 1.5*pi ) + 1)
  return(sea)
}

gauss_season <- function(day) {
  k <- 60
  sea <- exp(-(day-ndayyear/2)^2/(2*k^2))
}

ndayyear <- 365
days <- seq(1,ndayyear,1)
# ninorg <- sapply( days, FUN=sin_season )
ninorg <- sapply( days, FUN=sin_season )

# shift_range <- seq(1,ndayyear,1)
shift_range <- 0
tba <- rep(NA,length(shift_range))
nup <- rep(NA,length(shift_range))
psi <- array(NA,dim=c(ndayyear,length(shift_range)))

plot( days, ninorg, type="l" )

for (idx in 1:length(shift_range)){

  croot <- sapply( days-shift_range[idx], FUN=sin_season )
  # croot <- rep(0.5,ndayyear)

  ## total below-ground allocation
  y <- 0.25
  r <- 0.1
  e <- 0.2
  k <- 1/ndayyear
  Q <- 1

  resp_grow <- rep(0,ndayyear)
  resp_main <- rep(0,ndayyear)
  export    <- rep(0,ndayyear)
  growth    <- rep(0,ndayyear)

  growth[1] <- croot[1]
  resp_grow[1] <- y * growth[1]
  resp_main[1] <- r * croot[1]
  export[1] <- e * croot[1]

  for (i in 2:ndayyear){

    ## growth
    growth[i] <- croot[i] - (croot[i-1] - k * croot[i-1])

    ## growth respiration
    resp_grow[i] <- y * growth[i]

    ## maintenance respiration
    resp_main[i] <- r * croot[i]

    ## export (exudation)
    export[i] <- e * croot[i]

  }

  ## total below-ground allocation
  tba[idx] <- sum( growth + resp_grow + resp_main + export )

  ## N uptake
  nup[idx] <- sum( Q * croot * ninorg )

  ## cost of N uptake
  psi[,idx] <- 1 / ( Q * ninorg )


  # Kc <- 2
  # Kn <- 1 
  # Cacq <- Kc * croot
  # Nacq <- Kn * Kc * Cacq * ninorg


# plot(  days, ninorg, type="l", ylim=c(0,1) )
# lines( days, croot, col="red" )
}

print( paste( "maximum TBA     ", max(tba) ) )
print( paste( "maximum N uptake", max(nup) ) )


# plot(  shift_range, nup, type="l" )
# lines( shift_range, tba, col="red")

# plot( 1:ndayyear, psi[,1], type="l", ylim=c(0,1000) )
# abline( h=10, col="red" )

