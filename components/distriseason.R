sin_season <- function(day){
  sea <- 1/2* (sin( (2* pi* day/365) + 1.5*pi ) + 1)
  return(sea)
}

gauss_season <- function(day) {
  k <- 60
  sea <- exp(-(day-ndayyear/2)^2/(2*k^2))
}

y <- 0.25
r <- 0.1
e <- 0.2
k <- 1/ndayyear
Q <- 1

ndayyear <- 365
days <- seq(1,ndayyear,1)
ninorg <- sapply( days, FUN=sin_season )
guess  <- ninorg

## assume we put 100 units below ground - what amount of N do we get in return?
tba_trgt <- 100
scale_range <- seq( 1, 3, 0.001 )
# scale_range <- 1.81
tba <- rep(NA,length(scale_range))
nup <- rep(NA,length(scale_range))

for (idx in 1:length(scale_range)) {

  ## scale C root, fixing seasonality
  croot <- scale_range[idx] * guess

  ## get TBA for this scaling
  resp_grow <- rep(0,ndayyear)
  resp_main <- rep(0,ndayyear)
  export    <- rep(0,ndayyear)
  growth    <- rep(0,ndayyear)

  growth[1] <- croot[1]
  resp_grow[1] <- y * growth[1]
  resp_main[1] <- r * croot[1]
  export[1] <- e * croot[1]

  for (jdx in 2:ndayyear){

    ## growth
    growth[jdx] <- croot[jdx] - (croot[jdx-1] - k * croot[jdx-1])

    ## growth respjdxratjdxon
    resp_grow[jdx] <- y * growth[jdx]

    ## majdxntenance respjdxratjdxon
    resp_main[jdx] <- r * croot[jdx]

    ## export (exudatjdxon)
    export[jdx] <- e * croot[jdx]

  }

  ## total below-ground allocation
  tba[idx] <- sum( growth + resp_grow + resp_main + export )

  ## N uptake
  nup[idx] <- sum( Q * e * croot * ninorg)

}

scale_trgt <- scale_range[ which.min( abs(tba-tba_trgt) ) ]
scale_pred <- tba_trgt / ( guess[1] * ( (1+y)*k + e + r ) * sum(guess/guess[1]) )

nup_pred <- sum( Q * e * guess * scale_pred * ninorg ) 

print(paste("scale_trgt",scale_trgt))
print(paste("scale_pred",scale_pred))
print(paste("TBA       ",tba[ which.min( abs(tba-tba_trgt) ) ]))
print(paste("NUP       ",nup[ which.min( abs(tba-tba_trgt) ) ]))
print(paste("NUP pred  ",nup_pred ))

plot( scale_range, tba, type="l" )
lines( scale_range, nup, col="red" )


