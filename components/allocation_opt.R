## This is my implementation of the allocation optimisation by Makela et al., 2008. 
## All parameter values and equations are from that publication.

options(digits=8)
library(polynom) 

check.sanity <- TRUE
do.pine      <- TRUE
do.spruce    <- FALSE
verbose      <- FALSE

## set the independent variable: representing N availability (10 kgN/t-DW/yr = 10 gN/kg-DW/yr = 20 gN/kgC/yr = 0.02 gN/gC/yr)
#sigma_rM <- 0.02

if (do.pine||do.spruce) {

  if (do.pine){
    ## Parameters for Pinus sylvestris
    K_r       <- 2000
    K_f       <- 2500
    T_f       <- 3.3
    T_r       <- 1.25
    T_w       <- 40
    Y         <- 1.54
    r_m       <- 16
    sigma_fM0 <- 8.0
    n_r       <- 1
    n_w       <- 0.07
    f         <- 0.3
    alpha_w   <- 0.8
    c_H       <- 2800
    r_0       <- 0.009
    r_ref     <- 0.002
  } else if (do.spruce){
    ## Parameters for Picea abies
    K_r       <- 2000
    K_f       <- 8000
    T_f       <- 8
    T_r       <- 1.25
    T_w       <- 33.3
    Y         <- 1.54
    r_m       <- 16
    sigma_fM0 <- 4.0
    n_r       <- 1
    n_w       <- 0.07
    f         <- 0.3
    alpha_w   <- 0.4
    c_H       <- 3400
    r_0       <- 0.008
    r_ref     <- 0.002
  }

  ## light-saturated foliage-specific rate of photosynthesis (eq. 10); [sigma_fM] = kgC/kgDW/yr
  sigma_fM <- function( r_p, sigma_fM0=sigma_fM0, r_ref=r_ref ) {
    sigma_fM <- sigma_fM0 * r_p / ( r_p + r_ref)
  }

  ## actual (productive) photosynthetic N concentration (eq. 11)
  r_p <- function( r_f, r_0=r_0 ) {
    r_p <- max( r_f - r_0, 0.0 )  
  }

  ## Set range of foliar N:C ratio (1.2 kgN/kgDW)
  #r_f <- 0.024 
  r_f0 <- 0.012
  n_r_f0 <- 100
  r_f_range <- seq( from=r_f0-r_f0*0.5, to=r_f0+r_f0*0.5 , by=((r_f0+r_f0*0.5)-(r_f0-r_f0*0.5))/n_r_f0 )
  len_r_f <- length(r_f_range)

  ## Set range of N availability (=sigma_rM)
  n_sigma_rM <- 100
  sigma_rM0 <- 10/1000  # convert from kgN/t-DW to kgN/kg-DW
  sigma_rM1 <- 50/1000
  sigma_rM_range <- seq( from=sigma_rM0, to=sigma_rM1, by=(sigma_rM1-sigma_rM0)/n_sigma_rM )
  len_sigma_rM <- length(sigma_rM_range)

  ## initialise output arrays
  W_f_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  W_r_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  W_w_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  G_out       <- array(NA,dim=c(len_r_f,len_sigma_rM))
  P_out       <- array(NA,dim=c(len_r_f,len_sigma_rM))
  R_m_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  U_out       <- array(NA,dim=c(len_r_f,len_sigma_rM))
  phi_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  r_f_out     <- array(NA,dim=c(len_r_f,len_sigma_rM))
  val_sigma_fM_out <- array(NA,dim=c(len_r_f,len_sigma_rM))
  if (check.sanity) {
    G_alt_out   <- array(NA,dim=c(len_r_f,len_sigma_rM))
    U_alt_out   <- array(NA,dim=c(len_r_f,len_sigma_rM))
    W_f_alt_out <- array(NA,dim=c(len_r_f,len_sigma_rM))  
  }

  for ( jdx in 1:len_sigma_rM ){
    ## ----------------------------------------------
    ## Vary N availability
    ## ----------------------------------------------
    sigma_rM <- sigma_rM_range[jdx]

    idx <- 0
    for (r_f in r_f_range) {
      idx <- idx + 1

      ## ----------------------------------------------
      ## Vary foliage C:N ratio
      ## ----------------------------------------------
      r_f_out[idx,jdx] <- r_f

      ## get value of sigma_fM (function)
      val_sigma_fM <- sigma_fM( r_p( r_f, r_0 ), sigma_fM0, r_ref )

      ## define coefficients of cubic function (eq. S10)
      beta1 <- Y * val_sigma_fM * K_f / ( 1/T_r + r_m * r_f * n_r )
      beta2 <- (1/T_f + r_f * ( alpha_w * c_H / T_w + Y * r_m * (1 + n_w * alpha_w * c_H * r_f ))) / (1/T_r + Y * r_m * r_f * n_r )
      beta3 <- sigma_rM * K_r / (r_f * (1-f) * n_r / T_r)
      beta4 <- (((1-f)/T_f)+((1-f) * n_w * alpha_w * c_H * r_f / T_w )) / ((1-f) * n_r / T_r)

      a1 <- ( beta1 - beta3 + K_r - K_f * (beta2+beta4)) / (- K_f)
      a2 <- (beta1 * beta4 - beta2 * beta3 + K_r * (beta2+beta4) - K_f * beta2 * beta4) / (-K_f)
      a3 <- K_r * beta2 * beta4 / (-K_f)

      ## use R polynom library to solve, this yields three solutions for phi
      poly <- polynomial( c( a3, a2, a1, 1 ) )
      out <- solve(poly)

      sols <- 0
      G <- 0
      for (i in out){
        ## ----------------------------------------------
        ## loop over the three solutions for phi
        ## ----------------------------------------------
        if (Im(i)==0 && Re(i)>0){

          ## i is solution in positive real space
          phi <- Re(i)
          sols <- sols + 1

          ## Get leaf, wood, and root mass with phi and r_f using eq. S3
          W_f <- beta1 / (beta2+phi) - K_f
          W_r <- phi * W_f
          W_w <- alpha_w * c_H * r_f * W_f

          ## Get productivity from W_f, phi, and r_f (= NPP)
          G_check <- W_f/T_f + W_r/T_r + W_w/T_w  # eq. S1

          if (G_check>G) {
            ## ----------------------------------------------
            ## This solution of phi is superior to the previous one w.r.t. growth
            ## ----------------------------------------------
            G <- G_check

            ## N uptake use eq. 6
            U <- r_f * W_f * ((1-f)/T_f + (1-f)*n_r*phi/T_r + (1-f)*n_w*alpha_w*c_H*r_f/T_w)

            ## GPP, eq. 6
            if ((W_f + K_f)!=0) {
              P <- val_sigma_fM * W_f * K_f / (W_f + K_f)
            } else {
              print('invalid solution: W_f + K_f = 0')
              stop
            }

            ## Maintenance respiration, eq. 5
            R_m <- r_m*r_f*W_f + r_m*r_f*n_r*W_r + r_m*n_w*r_f*W_w  # eq. 5

            ## write to output
            W_f_out[idx,jdx] <- W_f
            W_r_out[idx,jdx] <- W_r
            W_w_out[idx,jdx] <- W_w
            G_out[idx,jdx] <- G
            P_out[idx,jdx] <- P
            R_m_out[idx,jdx] <- R_m
            U_out[idx,jdx] <- U
            phi_out[idx,jdx] <- phi
            r_f_out[idx,jdx] <- r_f
            val_sigma_fM_out[idx,jdx] <- val_sigma_fM

            if (check.sanity) { 

              ## Alternatively, use eq. S7 and test if result is identical
              W_f_alt <- beta3 / (beta4+phi) - K_r / phi
              if (abs((W_f-W_f_alt)/W_f)>1e-9) { 
                print(paste('W_f     ', W_f))
                print(paste('W_f_alt ', W_f_alt))
                print(paste('diff    ', W_f_alt - W_f))
                print('S3 and S7 do not give identical results')
              }        

              ## Alternatively, use eq. 4 and test if result is identical
              G_alt <- Y * ( P - R_m )
              G_alt_out[idx,jdx] <- G_alt
              if (verbose) {
                if (abs((G-G_alt)/G)>1e-12) { 
                  print('not identical G calculated') 
                  print(paste('G:     ', G))
                  print(paste('G_alt: ', G_alt))
                } 
              }

              ## For N uptake, alternativaly use eq. 9
              U_alt <- sigma_rM*W_r*K_r / (W_r+K_r)
              U_alt_out[idx,jdx] <- U_alt
              if (verbose){
                if (abs((U-U_alt)/U)>1e-12) { 
                  print('not identical U calculated') 
                  print(paste('U:     ', U))
                  print(paste('U_alt: ', U_alt))
                } 
              }


            }


          }

        }
      }
      # print(paste('number of solutions: ', sols))
    }
  }

  ## for each N availability, find foliar-C:N index for which G is max
  idx.max <- rep(NA,len_sigma_rM)
  for (jdx in 1:len_sigma_rM){
    idx.max[jdx] <- which.max(G_out[,jdx])
  }

  ## Reduce arrays to maximum in r_f (1st dimension of *_out arrays)
  W_f_max   <- rep(NA,len_sigma_rM)
  W_r_max   <- rep(NA,len_sigma_rM)
  W_w_max   <- rep(NA,len_sigma_rM)
  G_max     <- rep(NA,len_sigma_rM) # = NPP
  P_max     <- rep(NA,len_sigma_rM) # = GPP
  R_m_max   <- rep(NA,len_sigma_rM)
  U_max     <- rep(NA,len_sigma_rM)
  phi_max   <- rep(NA,len_sigma_rM)
  r_f_max   <- rep(NA,len_sigma_rM)
  val_sigma_fM_max <- rep(NA,len_sigma_rM)
  if (check.sanity) {
    G_alt_max <- rep(NA,len_sigma_rM)
    U_alt_max <- rep(NA,len_sigma_rM)
  }
  for (jdx in 1:len_sigma_rM){
    W_f_max[jdx]   <- W_f_out[idx.max[jdx],jdx]
    W_r_max[jdx]   <- W_r_out[idx.max[jdx],jdx]
    W_w_max[jdx]   <- W_w_out[idx.max[jdx],jdx]
    G_max  [jdx]   <- G_out[idx.max[jdx],jdx]
    P_max  [jdx]   <- P_out[idx.max[jdx],jdx]
    R_m_max[jdx] <- R_m_out[idx.max[jdx],jdx]
    U_max  [jdx]   <- U_out[idx.max[jdx],jdx]
    phi_max[jdx]   <- phi_out[idx.max[jdx],jdx]
    r_f_max[jdx]   <- r_f_out[idx.max[jdx],jdx]
    val_sigma_fM_max[jdx] <- val_sigma_fM_out[idx.max[jdx],jdx]
    if (check.sanity) {
      G_alt_max[jdx] <- G_alt_out[idx.max[jdx],jdx]
      U_alt_max[jdx] <- U_alt_out[idx.max[jdx],jdx]
    }
  }
}

## Visualisations
par(las=1)

## Growth vs. leaf C:N, for a given N availability
# plot( r_f_range, G_out[,50], type="l", xlab="foliar C:N [kg-DW/kgN]", ylab="grwoth [kg-DW/yr]" )
# if (check.sanity) {lines( r_f_range, G_alt_out[,50], col="red" )}
# points( r_f_range[which.max(G_out[,50])], G_out[which.max(G_out[,50]),50] )
# points( r_f_range[which.max(G_alt_out[,50])], G_alt_out[which.max(G_alt_out[,50]),50], col="red" )

## N uptake vs. leaf C:N, for a given N availability
# plot( r_f_range, U_out[,50], type="l", xlab="foliar C:N [kg-DW/kgN]", ylab="N uptake [kg-N/yr]" )
# lines( r_f_range, U_alt_out[,50], col="red" )

## Fig. 1a: biomass vs. N availability
# plot( sigma_rM_range*1000, W_f_max/1000, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="biomass [t-DW/ha]", ylim=c(0,18) )
# lines( sigma_rM_range*1000, W_r_max/1000, col="red" )
# title("Fig.1a")

## Fig. 1b: foliage:fine roots vs. N availability
# plot( sigma_rM_range*1000, 1/phi_max, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="Foliage:fine roots", ylim=c(0,16) )
# title("Fig.1b")

## Fig. 2a: foliar C:N vs. N availability
# plot( sigma_rM_range*1000, r_f_max*100, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="Foliar [N] (%)", ylim=c(0.8,1.6) )
# title("Fig.2a")

## Fig. 2b: N uptake vs. N availability
# plot( sigma_rM_range*1000, U_max, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="Foliar [N] (%)", ylim=c(0,55) )
# title("Fig.2b")

## Fig. 3a: GPP/NPP vs. N availability
# plot( sigma_rM_range*1000, P_max/1000, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="Production (tC/ha/yr)", ylim=c(0,16) )
# lines(sigma_rM_range*1000, G_max/1000, col="red" )
# lines(sigma_rM_range*1000, G_alt_max/1000, col="red", lty=2 )
# title("Fig.3a")

## Fig. 3b: NPP:GPP ratio
# plot( sigma_rM_range*1000, G_alt_max/P_max, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="NPP:GPP" )
# title("Fig.3b")
# # ylim=c(0.2,0.6)

## Fig. 3c: Photosynthesis (leaf-specific)
plot( sigma_rM_range*1000, val_sigma_fM_max, type="l", xlab="N availability [kg N / t fine root / yr]", ylab="Photosynthesis (tC/tDW/yr)" )
title("Fig.3c")




