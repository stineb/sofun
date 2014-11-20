## This is my implementation of the allocation optimisation by Makela et al., 2008. 
## All parameter values and equations are from that publication.

options(digits=16)
library(polynom) 

## set the independent variable: representing N availability (10 kgN/t-DW/yr = 10 gN/kg-DW/yr = 20 gN/kgC/yr = 0.02 gN/gC/yr)
#sigma_rM <- 0.02
sigma_rM <- 10

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


## light-saturated foliage-specific rate of photosynthesis (eq. 10)
sigma_fM <- function( r_p, sigma_fM0=sigma_fM0, r_ref=r_ref ) {
  sigma_fM <- sigma_fM0 * r_p / ( r_p + r_ref)
}

## actual (productive) photosynthetic N concentration (eq. 11)
r_p <- function( r_f, r_0=r_0 ) {
  r_p <- max( r_f - r_0, 0.0 )  
}

## Set foliar N:C ratio (1.2 kgN/kgDW)
#r_f <- 0.024 
r_f <- 0.012

val_sigma_fM <- sigma_fM( r_p( r_f, r_0 ), sigma_fM0, r_ref )

## define coefficients of cubic function (eq. S10)
beta1 <- Y * val_sigma_fM * K_f / ( 1/T_r + r_m * r_f * n_r )
beta2 <- (1/T_f + r_f * ( alpha_w * c_H / T_w + Y * r_m * (1 + n_w * alpha_w * c_H * r_f ))) / (1/T_r + Y * r_m * r_f * n_r )
beta3 <- sigma_rM * K_r / (r_f * (1-f) * n_r / T_r)
beta4 <- (((1-f)/T_f)+((1-f) * n_w * alpha_w * c_H * r_f / T_w )) / ((1-f) * n_r / T_r)

a1 <- ( beta1 - beta3 + K_r - K_f * (beta2+beta4)) / (- K_f)
a2 <- (beta1 * beta4 - beta2 * beta3 + K_r * (beta2+beta4) - K_f * beta2 * beta4) / (-K_f)
a3 <- K_r * beta2 * beta4 / (-K_f)

## use R polynom library to solve
poly <- polynomial( c( a3, a2, a1, 1 ) )
out <- solve(poly)

for (i in out){
  if (Im(i)==0 && Re(i)>0){
    ## i is solution in real space >0
    phi <- Re(i)
    print(paste("solution in positive real space: ",phi))

    ## Get leaf, wood, and root mass with phi and r_f using eq. S3
    W_f <- beta1 / (beta2+phi) - K_f
    W_r <- phi * W_f
    W_w <- alpha_w * c_H * r_f * W_f

    ## Alternatively, use eq. S7 and test if result is identical
    W_f_alt <- beta3 / (beta4+phi) - K_r / phi
    if (abs((W_f-W_f_alt)/W_f)>1e-12) { 
      print(paste('W_f     ', W_f))
      print(paste('W_f_alt ', W_f_alt))
      print(paste('diff    ', W_f_alt - W_f))
      print('S3 and S7 do not give identical results')
    }

    ## Get productivity from W_f, phi, and r_f
    G <- W_f/T_f + W_r/T_r + W_w/T_w  # eq. S1

    ## Alternatively, use eq. 4 and test if result is identical
    P <- val_sigma_fM * W_f * K_f / (W_f + K_f)
    R_m <- r_m*r_f*W_f + r_m*r_f*n_r*W_r + r_m*n_w*r_f*W_w
    G_alt <- Y * ( P - R_m )
    if (abs((G-G_alt)/G)>1e-12) { 
      print('not identical G calculated') 
      print(paste('G:     ', G))
      print(paste('G_alt: ', G_alt))
    }

  }
}

