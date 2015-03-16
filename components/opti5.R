##----------------------------------------------------------------
## NUMERICAL APPROACH
## vary alpha and calculate implicit growth/N-uptake ratio
## this should be equal to the prescribed C:N ratio of new growth
## interpret this as finding the root of the function eval_cton
##----------------------------------------------------------------

source('pmodel.R')
source('/alphadata01/bstocker/utilities/daily2monthly.R')

##----------------------------------------------------------------
## PARAMETERS
##----------------------------------------------------------------

## N availability over season
sin_season <- function(day){
  sea <- 0.5 * ( sin( (2 * pi * day / 365 ) + 1.5 * pi ) + 1 )
  return(sea)
}
ndayyear <- 365
days <- 1:ndayyear
nmonth <- 12

## Input: inorganic N availability
## xxx as it stands now, it's extremely sensitive to inorganic N 
## availability. Linear N uptake function may not be appropriate.
Nscale <- 0.09
ninorg <- Nscale * sapply( days, FUN=sin_season )

## Assumption: root phenology is congruent with 
## seasonality of inorganic N availability
root_season <- sapply( days, FUN=sin_season )

## Assumption: 
## evergreen: leaf phenology is uniform over the year
##            => N in cell walls and N in Rubisco are uniform
leaf_season <- rep( 1, ndayyear )
# leaf_season <- sapply( days, FUN=sin_season )

## Model parameters, from Li et al., 2014 (where given)
params <- list( 
              ndayyear = 365,
              nmonth   = 12,
              y        = 0.25,   # changed from 0.4
              r_root   = 0.005,  # yields ~0.913 
              r_leaf   = 0.0,    # Rd is directly substracted from GPP
              exu      = 0.05,
              k_root   = 1/(1.04*365),
              k_leaf   = 1/365,  # xxx changed from 1/(4.0*365)
              psi      = 0.1,
              kbeer    = 0.5,
              sla      = 0.0083,  # 0.0083: Vile et al., 2005 (Ann. Botany; Table 2, Trees)     xxx changed from 0.0014
              rntoc    = 1/50 
              )

## test
calc_sla <- function(k_leaf){
  long_leaf <- 1/(365*k_leaf)
  sla <- 2.0 - 4.0 * exp( 6.150 - 0.460 * log(k_leaf) * 12.00 )
  return(sla)  
}


##----------------------------------------------------------------
## ENVIRONMENTAL CONDITIONS
##----------------------------------------------------------------

## Set fixed boundary conditions
elv    <- 450.0
co2    <- 376.0
cpalpha <- rep(1.26,nmonth)

## PPFD FROM SOFUN (STASH implementation)
## read PPFD calculated by STASH (sofun implementation, '*.m.qm.out')
filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.m.qm.out"
df.ppfd <- read.table( filnam, col.names=c("year","ppfd") )

## take subset of year 2000
istart <- which.min( abs(df.ppfd$year-2002.0) )
df.ppfd <- df.ppfd[ istart:(istart+nmonth-1), ]
mppfd  <- df.ppfd$ppfd

## VPD
filnam <- "/alphadata01/bstocker/sofun/trunk/components/mvpd_CH-Oe1_2002.txt"
df.vpd <- read.csv(filnam)
mvpd   <- df.vpd$vpd

## TEMPERATURE
df.temp <- read.csv( "/alphadata01/bstocker/sofun/trunk/components/mtemp_CH-Oe1_2002.txt" )
mtemp <- df.temp$temp


##----------------------------------------------------------------
## FUNCTION DEFINITIONS
##----------------------------------------------------------------

## Net light use efficiency: (gpp - rd) per unit light absorbed
## Returns a vector with monthly values. Update this vector every year.
calc_mluenet <- function( co2, mtemp, cpalpha, mvpd, elv ){

  mluenet <- rep( NA, nmonth )
  for (moy in 1:nmonth){
    mluenet[moy] <- pmodel( fpar=NA, ppfd=NA, co2, mtemp[moy], cpalpha[moy], mvpd[moy], elv )$luenet
  }
  return(mluenet)

}

## Beer's Law
calc_fpar <- function( lai, params ){
  with( params,{
    fpar <- ( 1 - exp( -kbeer * lai ) )
    return(fpar)
  })
}


## Returns annual "net" GPP (= GPP - Rd)
calc_agpp_net <- function( alpha, leaf_season, mluenet, mppfd, params ){

  with( params, {

    ## get monthly vector of fraction of absorbed PAR
    mleaf_season <- daily2monthly( leaf_season, method="mean")
    mfpar <- calc_fpar( sla * alpha * mleaf_season, params )

    ## calculate monthly gpp vector, convert from mol/m2/month to gC/m2/month
    mgpp <- mppfd * mfpar * mluenet * 12.0107

    ## annual GPP as sum of monthly GPP
    agpp <- sum(mgpp)

    return( agpp )    

  })
}


# calc_maxcleaf <- function( alpha, leaf_season ){
#   maxcleaf <- max( leaf_season * alpha )
#   return( maxcleaf )
# }

## total above-ground allocation
calc_taa <- function( alpha, par, params ) {
  with( as.list(params), {

    ## deciduous tree (no continuous turnover)
    taa <- ( 1.0 + y ) * alpha

    # ## continuous turnover/replacement
    # taa <- sum(leaf_season) * alpha * ( (1+y) * k_leaf ) 

    return( taa )    

  })
}


## beta as a function of TBA
beta_of_tba <- function( tba, root_season, params ){
  # Assumes that
  # TBA = sum( beta * root_sason(day) * ( (1+y) * k_root + r_root + exu ) )
  with( as.list(params), {
    beta <- tba / ( sum(root_season) * ( ( 1.0 + y ) * k_root + r_root + exu ) )
    return( beta )
  })
}


## daily C exudation as a function of N uptake for a given initial N
calc_cexu <- function( nup, n0, params ){
  with( as.list(params),{
    cexu <- 1/psi * (log( n0 ) - log( n0 - nup ))
    return(cexu)
  })
}


## daily N uptake as a function of C exuded and initial N content
calc_dnup <- function( cexu, n0, params ){
  # This follows from FUN approach with
  # dCexu/dNup = K / (N0 - Nup); K=1/psi
  with( as.list(params),{
    # print(paste("psi",psi))
    # print(paste("cexu",cexu))
    nup <- n0 * ( 1.0 - exp( - psi * cexu ) )
    # print(paste("fraction N taken up:", nup/n0 ))`
    return(nup)
  })
}


## annual N uptake as a sum over daily N uptake
anup <- function( beta, root_season, ninorg, params, method=NA ){
  with( as.list(params),{

    if (is.na(method)){
      ## increasing cost as limited by available N (like FUN)
      nup <- 0
      for (doy in 1:ndayyear){
        cexu <- beta * root_season[doy] * exu
        nup <- nup + calc_dnup( cexu, ninorg[doy], params )
        # print(paste("fraction N taken up:", nup/ninorg[doy] ))
      }  

    } else if (method=="approx") {
      ## linear increase of N uptake with C exudation
      nup <- beta * psi * exu * sum( root_season * ninorg )
    }

    return( nup )
  })
}

# nup_range <- seq( 0, 0.99, 0.01 )
# cex_range <- sapply( nup_range, FUN = function(x) calc_cexu(x, 1, params))
# nup_range_result <- sapply( cex_range, FUN = function(x) calc_dnup(x, 1, params) )

# # plot( nup_range, cex_range, type="l")
# plot( cex_range, nup_range_result, type="l" )


aanreq <- function( alpha, leaf_season, params ) {
  with( as.list(params),{

    ## deciduous tree (no continuous turnover)
    nreq <- 

    ## evergreen tree (continuous turnover)
    nreq <- rntoc * ( sum(leaf_season) * alpha * k_leaf )

    return( nreq )
  })
}

# abnreq <- function( beta, root_season, params ) {
#   with( as.list(params),{
#     nreq <- rntoc * ( sum(root_season) * beta * k_root )
#     return( nreq )
#   })
# }

# calc_rd <- function( gpp ){
#   rd <- 0.1 * gpp
#   return(rd)
# }


# calc_aresp <- function( alpha, leaf_season, params ){
#   with( as.list(params), {
#     aresp <- sum(leaf_season) * alpha * ( y * k_leaf ) 
#     return(aresp)
#   })
# }

# calc_bresp <- function( beta, root_season, params ){
#   with( as.list(params), {
#     bresp <- sum(root_season) * beta  * ( y * k_root + r_root )
#     return(bresp)
#   })
# }

# calc_exu <- function( beta, root_season, params ){
#   with( as.list(params), {
#     exu <- sum(root_season) * beta  * exu
#     return(exu)
#   })
# }

# calc_npp <- function( alpha, beta, P0, leaf_season, root_season, params){

#   with( as.list(params), {
    
#     ## GPP as afunction of leaf mass (alpha)
#     gpp <- calc_agpp( alpha, ppfd=ppfd, co2=co2, tc=tc, cpalpha=cpalpha, vpd=vpd, elv=elv, params=params, method="full" )

#     ## subtract dark respiration directly
#     rd  <- calc_rd( gpp )

#     ## substract all respiration terms (except exudation)
#     aresp <- calc_aresp( alpha, leaf_season, params )
#     bresp <- calc_bresp( beta, root_season, params )
#     npp   <- gpp - aresp - bresp
#     return( npp )

#   })
# }

# eval_beta <- function( alpha, leaf_season, root_season, mluenet, mppfd, params ){
#   ## given alpha, we can calculate how much C is fixed,
#   ## total above-ground allocation, total below-ground 
#   ## allocation, from which we can derive beta
#   ## all on an annual basis.
#   with( as.list(params), {

#     ## Annual net GPP (= GPP - Rd) as afunction of leaf mass (alpha)
#     gpp <- calc_agpp_net( alpha, leaf_season, mluenet, mppfd, params )

#     ## total above-ground allocation, TAA
#     taa <- calc_taa( alpha, leaf_season, params )

#     ## total below-ground allocation, TBA
#     tba <- gpp - taa

#     ## beta as a function of TBA
#     beta <- beta_of_tba( tba, root_season, params )
#     return(beta)
    
#   })
# }


# eval_imbalance <- function( alpha, leaf_season, root_season, ninorg, mluenet, params ){
#   ## ///////////////////////////////////////////////
#   ## Evaluates the imbalance between N requirement 
#   ## and N uptake, given leaf mass scalar
#   ## -----------------------------------------------

#   with( as.list(params), {

#     ## 1. GPP minus dark respiration as afunction of leaf mass (alpha)
#     gpp_net <- calc_agpp_net( alpha, leaf_season, mluenet, mppfd, params )

#     ## 3. total above-ground allocation, TAA
#     taa <- calc_taa( alpha, leaf_season, params )

#     ## 4. total below-ground allocation, TBA
#     tba <- gpp_net - taa

#     ## 5. beta as a function of TBA
#     beta <- beta_of_tba( tba, root_season, params )

#     # ## derive beta from this alpha (function executes steps 1-4)
#     # ## Steps 1-4 may also be calculated on a sheet of paper!
#     # beta <- eval_beta( alpha, P0, leaf_season, root_season, params )

#     ## 6. N uptake, given root mass (beta)
#     nup <- anup( beta, root_season, ninorg, params )

#     ## 7. N requirement for growth
#     nreq <- aanreq( alpha, leaf_season, params ) + abnreq(  beta, root_season, params )

#     ## 8. Imbalance between N uptake and N requirement
#     eval <- nreq - nup
#     return(eval)

#   })
# }

# ## fine-root specific respiration rate based on seasonal maximum
# calc_r_root_impl <- function( beta, root_season, params ){
#   with( as.list(params),{
#     r_root_impl <- sum( beta * root_season * r_root ) / max( beta * root_season )
#     return(r_root_impl)
#   })
# }

##----------------------------------------------------------------
## "MAIN" PROGRAM
##----------------------------------------------------------------

## get monthly vector of LUE
mluenet <- calc_mluenet( co2, mtemp, cpalpha, mvpd, elv )


# # out.alpha_root  <- uniroot( function(x) eval_imbalance( x, P0, leaf_season, root_season, ninorg, params), c(0.01,500) )
# # alpha_root <- out.alpha_root$root
# # beta_root  <- eval_beta( alpha_root, P0, leaf_season, root_season, params  )
# # print( paste("leaf:root ratio", alpha_root * max(leaf_season) / beta_root * max(root_season) ) )

alpha_range   <- seq( 1, 5000, 10 )
lai_range     <- alpha_range * params$sla
# cleaf_range   <- sapply( alpha_range, FUN = function(x) calc_maxcleaf( x, leaf_season ) )
# imbal_range   <- sapply( alpha_range, FUN = function(x) eval_imbalance( x, P0, leaf_season, root_season, ninorg, params ) )
gpp_range     <- sapply( alpha_range, FUN = function(x) calc_agpp_net( x, leaf_season, mluenet, mppfd, params ) )
taa_range     <- sapply( alpha_range, FUN = function(x) calc_taa( x, leaf_season, params ) )
tba_range     <- gpp_range - taa_range
beta_range    <- sapply( tba_range, FUN = function(x) beta_of_tba( x, root_season, params ) )
alpha_of_tba0 <- alpha_range[ which.min(abs(tba_range)) ]
lai_of_tba0   <- params$sla * alpha_of_tba0
      
beta_range_pos  <- beta_range[ beta_range >= 0 ]
alpha_range_pos <- alpha_range[ beta_range >= 0 ]

nup_range_pos         <- sapply( beta_range_pos,  FUN = function(x) anup( x, root_season, ninorg, params ) )
nup_range_pos_approx  <- sapply( beta_range_pos,  FUN = function(x) anup( x, root_season, ninorg, params, method="approx" ) )

# beta_range_ind <- seq(1, 1000, 10)
# nup_range_ind  <- sapply( beta_range_ind,  FUN = function(x) anup( x, root_season, ninorg, params ) )


# nreq_range  <- sapply( alpha_range, FUN = function(x) aanreq( x, leaf_season, params ) ) + sapply( beta_range, FUN = function(x) abnreq( x, root_season, params ) )
# aresp_range <- sapply( alpha_range, FUN = function(x) calc_aresp( x, leaf_season, params ) )
# bresp_range <- sapply( beta_range,  FUN = function(x) calc_bresp( x, root_season, params ) )
# npp_range   <- gpp_range - aresp_range - bresp_range
# exu_range   <- sapply( beta_range,  FUN = function(x) calc_exu( x, root_season, params ) )
# # r_root_impl_range <- sapply( beta_range, FUN = function(x) calc_r_root_impl( x, root_season, params ) )

##----------------------------------------------------------------
## PLOTS
##----------------------------------------------------------------
par( las=1, mar=c(4,4,4,4) )

# # plot( alpha_range, imbal_range, type="l" )
# plot( alpha_range, cleaf_range, type="l" )

## Above-ground balance
plot( alpha_range, gpp_range, type="l", ylab="annual flux (gC/m2/yr)", xlab="leaf mass (gC/m2)" )
lines( alpha_range, tba_range, col="blue" )
lines( alpha_range, taa_range, lty=2, col="red")
abline( h=0, col="grey70")
abline( v=alpha_of_tba0, col="grey70")
legend( "bottomright", c("Net GPP", "TBA", "TAA"), lty=c(1,1,2), bty="n", col=c("black","blue","red") )
par(new=TRUE)
plot( lai_range, gpp_range, type="n", xlab="", ylab="", axes=F )
axis( 3, col="green")
mtext( "LAI", 3, col="green", line=3 )
# par(new=TRUE)
# plot( alpha_range, beta_range, type="n", xlab="", ylab="", col="brown", axes=F, ylim=c(0,355))
# axis( 4, col="blue" )
# mtext( "root mass (gC/m2)", col="blue", side=4, line=3, las=1 )

## Below-ground balance (as a function of alpha)
plot( alpha_range_pos, beta_range_pos, type="l", col="brown", xlab="leaf mass (gC/m2)", ylim=c(0,max(beta_range)), xlim=range(alpha_range) )
axis( 2, col="brown" )
abline( v=alpha_of_tba0, col="grey70")
par(new=TRUE)
plot( alpha_range_pos, nup_range_pos, col="blue", type="l", xlab="", ylab="", axes=FALSE, ylim=c(0,max(nup_range_pos_approx)), xlim=range(alpha_range) )
lines( alpha_range_pos, nup_range_pos_approx, col="blue", lty=2 )
axis( 4, col="blue" )
legend( "bottomright", c("root mass", "annual N uptake"), lty=c(1,1,2), bty="n", col=c("brown","blue") )

## Below-ground balance (as a function of beta)
plot( beta_range_pos, nup_range_pos, type="l", col="blue" )
lines( beta_range_pos, nup_range_pos_approx, lty=2, col="blue" )

# plot( alpha_range, taa_range, type="l" )
# plot( alpha_range, tba_range, type="l" )
# plot( alpha_range, beta_range, type="l" )
# # plot( beta_range, nup_range, type="l" )
# # plot( alpha_range, npp_range, type="l" )
# # plot( alpha_range, exu_range/npp_range, type="l" )

# # # pdf("bmodel_nreq_vs_nup.pdf")
# # plot( alpha_range, nreq_range, type="l" )
# # lines( alpha_range, nup_range, col="red" )
# # abline( v=alpha_root, col='grey')
# # legend( "topleft", c("N required", "N up"), lty=1, bty="n", col=c("black","red") )
# # dev.off()

# # plot( beta_range, nreq_range, type="l" )
# # lines( beta_range, nup_range, col="red" )
# # abline( v=beta_root, col='grey')
# # legend( "topleft", c("N required", "N up"), lty=1, bty="n", col=c("black","red") )

