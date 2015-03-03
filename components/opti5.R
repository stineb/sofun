##----------------------------------------------------------------
## NUMERICAL APPROACH
## vary alpha and calculate implicit growth/N-uptake ratio
## this should be equal to the prescribed C:N ratio of new growth
## interpret this as finding the root of the function eval_cton
##----------------------------------------------------------------

## N availability over season
sin_season <- function(day){
  sea <- 0.5 * ( sin( (2 * pi * day / 365 ) + 1.5 * pi ) + 1 )
  return(sea)
}
ndayyear <- 365
days <- 1:ndayyear

## Input: inorganic N availability
## xxx as it stands now, it's extremely sensitive to inorganic N 
## availability. Linear N uptake function may not be appropriate.
Nscale <- 0.09
ninorg <- Nscale * sapply( days, FUN=sin_season )

## Assumption: root phenology is congruent with 
## seasonality of inorganic N availability
root_season <- sapply( days, FUN=sin_season )

## Assumption: leaf phenology is uniform over the year
leaf_season <- rep( 1, ndayyear )
# leaf_season <- sapply( days, FUN=sin_season )

## Model parameters, from Li et al., 2014 (where given)
params <- c( 
            ndayyear = 365,
            y        = 0.4,
            r_root   = 0.005,  # yields ~0.913 
            r_leaf   = 0.0,    # Rd is directly substracted from GPP
            exu      = 0.05,
            k_root   = 1/(1.04*365),
            k_leaf   = 1/(4.0*365),
            psi      = 0.1,
            kbeer    = 0.5,
            sla      = 0.0014,
            rntoc    = 1/50 
            )

P0 <- 30 # photosynthesis rate

calc_agpp <- function( alpha, P0, leaf_season, params ){
  ## annual GPP
  with( as.list(params), {
    gpp <-  sum( P0 * ( 1 - exp( -kbeer * sla * alpha * leaf_season ) ) )
    return( gpp )    
  })
}

calc_maxcleaf <- function( alpha, leaf_season ){
  maxcleaf <- max( leaf_season * alpha )
  return( maxcleaf )
}

calc_taa <- function( alpha, par, params ) {
  ## total above-ground allocation
  with( as.list(params), {
    taa <- sum(leaf_season) * alpha * ( (1+y) * k_leaf + r_leaf ) 
    return( taa )    
  })
}

beta_of_tba <- function( tba, root_season, params ){
  ## beta as a function of TBA
  with( as.list(params), {
    beta <- tba / ( sum(root_season) * ( (1+y) * k_root + r_root + exu ) )
    return( beta )
  })
}

calc_cexu <- function( nup, n0, params ){
  with( as.list(params),{
    cexu <- 1/psi * (log( n0 ) - log( n0 - nup ))
    return(cexu)
  })
}

calc_dnup <- function( cexu, n0, params ){
  with( as.list(params),{
    # print(paste("psi",psi))
    # print(paste("cexu",cexu))
    nup <- n0 - exp( log(n0) - psi * cexu )
    # print(paste("fraction N taken up:", nup/n0 ))`
    return(nup)
  })
}

anup <- function( beta, root_season, ninorg, params ){
  with( as.list(params),{

    # ## increasing cost (like FUN)
    # nup <- 0
    # for (doy in 1:ndayyear){
    #   cexu <- beta * root_season[doy] * exu
    #   nup <- nup + calc_dnup( cexu, ninorg[doy], params )
    #   # print(paste("fraction N taken up:", nup/ninorg[doy] ))
    # }

    # ## linear increase of N uptake with C exudation
    # nup <- 0
    # for (doy in 1:ndayyear){
    #   dnup <- beta * root_season[doy] * exu * psi * ninorg[doy]
    #   nup  <- nup + dnup
    #   # print(paste("fraction N taken up:", dnup/ninorg[doy] ))
    # }

    nup <- beta * psi * exu * sum( root_season * ninorg )
    # print(paste("fraction N taken up:", nup/ninorg[doy] ))

    return( nup )
  })
}

# nup_range <- seq( 0, 0.99, 0.01 )
# cex_range <- sapply( nup_range, FUN = function(x) calc_cexu(x, 1, params))
# nup_range_result <- sapply( cex_range, FUN = function(x) calc_dnup(x, 2, params) )

# plot( nup_range, cex_range, type="l")
# plot( cex_range, nup_range_result, type="l" )


aanreq <- function( alpha, leaf_season, params ) {
  with( as.list(params),{
    nreq <- rntoc * ( sum(leaf_season) * alpha * k_leaf )
    return( nreq )
  })
}

abnreq <- function( beta, root_season, params ) {
  with( as.list(params),{
    nreq <- rntoc * ( sum(root_season) * beta * k_root )
    return( nreq )
  })
}

calc_rd <- function( gpp ){
  rd <- 0.1 * gpp
  return(rd)
}

eval_beta <- function( alpha, P0, leaf_season, root_season, params ){
  ## given alpha, we can calculate how much C is fixed,
  ## total above-ground allocation, total below-ground 
  ## allocation, from which we can derive beta
  with( as.list(params), {

    ## GPP as afunction of leaf mass (alpha)
    gpp <- calc_agpp( alpha, P0, leaf_season, params )

    ## subtract dark respiration directly
    rd  <- calc_rd( gpp )
    gpp <- gpp - rd

    ## total above-ground allocation, TAA
    taa <- calc_taa( alpha, leaf_season, params )

    ## total below-ground allocation, TBA
    tba <- gpp - taa

    ## beta as a function of TBA
    beta <- beta_of_tba( tba, root_season, params )
    return(beta)
    
  })
}

calc_aresp <- function( alpha, leaf_season, params ){
  with( as.list(params), {
    aresp <- sum(leaf_season) * alpha * ( y * k_leaf + r_leaf ) 
    return(aresp)
  })
}

calc_bresp <- function( beta, root_season, params ){
  with( as.list(params), {
    bresp <- sum(root_season) * beta  * ( y * k_root + r_root )
    return(bresp)
  })
}

calc_exu <- function( beta, root_season, params ){
  with( as.list(params), {
    exu <- sum(root_season) * beta  * exu
    return(exu)
  })
}

calc_npp <- function( alpha, beta, P0, leaf_season, root_season, params){

  with( as.list(params), {
    
    ## GPP as afunction of leaf mass (alpha)
    gpp <- calc_agpp( alpha, P0, leaf_season, params )

    ## subtract dark respiration directly
    rd  <- calc_rd( gpp )

    ## substract all respiration terms (except exudation)
    aresp <- calc_aresp( alpha, leaf_season, params )
    bresp <- calc_bresp( beta, root_season, params )
    npp   <- gpp - aresp - bresp
    return( npp )

  })
}

eval_imbalance <- function( alpha, P0, leaf_season, root_season, ninorg, params ){
  ## ///////////////////////////////////////////////
  ## Evaluates the imbalance between N requirement 
  ## and N uptake, given leaf mass scalar
  ## -----------------------------------------------

  with( as.list(params), {

    ## 1. GPP as afunction of leaf mass (alpha)
    gpp <- calc_agpp( alpha, P0, leaf_season, params )

    ## 2. subtract dark respiration directly
    rd  <- calc_rd( gpp )

    ## 3. total above-ground allocation, TAA
    taa <- calc_taa( alpha, leaf_season, params )

    ## 4. total below-ground allocation, TBA
    tba <- gpp - taa

    ## 5. beta as a function of TBA
    beta <- beta_of_tba( tba, root_season, params )

    # ## derive beta from this alpha (function executes steps 1-4)
    # ## Steps 1-4 may also be calculated on a sheet of paper!
    # beta <- eval_beta( alpha, P0, leaf_season, root_season, params )

    ## 6. N uptake, given root mass (beta)
    nup <- anup( beta, root_season, ninorg, params )

    ## 7. N requirement for growth
    nreq <- aanreq( alpha, leaf_season, params ) + abnreq(  beta, root_season, params )

    ## 8. Imbalance between N uptake and N requirement
    eval <- nreq - nup
    return(eval)

  })
}

## fine-root specific respiration rate based on seasonal maximum
calc_r_root_impl <- function( beta, root_season, params ){
  with( as.list(params),{
    r_root_impl <- sum( beta * root_season * r_root ) / max( beta * root_season )
    return(r_root_impl)
  })
}

# out.alpha_root  <- uniroot( function(x) eval_imbalance( x, P0, leaf_season, root_season, ninorg, params), c(0.01,500) )
# alpha_root <- out.alpha_root$root
# beta_root  <- eval_beta( alpha_root, P0, leaf_season, root_season, params  )
# print( paste("leaf:root ratio", alpha_root * max(leaf_season) / beta_root * max(root_season) ) )

alpha_range <- seq( 0, 500, 10 )
cleaf_range <- sapply( alpha_range, FUN = function(x) calc_maxcleaf( x, leaf_season ) )
imbal_range <- sapply( alpha_range, FUN = function(x) eval_imbalance( x, P0, leaf_season, root_season, ninorg, params ) )
gpp_range   <- sapply( alpha_range, FUN = function(x) calc_agpp( x, P0, leaf_season, params ) )
taa_range   <- sapply( alpha_range, FUN = function(x) calc_taa(x, leaf_season, params ) )

beta_range  <- sapply( alpha_range, FUN = function(x) eval_beta( x, P0, leaf_season, root_season, params ) )
beta_range[ beta_range < 0 ]  <- NA
alpha_range[ beta_range < 0 ] <- NA

nup_range   <- sapply( beta_range,  FUN = function(x) anup( x, root_season, ninorg, params ) )
nreq_range  <- sapply( alpha_range, FUN = function(x) aanreq( x, leaf_season, params ) ) + sapply( beta_range, FUN = function(x) abnreq( x, root_season, params ) )
aresp_range <- sapply( alpha_range, FUN = function(x) calc_aresp( x, leaf_season, params ) )
bresp_range <- sapply( beta_range,  FUN = function(x) calc_bresp( x, root_season, params ) )
npp_range   <- gpp_range - aresp_range - bresp_range
exu_range   <- sapply( beta_range,  FUN = function(x) calc_exu( x, root_season, params ) )
# r_root_impl_range <- sapply( beta_range, FUN = function(x) calc_r_root_impl( x, root_season, params ) )

# plot( alpha_range, imbal_range, type="l" )
plot( alpha_range, cleaf_range, type="l" )
# plot( alpha_range, gpp_range, type="l" )
# plot( alpha_range, taa_range, type="l" )
# plot( alpha_range, beta_range, type="l" )
# plot( beta_range, nup_range, type="l" )
# plot( alpha_range, npp_range, type="l" )
# plot( alpha_range, exu_range/npp_range, type="l" )

# # pdf("bmodel_nreq_vs_nup.pdf")
# plot( alpha_range, nreq_range, type="l" )
# lines( alpha_range, nup_range, col="red" )
# abline( v=alpha_root, col='grey')
# legend( "topleft", c("N required", "N up"), lty=1, bty="n", col=c("black","red") )
# dev.off()

# plot( beta_range, nreq_range, type="l" )
# lines( beta_range, nup_range, col="red" )
# abline( v=beta_root, col='grey')
# legend( "topleft", c("N required", "N up"), lty=1, bty="n", col=c("black","red") )

