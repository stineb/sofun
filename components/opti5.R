##----------------------------------------------------------------
## NUMERICAL APPROACH
## vary alpha and calculate implicit growth/N-uptake ratio
## this should be equal to the prescribed C:N ratio of new growth
## interpret this as finding the root of the function eval_cton
##----------------------------------------------------------------

source('/alphadata01/bstocker/sofun/trunk/components/pmodel.R')
source('/alphadata01/bstocker/sofun/trunk/components/sofun.R')
source('/alphadata01/bstocker/utilities/daily2monthly.R')


## test
calc_sla <- function(k_leaf){
  long_leaf <- 1/(365*k_leaf)
  sla <- 2.0 - 4.0 * exp( 6.150 - 0.460 * log(k_leaf) * 12.00 )
  return(sla)  
}


##----------------------------------------------------------------
## ENVIRONMENTAL CONDITIONS
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
Nscale <- 0.1
ninorg <- Nscale * sapply( days, FUN=sin_season )
ntrans <- 10.0  # free N taken up by transpiration stream

## Assumption: root phenology is congruent with 
## seasonality of inorganic N availability
root_season <- sapply( days, FUN=sin_season )

## Assumption: 
## evergreen: leaf phenology is uniform over the year
##            => N in cell walls and N in Rubisco are uniform
leaf_season <- rep( 1, ndayyear )
# leaf_season <- sapply( days, FUN=sin_season )
if (max(leaf_season)>1.0) { break }

## Set fixed boundary conditions
elv    <- 450.0
co2    <- 376.0
cpalpha <- rep(1.26,nmonth)

## PPFD FROM SOFUN (STASH implementation)
## read PPFD calculated by STASH (sofun implementation, '*.m.qm.out')
filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.m.qm.out"
df.ppfd <- read.table( filnam, col.names=c("year","ppfd") )

## take subset of year 2002
istart <- which.min( abs(df.ppfd$year-2002.0) )
df.ppfd <- df.ppfd[ istart:(istart+nmonth-1), ]
mppfd  <- df.ppfd$ppfd

## DAY LENGTH 
filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.d.ds.out"
df.dl <- read.table( filnam, col.names=c("year","dayl_h" ) )

## take subset of year 2002
istart <- which.min( abs(df.dl$year-2002.0) )
df.dl <- df.dl[ istart:(istart+ndayyear-1), ]
ddl  <- df.dl$dayl_h

## VPD
filnam <- "/alphadata01/bstocker/sofun/trunk/components/mvpd_CH-Oe1_2002.txt"
df.vpd <- read.csv(filnam)
mvpd   <- df.vpd$vpd

## TEMPERATURE
df.temp <- read.csv( "/alphadata01/bstocker/sofun/trunk/components/mtemp_CH-Oe1_2002.txt" )
mtemp <- df.temp$temp


# ## fine-root specific respiration rate based on seasonal maximum
# calc_r_root_impl <- function( beta, root_season, params ){
#   with( params,{
#     r_root_impl <- sum( beta * root_season * r_root ) / max( beta * root_season )
#     return(r_root_impl)
#   })
# }


##----------------------------------------------------------------
## "MAIN" PROGRAM
##----------------------------------------------------------------

## Run P-model for each month 

mluenet <- rep( NA, nmonth )
mlue    <- rep( NA, nmonth )
mnapar  <- rep( NA, nmonth )
for (moy in 1:nmonth){

  ## Execute P-model
  out <- pmodel( fpar=NA, ppfd=NA, co2, mtemp[moy], cpalpha[moy], mvpd[moy], elv )

  ## Light use efficiency: (gpp - rd) per unit light absorbed
  mlue[moy] <- out$lue

  ## Net light use efficiency: (gpp - rd) per unit light absorbed
  mluenet[moy] <- out$luenet

  ## conversion factor to get from APAR to Rubisco-N
  mnapar[moy]  <- out$n_apar

}


## Calculate Rubisco-N
daysecs   <- df.dl$dayl_h * 60 * 60  # number of daylight seconds in a year
monsecs   <- daily2monthly( daysecs, method="sum" )
meanmppfd <- mppfd / monsecs         # mol m-2 s-1



## VARY LAI (alpha)
alpha_range   <- seq( 0.01, 30, 0.1 )

gpp_net_range  <- sapply( alpha_range, FUN = function(x) calc_agpp( x, leaf_season, mluenet, mppfd, params ) )
gpp_range      <- sapply( alpha_range, FUN = function(x) calc_agpp( x, leaf_season, mlue, mppfd, params ) )
nr_canop_range <- sapply( alpha_range, FUN = function(x) calc_nr_canop( x, leaf_season, meanmppfd, mnapar, params ) )
nr_leaf_range  <- sapply( alpha_range, FUN = function(x) calc_nr_leaf( x, leaf_season, meanmppfd, mnapar, params ) )

ncw_canop_range <- sapply( alpha_range, FUN = function(x) calc_ncw_canop( x, leaf_season, meanmppfd, mnapar, params ) )
ncw_canop_range_nomin <- sapply( alpha_range, FUN = function(x) calc_ncw_canop_nomin( x, leaf_season, meanmppfd, mnapar, params ) )
ncw_leaf_range  <- sapply( alpha_range, FUN = function(x) calc_ncw_leaf( x, leaf_season, meanmppfd, mnapar, params ) )

c_leaf_range   <- ncw_leaf_range  * mol_weight_n * params$r_cton_leaf # g C m-2
c_canop_range  <- ncw_canop_range * mol_weight_n * params$r_cton_leaf # g C m-2
# c_leaf_range   <- (ncw_leaf_range  + nr_leaf_range  ) * mol_weight_n * params$r_cton_leaf # g C m-2
# c_canop_range  <- (ncw_canop_range + nr_canop_range ) * mol_weight_n * params$r_cton_leaf # g C m-2

lma_range      <- 2.0 * c_leaf_range

taa_range      <- sapply( c_canop_range, FUN = function(x) calc_taa( x, leaf_season, params ) )
tba_range      <- gpp_net_range - taa_range
beta_range     <- sapply( tba_range, FUN = function(x) beta_of_tba( x, root_season, params ) )

aresp_range <- sapply( c_canop_range, FUN = function(x) calc_aresp( x, leaf_season, params ) )
bresp_range <- sapply( beta_range,  FUN = function(x) calc_bresp( x, root_season, params ) )
npp_range   <- gpp_net_range - aresp_range - bresp_range
exu_range   <- sapply( beta_range,  FUN = function(x) calc_exu( x, root_season, params ) )

root2shoot_range <- beta_range / c_leaf_range

nup_range         <- sapply( beta_range,  FUN = function(x) anup( x, root_season, ninorg, params ) ) + ntrans
nup_range_approx  <- sapply( beta_range,  FUN = function(x) anup( x, root_season, ninorg, params, method="approx" ) ) + ntrans

## balance in g N m-2
aanreq_range  <- (nr_canop_range + ncw_canop_range) * mol_weight_n
aanreq_range_nomin <- (nr_canop_range + ncw_canop_range_nomin) * mol_weight_n
abnreq_range  <- sapply( beta_range, FUN = function(x) abnreq( x, root_season, params ) )
anreq_range   <- aanreq_range + abnreq_range
anreq_range_nomin  <- aanreq_range_nomin + abnreq_range

imbal_range   <- sapply( alpha_range, FUN = function(x) eval_imbalance( x, leaf_season, root_season, ninorg, ntrans, mluenet, mppfd, meanmppfd, mnapar, params ) )
out.alpha_root  <- uniroot( function(x) eval_imbalance( x, leaf_season, root_season, ninorg, ntrans, mluenet, mppfd, meanmppfd, mnapar, params ), range(alpha_range) )
alpha_root      <- out.alpha_root$root


##----------------------------------------------------------------
## PLOTS
##----------------------------------------------------------------
data <- read.csv(file="/alphadata01/bstocker/data/hikosaka/TableS1_reformatted.csv",header = TRUE)

par( las=1, mar=c(4,4,4,4) )

## Primary productivity
plot(  alpha_range, gpp_range, type="l", ylab="annual flux (gC/m2/yr)", xlab="LAI" )
lines( alpha_range, gpp_net_range, col="magenta" )
lines( alpha_range, npp_range, col="green")
lines( alpha_range, exu_range, col="blue" )
abline( h=0, col="grey70")
abline( v=alpha_root, col="magenta ")
legend( "bottomright", c("GPP","GPP - Rd", "NPP", "EXU"), lty=1, bty="n", col=c("black","magenta","green","blue") )

## Exudation as a fraction of NPP
plot( alpha_range, exu_range/npp_range, type="l")
abline( v=alpha_root, col="magenta ")

## Above-ground balance
plot(  alpha_range, gpp_range, type="l", lty=2, ylab="annual flux (gC/m2/yr)", xlab="LAI" )
lines( alpha_range, gpp_net_range )
lines( alpha_range, taa_range, col="red")
lines( alpha_range, tba_range, col="blue" )
abline( h=0, col="grey70")
abline( v=alpha_root, col="magenta ")
legend( "bottomright", c("GPP","GPP - Rd", "TBA", "TAA"), lty=c(2,1,1,1), bty="n", col=c("black","black","blue","red") )

## Below-ground balance (as a function of alpha)
plot( alpha_range, beta_range, type="l", col="brown", xlab="leaf mass (gC/m2)", ylim=c(0,max(beta_range)), xlim=range(alpha_range) )
par(new=TRUE)
plot( alpha_range, nup_range, col="blue", type="l", xlab="", ylab="", axes=FALSE, ylim=c(0,max(nup_range_approx)), xlim=range(alpha_range) )
lines( alpha_range, nup_range_approx, col="blue", lty=2 )
axis( 4, col="blue" )
abline( v=alpha_root, col="magenta ")
legend( "bottomright", c("root mass", "annual N uptake","annual N uptake, approx."), lty=c(1,1,2), bty="n", col=c("brown","blue","blue") )

plot( alpha_range, root2shoot_range, type="l")
abline( v=alpha_root, col="magenta ")

# ## Below-ground balance (as a function of beta)
# plot(  beta_range, nup_range, type="l", col="blue" )
# lines( beta_range, nup_range_approx, lty=3, col="blue" )
# legend( "bottomright", c("annual N uptake","annual N uptake, approx."), lty=c(1,2), bty="n", col=c("blue","blue") )

## Below-ground balance (as a function of alpha)
# plot(  alpha_range, nup_range_approx, type="l", lty=2, col="blue", ylab="annual N uptake (gN/m2/yr)", xlab="LAI", ylim=c(0,max(anreq_range))  )
plot(  alpha_range, nup_range, type="l", col="blue", ylab="annual N uptake (gN/m2/yr)", xlab="LAI", ylim=c(0,max(anreq_range))  )
lines( alpha_range, anreq_range, col="red" )
# lines( alpha_range, imbal_range, col="magenta" )
abline( v=alpha_root, col="magenta ")
abline( h=0.0, col="grey70" )
# lines( alpha_range, anreq_range_nomin, col="red", lty=2 )
legend( "bottomright", c("annual N uptake","annual N uptake, approx.","N required, n_cw ~ Vcmax25 + LAI x const.","N required, n_cw ~ Vcmax25"), lty=c(1,2,1,2), bty="n", col=c("blue","blue","red","red") )

# ## relationship between N per leaf area and Rubisco per leaf area (y-intersect is a measure for minimum content)
# Rarea <- data$RNF * data$Narea
# reg <- lm( data$Narea ~ Rarea )
# plot(  data$Narea ~ Rarea, ylim=c(0,max(data$Narea, na.rm=TRUE)), xlim=c(0,max(Rarea, na.rm=TRUE)), xlab="rNarea (mol Rubisco N m-2)", ylab="Narea (mol N m-2)" )
# lines( nr_leaf_range, nr_leaf_range + ncw_leaf_range, col="red" )
# text( 0, 0.03, paste("y-intersect:", reg$coefficients[1] ), adj=c(0,0))
# text( 0, 0.01, paste("slope:      ", reg$coefficients[2] ), adj=c(0,0))
# abline( reg )


## relationship between LMA and Rubisco per leaf area (y-intersect is a measure for minimum content)
Rarea <- data$RNF * data$Narea
reg <- lm( data$LMA ~ Rarea )
plot(  data$LMA ~ Rarea, ylim=c(0,max(data$LMA, na.rm=TRUE)), xlim=c(0,max(Rarea, na.rm=TRUE)), xlab="rLMA (mol Rubisco N m-2)", ylab="LMA (g m-2)" )
points( nr_leaf_range, lma_range, col="red" )
text( 0, 0.03, paste("y-intersect:", reg$coefficients[1] ), adj=c(0,0))
text( 0, 0.01, paste("slope:      ", reg$coefficients[2] ), adj=c(0,0))
abline( reg )



# ## Canopy N as a function of LAI
# plot( alpha_range, nr_canop_range, type="l", xlab="LAI", ylab="Nr canopy (g N m-2)" )

## Leaf Rubisco-N as a function of LAI, overlay data: RNF * Narea (Rubisco N fraction * Narea = Rubisco N m-2)
plot( alpha_range, nr_leaf_range, type="l", xlab="LAI", ylab="Nr leaf (g N m-2)" )
abline( h= (data$RNF * data$Narea ), col=rgb(0,0,0,0.3) )

## Leaf cell wall-N as a function of LAI
plot( alpha_range, ncw_leaf_range, type="l", xlab="LAI", ylab="Ncw leaf (g N m-2)" )
# ## overlay data: WN (cell wall N per area, N m-2) xxx DOES NOT WORK PROPERLY: THERE MUST BE AN ADDITIONAL COMPONENT OF LEAF N XXX
# abline( h= (data$WN ), col=rgb(0,0,0,0.3) )
## overlay data: Narea - data$RNF * data$Narea (this is what n_cw really represents)
abline( h = data$Narea - data$RNF * data$Narea, col=rgb(0,0,0,0.3) )

## Leaf-N as a function of LAI, overlay data: WN (cell wall N per area, N m-2)
plot( alpha_range, ncw_leaf_range+nr_leaf_range, type="l", xlab="LAI", ylab="Narea (g N m-2)", ylim=range(data$Narea,ncw_leaf_range+nr_leaf_range) )
abline( h= (data$Narea ), col=rgb(0,0,0,0.3) )

## LMA as a function of LAI, overlay data: LMA (g m-2)
plot( alpha_range, lma_range, type="l", xlab="LAI", ylab="LMA (gC m-2)", ylim=c(range(data$LMA, lma_range)))
abline( h= (data$LMA ), col=rgb(0,0,0,0.3) )


# ## Canopy C as a function of LAI
# plot( alpha_range, c_canop_range, type="l", xlab="LAI", ylab="canopy C (g C m-2)"  )


# plot( alpha_range, taa_range, type="l" )
# plot( alpha_range, tba_range, type="l" )
# plot( alpha_range, beta_range, type="l" )
# # plot( beta_range, nup_range, type="l" )
# plot( alpha_range, aresp_range, type="l" )
# plot( alpha_range, bresp_range, type="l" )
# plot( alpha_range, aresp_range + bresp_range, type="l" )

# plot( alpha_range, npp_range, type="l" )
# plot( alpha_range, exu_range/npp_range, type="l" )

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


