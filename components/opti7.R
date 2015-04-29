##----------------------------------------------------------------
## DAILY GROWTH/ALLOCATION
## No prescribed phenology. Leaf/root allocation is optimised to
## satisfy C and N budgets daily.
##----------------------------------------------------------------

source('/alphadata01/bstocker/sofun/trunk/components/pmodel.R')
source('/alphadata01/bstocker/sofun/trunk/components/sofun.R')
source('/alphadata01/bstocker/utilities/daily2monthly.R')

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
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

## Input: inorganic N availability
## xxx as it stands now, it's extremely sensitive to inorganic N 
## availability. Linear N uptake function may not be appropriate.
Nscale <- 0.1
ninorg <- Nscale * sapply( days, FUN=sin_season )
ntrans <- 5.0  # free N taken up by transpiration stream

## Assumption: root phenology is congruent with 
## seasonality of inorganic N availability
root_season <- sapply( days, FUN=sin_season )

## Set fixed boundary conditions
elv    <- 450.0
co2    <- 376.0
cpalpha <- rep(1.26,nmonth)

## PPFD FROM SOFUN (STASH implementation)
## read PPFD calculated by STASH (sofun implementation, '*.m.qm.out')
filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.m.qm.out"
df.ppfd <- read.table( filnam, col.names=c("year","ppfd") )
istart <- which.min( abs(df.ppfd$year-2002.0) )
df.ppfd <- df.ppfd[ istart:(istart+nmonth-1), ]
mppfd  <- df.ppfd$ppfd

filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.d.qn.out"
df.ppfd <- read.table( filnam, col.names=c("year","ppfd") )
istart <- which.min( abs(df.ppfd$year-2002.0) )
df.ppfd <- df.ppfd[ istart:(istart+ndayyear-1), ]
dppfd  <- df.ppfd$ppfd

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

##----------------------------------------------------------------
## PARAMETERS
##----------------------------------------------------------------
## Model parameters, from Li et al., 2014 (where given)
params <- list( 
              ndayyear = 365,
              nmonth   = 12,
              y        = 0.6,   # yield factor in Li et al., 2014
              # r_root   = 0.005,  # yields ~0.913 
              # r_root   = 0.05, nice intersect
              r_root   = 0.005,
              r_leaf   = 0.0,    # Rd is directly substracted from GPP
              # exu      = 0.0001,
              # exu      = 0.0004, nice intersect
              exu      = 0.01,
              k_root   = 1/(1.04*365),
              k_leaf   = 1/365,  # xxx changed from 1/(4.0*365)
              psi      = 100.0,
              kbeer    = 0.5,
              c_content_of_biomass = 0.46, # McMurtrie & Dewar, 2011
              # lma      = 30,     # g C m-2; is mean of Table S1 LMA in Hikosaka& ... divided by 2 (conversion from biomass to C)
              # lma      = 60,     # g C m-2; SLA = 8.3 m2 kg-1, Vile et al., 2005 (Ann. Botany; Table 2, Trees), 0.5 g C / g biomss 
              sla      = 1/30,  # 0.0083: Vile et al., 2005 (Ann. Botany; Table 2, Trees)     xxx changed from 0.0014
              r_cton_leaf = 19,
              r_cton_root = 50,
              r_ntoc_leaf = 1/19,
              r_ntoc_root = 1/50,
              c_molmass= 12.0107,  # g C / mol C
              kphio    = 0.093,
              mol_weight_n = 14.0067,  # molecular weight of N, (g N)(mol N)-1
              r_n_cw_v = 1.23223,  #0.32  # 0.32 is mean in Hikosaka data as WNF/RNF
              ncw_min = 0.056,  # from Hikosaka analysis
              mol_weight_rubisco = 5.5e5,    # molecular weight of Rubisco, (g R)(mol R)-1
              n_conc_rubisco     = 1.14e-2,  # N concentration in rubisco, (mol N)(g R)-1
              cat_turnover_per_site = 2.33,  # catalytic turnover rate per site at 25 deg C, (mol CO2)(mol R sites)-1; use 2.33 instead of (3.5) as not all Rubisco is active (see Harrison et al., 2009)  
              cat_sites_per_mol_R   = 8.0,   # number of catalytic sites per mol R, (mol R sites)(mol R)-1
              n_v =  5.5e5 * 1.14e-2 / ( 2.33 * 8.0 ) # mol_weight_rubisco * n_conc_rubisco / ( cat_turnover_per_site * cat_sites_per_mol_R )
              )

##----------------------------------------------------------------
## FUNCTIONS
##----------------------------------------------------------------

calc_fapar <- function( lai, params ){
  ## Beer's Law
  with( params,{
    fapar <- ( 1.0 - exp( -kbeer * lai ) )
    return(fapar)
  })
}

calc_dgpp <- function( cleaf, mluenet, dppfd, params ){

  ## Returns annual "net" GPP (= GPP - Rd) as a function of LAI (alpha)
  with( params, {

    lai    <- cleaf * sla
    fapar  <- calc_fapar( lai, params )

    ## calculate monthly gpp vector, convert from mol/m2/month to gC/m2/month
    dgpp_net <- sum( dppfd * fapar * mluenet ) * c_molmass 

    return(dgpp_net)  
  })
}

calc_dnup <- function( croot, n0, params ){
  # This follows from FUN approach with
  # dCexu/dNup = K / (N0 - Nup); K=1/psi
  with( params,{
    # print(paste("psi",psi))
    # print(paste("cexu",cexu))
    cexu <- exu * croot
    nup <- n0 * ( 1.0 - exp( - psi * cexu ) )
    # print(paste("fraction N taken up:", nup/n0 ))`
    return(nup)
  })
}



# calc_ra <- function( cleaf, croot, params ){
#   with( params, {

#     ra <- 

#   })
# }


eval_imbalance <- function( dcleaf, cleaf, nleaf, croot, nroot, clabl, nlabl, mluenet, dppfd, ninorg, params ){
  ## /////////////////////////////////////////////////////////
  ## Evaluates C:N ratio of new assimilation after allocation 
  ## versus whole-plant C:N ratio after allocation. Optimal 
  ## allocation is where the two are equal. 
  ## ---------------------------------------------------------

  with( params, {

    ## Allocate
    dnleaf <- dcleaf * params$r_ntoc_leaf
    clabl  <- clabl - 1.0 / params$y * dcleaf
    nlabl  <- nlabl - dnleaf
    cleaf  <- cleaf + dcleaf
    nleaf  <- nleaf + dnleaf
    
    dcroot <- min( params$y * clabl, params$r_cton_root * nlabl )
    dnroot <- dcroot * params$r_ntoc_root
    clabl  <- clabl - 1.0 / params$y * dcroot
    nlabl  <- nlabl - dnroot
    croot  <- croot + dcroot
    nroot  <- nroot + dnroot

    ## Allocation and decay
    cleaf <- cleaf * (1.0 - k_leaf)
    croot <- croot * (1.0 - k_root)
    nleaf <- nleaf * (1.0 - k_leaf)
    nroot <- nroot * (1.0 - k_root)

    ## Calculate next day's C and N return after assumed allocation (tissue turnover happens before!)
    dc <- calc_dgpp( cleaf, mluenet, dppfd, params ) - croot * (r_root + exu)
    dn <- calc_dnup( croot, ninorg, params )

    ## Evaluation quantity is the difference between the 
    ## C:N ratio of new assimilates and the C:N ratio 
    ## of the whole plant after allocation.
    # eval <- (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )

    ## or should it be: 
    eval <- y * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )

    return( eval )

  })
}

eval_const_cost <- function( croot, cost, n0, params){
  with( params, {

    eval <- cost * n0 * ( 1.0 - exp( - psi * exu * croot ) ) - exu * croot

    return(eval)  
  })
}

eval_cost <- function( croot, n0, params){
  with( params, {
    cost <- ( exu * croot ) / ( n0 * (1.0 - exp( - psi * exu * croot ) ) )
    return(cost)  
  })
}

##----------------------------------------------------------------
## "MAIN" PROGRAM
##----------------------------------------------------------------

## Run P-model for each month 

mluenet <- rep( NA, nmonth )
mlue    <- rep( NA, nmonth )
mnapar  <- rep( NA, nmonth )
factor25<- rep( NA, nmonth )

for (moy in 1:nmonth){

  ## Execute P-model
  out <- pmodel( fpar=NA, ppfd=NA, co2, mtemp[moy], cpalpha[moy], mvpd[moy], elv )

  ## Light use efficiency: (gpp - rd) per unit light absorbed
  mlue[moy] <- out$lue

  ## Net light use efficiency: (gpp - rd) per unit light absorbed
  mluenet[moy] <- out$luenet

  ## conversion factor to get from APAR to Rubisco-N
  mnapar[moy]  <- out$n_apar

  ## factor to convert from 25 deg-normalised to ambient T
  factor25[moy] <- out$factor25_vcmax

}

## Calculate Rubisco-N
daysecs   <- df.dl$dayl_h * 60 * 60  # number of daylight seconds in a year
monsecs   <- daily2monthly( daysecs, method="sum" )
meanmppfd <- mppfd / monsecs         # mol m-2 s-1


##----------------------------------------------------------------
## LOOP THROUGH DAYS
##----------------------------------------------------------------
cleaf0 <- 0.1
croot0 <- 0.1
nlabl0 <- 0.0
clabl0 <- 0.0

out_cton_labl <- rep( NA, ndayyear )
out_cleaf     <- rep( NA, ndayyear )
out_nleaf     <- rep( NA, ndayyear )
out_clabl     <- rep( NA, ndayyear )
out_nlabl     <- rep( NA, ndayyear )
out_croot     <- rep( NA, ndayyear )
out_nroot     <- rep( NA, ndayyear )
out_dcleaf    <- rep( NA, ndayyear )
out_dnleaf    <- rep( NA, ndayyear )
out_dcroot    <- rep( NA, ndayyear )
out_dnroot    <- rep( NA, ndayyear )
out_lai       <- rep( NA, ndayyear )
out_ncost     <- rep( NA, ndayyear )
out_nup       <- rep( NA, ndayyear )
out_dclabl    <- rep( NA, ndayyear )
out_gpp       <- rep( NA, ndayyear )

r_cton_leaf <- 19
r_cton_root <- 50

cleaf <- cleaf0
nleaf <- cleaf / r_cton_leaf
croot <- croot0
nroot <- croot / r_cton_root
nlabl <- nlabl0
clabl <- clabl0

mess <- "all ok"
const_cost <- FALSE
deactivate <- FALSE

doy <- 0
for (moy in 1:nmonth){
  for (dm in 1:ndaymonth[moy]){
    doy <- doy + 1 
    doy <- min( 364, doy )

    # ## Two modes of optimisation: 
    # ## 1. balanced C:N. During early season, maximises growth
    # ## 2. constant C cost of N uptake. During late season
    # if (doy>1 && ninorg[doy]<ninorg[doy-1]) {
    #   ## Optimisation by constant C cost of N uptake
    #   const_cost <- TRUE
    # }

    ## Continuous root turnover
    croot <- croot * ( 1.0 - params$k_root )
    nroot <- nroot * ( 1.0 - params$k_root )            

    ## Gross primary production minus leaf respiration
    gpp_net <- calc_dgpp( cleaf, mluenet[moy], dppfd[doy], params )
    out_gpp[doy] <- gpp_net
    # print(paste("gpp",gpp_net))

    # Root maintenance respiration
    rm_root <- croot * params$r_root
    # print(paste("rm_root",rm_root))

    ## Root exudation
    cexu <- croot * params$exu
    # print(paste("rm_root+cexu",rm_root+cexu))

    ## # Root growth respiration of root growth to maintain root mass
    ## rg_root <- params$y * params$k_root * croot 
    # print(paste("rg_root",rg_root))

    ## Add remainder C to labile pool
    dclabl <- gpp_net - cexu - rm_root
    clabl  <- clabl + dclabl
    out_dclabl[doy] <- dclabl

    if ((gpp_net - cexu - rm_root)<0.0) { mess <- "net C assimilation neg."; break }

    ## Nitrogen uptake
    nup <- calc_dnup( croot, ninorg[doy], params )
    # print(paste("Nup",nup))

    out_nup[doy] <- nup

    out_ncost[doy] <- cexu / nup

    print(paste("croot",croot))
    print(paste("cexu",cexu))
    print(paste("nup ",nup))
    print(paste("cost",out_ncost[doy]))

    ## Two modes of optimisation: 
    ## 1. balanced C:N. During early season, maximises growth
    ## 2. constant C cost of N uptake. During late season
    if (doy>30 && out_ncost[doy] > 50 && !const_cost ) {
      ## Optimisation by constant C cost of N uptake
      const_cost <- TRUE
    }


    ## Nitrogen addition to labile pool
    nlabl <- nlabl + nup

    if (nlabl<0) { mess <- "nlabl neg."; break }
    if (clabl<0) { mess <- "clabl neg."; break }

    print(paste("day of yr", doy))
    # print(paste("C labl. a", clabl))
    # print(paste("N labl.  ", nlabl))
    # print(paste("C:N labl.", clabl/nlabl))

    out_cton_labl[doy] <- clabl/nlabl

    ## Maximum is the lower of all labile C and the C to be matched by all labile N,
    ## discounted by the yield factor.
    max_dcleaf_n_constraint <- nlabl * r_cton_leaf 
    max_dcroot_n_constraint <- nlabl * r_cton_root
    max_dcleaf <- min( params$y * clabl, max_dcleaf_n_constraint )
    max_dcroot <- min( params$y * clabl, max_dcroot_n_constraint )

    # print(paste("C in labile pool      ", clabl ) )
    # print(paste("C avl. for growth     ", params$y * clabl ) )
    # print(paste("C for roots, avl. by N", nlabl * r_cton_root ) )
    # print(paste("C for leafs, avl. by N", nlabl * r_cton_leaf ) )
    # print(paste("max_dcleaf            ", max_dcleaf))

    # print(paste("clabl              ",clabl))
    # print(paste("nlabl * r_cton_leaf",nlabl * r_cton_leaf))

    if (const_cost) {

      ## Optimisation by constant C cost of N uptake
      # print(paste("croot",croot))
      if (ninorg[doy+1]>0.0){

        if ( croot==0.0 || eval_cost( 1e-12, ninorg[doy+1], params ) > out_ncost[doy] ) {

          croot_tgt <- 0.0

        } else {

          # ## test
          # croot_range <- seq( 0.0, croot, by=croot/100)
          # eval_range  <- sapply( croot_range, FUN = function(x) eval_cost( x, ninorg[doy+1], params ) )

          out.root  <- uniroot( function(x) eval_const_cost( x, out_ncost[doy], ninorg[doy+1], params ), interval=c(1e-12,croot) )        
          croot_tgt <- out.root$root

        }

      } else {

        croot_tgt <- 0.0

      }

      print(paste("cost today:   ",out_ncost[doy]))
      print(paste("cost tomorrow:",croot_tgt*params$exu/calc_dnup( croot_tgt, ninorg[doy+1],params)))
      print(paste("croot_tgt",croot_tgt))

      ## Find root C increment to satisfy targeted root mass after allocation and continuous decay
      dcroot <- croot_tgt / (1.0 - params$k_root ) - croot
      print(paste("dcroot",dcroot))

      if ( dcroot <= 0.0 ){
        kill_croot <- (-1.0) * dcroot
        deactivate <- TRUE
        dcroot <- 0.0
        print("deactivate")
      }

      if (deactivate){

        ## Deactivate root mass, use all assimilates for leaf growth
        croot  <- croot - kill_croot 
        nroot  <- croot * params$r_ntoc_root

        ## Stop leaf growth

        ## Use all assimilates for leaf growth
        dcleaf <- max_dcleaf
        dnleaf <- dcleaf * params$r_ntoc_leaf
        clabl  <- clabl - 1.0 / params$y * dcleaf
        nlabl  <- nlabl - dnleaf
        cleaf  <- cleaf + dcleaf
        nleaf  <- nleaf + dnleaf

      } else {

        ## Allocate based on 'dcroot' determined above
        dnroot <- dcroot * params$r_ntoc_root
        clabl  <- clabl - 1.0 / params$y * dcroot
        nlabl  <- nlabl - dnroot
        croot  <- croot + dcroot
        nroot  <- nroot + dnroot
        
        dcleaf <- min( params$y * clabl, params$r_cton_leaf * nlabl )

        dnleaf <- dcleaf * params$r_ntoc_leaf
        clabl  <- clabl - 1.0 / params$y * dcleaf
        nlabl  <- nlabl - dnleaf
        cleaf  <- cleaf + dcleaf
        nleaf  <- nleaf + dnleaf

      }

    } else {

      ## Optimisation by balanced growth
      ## Test I: Evaluate balance if all is put to roots.
      ## If C:N ratio of return is still greater than whole-plant C:N ratio, then put all to roots.
      findroot <- TRUE
      eval_allroots  <- eval_imbalance( 0.0, cleaf, nleaf, croot, nroot, clabl, nlabl, mluenet[moy], dppfd[doy+1], ninorg[doy+1], params )
      if (eval_allroots > 0.0) { dcleaf <- 0.0; findroot <- FALSE }

      ## Test II: Evaluate balance if all is put to leaves.
      ## If C:N ratio of return is still lower than whole-plant C:N ratio, then put all to leaves.
      eval_allleaves <- eval_imbalance( max_dcleaf, cleaf, nleaf, croot, nroot, clabl, nlabl, mluenet[moy], dppfd[doy+1], ninorg[doy+1], params )
      if (eval_allleaves < 0.0) { dcleaf <- max_dcleaf; findroot <- FALSE}

      if (findroot) {

        ## Find root
        print("finding root")
        out.root <- NA
        try ( 
          out.root <- uniroot( function(x) eval_imbalance( x, cleaf, nleaf, croot, nroot, clabl, nlabl, mluenet[moy], dppfd[doy], ninorg[doy], params ), interval=c(0,max_dcleaf) )
          )        
        if( is.na(out.root) ){
          dcleaf <- 0.0
        } else { 
          dcleaf <- out.root$root
        }
      }

      ## Allocate based on 'dcleaf' determined above
      dnleaf <- dcleaf * params$r_ntoc_leaf
      clabl  <- clabl - 1.0 / params$y * dcleaf
      nlabl  <- nlabl - dnleaf
      cleaf  <- cleaf + dcleaf
      nleaf  <- nleaf + dnleaf
      
      dcroot <- min( params$y * clabl, params$r_cton_root * nlabl )

      # print(paste("Allocation decision, dcroot:",dcroot))
      dnroot <- dcroot * params$r_ntoc_root
      clabl  <- clabl - 1.0 / params$y * dcroot
      nlabl  <- nlabl - dnroot
      croot  <- croot + dcroot
      nroot  <- nroot + dnroot

    }

    # print(paste("Allocation decision, dcleaf:",dcleaf))

    if ( clabl < -1e-12 ){
      print(paste("problem: neg. clabl:",clabl))
    }
    if ( nlabl < -1e-12 ){
      print(paste("problem: neg. nlabl",nlabl))
    }

    
    ## Write to output
    out_cleaf[doy]  <- cleaf
    out_nleaf[doy]  <- nleaf
    out_croot[doy]  <- croot
    out_nroot[doy]  <- nroot
    out_dcleaf[doy] <- dcleaf
    out_dcroot[doy] <- dcroot
    out_dnleaf[doy] <- dnleaf
    out_dnroot[doy] <- dnroot
    out_lai[doy]    <- cleaf * params$sla
    out_clabl[doy]  <- clabl
    out_nlabl[doy]  <- nlabl

  }
}

print(mess)

# plot( 1:doy, out_cton_labl, type="l", ylim=c(-80,80) )

# pdf( "cmass_vs_doy.pdf", width=6, height=5 )
plot( 1:doy, out_croot[1:doy], type="l", ylab="C mass (gC/m2)", xlab="DOY" )
lines( 1:doy, out_cleaf[1:doy], type="l", col="red" )
lines( 1:doy, out_clabl[1:doy], col="blue" )
legend( "topleft", c("root C","leaf C", "labile C"), lty=1, bty="n", col=c("black","red","blue") )
# dev.off()

plot( 1:doy, out_clabl[1:doy], type="l", ylim=range(c( out_clabl[1:doy], out_nlabl[1:doy]),na.rm=T) )
lines( 1:doy, out_nlabl[1:doy], col="red" )

plot( 1:doy, out_dcroot[1:doy], type="l" )
lines( 1:doy, out_dcleaf[1:doy], col="red" )

plot( 1:doy, out_cleaf[1:doy]/out_nleaf[1:doy], type="l" )
plot( 1:doy, out_croot[1:doy]/out_nroot[1:doy], type="l"  )

# pdf( "lai_vs_doy.pdf", width=6, height=5 )
plot( 1:doy, out_lai[1:doy], type="l", ylab="LAI", xlab="DOY")
# dev.off()

# pdf( "cost_of_n_vs_doy.pdf", width=6, height=5 )
plot( 1:doy, out_ncost[1:doy], type="l", ylim=c(0,200), ylab="C cost per N uptake (gC/gN)", xlab="DOY")
# dev.off()

# pdf( "nup_vs_doy.pdf", width=6, height=5 )
plot( 1:doy, out_nup[1:doy], type="l", ylab="N uptake (gN/day)", xlab="DOY")
# dev.off()

plot( 1:doy, out_dclabl[1:doy], type="l", ylab="C balance (gC/m2/day)", xlab="DOY")
plot( 1:doy, out_gpp[1:doy], type="l", ylab="GPP (gC/m2/day)", xlab="DOY")
