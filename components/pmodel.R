# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Contains P model functions adopted from gepisat
# ////////////////////////////////////////////////////////////////////////

kc <- 0.41          # Jmax cost coefficient
kphio <- 0.093      # quantum efficiency (Long et al., 1993)
kPo <- 101325.0     # standard atmosphere, Pa (Allen, 1973)
kTo <- 25.0         # base temperature, deg C (Prentice, unpublished)

sla <- 0.0014       # specific leaf area (m2/gC)

ncw <- 0.0          # N:C ratio in cell walls, working hypothesis: leaf N is solely determined by Vcmax25
nv  <- 1.0/40.96    # gN µmol-1 s-1. Value 40.96 is 'sv' in Table 2 in Kattge et al., 2009, GCB, C3 herbaceous

#-----------------------------------------------------------------------
# Email from Tyler (10.3.2015):
# I was estimating values of β based on the Wang Han approximation equation 
# of χ using both the simplified and "more precise" expressions for χ and ξ 
# (Prentice et al., 2014, Ecology Letters).  After examination, Colin and I 
# noticed that the value of β is not significantly influenced by the 
# expressions for χ and ξ. Since then, Colin has theorised the use of a 
# "ground state" universal value of β, which is derived from the Wang Han 
# equation at sea level (i.e., z = 0 m and Patm = 101325 Pa), standard temp-
# erature (i.e., Tair = 25 deg C) and a non-influencial VPD (i.e., 
# D = 1000 Pa). Based on these climatological values, the following were 
# calculated:
#   a. Γ* = 4.220 Pa
#   b. K = 70.842 Pa
#   c. η* = 1.0
#   d. χ = 0.767
#   e. β = 244.033
# Results from modelled versus "observed" monthly GPP based on the universal 
# value of β are promising. Colin and I are currently in the works on the next 
# set of improvements, which, as I far as I know, will be based on this uni-
# versal value of β.
#-----------------------------------------------------------------------
beta <- 244.033

params <- list( kc, kphio, kPo, kTo, beta )


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#  Function Definitions
# ////////////////////////////////////////////////////////////////////////
pmodel <- function( fpar, ppfd, co2, tc, cpalpha, vpd, elv, params, method="full" ){
  #-----------------------------------------------------------------------
  # Input:    - fpar (unitless)    : monthly fraction of absorbed photosynthetically active radiation
  #           - ppfd (mol/m2)      : monthly photon flux density
  #           - co2 (ppm)          : atmospheric CO2 concentration
  #           - tc (deg C)       : monthly air temperature
  #           - cpalpha (unitless, within [0,1.26]) : monthly Cramer-Prentice-alpha
  #           - vpd (Pa)           : mean monthly vapor pressure -- CRU data is in hPa
  #           - elv (m)            : elevation above sea-level
  # Output:   gpp (gC/m2/month)    : gross primary production
  #-----------------------------------------------------------------------

  with( params, {

    ## absorbed photosynthetically active radiation (mol/m2)
    iabs <- fpar * ppfd

    ## atmospheric pressure as a function of elevation (Pa)
    patm <- calc_patm( elv )

    ## ambient CO2 partial pression (Pa)
    ca   <- co2_to_ca( co2, patm )

    ## photorespiratory compensation point - Gamma-star (Pa)
    gs   <- calc_gstar_gepisat( tc )

    ## function of alpha to reduce GPP in strongly water-stressed months (unitless)
    fa   <- calc_fa( cpalpha )


    if (method=="approx"){
      ##-----------------------------------------------------------------------
      ## A. APPROXIMATIVE METHOD
      ##-----------------------------------------------------------------------

      lue.out <- lue_approx( tc, vpd, elv, ca, gs )

    } else {
      ##-----------------------------------------------------------------------
      ## B. THEORETICAL METHOD
      ##-----------------------------------------------------------------------

      ## Michaelis-Menten coef. (Pa)
      kmm  <- calc_k( tc, patm )

      ## viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
      ns      <- viscosity_h2o( tc, patm )  # Pa s 
      ns25    <- viscosity_h2o( kTo, kPo )  # Pa s 
      ns_star <- ns/ns25  # (unitless)

      if (method=="simpl") {

        ## B.1 SIMPLIFIED FORMULATION 
        ##-----------------------------------------------------------------------
        lue.out <- lue_vpd_simpl( kmm, gs, ns, ca, vpd, params  )

      } else if (method=="full"){

        ## B.2 FULL FORMULATION
        ##-----------------------------------------------------------------------
        lue.out <- lue_vpd_full( kmm, gs, ns_star, ca, vpd, params  )

      }

    }

    ## LUE-functions return m, n, and chi
    m <- lue.out$m
    n <- lue.out$n

    ##-----------------------------------------------------------------------
    ## Calculate function return variables
    ##-----------------------------------------------------------------------

    ## GPP per unit ground area is the product of the intrinsic quantum 
    ## efficiency, the absorbed PAR, the function of alpha (drought-reduction),
    ## and 'm'
    gpp <- iabs * kphio * fa * m 

    ## Vcmax per unit ground area is the product of the intrinsic quantum 
    ## efficiency, the absorbed PAR, and 'n'
    vcmax <- iabs * kphio * n

    ## Dark respiration
    rd <- 0.015 * vcmax

    ## Vcmax25 (vcmax normalized to 25 deg C)
    vcmax25 <- calc_vcmax25( vcmax, tc )

    ## Narea (nitrogen content per unit leaf area, gN/m2-leaf)
    narea <- ncw * 1.0/sla + nv * vcmax25


    out <- list( gpp=gpp, vcmax=vcmax, rd=rd, narea=narea )
    return( out )
  })
}


calc_leaf_N <- function( lai, narea ){

  nleaf <- lai * narea
  return(nleaf)
}


# calc_gpp_gepisat <- function( fpar, ppfd, vpd, alpha, tair, co2, patm, elv, nxtgn=FALSE ){
#   #-----------------------------------------------------------------------
#   # Input:    - float, monthly GPP, mol/m2 (gpp)
#   #           - float, associated GPP error, mol/m2 (gpp_err)
#   #           - float, monthly FAPAR (fpar)
#   #           - float, monthly PPFD, mol/m2 (ppfd)
#   #           - float, monthly VPD, Pa (vpd)
#   #           - float, monthly CPA (alpha)
#   #           - float, monthly air temp, degC (tmp)
#   #           - float, annual atm. CO2, ppm (co2)
#   #           - float, monthly atm. pressure, Pa (patm)
#   #           - float, elevation, m (elv)
#   # Output:   None.
#   # Features: Appends a set of monthly values to the value dictionary
#   # Depends:  - calc_gstar
#   #           - calc_k
#   #           - nxgn
#   #           - viscosity_h2o
#   #           - kPo
#   #-----------------------------------------------------------------------

#   iabs <- fpar * ppfd                    # mol/m2, abs. PPFD
#   ca   <- co2_to_ca( co2, patm )         # Pa, atms. CO2
#   gs   <- calc_gstar_gepisat( tair )     # Pa, photores. comp. point

#   if (nxtgn) {

#     ## Formulation based on theoretical model
#     k    <- calc_k( tair, patm )           # Pa, Michaelis-Menten coef.
#     ns   <- viscosity_h2o( tair, patm ) / n25 # Pa s, viscosity
#     fa   <- ( alpha / 1.26 )^(0.25)       # unitless, func. of alpha

#     out.beta <- beta_estimate(ca, vpd, k, gs, ns, tair, elv)
#     beta1 <- out.beta[1]
#     beta2 <- out.beta[2]

#     # m <- calc_m_nxtgn(ca, gs, vpd, k, ns, fa, beta1)
#     m <- calc_m_nxtgn(ca, gs, vpd, k, ns, fa, beta2)

#   } else {

#     ## use approximation function for chi
#     atilde <- 0.8

#     ## get chi using the approximative function based on Wang Han's equation
#     chi <- lue_approx( tair, vpd, elv, ca, gs )$chi

#     ## substitution
#     gamma <- gs / ca
#     m <- atilde * ( chi - gamma ) / ( chi + 2.0 * gamma )

#   }

#   ## GPP is the productof the intrinsic quantum efficiency, the absorbed PAR, and the light use efficiency
#   gpp <- iabs * kphio * m
#   return( gpp )
  
# }


calc_m_nxtgn <- function(ca, gs, d, k, ns, beta){
  #-----------------------------------------------------------------------
  # Input:    - float, 'ca' : Pa, atmospheric CO2
  #           - float, 'Gs' : Pa, photores. comp. point (Gamma-star)
  #           - float, 'D'  : Pa, vapor pressure deficit
  #           - float, 'K'  : Pa, Michaelis-Menten coeff.
  #           - float, 'ns' : Pa s, viscosity of water
  #           - float, beta parameter (beta)
  # Output:   float, estimate of GPP (gpp)
  # Features: Returns an estimate of GPP based on the next-generation light 
  #           and water use efficiency model. This corresponds to Eq.6 in 
  #           Simplifying_m.pdf with b/(a*vpd) = beta
  # Depends:  - kc
  #           - kphio
  #-----------------------------------------------------------------------

  # Define variable substitutes:
  ns   <- ns * 1.0e3 # convert from Pa s to mPa s
  vdcg <- ca - gs
  vacg <- ca + 2.0 * gs
  vbkg <- beta * (k + gs)

  # Check for negatives:
  if (vbkg > 0){
    vsr <- sqrt(1.6*ns*d/(vbkg))

    # Based on the m' formulation (see Regressing_LUE.pdf)
    m <- vdcg/(vacg + 3.0*gs*vsr)

    mpi <- m^2 - kc^(2.0/3.0)*(m^(4.0/3.0))

    # Check for negatives:
    if (mpi > 0){
      mp <- sqrt(mpi)
      # gpp <- kphio*iabs*fa*mp
      # mp <- mp * fa
    }
  }

  return( mp )

}

lue_approx <- function( temp, vpd, elv, ca, gs ){
  #-----------------------------------------------------------------------
  # Input:    - float, 'temp' : deg C, air temperature
  #           - float, 'vpd'  : Pa, vapour pressure deficit
  #           - float, 'elv'  : m, elevation above sea level
  #           - float, 'ca'   : Pa, ambient CO2 partial pressure
  #           - float, 'gs'   : Pa, photores. comp. point (Gamma-star)
  # Output:   list: 'm' (unitless), 'chi' (unitless)
  # Features: Returns list containing light use efficiency (m) and ci/ci ratio (chi)
  ##          based on the approximation of the theoretical relationships
  ##          of chi with temp, vpd, and elevation
  #-----------------------------------------------------------------------

  ## Wang-Han Equation
  whe <- exp( 
    1.19 
    + 0.0545 * ( temp - 25.0 )
    - 0.5 * log( 1e-3 * vpd )    # convert vpd from Pa to kPa 
    - 0.0815 * ( 1e-3 * elv )    # convert elv from m to km
    )

  ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  chi <- whe / ( 1.0 + whe )

  ##  m
  gamma <- gs / ca
  m <- (chi - gamma) / (chi + 2 * gamma)

  out <- list( chi=chi, m=m, n=NA )
  return(out)
}


lue_vpd_simpl <- function( kmm, gs, ns_star, ca, vpd, params ){
  #-----------------------------------------------------------------------
  # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
  #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
  #           - float, 'vpd' : Pa, vapor pressure deficit
  # Output:   float, ratio of ci/ca (chi)
  # Features: Returns an estimate of leaf internal to ambient CO2
  #           partial pressure following the "simple formulation".
  # Depends:  - kc
  #           - ns
  #           - vpd
  #-----------------------------------------------------------------------
  with( params, {

    ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  <- sqrt( beta * kmm / (1.6 * ns_star))
    chi <- xi / (xi + sqrt(vpd))

    ## light use efficiency (m)
    ## consistent with this, directly return light-use-efficiency (m)
    m <- ( xi * (ca - gs) - gs * sqrt( vpd ) ) / ( xi * (ca + 2.0 * gs) + 2.0 * gs * sqrt( vpd ) )

    ## n 
    gamma <- gs / ca
    kappa <- kmm / ca
    n <- (chi + kappa) / (chi + 2 * gamma)

    out <- list( chi=chi, m=m, n=n )
    return(out)
  })
}


lue_vpd_full <- function( kmm, gs, ns_star, ca, vpd, params ){
  #-----------------------------------------------------------------------
  # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
  #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
  #           - float, 'vpd' : Pa, vapor pressure deficit
  # Output:   float, ratio of ci/ca (chi)
  # Features: Returns an estimate of leaf internal to ambient CO2
  #           partial pressure following the "simple formulation".
  # Depends:  - kc
  #           - ns
  #           - vpd
  #-----------------------------------------------------------------------
  with( params, {

    ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  <- sqrt( (beta * ( kmm + gs ) ) / ( 1.6 * ns_star ) )
    chi <- gs / ca + ( 1.0 - gs / ca ) * xi / ( xi + sqrt(vpd) )

    ## consistent with this, directly return light-use-efficiency (m)
    ## see Eq. 13 in 'Simplifying_LUE.pdf'

    ## light use efficiency (m)
    # m <- (ca - gs)/(ca + 2.0 * gs + 3.0 * gs * sqrt( (1.6 * vpd) / (beta * (K + gs) / ns_star ) ) )

    # Define variable substitutes:
    vdcg <- ca - gs
    vacg <- ca + 2.0 * gs
    vbkg <- beta * (kmm + gs)

    # Check for negatives:
    if (vbkg > 0){
      vsr <- sqrt( 1.6 * ns_star * vpd / vbkg )

      # Based on the m' formulation (see Regressing_LUE.pdf)
      m <- vdcg / ( vacg + 3.0 * gs * vsr )
    }

    ## n 
    gamma <- gs / ca
    kappa <- kmm / ca
    n <- (chi + kappa) / (chi + 2 * gamma)


    out <- list( chi=chi, m=m, n=n )
    return(out)
  })
}


mprime <- function( m, params ){
  with( params,{
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    mpi <- m^2 - kc^(2.0/3.0) * (m^(4.0/3.0))

    # Check for negatives:
    if (mpi > 0){ mp <- sqrt(mpi) }
    return(mpi)    
  })
}


calc_fa <- function( cpalpha ){
  #-----------------------------------------------------------------------
  # Input:  cpalpha (unitless, within [0,1.26]): monthly Cramer-Prentice-alpha
  # Output: fa (unitless, within [0,1]): function of alpha to reduce GPP 
  #                                      in strongly water-stressed months
  #-----------------------------------------------------------------------
  fa <- ( cpalpha / 1.26 )^(0.25)
  return(fa)
}


co2_to_ca <- function( co2, patm ){
  #-----------------------------------------------------------------------
  # Input:    - float, annual atm. CO2, ppm (co2)
  #           - float, monthly atm. pressure, Pa (patm)
  # Output:   - ca in units of Pa
  # Features: Converts ca (ambient CO2) from ppm to Pa.
  #-----------------------------------------------------------------------
  ca   <- ( 1.e-6 ) * co2 * patm         # Pa, atms. CO2
  return( ca )
}


ca_to_co2 <- function( ca, patm ){
  #-----------------------------------------------------------------------
  # Input:    - float, ambient CO2, Pa (ca)
  #           - float, monthly atm. pressure, Pa (patm)
  # Output:   - co2 in units of Pa
  # Features: Converts ca (ambient CO2) from Pa to ppm.
  #-----------------------------------------------------------------------
  co2   <- ca * ( 1.e6 ) / patm
  return( co2 )
}




# temp_range <- seq( 0, 50, 1)
# gstar_gepisat <- sapply( temp_range, FUN = calc_gstar_gepisat)
# gstar_wh <- sapply( temp_range, FUN = calc_gstar_wh)
# gstar_colin <- sapply( temp_range, FUN = calc_gstar_colin)

# plot( temp_range, gstar_gepisat, type="l"  )
# lines( temp_range, gstar_wh, col="red" )
# lines( temp_range, gstar_colin, col="blue" )

beta_estimate <- function(my_ca, my_d, my_k, my_gs, my_ns, my_t, my_z) {
  #-----------------------------------------------------------------------
  # Input:    - float, atmospheric CO2 concentration, Pa (my_ca)
  #           - float, vapor pressure deficit, Pa (my_d)
  #           - float, Michaelis-Menten coeff, Pa (my_k)
  #           - float, photorespiratory comp point, Pa (my_gs)
  #           - float, viscosity, unitless (my_ns)
  #           - float, air temperature, deg C (my_t)
  #           - float, elevation, m (my_z)
  # Output:   tuple 
  #           - float, predicted beta from simple expression (beta_p1)
  #           - float, predicted beta from 
  # Features: Returns an estimate for beta based on the Wang Han equation.
  #           See also 'Estimation_of_beta.pdf'
  #-----------------------------------------------------------------------

  chi <- lue_approx( my_t, my_d, my_z, my_ca, my_gs )$chi

  # beta following Estimation_of_beta.pdf, p.5 (Method 2, the simple expression)
  beta_p1 <- 1.6*my_ns*my_d*(chi^2)
  beta_p1 <- beta_p1 / (1.0 - chi)^2
  beta_p1 <- beta_p1 / my_k

  # beta following Estimation_of_beta.pdf, p.4
  beta_p2 <- 1.6*my_ns*my_d
  beta_p2 <- beta_p2 / (my_k + my_gs)
  beta_p2 <- beta_p2 * (chi*my_ca - my_gs)^2
  beta_p2 <- beta_p2 / (my_ca^2)
  beta_p2 <- beta_p2 / (chi - 1.0)^2

  return (beta_p1, beta_p2)
}


calc_vpd <- function( temp, vap, tmin=NA, tmax=NA ){
  #-----------------------------------------------------------------------
  # Input:    - mean monthly temperature, deg C (temp)
  #           - mean monthly vapor pressure, hPa (vap) -- CRU data is in hPa
  #           - (optional) mean monthly min daily air temp, deg C (tmin)
  #           - (optional) mean monthly max daily air temp, deg C (tmax)
  # Output:   mean monthly vapor pressure deficit, Pa (vpd)
  # Features: Returns mean monthly vapor pressure deficit
  # Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure 
  #           Calculation Methods, in Evaporation and Evapotranspiration: 
  #           Measurements and Estimations, Springer, London.
  #             vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
  #             where:
  #                 temp = average daily air temperature, deg C
  #                 vap  = actual vapor pressure, kPa
  #-----------------------------------------------------------------------
  if ( !is.na(tmin) && !is.na(tmax) ){
    temp <- 0.5 * (tmin + tmax)
  }

  ## calculate VPD in units of kPa
  vpd <- ( 0.611 * exp( (17.27 * temp)/(temp + 237.3) ) - 0.10 * vap )    

  ## convert to Pa
  vpd <- vpd * 1000

  return( vpd )
}


calc_k <- function(tc, patm) {
  #-----------------------------------------------------------------------
  # Input:    - float, air temperature, deg C (temp)
  #           - float, atmospheric pressure, Pa (patm)
  # Output:   float, Pa (mmk)
  # Features: Returns the temperature & pressure dependent Michaelis-Menten
  #           coefficient, K (Pa).
  # Ref:      Bernacchi et al. (2001), Improved temperature response 
  #           functions for models of Rubisco-limited photosynthesis, 
  #           Plant, Cell and Environment, 24, 253--259.
  #-----------------------------------------------------------------------

  kc25 <- 39.97      # Pa, assuming 25 deg C & 98.716 kPa
  ko25 <- 2.748e4    # Pa, assuming 25 deg C & 98.716 kPa
  dhac <- 79430      # J/mol
  dhao <- 36380      # J/mol
  kR   <- 8.3145     # J/mol/K
  kco  <- 2.09476e5  # ppm, US Standard Atmosphere

  kc <- kc25 * exp( dhac * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 
  ko <- ko25 * exp( dhao * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 

  po <- kco * (1e-6) * patm # O2 partial pressure
  k  <- kc * (1.0 + po/ko)

  return(k)
}


calc_k_colin <- function( tc ) {
  #-----------------------------------------------------------------------
  # Input:    - float, air temperature, deg C (tc)
  # Output:   float, Pa (mmk)
  # Features: Returns the temperature & pressure dependent Michaelis-Menten
  #           coefficient, K (Pa).
  # Ref:      Colin's documents
  #-----------------------------------------------------------------------

  ## conversion to temperature in Kelvin
  tk <- tc + 273.15

  kc25 <- 41.03      # Pa, assuming 25 deg C & 98.716 kPa
  ko25 <- 28210      # Pa
  dhac <- 79430      # J/mol
  dhao <- 36380      # J/mol
  kR   <- 8.3145     # J/mol/K
  kco  <- 2.09476e5  # ppm, US Standard Atmosphere

  kc <- kc25 * exp( dhac / kR * (1.0/298.15 - 1.0/tk) )
  ko <- ko25 * exp( dhao / kR * (1.0/298.15 - 1.0/tk) )

  po <- kco * (1e-6) * patm # O2 partial pressure
  k  <- kc * (1.0 + po/ko )

  return(k)
}

calc_gstar_gepisat <- function( tc ) {
  #-----------------------------------------------------------------------
  # Input:    float, air temperature, degrees C (tc)
  # Output:   float, gamma-star, Pa (gs)
  # Features: Returns the temperature-dependent photorespiratory 
  #           compensation point, Gamma star (Pascals), based on constants 
  #           derived from Bernacchi et al. (2001) study.
  # Ref:      Bernacchi et al. (2001), Improved temperature response 
  #           functions for models of Rubisco-limited photosynthesis, 
  #           Plant, Cell and Environment, 24, 253--259.
  #-----------------------------------------------------------------------

  gs25 <- 4.220    # Pa, assuming 25 deg C & 98.716 kPa)
  dha  <- 37830    # J/mol
  kR   <- 8.3145   # J/mol/K

  gs <- gs25 * exp( dha * ( tc - 25.0 ) / ( 298.15 * kR * ( tc + 273.15 ) ) )

  return( gs )
}


calc_gstar_wh <- function( temp ){
  #-----------------------------------------------------------------------
  ## Returns the temperature-dependent photorespiratory compensation point, 
  ## Gamma star (Pascals), based on constants derived from 
  ## Bernacchi et al. (2001) study.
  #-----------------------------------------------------------------------

  ## wang han's values
  gs25 <- 42.75  # corresponds to what Colin gave me
  k    <- 0.0512

  gs <- gs25 * exp( k * ( temp - 25 ) )

  return( gs )

}

calc_gstar_colin <- function( tc ){
  #-----------------------------------------------------------------------
  # Input:    tc: air temperature (degrees C)
  # Output:   gs: gamma-star (Pa)
  # Features: Returns the temperature-dependent photorespiratory 
  #           compensation point, Gamma star (Pascals), based on constants 
  #           derived from Bernacchi et al. (2001) study.
  # Ref:      Colin's document
  #-----------------------------------------------------------------------

  ## conversion to temperature in Kelvin
  tk <- tc + 273.15

  gs25 <- 4.275    # corresponds to what Colin gave me
  kR   <- 8.3145   # J/mol/K
  dha  <- 37830    # J/mol

  gs <- gs25 * exp( ( dha / kR ) * ( 1/298.15 - 1.0/tk ) )
  
  return( gs )
}


calc_vcmax25 <- function( vcmax, tc ){
  #-----------------------------------------------------------------------
  # Input:    - gcmax  : Vcmax at a given temperature tc 
  #           - tc     : air temperature (degrees C)
  # Output:   vcmax25  : Vcmax at 25 deg C
  # Features: Returns the temperature-corrected Vcmax at 25 deg C
  # Ref:      Analogue function like 'calc_gstar_gepisat'
  #-----------------------------------------------------------------------

  dhav <- 65330    # J/mol
  kR   <- 8.3145   # J/mol/K

  vcmax25 <- vcmax * exp( -dhav * ( tc - 25.0 ) / ( 298.15 * kR * ( tc + 273.15 ) ) )
  return( vcmax25 )
}


calc_vcmax25_colin <- function( vcmax, tc ){
  #-----------------------------------------------------------------------
  # Input:    - gcmax  : Vcmax at a given temperature tc 
  #           - tc     : air temperature (degrees C)
  # Output:   vcmax25  : Vcmax at 25 deg C
  # Features: Returns the temperature-corrected Vcmax at 25 deg C
  # Ref:      Colin's document
  #-----------------------------------------------------------------------

  ## conversion to temperature in Kelvin
  tk <- tc + 273.15

  dhav <- 65330    # J/mol
  kR   <- 8.3145   # J/mol/K

  vcmax25 <- vcmax * exp( -dhav/kR * (1/298.15 - 1/tk) )
  return( vcmax25 )
}


calc_patm <- function( elv ){
  #-----------------------------------------------------------------------
  # Input:    - elevation, m (elv)
  # Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
  # Features: Returns the atmospheric pressure as a function of elevation
  #           and standard atmosphere (1013.25 hPa)
  # Depends:  - connect_sql
  #           - flux_to_grid
  #           - get_data_point
  #           - get_msvidx
  # Ref:      Allen et al. (1998)
  #-----------------------------------------------------------------------

  # Define constants:
  kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo <- 298.15   # base temperature, K (Prentice, unpublished)
  kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
  kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)

  # Convert elevation to pressure, Pa:
  patm <- kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
  
  return (patm)
}


density_h2o <- function( tc, p ){
  #-----------------------------------------------------------------------
  # Input:    - float, air temperature (tc), degrees C
  #           - float, atmospheric pressure (p), Pa
  # Output:   float, density of water, kg/m^3
  # Features: Calculates density of water at a given temperature and 
  #           pressure using the Tumlirz Equation
  # Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of 
  #           pure water and sea water, Tech. Rept., Marine Physical 
  #           Laboratory, San Diego, CA.
  #-----------------------------------------------------------------------

  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 + 
          21.55053*tc + 
        -0.4695911*tc*tc + 
     (3.096363e-3)*tc*tc*tc + 
    -(7.341182e-6)*tc*tc*tc*tc

  # Calculate po, bar
  po <- 5918.499 + 
           58.05267*tc + 
         -1.1253317*tc*tc + 
     (6.6123869e-3)*tc*tc*tc + 
    -(1.4661625e-5)*tc*tc*tc*tc

  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
     (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
     (9.829576e-9)*tc*tc*tc*tc +
   -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
   -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
     (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
   -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p
  
  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)

  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)

  return(rho)
}


viscosity_h2o <- function( tc, p ) {
  #-----------------------------------------------------------------------
  # Input:    - float, ambient temperature (tc), degrees C
  #           - float, ambient pressure (p), Pa
  # Return:   float, viscosity of water (mu), Pa s
  # Features: Calculates viscosity of water at a given temperature and 
  #           pressure.
  # Depends:  density_h2o
  # Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
  #           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
  #           international formulation for the viscosity of H2O, J. Phys. 
  #           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
  #-----------------------------------------------------------------------

  # Define reference temperature, density, and pressure values:
  tk_ast  <- 647.096    # Kelvin
  rho_ast <- 322.0      # kg/m^3
  mu_ast  <- 1e-6       # Pa s

  # Get the density of water, kg/m^3
  rho <- density_h2o(tc, p)

  # Calculate dimensionless parameters:
  tbar <- (tc + 273.15)/tk_ast
  tbarx <- tbar^(0.5)
  tbar2 <- tbar^2
  tbar3 <- tbar^3
  rbar <- rho/rho_ast

  # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
  mu0 <- 1e2*tbarx/mu0

  # Create Table 3, Huber et al. (2009):
  h_array <- array(0.0, dim=c(7,6))
  h_array[1,] <- c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)  # hj0
  h_array[2,] <- c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573) # hj1
  h_array[3,] <- c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0) # hj2
  h_array[4,] <- c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0) # hj3
  h_array[5,] <- c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0) # hj4
  h_array[6,] <- c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0) # hj5
  h_array[7,] <- c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264) # hj6
  # print("h_array")
  # print(dim(h_array))
  # print(h_array)

  # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  mu1 <- 0.0
  ctbar <- (1.0/tbar) - 1.0
  # print(paste("ctbar",ctbar))
  # for i in xrange(6):
  for (i in 1:6){
    coef1 <- ctbar^(i-1)
    # print(paste("i, coef1", i, coef1))
    coef2 <- 0.0
    for (j in 1:7){
      coef2 <- coef2 + h_array[j,i] * (rbar - 1.0)^(j-1)
    }
    mu1 <- mu1 + coef1 * coef2    
  }
  mu1 <- exp( rbar * mu1 )
  # print(paste("mu1",mu1))

  # Calculate mu_bar (Eq. 2, Huber et al., 2009)
  #   assumes mu2 <- 1
  mu_bar <- mu0 * mu1

  # Calculate mu (Eq. 1, Huber et al., 2009)
  mu <- mu_bar * mu_ast    # Pa s

  return( mu )
}


viscosity_h2o_vogel <- function( tc ) {
  #-----------------------------------------------------------------------
  # Input:    - float, ambient temperature (tc), degrees C
  # Return:   float, viscosity of water (mu), Pa s
  # Features: Calculates viscosity of water at a given temperature and 
  #           pressure.
  # Depends:  density_h2o
  # Ref:      Vogel ...
  #-----------------------------------------------------------------------

  tk <- tc + 272.15

  a <- -3.7188
  b <- 578.919
  c <- 137.546

  visc <- 1e-3 * exp(a + b/(tk - c))

  # visc <- 10^(-3) * exp( -3.719+580 / (-138 + tk) )

  return( visc )
}

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
##  Test: Calculate GPP for monthly input data for CH-Oe1
##  year: 2000, elv: 450, lon: 7.73, lat: 47.3
## ////////////////////////////////////////////////////////////////////////
source("/alphadata01/bstocker/utilities/daily2monthly.R")
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
doy_vec <- 

nmonth <- 12

## GET MONTHLY CLIMATE
## temperature, convert from daily to monthly values
if (file.exists("/alphadata01/bstocker/sofun/trunk/components/mtemp_CH-Oe1_2002.txt")){

    ## read monthly file
    df.mtemp <- read.csv( "/alphadata01/bstocker/sofun/trunk/components/mtemp_CH-Oe1_2002.txt" )
    mtemp <- df.mtemp$temp

} else {

  ## read climate from year 2002 for station CH-Oe1
  dtemp <- read.table( "/alphadata01/bstocker/sofun/trunk/input/dtemp_CH-Oe1_2002.txt" )

  ## use function from our 'utilities' repository
  mtemp <- daily2monthly( dtemp$V1, "mean" )

  ## write monthly data to file
  df.mtemp <- data.frame( month=1:nmonth, temp=mtemp )
  write.csv( df.mtemp, file="/alphadata01/bstocker/sofun/trunk/components/mtemp_CH-Oe1_2002.txt", row.names=FALSE )

}

## precipitation, convert from daily to monthly values
if (file.exists("/alphadata01/bstocker/sofun/trunk/components/mprec_CH-Oe1_2002.txt")){

    ## read monthly file
    df.mprec <- read.csv( "/alphadata01/bstocker/sofun/trunk/components/mprec_CH-Oe1_2002.txt" )
    mprec <- df.mprec$prec

} else {

  ## read climate from year 2002 for station CH-Oe1
  dprec <- read.table( "/alphadata01/bstocker/sofun/trunk/input/dprec_CH-Oe1_2002.txt" )

  ## use function from our 'utilities' repository
  mprec <- daily2monthly( dprec$V1, "sum" )

  ## write monthly data to file
  df.mprec <- data.frame( month=1:nmonth, prec=mprec )
  write.csv( df.mprec, file="/alphadata01/bstocker/sofun/trunk/components/mprec_CH-Oe1_2002.txt", row.names=FALSE )

}  

## sunshine fraction
mfsun <- read.table("/alphadata01/bstocker/sofun/trunk/input/mfsun_CH-Oe1_2002.txt")$V1
filnam <- "/alphadata01/bstocker/sofun/trunk/components/mfsun2_CH-Oe1_2002.txt"
if (!file.exists(filnam)) {
  df.mfsun <- data.frame( month=1:nmonth, fsun=mfsun )
  write.csv( df.mfsun, file=filnam, row.names=FALSE )
}

## vapour pressure (mvapr, hPa), vapour pressure deficit (mvpd, kPa)
mvapr <- read.table("/alphadata01/bstocker/sofun/trunk/input/mvapr_CH-Oe1_2002.txt")$V1
filnam <- "/alphadata01/bstocker/sofun/trunk/components/mvapr2_CH-Oe1_2002.txt"
if (!file.exists(filnam)) {
  df.mvapr <- data.frame( month=1:nmonth, fsun=mvapr )
  write.csv( df.mvapr, file=filnam, row.names=FALSE )
}
mvpd <- rep(NA,nmonth)
for (moy in 1:nmonth){
  mvpd[moy] <- calc_vpd( mtemp[moy], mvapr[moy] )
}


## WANG HAN'S STATION DATA
## read input data for GPP calculation
indata <- read.csv( "/alphadata01/bstocker/sofun/trunk/input_raw/Wang_Han_beta_estimates.csv", header=TRUE )

## get day of year (doy) and year
indata$doy  <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$yday+1
indata$year <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$year+1900
indata$moy  <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$mon+1

## take sub-set of this for station CH-Oe1
indata <- indata[ indata$station == "CH-Oe1", ]
indata <- indata[ indata$year == 2002, ]

## PPFD FROM SOFUN (STASH implementation)
## read PPFD calculated by STASH (sofun implementation, '*.m.qm.out')
filnam <- "/alphadata01/bstocker/sofun/trunk/output/CH-Oe1_2002.m.qm.out"
df.ppfd <- read.table( filnam, col.names=c("year","ppfd") )
if (!file.exists(filnam)) {
  write.csv( df.ppfd, file=filnam, row.names=FALSE )
}

## take subset of year 2000
istart <- which.min( abs(df.ppfd$year-2000.0) )
df.ppfd <- df.ppfd[ istart:(istart+nmonth-1), ]


##------------------------------------------------------------
## CALCULATE STUFF for each month with:
##    - fapar = 1
##    - ppfd from STASH (sofun implementation), see input file
##    - alpha = 1.26
##    - co2 = 376 ppm
##    - elv = 450 m
##------------------------------------------------------------
elv    <- 450.0
patm   <- calc_patm(elv)
co2    <- 376.0
ca     <- co2_to_ca( co2, patm )         # Pa, atms. CO2
beta   <- 244
visc25 <- viscosity_h2o( kTo, kPo )      # Pa s


## CALCULATE GPP
mgpp <- rep( NA, nmonth )
for ( moy in 1:nmonth ){
  mgpp[moy] <- pmodel( fpar=1.0, ppfd=df.ppfd$ppfd[moy], co2=co2, tc=mtemp[moy], cpalpha=1.0, vpd=mvpd[moy], elv=elv, params=params, method="full" )$gpp
}

## CALCULATE K
## Michaelis-Menten coefficient (Pa)
kmm       <- sapply( mtemp, FUN = function(x) calc_k(x, patm) )
kmm_colin <- sapply( mtemp, FUN = calc_k_colin )

## CALCULATE VISCOSITY
visc25    <- viscosity_h2o( kTo, kPo )
visc      <- mapply( viscosity_h2o, mtemp, patm )
visc_star <- visc/visc25
visc_vogel<- sapply( mtemp, FUN = viscosity_h2o_vogel )
visc_out  <- visc * 1e3  # convert from Pa s to mPa s
visc_vogel_out <- visc_vogel * 1e3  # convert from Pa s to mPa s

## CALCULATE GAMMA-STAR
gstar       <- sapply( mtemp, FUN = calc_gstar_gepisat )
gstar_colin <- sapply( mtemp, FUN = calc_gstar_colin )

## CALCULATE CHI 
chi_wh     <- rep( NA, nmonth )
chi_simpl  <- rep( NA, nmonth )
chi_full   <- rep( NA, nmonth )
for (moy in 1:nmonth){
  chi_wh[moy]    <- lue_approx( mtemp[moy], mvpd[moy], elv, ca, gstar )$chi
  chi_simpl[moy] <- lue_vpd_simpl( kmm[moy], gstar[moy], visc_star[moy], ca, mvpd[moy], params )$chi
  chi_full[moy]  <- lue_vpd_full(  kmm[moy], gstar[moy], visc_star[moy], ca, mvpd[moy], params )$chi
}


## SAVE DATA FOR USE BY TYLER AND WANG-HAN
benioutput <- data.frame( year=rep(2002,nmonth), moy=1:nmonth, kpp_Pa=kmm
  , visc_Pa_s=visc, gstar_Pa=gstar, chi_WH_method=chi_wh, chi_simpl=chi_simpl, chi_full=chi_full 
  )
save( benioutput, file="benioutput.Rdata" )
write.csv( benioutput, file="benioutput.txt", row.names=FALSE )

## COMPARE DATA
## temperature: ok
# plot( 1:nmonth, mtemp, type="l", col="red" )
# for (idx in 1:dim(indata)[1]){
#   points( indata$moy[idx], indata$Tair_degC[idx] )
# }

## VPD xxx not ok xxx: is WH's data really in kPa (not Pa)?
pdf("vpd_comparison.pdf")
plot( 1:nmonth, mvpd, type="l", col="red", xlab="MOY"  )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$D_kPa[idx] )
}
dev.off()

## K - Michaelis-Menten coefficient: ok
pdf("Kc_comparison.pdf")
plot( 1:nmonth, kmm, type="l", col="red", xlab="MOY"  )
lines(1:nmonth, kmm_colin, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$K_Pa[idx] )
}
legend( "topleft", c("GePiSaT method", "Colin method"), lty=1, bty="n", col=c("red","green") )
dev.off()

## Viscosity: systematically lower values by my code, but tested to be identical with python implementation in gepisat
pdf("viscosity_comparison.pdf")
plot( 1:nmonth, visc_out, type="l", col="red", xlab="MOY", ylab="viscosity [mPa s]" )
lines(1:nmonth, visc_vogel_out, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$ns[idx] )
}
legend( "bottomleft", c("Huber method", "Vogel method"), lty=1, bty="n", col=c("red","green") )
dev.off()

## Gamma-star: ok
pdf("gammastar_comparison.pdf")
plot( 1:nmonth, gstar, type="l", col="red", xlab="MOY"  )
lines(1:nmonth, gstar_colin, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$Gs_Pa[idx] )
}
legend( "topleft", c("GePiSaT method", "Colin method"), lty=1, bty="n", col=c("red","green") )
dev.off()

## CHI - WH method
pdf("chi_comparison.pdf")
plot( 1:nmonth, chi_wh, type="l", col="red", xlab="MOY"  )
lines( 1:nmonth, chi_full, col="blue" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$chi_wang_han[idx] )
  points( indata$moy[idx], indata$chi_vpd[idx], col="green" )
}
legend( "topleft", c("Wang-Han method", "theoretical full method"), lty=1, bty="n", col=c("red","blue") )
dev.off()

## GPP 
plot( 1:nmonth, mgpp, type="l", col="red", xlab="MOY"  )
