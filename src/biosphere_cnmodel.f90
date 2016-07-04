subroutine biosphere( c_uptake )
  !////////////////////////////////////////////////////////////////
  ! This subroutine acts as a "wrapper" for modules simulating 
  ! different processes in ecosystem C and N cycling. This wrapper
  ! is specific for each level of model integration (see compilation 
  ! profiles in Makefile). E.g., allocation module is not used in
  ! 'pmodel' setup (using prescribed fAPAR and simulating only 
  ! photosynthesis). 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_interface
  use md_params_core
  use md_params_siml
  use md_params_site
  use md_plant, only: initdaily_plant, initglobal_plant, initpft, getout_daily_plant, &
    getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant, getout_annual_plant
  use md_soiltemp, only: soiltemp, initoutput_soiltemp, initio_soiltemp, getout_daily_soiltemp, &
    writeout_ascii_soiltemp
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: waterbal, getsolar_alldays, initdaily_waterbal, initglobal_waterbal, &
    initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, &
    writeout_ascii_waterbal
  use md_phenology, only: gettempphenology, getpar_modl_phenology
  use md_gpp, only: getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, &
    getout_daily_gpp, getout_annual_gpp, writeout_ascii_gpp
  use md_npp, only: npp
  use md_turnover, only: turnover
  use md_vegdynamics, only: vegdynamics
  use md_littersom, only: getpar_modl_littersom, initio_littersom, initoutput_littersom, &
    getout_annual_littersom, writeout_ascii_littersom, littersom, initdaily_littersom, initglobal_littersom
  use md_ntransform, only: pninorg, ntransform, getpar_modl_ntransform, initglobal_ntransform, &
    initdaily_ntransform, initio_ntransform, initoutput_ntransform, getout_daily_ntransform, &
    writeout_ascii_ntransform
  use md_nuptake, only: nuptake, getpar_modl_nuptake, initdaily_nuptake, initio_nuptake, &
    initoutput_nuptake, getout_daily_nuptake, writeout_ascii_nuptake
  use md_allocation, only: allocation_daily, initio_allocation, initoutput_allocation, &
    getout_daily_allocation, writeout_ascii_allocation
  use md_landuse, only: grharvest, initoutput_landuse, initio_landuse, getout_annual_landuse, &
    writeout_ascii_landuse, initglobal_landuse, init_mharv

  ! Verbose: uncomment these lines if 'verbose' is set to true
  use md_plant
  use md_gpp
  use md_classdefs
  use md_landuse
  use md_littersom
  use md_ntransform
  use md_waterbal

  implicit none

  ! return variable
  real, intent(out) :: c_uptake   ! annual net global C uptake by biosphere (gC/yr)

  ! local variables
  integer :: dm, moy, jpngr, day, usemoy, usedoy

  ! Variables used for verbose mode and mass conservation test 'baltest'
  logical, parameter :: verbose = .false.  ! set to true to activate verbose mode
  logical, parameter :: baltest = .false.  ! set to true to do mass conservation test
  real            :: cbal1, cbal2
  type( orgpool ) :: orgtmp1, orgtmp2, orgbal1, orgbal2
  real :: eps = 9.999e-8                   ! numerical imprecision allowed in mass conservation tests

  !----------------------------------------------------------------
  ! INITIALISATIONS
  ! at the start of the simulation
  ! ----------------------------------------------------------------
  if (interface%steering%init) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'getting model parameters ...'
    call getpar_modl_plant()
    call getpar_modl_waterbal()
    call getpar_modl_gpp()
    call getpar_modl_phenology()
    call getpar_modl_littersom()
    call getpar_modl_ntransform()
    call getpar_modl_nuptake()
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'initialising variables ...'
    call initglobal_plant()
    call initglobal_waterbal()
    call initglobal_littersom()
    call initglobal_ntransform()
    call initglobal_landuse()
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'initialising IO ...'
    call initio_waterbal()
    call initio_soiltemp()
    call initio_gpp()
    call initio_plant()
    call initio_littersom()
    call initio_ntransform()
    call initio_nuptake()
    call initio_allocation()
    call initio_landuse()
    if (verbose) write(0,*) '... done'

  endif 

  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  if (verbose) write(0,*) 'initialising output variables ...'
  call initoutput_waterbal()
  call initoutput_soiltemp()
  call initoutput_gpp()
  call initoutput_plant()
  call initoutput_littersom()
  call initoutput_ntransform()
  call initoutput_nuptake()
  call initoutput_allocation()
  call initoutput_landuse()
  if (verbose) write(0,*) '... done'

  !----------------------------------------------------------------
  ! LOOP THROUGH GRIDCELLS
  ! So far, 'maxgrid' (number of gridcells) is 1. Multiple-cell (spatial)
  ! simulations are not possible but the structure of the program is
  ! designed for a extending it straight-forward.
  !----------------------------------------------------------------
  do jpngr=1,maxgrid

    !----------------------------------------------------------------
    ! Get radiation based on daily temperature, sunshine fraction, and 
    ! elevation.
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a daily loop within 'getsolar'!
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'calling getsolar_alldays() ...'
    call getsolar_alldays( &
      interface%grid(jpngr)%lat, & 
      interface%grid(jpngr)%elv, & 
      interface%climate(jpngr)%dfsun(:) & 
      )
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! Get monthly light use efficiency, and Rd per unit of light absorbed.
    ! Photosynthetic parameters acclimate at monthly time scale.
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'calling getlue() ...'
    call getlue( &
      jpngr, & 
      interface%pco2, & 
      interface%climate(jpngr)%dtemp(:), & 
      interface%climate(jpngr)%dvpd(:), & 
      interface%grid(jpngr)%elv & 
      )
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! Get phenology variables (temperature-drivenå)
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'calling gettempphenology() ...'
    call gettempphenology( jpngr, interface%climate(jpngr)%dtemp(:) )
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    day=0
    do moy=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
      do dm=1,ndaymonth(moy)
        day=day+1

        if (verbose) write(0,*) '----------------------'
        if (verbose) write(0,*) 'YEAR, DAY ', interface%steering%year, day
        if (verbose) write(0,*) '----------------------'

        !----------------------------------------------------------------
        ! initialise daily updated variables 
        !----------------------------------------------------------------
        if (day==1) call init_mharv(jpngr) ! xxx try
        call initdaily_plant()
        call initdaily_waterbal()
        call initdaily_gpp()
        call initdaily_plant()
        call initdaily_littersom()
        call initdaily_ntransform()
        call initdaily_nuptake()

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling waterbal() ... '
        call waterbal( &
          jpngr, day, & 
          interface%grid(jpngr)%lat, & 
          interface%grid(jpngr)%elv, & 
          interface%climate(jpngr)%dprec(day), & 
          interface%climate(jpngr)%dtemp(day), & 
          interface%climate(jpngr)%dfsun(day)  &
          )
        if (verbose) write(0,*) '... done'

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling soiltemp() ... '
        call soiltemp( &
                      jpngr, & 
                      moy, & 
                      day, & 
                      interface%climate(jpngr)%dtemp(:) &
                      )
        if (verbose) write(0,*) '... done'

        !----------------------------------------------------------------
        ! update canopy and stand variables and simulate daily 
        ! establishment / sprouting
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling vegdynamics() ... '
        call vegdynamics( jpngr, day ) 
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! calculate GPP
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling gpp() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              fapar = ', canopy(:)%fapar_ind
        if (verbose) write(0,*) '              dppfd = ', solar%dppfd(day)
        if (verbose) write(0,*) '              lue   = ', mlue(1,moy)
        if (verbose) write(0,*) '              temp  = ', interface%climate(jpngr)%dtemp(day)
        if (verbose) write(0,*) '              CPA   = ', evap(1)%cpa
        !----------------------------------------------------------------
        call gpp( &
          jpngr, day, moy, & 
          interface%climate(jpngr)%dtemp(day), & 
          interface%mfapar_field(moy,jpngr) & 
          )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              dgpp  = ', dgpp(:)
        if (verbose) write(0,*) '              drd   = ', drd(:)
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling npp() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) orgtmp1 =  plabl(1,jpngr)
        !----------------------------------------------------------------
        call npp( jpngr, interface%climate(jpngr)%dtemp(day) )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              dnpp  = ', dnpp(:)
        if (verbose) write(0,*) '              dcex  = ', dcex(:)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              dlabl = ', orgminus( plabl(1,jpngr), orgtmp1 )
        if (verbose) write(0,*) '    --- balance: '
        if (verbose) cbal1 = dgpp(1) - dnpp(1)%c12 - drleaf(1) - drroot(1)
        if (verbose) cbal2 = dgpp(1) - ( plabl(1,jpngr)%c%c12 - orgtmp1%c%c12 ) - dcex(1) - drleaf(1) - drroot(1)
        if (verbose) write(0,*) '        gpp - npp - ra_maint          = ', cbal1
        if (verbose) write(0,*) '        gpp - dlabl - dcex - ra_maint = ', cbal2
        if (verbose .and. baltest .and. abs(cbal1)>eps) stop 'balance 1 not satisfied'
        if (verbose .and. baltest .and. abs(cbal2)>eps) stop 'balance 2 not satisfied'
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! calculate N acquisition as a function of C exudation
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling nuptake() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              dcex  = ', dcex(:)
        !----------------------------------------------------------------
        call nuptake( jpngr )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              dnup  = ', dnup(:)
        if (verbose) write(0,*) '... done'

        ! if (break_after_alloc) stop 'check quantities after allocation'

        !/////////////////////////////////////////////////////////////////
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling turnover() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              plitt = ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) orgtmp1 = orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) )
        if (verbose) orgtmp2 = orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        !----------------------------------------------------------------
        call turnover( jpngr, day )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              plitt = ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) write(0,*) '   --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus(   orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) ),   orgtmp2   ), orgminus(   orgtmp1,   orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) )   ) )
        if (verbose) write(0,*) '       dlitt - dplant                = ', orgbal1
        if (verbose .and. baltest .and. abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        if (verbose .and. baltest .and. abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! grass / crop harvest
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling grharvest() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              mharv = ', mharv(:,jpngr)
        if (verbose) orgtmp1 =  orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) )
        if (verbose) orgtmp2 =  mharv(1,jpngr)
        !----------------------------------------------------------------
        call grharvest( jpngr, day )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              mharv = ', mharv(:,jpngr)
        if (verbose) write(0,*) '    --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus( orgtmp1, orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) ) ), orgminus( mharv(1,jpngr), orgtmp2 ) )
        if (verbose) write(0,*) '        dharv - dplant  = ', orgbal1
        if (verbose .and. baltest .and. abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        if (verbose .and. baltest .and. abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! litter and soil decomposition and N mineralisation
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling littersom() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              plitt_af =  ', plitt_af(1,jpngr)
        if (verbose) write(0,*) '              plitt_as =  ', plitt_as(1,jpngr)
        if (verbose) write(0,*) '              plitt_bg =  ', plitt_bg(1,jpngr)
        if (verbose) write(0,*) '              plitt tot=  ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) write(0,*) '              psoil_fs  = ', psoil_fs(1,jpngr)
        if (verbose) write(0,*) '              psoil_sl  = ', psoil_sl(1,jpngr)
        if (verbose) write(0,*) '              psoil tot = ', orgplus( psoil_fs(1,jpngr), psoil_sl(1,jpngr) )
        if (verbose) write(0,*) '              pexud     = ', pexud(1,jpngr)
        if (verbose) write(0,*) '              pninorg=    ', pninorg(1,jpngr)
        if (verbose) write(0,*) '              drhet  =    ', drhet(1)
        if (verbose) orgtmp1 = orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr), psoil_fs(1,jpngr), psoil_sl(1,jpngr) )
        if (verbose) orgtmp2 = orgpool( drhet(1), pninorg(1,jpngr) )
        !----------------------------------------------------------------
        call littersom( jpngr, day, interface%climate(jpngr)%dtemp(day) )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              plitt  = ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) write(0,*) '              psoil  = ', orgplus( psoil_fs(1,jpngr), psoil_sl(1,jpngr) )
        if (verbose) write(0,*) '              pninorg= ', pninorg(1,jpngr)
        if (verbose) write(0,*) '              drhet  = ', drhet(1)
        if (verbose) write(0,*) '   --- balance: '
        if (verbose) write(0,*) '       d( litt + soil ) - d(drhet,ninorg) = ', orgminus( &
                                                                                          orgminus( &
                                                                                            orgtmp1, &
                                                                                            orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr), psoil_fs(1,jpngr), psoil_sl(1,jpngr) ) &
                                                                                            ), &
                                                                                          orgminus( &
                                                                                            orgpool( drhet(1), pninorg(1,jpngr) ), &
                                                                                            orgtmp2 ) &
                                                                                          )
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! inorganic soil N dynamics
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling ntransform() ... '
        !----------------------------------------------------------------
        call ntransform( dm, moy, jpngr, interface%ninput_field(jpngr)%dtot(day), sum(interface%climate(jpngr)%dprec(:)) )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '... done'

        !/////////////////////////////////////////////////////////////////
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling allocation() ... '
        if (verbose) write(0,*) '              with state variables:'
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr)
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr)
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              drgrow= ', drgrow(:)
        if (verbose) orgtmp1 =  plabl(1,jpngr)
        if (verbose) orgtmp2 =  orgplus( pleaf(1,jpngr), proot(1,jpngr) )
        !----------------------------------------------------------------
        if (dm==ndaymonth(moy)) then
          usemoy = moy + 1
          if (usemoy==13) usemoy = 1
        else
          usemoy = moy
        end if
        if (day==ndayyear) then
          usedoy = 1
        else
          usedoy = day + 1
        end if
        !----------------------------------------------------------------
        call allocation_daily( jpngr, usedoy, usemoy, interface%climate(jpngr)%dtemp(usedoy) )
        !----------------------------------------------------------------
        if (verbose) write(0,*) '              ==> returned: '
        if (verbose) write(0,*) '              pleaf = ', pleaf(:,jpngr), ' C:N = ', cton( pleaf(1,jpngr), default=0.0 )
        if (verbose) write(0,*) '              proot = ', proot(:,jpngr), ' C:N = ', cton( proot(1,jpngr), default=0.0 )
        if (verbose) write(0,*) '              plabl = ', plabl(:,jpngr)
        if (verbose) write(0,*) '              drgrow= ', drgrow(:)
        if (verbose) write(0,*) '   --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus( orgminus( orgtmp1, plabl(1,jpngr) ), orgpool( carbon(drgrow(1)), nitrogen(0.0) ) ), orgminus( orgplus( pleaf(1,jpngr), proot(1,jpngr) ), orgtmp2 ) )        
        if (verbose) write(0,*) '       dlabl - drgrow - dleaf - droot=', orgbal1
        if (verbose .and. baltest .and. abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        if (verbose) write(0,*) '... done'


        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling getout_daily_*() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        call getout_daily_gpp( jpngr, moy, day )
        call getout_daily_plant( jpngr, moy, day )
        call getout_daily_ntransform( jpngr, moy, day )
        call getout_daily_nuptake( jpngr, moy, day )
        call getout_daily_allocation( jpngr, moy, day )
        if (verbose) write(0,*) '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! collect annually updated output variables
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'calling getout_annual() ... '
    call getout_annual_gpp( jpngr )
    call getout_annual_plant( jpngr )
    call getout_annual_littersom( jpngr )
    call getout_annual_landuse( jpngr )
    if (verbose) write(0,*) '... done'

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    if (verbose) write(0,*) 'calling writeout_ascii_() ... '
    call writeout_ascii_waterbal( interface%steering%year )
    call writeout_ascii_soiltemp( interface%steering%year )
    call writeout_ascii_gpp( interface%steering%year )
    call writeout_ascii_plant( interface%steering%year )
    call writeout_ascii_ntransform( interface%steering%year )
    call writeout_ascii_nuptake( interface%steering%year )
    call writeout_ascii_allocation( interface%steering%year )
    call writeout_ascii_littersom( interface%steering%year )
    call writeout_ascii_landuse( interface%steering%year )
    if (verbose) write(0,*) '... done'

  end do

  ! xxx insignificant
  c_uptake = 0.0

end subroutine biosphere

