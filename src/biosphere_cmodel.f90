subroutine biosphere( c_uptake )
  !////////////////////////////////////////////////////////////////
  ! subroutine BIOSPHERE calculates net ecosystem exchange (nee)
  ! in response to environmental boundary conditions (atmospheric 
  ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
  ! LPJ, also formulated as subroutine.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_interface
  use md_params_core
  use md_plant, only: initdaily_plant, initglobal_plant, initpft, getout_daily_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant, getout_annual_plant
  use md_soiltemp, only: soiltemp, initoutput_soiltemp, initio_soiltemp, getout_daily_soiltemp, writeout_ascii_soiltemp
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: waterbal, getsolar_alldays, initdaily_waterbal, initglobal_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_phenology, only: gettempphenology, getpar_modl_phenology
  use md_gpp, only: getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, getout_daily_gpp, writeout_ascii_gpp
  use md_npp, only: npp
  use md_turnover, only: turnover
  use md_allocation, only: allocation_daily
  use md_vegdynamics, only: vegdynamics
  use md_landuse, only: grharvest, initoutput_landuse, initio_landuse, getout_annual_landuse, &
    writeout_ascii_landuse, initglobal_landuse, init_mharv
  
  ! xxx verbose
  use md_plant
  use md_gpp
  use md_classdefs
  use md_landuse

  implicit none

  ! return variable
  real, intent(out) :: c_uptake   ! annual net global C uptake by biosphere (gC/yr)

  ! local variables
  integer :: dm, moy, jpngr, day

  ! xxx verbose
  logical, parameter :: verbose = .true.
  real            :: cbal1, cbal2
  type( orgpool ) :: orgtmp1, orgtmp2, orgbal1, orgbal2
  real :: eps = 9.999e-11

  ! ! XXX PMODEL_TEST
  ! print*, 'WARNING: FAPAR = 1.00 USED IN PMODEL'

  !----------------------------------------------------------------
  ! INITIALISATIONS
  !----------------------------------------------------------------
  if (interface%steering%init) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    ! print*,'getting model parameters ...'
    call getpar_modl_plant()
    call getpar_modl_waterbal()
    call getpar_modl_gpp()
    call getpar_modl_phenology()
    ! print*,'... done'

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    ! print*,'initialising variables ...'
    call initglobal_plant()
    call initglobal_waterbal()
    call initglobal_landuse()
    ! print*,'... done'

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    ! print*,'initialising IO ...'
    call initio_waterbal()
    call initio_soiltemp()
    call initio_gpp()
    call initio_plant()
    call initio_landuse()
    ! print*,'... done'

  endif 

  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  ! print*,'initialising output variables ...'
  call initoutput_waterbal()
  call initoutput_soiltemp()
  call initoutput_gpp()
  call initoutput_plant()
  call initoutput_landuse()
  ! print*,'... done'

  !----------------------------------------------------------------
  ! LOOP THROUGH GRIDCELLS
  !----------------------------------------------------------------
  do jpngr=1,maxgrid

    !----------------------------------------------------------------
    ! Get radiation based on daily temperature, sunshine fraction, and 
    ! elevation.
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a daily loop within 'getsolar'!
    !----------------------------------------------------------------
    call getsolar_alldays( &
      interface%grid(jpngr)%lat, & 
      interface%grid(jpngr)%elv, & 
      interface%climate(jpngr)%dfsun(:) & 
      )

    !----------------------------------------------------------------
    ! Get monthly light use efficiency, and Rd per unit of light absorbed
    ! Photosynthetic parameters acclimate at monthly time scale
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a monthly loop within 'getlue'!
    !----------------------------------------------------------------
    call getlue( &
      jpngr, & 
      interface%pco2, & 
      interface%climate(jpngr)%dtemp(:), & 
      interface%climate(jpngr)%dvpd(:), & 
      interface%grid(jpngr)%elv & 
      )

    !----------------------------------------------------------------
    ! Get radiation based on daily temperature, sunshine fraction, and 
    ! elevation.
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a daily loop within 'getsolar'!
    !----------------------------------------------------------------
    call gettempphenology( jpngr, interface%climate(jpngr)%dtemp(:) )

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

        !----------------------------------------------------------------
        ! initialise daily updated variables 
        !----------------------------------------------------------------
        if (day==1) call init_mharv(jpngr) ! xxx try
        call initdaily_plant()
        call initdaily_waterbal()
        call initdaily_gpp()
        call initdaily_plant()

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        ! print*, 'calling waterbal() ... '
        call waterbal( &
          jpngr, day, & 
          interface%grid(jpngr)%lat, & 
          interface%grid(jpngr)%elv, & 
          interface%climate(jpngr)%dprec(day), & 
          interface%climate(jpngr)%dtemp(day), & 
          interface%climate(jpngr)%dfsun(day)  &
          )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        ! print*, 'calling soiltemp() ... '
        call soiltemp( &
                      jpngr, & 
                      moy, & 
                      day, & 
                      interface%climate(jpngr)%dtemp(:) &
                      )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! update canopy and stand variables and simulate daily 
        ! establishment / sprouting
        !----------------------------------------------------------------
        ! print*, 'calling vegdynamics() ... '
        call vegdynamics( jpngr, day ) 
        ! print*, '... done'

        !----------------------------------------------------------------
        ! calculate GPP
        !----------------------------------------------------------------
        if (verbose) print*, 'calling gpp() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              fapar = ', canopy(:)%fapar_ind
        ! write(0,*) 'WARNING IN BIOSPHERE: CAPPED DAILY TEMPERATURE AT 25 DEG C.'
        ! if (interface%climate(jpngr)%dtemp(day).gt.25.0) print*,'dtemp = ', interface%climate(jpngr)%dtemp(day)
        ! call gpp( jpngr, day, moy, min( 25.0, interface%climate(jpngr)%dtemp(day) ), mfapar_field(moy,jpngr) )
        call gpp( &
          jpngr, day, moy, & 
          interface%climate(jpngr)%dtemp(day), & 
          interface%mfapar_field(moy,jpngr) & 
          )
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              dgpp  = ', dgpp(:)
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        if (verbose) print*, 'calling npp() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) orgtmp1 =  plabl(1,jpngr)
        !----------------------------------------------------------------
        call npp( jpngr, interface%climate(jpngr)%dtemp(day), day )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              dnpp  = ', dnpp(:)
        if (verbose) print*, '              dcex  = ', dcex(:)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '    --- balance: '
        if (verbose) cbal1 = dgpp(1) - dnpp(1)%c12 - drleaf(1) - drroot(1)
        if (verbose) cbal2 = dgpp(1) - ( plabl(1,jpngr)%c%c12 - orgtmp1%c%c12 ) - dcex(1) - drleaf(1) - drroot(1)
        if (verbose) print*, '        gpp - npp - ra_maint          = ', cbal1
        if (verbose) print*, '        gpp - dlabl - dcex - ra_maint = ', cbal2
        if (verbose.and.abs(cbal1)>eps) stop 'balance 1 not satisfied'
        if (verbose.and.abs(cbal2)>eps) stop 'balance 2 not satisfied'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        if (verbose) print*, 'calling turnover() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '              plitt = ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) orgtmp1  =  orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) )
        if (verbose) orgtmp2 =  orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        !----------------------------------------------------------------
        call turnover( jpngr, day )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '              plitt = ', orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) )
        if (verbose) print*, '   --- balance: '
        if (verbose) orgbal1 = orgminus( &
                                        orgminus( &
                                          orgplus( plitt_af(1,jpngr), plitt_as(1,jpngr), plitt_bg(1,jpngr) ), &
                                          orgtmp2 &
                                          ), &
                                        orgminus( &
                                          orgtmp1, &
                                          orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) ) &
                                          ) &
                                        )
        if (verbose) print*, '       dlitt - dplant                = ', orgbal1
        if (verbose.and.abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        if (verbose.and.abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! grass / crop harvest
        !----------------------------------------------------------------
        if (verbose) print*, 'calling grharvest() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '              mharv = ', mharv(:,jpngr)
        if (verbose) orgtmp1 =  orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) )
        if (verbose) orgtmp2 =  mharv(1,jpngr)
        !----------------------------------------------------------------
        call grharvest( jpngr, day )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '              mharv = ', mharv(:,jpngr)
        if (verbose) print*, '    --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus( orgtmp1, orgplus( pleaf(1,jpngr), proot(1,jpngr), plabl(1,jpngr) ) ), orgminus( mharv(1,jpngr), orgtmp2 ) )
        if (verbose) print*, '        dharv - dplant  = ', orgbal1
        if (verbose.and.abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        if (verbose.and.abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        if (verbose) print*, 'calling allocation() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) orgtmp1 =  plabl(1,jpngr)
        if (verbose) orgtmp2 =  orgplus( pleaf(1,jpngr), proot(1,jpngr) )
        !----------------------------------------------------------------
        call allocation_daily( jpngr, day, moy, dm )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr), ' C:N = ', cton( pleaf(1,jpngr) )
        if (verbose) print*, '              proot = ', proot(:,jpngr), ' C:N = ', cton( proot(1,jpngr) )
        if (verbose) print*, '              plabl = ', plabl(:,jpngr)
        if (verbose) print*, '   --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus( orgminus( orgtmp1, plabl(1,jpngr) ), orgpool( carbon(drgrow(1)), nitrogen(0.0) ) ), orgminus( orgplus( pleaf(1,jpngr), proot(1,jpngr) ), orgtmp2 ) )        
        ! dlabl  = orgminus( orgtmp1, plabl(1,jpngr) )
        ! drgrow = orgpool( carbon(drgrow(1)), nitrogen(0.0) )
        ! dplant = orgminus( orgplus( pleaf(1,jpngr), proot(1,jpngr) ), orgtmp2 )
        if (verbose) print*, '       dlabl - drgrow - dleaf - droot=', orgbal1
        if (verbose.and.abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
        ! print*, 'dlabl  =', orgminus( orgtmp, plabl(1,jpngr) )
        ! print*, 'drgrow =', drgrow(1), 0.0
        ! print*, 'dplant =', orgminus( orgplus( pleaf(1,jpngr), proot(1,jpngr) ), orgtmp2 )
        if (verbose) print*, '... done'

        ! !----------------------------------------------------------------
        ! ! litter and soil decomposition and N mineralisation
        ! !----------------------------------------------------------------
        ! print*, 'calling littersom() ... '
        ! call littersom( jpngr, day )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        ! print*, 'calling getout_daily_() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        call getout_daily_gpp( jpngr, moy, day )
        call getout_daily_plant( jpngr, moy, day )
        ! print*, '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! collect annual output
    !----------------------------------------------------------------
    ! print*, 'calling getout_annual() ... '
    call getout_annual_plant( jpngr )
    call getout_annual_landuse( jpngr )
    ! print*, '... done'

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    ! print*, 'calling writeout_ascii_() ... '
    call writeout_ascii_waterbal( interface%steering%year )
    call writeout_ascii_soiltemp( interface%steering%year )
    call writeout_ascii_gpp( interface%steering%year )
    call writeout_ascii_plant( interface%steering%year )
    call writeout_ascii_landuse( interface%steering%year )
    ! print*, '... done'

  end do

  ! if (interface%steering%year==5) stop 'end of year no. 5'

  ! xxx insignificant
  c_uptake = 0.0

end subroutine biosphere

