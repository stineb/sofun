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
    getout_daily_gpp, writeout_ascii_gpp
  use md_npp, only: npp
  use md_turnover, only: turnover
  use md_vegdynamics, only: vegdynamics
  use md_littersom, only: getpar_modl_littersom, initio_littersom, initoutput_littersom, &
    getout_annual_littersom, writeout_ascii_littersom, littersom, initdaily_littersom, initglobal_littersom
  use md_ntransform, only: pninorg, ntransform, getpar_modl_ntransform, init_global_ntransform, &
    initdaily_ntransform, initio_ntransform, initoutput_ntransform, getout_daily_ntransform, &
    writeout_ascii_ntransform
  use md_nuptake, only: nuptake, getpar_modl_nuptake, initdaily_nuptake, initio_nuptake, &
    initoutput_nuptake, getout_daily_nuptake, writeout_ascii_nuptake
  use md_allocation, only: allocation_daily, initio_allocation, initoutput_allocation, &
    getout_daily_allocation, writeout_ascii_allocation

  
  implicit none

  ! return variable
  real, intent(out) :: c_uptake   ! annual net global C uptake by biosphere (gC/yr)

  ! local variables
  integer :: dm, moy, jpngr, day

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
    call getpar_modl_littersom()
    call getpar_modl_ntransform()
    call getpar_modl_nuptake()
    ! print*,'... done'

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    ! print*,'initialising variables ...'
    call initglobal_plant()
    call initglobal_waterbal()
    call initglobal_littersom()
    call init_global_ntransform()

    ! print*,'... done'

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    ! print*,'initialising IO ...'
    call initio_waterbal()
    call initio_soiltemp()
    call initio_gpp()
    call initio_plant()
    call initio_littersom()
    call initio_ntransform()
    call initio_nuptake()
    call initio_allocation()
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
  call initoutput_littersom()
  call initoutput_ntransform()
  call initoutput_nuptake()
  call initoutput_allocation()
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
        ! print*, 'calling gpp() ... '
        ! write(0,*) 'WARNING IN BIOSPHERE: CAPPED DAILY TEMPERATURE AT 25 DEG C.'
        ! if (interface%climate(jpngr)%dtemp(day).gt.25.0) print*,'dtemp = ', interface%climate(jpngr)%dtemp(day)
        ! call gpp( jpngr, day, moy, min( 25.0, interface%climate(jpngr)%dtemp(day) ), mfapar_field(moy,jpngr) )
        call gpp( &
          jpngr, day, moy, & 
          interface%climate(jpngr)%dtemp(day), & 
          interface%mfapar_field(moy,jpngr) & 
          )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        ! print*, 'calling npp() ... '
        call npp( jpngr, interface%climate(jpngr)%dtemp(day), day )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! calculate N acquisition as a function of C exudation
        !----------------------------------------------------------------
        ! print*, 'calling npp() ... '
        call nuptake( jpngr )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        ! print*, 'calling turnover() ... '
        call turnover( jpngr, day )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! litter and soil decomposition and N mineralisation
        !----------------------------------------------------------------
        ! print*, 'calling littersom() ... '
        call littersom( jpngr, day )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! inorganic soil N dynamics
        !----------------------------------------------------------------
        ! print*, 'calling ntransform() ... '
        call ntransform( dm, moy, jpngr, interface%ndep_field(jpngr)%dtot(day), sum(interface%climate(jpngr)%dprec(:)) )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        ! print*, 'calling allocation() ... '
        call allocation_daily( jpngr, day, moy, dm, interface%climate(jpngr)%dtemp(day) )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        ! print*, 'calling getout_daily_*() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        call getout_daily_gpp( jpngr, moy, day )
        call getout_daily_plant( jpngr, moy, day )
        call getout_daily_ntransform( jpngr, moy, day )
        call getout_daily_nuptake( jpngr, moy, day )
        call getout_daily_allocation( jpngr, moy, day )
        ! print*, '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! collect annually updated output variables
    !----------------------------------------------------------------
    call getout_annual_plant( jpngr )
    call getout_annual_littersom( jpngr )

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    ! print*, 'calling writeout_ascii_() ... '
    call writeout_ascii_waterbal( interface%steering%year )
    call writeout_ascii_soiltemp( interface%steering%year )
    call writeout_ascii_gpp( interface%steering%year )
    call writeout_ascii_plant( interface%steering%year )
    call writeout_ascii_ntransform( interface%steering%year )
    call writeout_ascii_nuptake( interface%steering%year )
    call writeout_ascii_allocation( interface%steering%year )
    ! print*, '... done'

  end do

  ! xxx insignificant
  c_uptake = 0.0

end subroutine biosphere

