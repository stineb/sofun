module md_biosphere

  use md_params_core
  use md_classdefs
  use md_plant, only: plant_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant, initio_nc_plant, writeout_nc_plant
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, getsolar, initdaily_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_gpp, only: outtype_pmodel, getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, getout_daily_gpp, getout_annual_gpp, writeout_ascii_gpp
  use md_vegdynamics, only: vegdynamics
  use md_tile, only: tile_type, initglobal_tile
  use md_interface, only: getout_daily_forcing, initoutput_forcing, initio_forcing, initio_nc_forcing, writeout_ascii_forcing, writeout_nc_forcing

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  type( tile_type ),  allocatable, dimension(:,:) :: tile
  type( plant_type ), allocatable, dimension(:,:) :: plant

  type( solartype )                              :: solar
  type( outtype_pmodel ), dimension(npft,nmonth) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

contains

  function biosphere_annual() result( c_uptake )
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface, only: interface
  
    ! return variable
    real :: c_uptake   ! annual net global C uptake by biosphere (gC/yr)

    ! local variables
    integer :: dm, moy, jpngr, doy

    ! xxx verbose
    logical, parameter :: verbose = .false.

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (interface%steering%init) then

      !----------------------------------------------------------------
      ! GET MODEL PARAMETERS
      ! read model parameters that may be varied for optimisation
      !----------------------------------------------------------------
      if (verbose) print*, 'getpar_modl() ...'
      call getpar_modl_plant()
      call getpar_modl_waterbal()
      call getpar_modl_gpp()
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      if (verbose) print*, 'nitglobal_() ...'
      allocate( tile(  nlu,  size(interface%grid) ) )
      allocate( plant( npft, size(interface%grid) ) )

      call initglobal_tile(  tile(:,:),  size(interface%grid) )
      call initglobal_plant( plant(:,:), size(interface%grid) )
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Open ascii output files
      !----------------------------------------------------------------
      if (verbose) print*, 'initio_() ...'
      call initio_waterbal()
      call initio_gpp()
      call initio_plant()
      call initio_forcing()
      if (verbose) print*, '... done'

    endif 

    !----------------------------------------------------------------
    ! Open NetCDF output files (one for each year)
    !----------------------------------------------------------------
    if (verbose) print*, 'initio_nc_() ...'
    call initio_nc_forcing()
    call initio_nc_plant()
    if (verbose) print*, '... done'

    !----------------------------------------------------------------
    ! Initialise output variables for this year
    !----------------------------------------------------------------
    if (verbose) print*, 'initoutput_() ...'
    call initoutput_waterbal( size(interface%grid) )
    call initoutput_gpp(      size(interface%grid) )
    call initoutput_plant(    size(interface%grid) )
    call initoutput_forcing(  size(interface%grid) )
    if (verbose) print*, '... done'

    !----------------------------------------------------------------
    ! LOOP THROUGH GRIDCELLS
    !----------------------------------------------------------------
    gridcellloop: do jpngr=1,size(interface%grid)
      if (interface%grid(jpngr)%dogridcell) then

        write(0,*) '----------------------'
        write(0,*) 'YEAR, jpngr ', interface%steering%year, jpngr
        write(0,*) '----------------------'

        !----------------------------------------------------------------
        ! Get radiation based on daily temperature, sunshine fraction, and 
        ! elevation.
        ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! there is a daily loop within 'getsolar'!
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling getsolar() ... '
        if (verbose) write(0,*) '    with argument lat = ', interface%grid(jpngr)%lat
        if (verbose) write(0,*) '    with argument elv = ', interface%grid(jpngr)%elv
        if (verbose) write(0,*) '    with argument dfsun (ann. mean) = ', sum( interface%climate(jpngr)%dfsun(:) / ndayyear )
        solar = getsolar( &
                          interface%grid(jpngr)%lat, & 
                          interface%grid(jpngr)%elv, & 
                          interface%climate(jpngr)%dfsun(:), & 
                          interface%climate(jpngr)%dppfd(:)  & 
                          )
        if (verbose) write(0,*) '... done'

        !----------------------------------------------------------------
        ! Get monthly light use efficiency, and Rd per unit of light absorbed
        ! Photosynthetic parameters acclimate at monthly time scale
        ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! there is a monthly loop within 'getlue'!
        !----------------------------------------------------------------
        if (verbose) write(0,*) 'calling getlue() ... '
        if (verbose) write(0,*) '    with argument CO2  = ', interface%pco2
        if (verbose) write(0,*) '    with argument temp.= ', interface%climate(jpngr)%dtemp(:)
        if (verbose) write(0,*) '    with argument VPD  = ', interface%climate(jpngr)%dvpd(:)
        if (verbose) write(0,*) '    with argument elv. = ', interface%grid(jpngr)%elv
        out_pmodel(:,:) = getlue( &
                                  interface%pco2, & 
                                  interface%climate(jpngr)%dtemp(:), & 
                                  interface%climate(jpngr)%dvpd(:), & 
                                  interface%grid(jpngr)%elv & 
                                  )
        if (verbose) write(0,*) '... done'

        !----------------------------------------------------------------
        ! LOOP THROUGH MONTHS
        !----------------------------------------------------------------
        doy=0
        monthloop: do moy=1,nmonth

          !----------------------------------------------------------------
          ! LOOP THROUGH DAYS
          !----------------------------------------------------------------
          dayloop: do dm=1,ndaymonth(moy)
            doy=doy+1

            ! if (verbose) write(0,*) '----------------------'
            ! if (verbose) write(0,*) 'YEAR, Doy ', interface%steering%year, doy
            ! if (verbose) write(0,*) '----------------------'

            !----------------------------------------------------------------
            ! initialise daily updated variables 
            !----------------------------------------------------------------
            call initdaily_plant()
            call initdaily_waterbal()
            call initdaily_gpp()

            !----------------------------------------------------------------
            ! get soil moisture, and runoff
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling waterbal() ... '
            call waterbal( &
                            tile(:,jpngr)%soil%phy, &
                            doy, & 
                            interface%grid(jpngr)%lat, & 
                            interface%grid(jpngr)%elv, & 
                            interface%climate(jpngr)%dprec(doy), & 
                            interface%climate(jpngr)%dtemp(doy), & 
                            interface%climate(jpngr)%dfsun(doy), &
                            interface%climate(jpngr)%dnetrad(doy)&
                            )
            if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! update canopy and tile variables and simulate daily 
            ! establishment / sprouting
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling vegdynamics() ... '
            call vegdynamics( tile(:,jpngr), plant(:,jpngr), solar, out_pmodel(:,:), interface%dfapar_field(doy,jpngr) )
            if (verbose) write(0,*) '... done'

            !/////////////////////////////////////////////////////////////////
            ! calculate GPP
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling gpp() ... '
            ! print*,'acrown', plant(1,1)%acrown
            ! print*,'in biosphere: fapar', plant(1,1)%fapar_ind
            ! print*,'in biosphere: acrown', plant(1,1)%acrown
            call gpp( &
                      out_pmodel(:,moy), solar, plant(:,jpngr), doy, moy, &
                      interface%climate(jpngr)%dtemp(doy) &
                      )
            if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! collect from daily updated state variables for annual variables
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling getout_daily() ... '
            call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy )
            call getout_daily_gpp( out_pmodel(:,moy), jpngr, doy )
            call getout_daily_plant( plant(:,jpngr), jpngr, moy, doy )
            call getout_daily_forcing( jpngr, moy, doy )
            if (verbose) write(0,*) '... done'

          end do dayloop

        end do monthloop

        !----------------------------------------------------------------
        ! collect annual output
        !----------------------------------------------------------------
        call getout_annual_plant( plant(:,jpngr), jpngr )
        call getout_annual_gpp( jpngr )

      end if
    end do gridcellloop

    !----------------------------------------------------------------
    ! Write to ascii output
    !----------------------------------------------------------------
    call writeout_ascii_waterbal()
    call writeout_ascii_gpp()
    call writeout_ascii_plant()
    call writeout_ascii_forcing()

    !----------------------------------------------------------------
    ! Write to NetCDF output
    !----------------------------------------------------------------
    call writeout_nc_forcing()
    call writeout_nc_plant()

    ! xxx insignificant
    c_uptake = 0.0

  end function biosphere_annual

end module md_biosphere