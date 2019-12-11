module md_biosphere

  use md_params_core
  use md_classdefs
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, get_solar, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, initio_nc_waterbal, writeout_nc_waterbal, init_rlm_waterbal, get_rlm_waterbal, getrlm_daily_waterbal
  use md_tile, only: tile_type, tile_fluxes_type, initglobal_tile, initdaily_tile
  use md_interface, only: getout_daily_forcing, initoutput_forcing, initio_nc_forcing, writeout_nc_forcing

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  ! derived types from L1 modules
  type( tile_type ),         allocatable, dimension(:,:) :: tile
  type( tile_fluxes_type ),  allocatable, dimension(:)   :: tile_fluxes

  ! derived types from L2 modules
  type( solartype ) :: solar

contains

  function biosphere_annual() result( out_biosphere )
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface, only: interface, outtype_biosphere
  
    ! return variable
    type(outtype_biosphere) :: out_biosphere

    ! local variables
    integer :: dm, moy, jpngr, doy

    ! xxx debug
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
      call getpar_modl_waterbal()
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      if (verbose) print*, 'initglobal_() ...'
      allocate( tile(  nlu,  size(interface%grid) ) )
      allocate( tile_fluxes(  nlu ) )

      call initglobal_tile(  tile(:,:),  size(interface%grid) )

      if (verbose) print*, '... done'

    endif

    !----------------------------------------------------------------
    ! Open NetCDF output files (one for each year)
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*, 'initio_nc_() ...'
      call initio_nc_forcing()
      call initio_nc_waterbal()
      if (verbose) print*, '... done'
    end if
    
    !----------------------------------------------------------------
    ! Initialise output variables for this year
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*, 'initoutput_() ...'
      call initoutput_forcing(  size(interface%grid) )
      call initoutput_waterbal( size(interface%grid) )
      if (verbose) print*, '... done'
    end if

    ! additional initialisation for rolling annual mean calculations (also needed in calibration mode)
    call init_rlm_waterbal( size(interface%grid) )

    !----------------------------------------------------------------
    ! LOOP THROUGH GRIDCELLS
    !----------------------------------------------------------------
    if (verbose) print*,'looping through gridcells ...'
    gridcellloop: do jpngr=1,size(interface%grid)

      if (interface%grid(jpngr)%dogridcell) then

        if (verbose) print*,'----------------------'
        if (verbose) print*,'JPNGR: ', jpngr
        if (verbose) print*,'----------------------'

        !----------------------------------------------------------------
        ! Get radiation based on daily temperature, sunshine fraction, and 
        ! elevation.
        ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! there is a daily loop within 'get_solar'!
        !----------------------------------------------------------------
        if (verbose) print*,'calling get_solar() ... '
        if (verbose) print*,'    with argument lat = ', interface%grid(jpngr)%lat
        if (verbose) print*,'    with argument elv = ', interface%grid(jpngr)%elv
        if (verbose) print*,'    with argument dfsun (ann. mean) = ', sum( interface%climate(jpngr)%dfsun(:) / ndayyear )
        if (verbose) print*,'    with argument dppfd (ann. mean) = ', sum( interface%climate(jpngr)%dppfd(:) / ndayyear )
        solar = get_solar( &
                          interface%grid(jpngr)%lat, & 
                          interface%grid(jpngr)%elv, & 
                          interface%climate(jpngr)%dfsun(:), & 
                          interface%climate(jpngr)%dppfd(:)  & 
                          )
        if (verbose) print*,'... done'

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

            if (verbose) print*,'----------------------'
            if (verbose) print*,'YEAR, DOY ', interface%steering%year, doy
            if (verbose) print*,'----------------------'

            !----------------------------------------------------------------
            ! initialise daily updated variables 
            !----------------------------------------------------------------
            if (verbose) print*,'calling initdaily_() ...'
            call initdaily_tile( tile_fluxes(:) )
            if (verbose) print*,'... done.'

            !----------------------------------------------------------------
            ! get soil moisture, and runoff
            !----------------------------------------------------------------
            if (verbose) print*,'calling waterbal() ... '
            ! print*,'lon,lat,ilon,ilat,jpngr', interface%grid(jpngr)%lon, interface%grid(jpngr)%lat, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, jpngr
            call waterbal( &
                          tile(:,jpngr)%soil, &
                          tile_fluxes(:), &
                          plant_fluxes(:), &
                          doy, &
                          jpngr, & 
                          interface%grid(jpngr)%lat,             & 
                          interface%grid(jpngr)%elv,             & 
                          interface%climate(jpngr)%dprec(doy),   & 
                          interface%climate(jpngr)%dsnow(doy),   & 
                          interface%climate(jpngr)%dtemp(doy),   & 
                          interface%climate(jpngr)%dfsun(doy),   &
                          interface%climate(jpngr)%dnetrad(doy), &
                          interface%climate(jpngr)%dvpd(doy),    &
                          interface%vegcover(jpngr)%dfapar(doy)  &
                          )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! collect from daily updated state variables for annual variables
            !----------------------------------------------------------------
            if (.not.interface%params_siml%is_calib) then
              if (verbose) print*,'calling getout_daily() ... '
              call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy, tile_fluxes(:) )
              call getout_daily_forcing( jpngr, moy, doy )
              if (verbose) print*,'... done'
            end if

            call getrlm_daily_waterbal( jpngr, doy )

            !----------------------------------------------------------------
            ! populate function return variable
            !----------------------------------------------------------------
            ! if (npft>1) stop 'think about npft > 1'
            out_biosphere%fapar(doy)   = 0.0
            out_biosphere%gpp(doy)     = 0.0
            out_biosphere%transp(doy)  = 0.0

          end do dayloop

        end do monthloop

        ! !----------------------------------------------------------------
        ! ! collect annual output
        ! !----------------------------------------------------------------
        ! if (.not.interface%params_siml%is_calib) then
        !   if (verbose) print*,'calling getout_annual_() ... '
        !   if (verbose) print*,'... done'
        ! end if

      end if

    end do gridcellloop

    !----------------------------------------------------------------
    ! Write to NetCDF output
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*,'calling writeout_nc_() ... '
      call writeout_nc_forcing()
      call writeout_nc_waterbal()
      if (verbose) print*,'... done'
    end if

    if (interface%steering%finalize) then
      !----------------------------------------------------------------
      ! Finazlize run: deallocating memory
      !----------------------------------------------------------------
      deallocate( tile )
      deallocate( tile_fluxes )
    end if

    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end function biosphere_annual

end module md_biosphere
