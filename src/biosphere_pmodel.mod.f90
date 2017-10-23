module md_biosphere

  use md_params_core
  use md_classdefs
  use md_plant, only: plant_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant, initio_nc_plant, writeout_nc_plant
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, getsolar, initdaily_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal, initio_nc_waterbal, writeout_nc_waterbal
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
    logical, parameter :: splashtest = .false.
    integer, parameter :: lev_splashtest = 2
    integer, parameter :: testdoy = 8

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
      if (verbose) print*, 'initglobal_() ...'
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

    print*,'----------------------'
    print*,'YEAR ', interface%steering%year
    print*,'----------------------'

    !----------------------------------------------------------------
    ! Open NetCDF output files (one for each year)
    !----------------------------------------------------------------
    if (verbose) print*, 'initio_nc_() ...'
    call initio_nc_forcing()
    call initio_nc_plant()
    call initio_nc_waterbal()
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
    print*,'looping through gridcells ...'
    gridcellloop: do jpngr=1,size(interface%grid)
    ! gridcellloop: do jpngr=48790,48790   ! negative PET in january 2010
      if (interface%grid(jpngr)%dogridcell) then

        if (verbose) print*,'----------------------'
        if (verbose) print*,'jpngr ', jpngr
        if (verbose) print*,'----------------------'

        !----------------------------------------------------------------
        ! Get radiation based on daily temperature, sunshine fraction, and 
        ! elevation.
        ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! there is a daily loop within 'getsolar'!
        !----------------------------------------------------------------
        if (verbose) print*,'calling getsolar() ... '
        if (verbose) print*,'    with argument lat = ', interface%grid(jpngr)%lat
        if (verbose) print*,'    with argument elv = ', interface%grid(jpngr)%elv
        if (verbose) print*,'    with argument dfsun (ann. mean) = ', sum( interface%climate(jpngr)%dfsun(:) / ndayyear )
        if (verbose) print*,'    with argument dppfd (ann. mean) = ', sum( interface%climate(jpngr)%dppfd(:) / ndayyear )
        if (splashtest) then
          ! for comparison with Python SPLASH
          interface%climate(jpngr)%dfsun(:) = 0.422999978
          interface%climate(jpngr)%dppfd(:) = dummy
          solar = getsolar( 54.75, 954.0, interface%climate(jpngr)%dfsun(:), interface%climate(jpngr)%dppfd(:), splashtest=splashtest, testdoy=testdoy )
          if (lev_splashtest==1) stop 'end of splash test level 1'
        else          
          solar = getsolar( &
                            interface%grid(jpngr)%lat, & 
                            interface%grid(jpngr)%elv, & 
                            interface%climate(jpngr)%dfsun(:), & 
                            interface%climate(jpngr)%dppfd(:),  & 
                            splashtest = splashtest, testdoy=testdoy &
                            )
        end if
        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! Get monthly light use efficiency, and Rd per unit of light absorbed
        ! Photosynthetic parameters acclimate at ~monthly time scale
        ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! there is a monthly loop within 'getlue'!
        !----------------------------------------------------------------
        if (verbose) print*,'calling getlue() ... '
        if (verbose) print*,'    with argument CO2  = ', interface%pco2
        if (verbose) print*,'    with argument temp.= ', interface%climate(jpngr)%dtemp(:)
        if (verbose) print*,'    with argument VPD  = ', interface%climate(jpngr)%dvpd(:)
        if (verbose) print*,'    with argument elv. = ', interface%grid(jpngr)%elv
        out_pmodel(:,:) = getlue( &
                                  interface%pco2, & 
                                  interface%climate(jpngr)%dtemp(:), & 
                                  interface%climate(jpngr)%dvpd(:), & 
                                  interface%grid(jpngr)%elv & 
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
            if (verbose) print*,'YEAR, Doy ', interface%steering%year, doy
            if (verbose) print*,'----------------------'

            !----------------------------------------------------------------
            ! initialise daily updated variables 
            !----------------------------------------------------------------
            call initdaily_plant()
            call initdaily_waterbal()
            call initdaily_gpp()

            !----------------------------------------------------------------
            ! get soil moisture, and runoff
            !----------------------------------------------------------------
            if (verbose) print*,'calling waterbal() ... '
            ! print*,'lon,lat,ilon,ilat,jpngr', interface%grid(jpngr)%lon, interface%grid(jpngr)%lat, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, jpngr
            call waterbal( &
                            tile(:,jpngr)%soil%phy, &
                            doy, & 
                            interface%grid(jpngr)%lat, & 
                            interface%grid(jpngr)%elv, & 
                            interface%climate(jpngr)%dprec(doy), & 
                            interface%climate(jpngr)%dtemp(doy), & 
                            interface%climate(jpngr)%dfsun(doy), &
                            interface%climate(jpngr)%dnetrad(doy), &
                            splashtest=splashtest, lev_splashtest=lev_splashtest, testdoy=testdoy &
                            )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! update canopy and tile variables and simulate daily 
            ! establishment / sprouting
            !----------------------------------------------------------------
            if (verbose) print*,'calling vegdynamics() ... '
            call vegdynamics( tile(:,jpngr), plant(:,jpngr), solar, out_pmodel(:,:), interface%dfapar_field(doy,jpngr) )
            if (verbose) print*,'... done'

            !/////////////////////////////////////////////////////////////////
            ! calculate GPP
            !----------------------------------------------------------------
            if (verbose) print*,'calling gpp() ... '
            ! print*,'acrown', plant(1,1)%acrown
            ! print*,'in biosphere: fapar', plant(1,1)%fapar_ind
            ! print*,'in biosphere: acrown', plant(1,1)%acrown
            call gpp( &
                      out_pmodel(:,moy), solar, plant(:,jpngr), doy, moy, &
                      interface%climate(jpngr)%dtemp(doy) &
                      )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! collect from daily updated state variables for annual variables
            !----------------------------------------------------------------
            if (verbose) print*,'calling getout_daily() ... '
            call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy )
            call getout_daily_gpp( out_pmodel(:,moy), jpngr, doy )
            call getout_daily_plant( plant(:,jpngr), jpngr, moy, doy )
            call getout_daily_forcing( jpngr, moy, doy )
            if (verbose) print*,'... done'

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
    call writeout_nc_waterbal()

    ! xxx insignificant
    c_uptake = 0.0

  end function biosphere_annual

end module md_biosphere