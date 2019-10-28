module md_biosphere

  use md_params_core
  use md_classdefs
  use md_plant, only: plant_type, plant_fluxes_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, get_solar, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, initio_nc_waterbal, writeout_nc_waterbal, init_rlm_waterbal, get_rlm_waterbal, getrlm_daily_waterbal
  use md_gpp, only: outtype_pmodel, getpar_modl_gpp, initoutput_gpp, getlue, gpp, getout_daily_gpp, getout_annual_gpp, initio_nc_gpp, writeout_nc_gpp
  use md_vegdynamics, only: vegdynamics
  use md_tile, only: tile_type, tile_fluxes_type, initglobal_tile, initdaily_tile
  use md_interface, only: getout_daily_forcing, initoutput_forcing, initio_nc_forcing, writeout_nc_forcing
  use md_soiltemp, only: getout_daily_soiltemp, soiltemp, initoutput_soiltemp

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  ! derived types from L1 modules
  type( tile_type ),         allocatable, dimension(:,:) :: tile
  type( tile_fluxes_type ),  allocatable, dimension(:)   :: tile_fluxes
  type( plant_type ),        allocatable, dimension(:,:) :: plant
  type( plant_fluxes_type ), allocatable, dimension(:)   :: plant_fluxes

  ! derived types from L2 modules
  type( solartype )                              :: solar
  type( outtype_pmodel ), dimension(npft,nmonth) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

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
    use md_sofunutils, only: daily2monthly
  
    ! return variable
    type(outtype_biosphere) :: out_biosphere

    ! local variables
    integer :: dm, moy, jpngr, doy
    real, dimension(nmonth) :: mtemp      ! monthly mean air temperature (deg C)
    real, dimension(nmonth) :: mvpd       ! monthly mean vapour pressure deficit (Pa)
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
      if (verbose) print*, 'initglobal_() ...'
      allocate( tile(  nlu,  size(interface%grid) ) )
      allocate( tile_fluxes(  nlu ) )
      allocate( plant( npft, size(interface%grid) ) )
      allocate( plant_fluxes( npft ) )

      call initglobal_plant( plant(:,:), size(interface%grid) )
      call initglobal_tile(  tile(:,:),  size(interface%grid) )
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Open ascii output files
      !----------------------------------------------------------------
      ! if (verbose) print*, 'initio_() ...'
      ! call initio_waterbal()
      ! call initio_gpp()
      ! call initio_plant()
      ! call initio_forcing()
      ! if (verbose) print*, '... done'

    endif 

    !----------------------------------------------------------------
    ! Open NetCDF output files (one for each year)
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*, 'initio_nc_() ...'
      call initio_nc_forcing()
      call initio_nc_gpp()
      call initio_nc_waterbal()
      if (verbose) print*, '... done'
    end if
    
    !----------------------------------------------------------------
    ! Initialise output variables for this year
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*, 'initoutput_() ...'
      call initoutput_waterbal( size(interface%grid) )
      call initoutput_gpp(      size(interface%grid) )
      call initoutput_plant(    size(interface%grid) )
      call initoutput_forcing(  size(interface%grid) )
      call initoutput_soiltemp( size(interface%grid) )
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
        ! Run P-model
        ! to get monthly light use efficiency, Rd, and Vcmax per unit of 
        ! light absorbed light.
        ! Photosynthetic parameters acclimate at ~monthly time scale
        !----------------------------------------------------------------
        if (verbose) print*,'calling getlue() ... '
        if (verbose) print*,'    with argument CO2  = ', interface%pco2
        if (verbose) print*,'    with argument temp.= ', interface%climate(jpngr)%dtemp(1:10)
        if (verbose) print*,'    with argument VPD  = ', interface%climate(jpngr)%dvpd(1:10)
        if (verbose) print*,'    with argument elv. = ', interface%grid(jpngr)%elv
        out_pmodel(:,:) = getlue( &
                                  interface%pco2, & 
                                  interface%climate(jpngr)%dtemp(:), & 
                                  interface%climate(jpngr)%dvpd(:), & 
                                  interface%grid(jpngr)%elv & 
                                  )
        ! ! xxx trevortest
        ! interface%climate(jpngr)%dtemp(:) = 20.0
        ! interface%climate(jpngr)%dvpd(:) = 1000.0
        ! out_pmodel(:,:) = getlue( &
        !                           400.0, & 
        !                           interface%climate(jpngr)%dtemp(:), & 
        !                           interface%climate(jpngr)%dvpd(:), & 
        !                           0.0 & 
        !                           )
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
            if (verbose) print*,'calling initdaily_() ...'
            call initdaily_plant( plant_fluxes(:) )
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
                          interface%grid(jpngr)%lat, & 
                          interface%grid(jpngr)%elv, & 
                          interface%climate(jpngr)%dprec(doy), & 
                          interface%climate(jpngr)%dsnow(doy), & 
                          interface%climate(jpngr)%dtemp(doy), & 
                          interface%climate(jpngr)%dfsun(doy), &
                          interface%climate(jpngr)%dnetrad(doy), &
                          interface%dfapar_field(doy,jpngr), &
                          interface%climate(jpngr)%dvpd(doy) &
                          )
            if (verbose) print*,'... done'

            ! !----------------------------------------------------------------
            ! ! calculate soil temperature
            ! !----------------------------------------------------------------
            ! if (verbose) print*, 'calling soiltemp() ... '
            ! call soiltemp(&
            !               tile(:,jpngr)%soil, &
            !               interface%climate(jpngr)%dtemp(:), &
            !               size(interface%grid), &
            !               interface%steering%init, &
            !               jpngr, & 
            !               moy, & 
            !               doy & 
            !               )
            ! if (verbose) print*, '... done'

            !----------------------------------------------------------------
            ! update canopy and tile variables and simulate daily 
            ! establishment / sprouting
            !----------------------------------------------------------------
            if (verbose) print*,'calling vegdynamics() ... '
            call vegdynamics( tile(:,jpngr), &
                              plant(:,jpngr), &
                              solar, &
                              out_pmodel(:,:), &
                              interface%dfapar_field(doy,jpngr), &
                              interface%fpc_grid(:,jpngr) &
                              )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! calculate GPP
            !----------------------------------------------------------------
            if (verbose) print*,'calling gpp() ... '
            call gpp( plant(:,jpngr), &
                      plant_fluxes(:), &
                      out_pmodel(:,moy), &
                      solar%dppfd(doy), &
                      solar%dayl(doy), &
                      solar%meanmppfd(moy), &
                      tile(:,jpngr)%soil%phy%wscal, &
                      tile(:,jpngr)%soil%phy%rlmalpha, &
                      doy, &
                      moy, &
                      interface%climate(jpngr)%dtemp(doy), &
                      interface%params_siml%soilmstress, &
                      interface%params_siml%tempstress &
                      )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! collect from daily updated state variables for annual variables
            !----------------------------------------------------------------
            if (.not.interface%params_siml%is_calib) then
              if (verbose) print*,'calling getout_daily() ... '
              call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy, tile_fluxes(:) )
              call getout_daily_gpp( out_pmodel(:,moy), plant_fluxes(:), jpngr, doy )
              call getout_daily_plant( plant(:,jpngr), plant_fluxes(:), jpngr, moy, doy )
              call getout_daily_forcing( jpngr, moy, doy )
              call getout_daily_soiltemp( jpngr, moy, doy, tile(:,jpngr)%soil%phy )
              if (verbose) print*,'... done'
            end if

            call getrlm_daily_waterbal( jpngr, doy )

            !----------------------------------------------------------------
            ! populate function return variable
            !----------------------------------------------------------------
            !if (npft>1) stop 'think about npft > 1'
            out_biosphere%fapar(doy)   = plant(1,jpngr)%fapar_ind
            out_biosphere%gpp(doy)     = plant_fluxes(1)%dgpp
            out_biosphere%transp(doy)  = plant_fluxes(1)%dtransp
            out_biosphere%latenth(doy) = plant_fluxes(1)%dlatenth

          end do dayloop

        end do monthloop

        !----------------------------------------------------------------
        ! collect annual output
        !----------------------------------------------------------------
        if (.not.interface%params_siml%is_calib) then
          if (verbose) print*,'calling getout_annual_() ... '
          call getout_annual_plant( plant(:,jpngr), jpngr )
          call getout_annual_gpp( jpngr )
          if (verbose) print*,'... done'
        end if

      end if
    end do gridcellloop

    !----------------------------------------------------------------
    ! Get rolling multi-year averages (needs to store entire arrays)
    !----------------------------------------------------------------
    call get_rlm_waterbal( tile(:,:)%soil%phy, interface%steering%init )

    !----------------------------------------------------------------
    ! Write to ascii output
    !----------------------------------------------------------------
    ! call writeout_ascii_waterbal()
    ! call writeout_ascii_gpp()
    ! call writeout_ascii_plant()
    ! call writeout_ascii_forcing()

    !----------------------------------------------------------------
    ! Write to NetCDF output
    !----------------------------------------------------------------
    if (.not.interface%params_siml%is_calib) then
      if (verbose) print*,'calling writeout_nc_() ... '
      call writeout_nc_forcing()
      call writeout_nc_gpp()
      call writeout_nc_waterbal()
      if (verbose) print*,'... done'
    end if


    if (interface%steering%finalize) then
      !----------------------------------------------------------------
      ! Finazlize run: deallocating memory
      !----------------------------------------------------------------
      deallocate( tile )
      deallocate( tile_fluxes )
      deallocate( plant )
      deallocate( plant_fluxes )
    end if

    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end function biosphere_annual

end module md_biosphere
