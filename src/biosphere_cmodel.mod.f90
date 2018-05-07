module md_biosphere

  use md_params_core
  use md_classdefs
  use md_plant, only: plant_type, plant_fluxes_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, getsolar, initdaily_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal, initio_nc_waterbal, writeout_nc_waterbal, get_rlm_waterbal
  use md_gpp, only: outtype_pmodel, getpar_modl_gpp, initio_gpp, initoutput_gpp, getlue, gpp, getout_daily_gpp, getout_annual_gpp, writeout_ascii_gpp, initio_nc_gpp, writeout_nc_gpp
  use md_vegdynamics, only: vegdynamics
  use md_tile, only: tile_type, initglobal_tile
  use md_interface, only: getout_daily_forcing, initoutput_forcing, initio_forcing, initio_nc_forcing, writeout_ascii_forcing, writeout_nc_forcing
  use md_phenology, only: temppheno_type, get_temppheno, getpar_modl_phenology
  use md_allocation, only: allocation_daily
  use md_soiltemp, only: getout_daily_soiltemp, soiltemp
  use md_turnover, only: turnover, initoutput_turnover
  use md_littersom, only: littersom
  use md_npp, only: npp

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  type( tile_type ),         allocatable, dimension(:,:) :: tile
  type( plant_type ),        allocatable, dimension(:,:) :: plant
  type( plant_fluxes_type ), allocatable, dimension(:)   :: plant_fluxes

  type( solartype )                              :: solar
  type( outtype_pmodel ), dimension(npft,nmonth) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

  type( temppheno_type ), dimension(ndayyear,npft) :: out_temppheno

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

    ! xxx debug
    logical, parameter :: verbose = .false.
    logical, parameter :: splashtest = .false.
    integer, parameter :: lev_splashtest = 2
    integer, parameter :: testdoy = 55
    real            :: cbal1, cbal2
    type( orgpool ) :: orgtmp1, orgtmp2, orgbal1, orgbal2

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
      call getpar_modl_phenology()
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      if (verbose) print*, 'initglobal_() ...'
      allocate( tile(  nlu,  size(interface%grid) ) )
      allocate( plant( npft, size(interface%grid) ) )
      allocate( plant_fluxes( npft ) )

      call initglobal_tile(  tile(:,:),  size(interface%grid) )
      call initglobal_plant( plant(:,:), size(interface%grid) )
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
    if (verbose) print*, 'initio_nc_() ...'
    call initio_nc_forcing()
    call initio_nc_gpp()
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
    call initoutput_turnover( size(interface%grid) )
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
          interface%climate(jpngr)%dfsun(:) = 0.562000036
          interface%climate(jpngr)%dppfd(:) = dummy
          solar = getsolar( 67.25, 87.0, interface%climate(jpngr)%dfsun(:), interface%climate(jpngr)%dppfd(:), splashtest=splashtest, testdoy=testdoy )
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
        ! Temperture-driven phenology is based simply on temperature.
        !----------------------------------------------------------------
        out_temppheno(:,:) = get_temppheno( jpngr, interface%climate(jpngr)%dtemp(:) )
        ! call get_temppheno( jpngr, interface%climate(jpngr)%dtemp(:) )

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
            call initdaily_plant( plant_fluxes(:) )
            call initdaily_waterbal()

            !----------------------------------------------------------------
            ! get soil moisture, and runoff
            !----------------------------------------------------------------
            if (verbose) print*,'calling waterbal() ... '
            ! print*,'lon,lat,ilon,ilat,jpngr', interface%grid(jpngr)%lon, interface%grid(jpngr)%lat, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, jpngr
            call waterbal( &
                            tile(:,jpngr)%soil%phy, &
                            doy, jpngr, & 
                            interface%grid(jpngr)%lat, & 
                            interface%grid(jpngr)%elv, & 
                            interface%soilparams(jpngr), &
                            interface%climate(jpngr)%dprec(doy), & 
                            interface%climate(jpngr)%dtemp(doy), & 
                            interface%climate(jpngr)%dfsun(doy), &
                            interface%climate(jpngr)%dnetrad(doy), &
                            splashtest=splashtest, lev_splashtest=lev_splashtest, testdoy=testdoy &
                            )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! calculate soil temperature
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling soiltemp() ... '
            call soiltemp(&
                          tile(:,jpngr)%soil%phy, &
                          interface%climate(jpngr)%dtemp(:), &
                          size(interface%grid), &
                          interface%steering%init, &
                          jpngr, & 
                          moy, & 
                          doy & 
                          )
            if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! update canopy and tile variables and simulate daily 
            ! establishment / sprouting
            !----------------------------------------------------------------
            if (verbose) print*,'calling vegdynamics() ... '
            call vegdynamics( plant(:,jpngr), out_temppheno(doy,:) )
            ! call vegdynamics( tile(:,jpngr), &
            !                   plant(:,jpngr), &
            !                   solar, &
            !                   out_pmodel(:,:), &
            !                   interface%dfapar_field(doy,jpngr), &
            !                   interface%fpc_grid(:,jpngr) &
            !                   )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! calculate GPP
            !----------------------------------------------------------------
            if (verbose) print*,'calling gpp() ... '
            ! print*,'acrown', plant(1,1)%acrown
            ! print*,'in biosphere: fapar', plant(1,1)%fapar_ind
            ! print*,'in biosphere: acrown', plant(1,1)%acrown
            call gpp( &
                      out_pmodel(:,moy), solar, plant(:,jpngr), &
                      plant_fluxes(:), &
                      tile(:,jpngr)%soil%phy, &
                      doy, moy, &
                      interface%climate(jpngr)%dtemp(doy), &
                      interface%params_siml%soilmstress &
                      )
            if (verbose) print*,'... done'

            !----------------------------------------------------------------
            ! substract autotrophic respiration to get NPP, remainder is added 
            ! to labile pool (plabl)
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling npp() ... '
            if (verbose) write(0,*) '              with state variables:'
            if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) orgtmp1 =  plant(1,jpngr)%plabl
            !----------------------------------------------------------------
            call npp( plant(:,jpngr), plant_fluxes(:), interface%climate(jpngr)%dtemp(doy) )
              ! jpngr, interface%climate(jpngr)%dtemp(day), day )
            !----------------------------------------------------------------
            if (verbose) write(0,*) '              ==> returned: '
            if (verbose) write(0,*) '              dnpp  = ', plant_fluxes(:)%dnpp
            if (verbose) write(0,*) '              dcex  = ', plant_fluxes(:)%dcex
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) write(0,*) '              dlabl = ', orgminus( plant(1,jpngr)%plabl, orgtmp1 )
            if (verbose) write(0,*) '    --- balance: '
            if (verbose) cbal1 = plant_fluxes(1)%dgpp - plant_fluxes(1)%dnpp%c12 - plant_fluxes(1)%drleaf - plant_fluxes(1)%drroot
            if (verbose) cbal2 = plant_fluxes(1)%dgpp - ( plant(1,jpngr)%plabl%c%c12 - orgtmp1%c%c12 ) - plant_fluxes(1)%dcex - plant_fluxes(1)%drleaf - plant_fluxes(1)%drroot
            if (verbose) write(0,*) '        gpp - npp - ra_maint          = ', cbal1
            if (verbose) write(0,*) '        gpp - dlabl - dcex - ra_maint = ', cbal2
            if (verbose.and.abs(cbal1)>eps) stop 'balance 1 not satisfied'
            if (verbose.and.abs(cbal2)>eps) stop 'balance 2 not satisfied'
            if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! leaf, sapwood, and fine-root turnover
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling turnover() ... '
            if (verbose) write(0,*) '              with state variables:'
            if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) write(0,*) '              plitt = ', orgplus( plant(1,jpngr)%plitt_af, plant(1,jpngr)%plitt_as, plant(1,jpngr)%plitt_bg )
            if (verbose) orgtmp1 = orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot, plant(1,jpngr)%plabl )
            if (verbose) orgtmp2 = orgplus( plant(1,jpngr)%plitt_af, plant(1,jpngr)%plitt_as, plant(1,jpngr)%plitt_bg )
            !----------------------------------------------------------------
            call turnover( plant(:,jpngr), out_pmodel(:,moy), solar, out_temppheno(doy,:) )
            !----------------------------------------------------------------
            if (verbose) write(0,*) '              ==> returned: '
            if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) write(0,*) '              plitt = ', orgplus( plant(1,jpngr)%plitt_af, plant(1,jpngr)%plitt_as, plant(1,jpngr)%plitt_bg )
            if (verbose) write(0,*) '   --- balance: '
            if (verbose) orgbal1 = orgminus( orgminus(   orgplus( plant(1,jpngr)%plitt_af, plant(1,jpngr)%plitt_as, plant(1,jpngr)%plitt_bg ),   orgtmp2   ), orgminus(   orgtmp1,   orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot, plant(1,jpngr)%plabl )   ) )
            if (verbose) write(0,*) '       dlitt - dplant                = ', orgbal1
            if (verbose .and. abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
            if (verbose .and. abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
            if (verbose) write(0,*) '... done'

            ! !----------------------------------------------------------------
            ! ! grass / crop harvest
            ! !----------------------------------------------------------------
            ! if (verbose) write(0,*) 'calling grharvest() ... '
            ! if (verbose) write(0,*) '              with state variables:'
            ! if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            ! if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            ! if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            ! if (verbose) write(0,*) '              mharv = ', mharv(:,jpngr)
            ! if (verbose) orgtmp1 =  orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot, plant(1,jpngr)%plabl )
            ! if (verbose) orgtmp2 =  mharv(1,jpngr)
            ! !----------------------------------------------------------------
            ! call grharvest( jpngr, day )
            ! !----------------------------------------------------------------
            ! if (verbose) write(0,*) '              ==> returned: '
            ! if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            ! if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            ! if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            ! if (verbose) write(0,*) '              mharv = ', mharv(:,jpngr)
            ! if (verbose) write(0,*) '    --- balance: '
            ! if (verbose) orgbal1 = orgminus( orgminus( orgtmp1, orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot, plant(1,jpngr)%plabl ) ), orgminus( mharv(1,jpngr), orgtmp2 ) )
            ! if (verbose) write(0,*) '        dharv - dplant  = ', orgbal1
            ! if (verbose.and.abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
            ! if (verbose.and.abs(orgbal1%n%n14)>eps) stop 'balance not satisfied for N'
            ! if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! allocation of labile pools to biomass
            !----------------------------------------------------------------
            if (verbose) write(0,*) 'calling allocation() ... '
            if (verbose) write(0,*) '              with state variables:'
            if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf
            if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) orgtmp1 =  plant(1,jpngr)%plabl
            if (verbose) orgtmp2 =  orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot )
            !----------------------------------------------------------------
            call allocation_daily( plant(:,jpngr), plant_fluxes(:), solar, out_pmodel(:,:), interface%climate(jpngr)%dtemp(doy) )
            !----------------------------------------------------------------
            if (verbose) write(0,*) '              ==> returned: '
            if (verbose) write(0,*) '              pleaf = ', plant(:,jpngr)%pleaf, ' C:N = ', cton( plant(1,jpngr)%pleaf )
            if (verbose) write(0,*) '              proot = ', plant(:,jpngr)%proot, ' C:N = ', cton( plant(1,jpngr)%proot )
            if (verbose) write(0,*) '              plabl = ', plant(:,jpngr)%plabl
            if (verbose) write(0,*) '   --- balance: '
            if (verbose) orgbal1 = orgminus( orgminus( orgminus( orgtmp1, plant(1,jpngr)%plabl ), orgpool( carbon(plant_fluxes(1)%drgrow), nitrogen(0.0) ) ), orgminus( orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot ), orgtmp2 ) )        
            ! dlabl  = orgminus( orgtmp1, plant(1,jpngr)%plabl )
            ! drgrow = orgpool( carbon(drgrow(1)), nitrogen(0.0) )
            ! dplant = orgminus( orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot ), orgtmp2 )
            if (verbose) write(0,*) '       dlabl - drgrow - dleaf - droot=', orgbal1
            if (verbose.and.abs(orgbal1%c%c12)>eps) stop 'balance not satisfied for C'
            ! print*, 'dlabl  =', orgminus( orgtmp, plant(1,jpngr)%plabl )
            ! print*, 'drgrow =', drgrow(1), 0.0
            ! print*, 'dplant =', orgminus( orgplus( plant(1,jpngr)%pleaf, plant(1,jpngr)%proot ), orgtmp2 )
            if (verbose) write(0,*) '... done'

            !----------------------------------------------------------------
            ! litter and soil decomposition and N mineralisation
            !----------------------------------------------------------------
            print*, 'calling littersom() ... '
            call littersom( plant(:,jpngr), tile(:,jpngr)%soil%phy, jpngr, doy, interface%climate(jpngr)%dtemp(doy) )
            print*, '... done'

            !----------------------------------------------------------------
            ! collect from daily updated state variables for annual variables
            !----------------------------------------------------------------
            if (verbose) print*,'calling getout_daily() ... '
            call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy )
            call getout_daily_gpp( out_pmodel(:,moy), plant_fluxes(:), jpngr, doy )
            call getout_daily_plant( plant_fluxes(:), jpngr, moy, doy )
            call getout_daily_forcing( jpngr, moy, doy )
            call getout_daily_soiltemp( jpngr, moy, doy, tile(:,jpngr)%soil%phy )
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
    print*,'... gridcellloop finished.'

    !----------------------------------------------------------------
    ! Get rolling multi-year averages (needs to store entire arrays)
    !----------------------------------------------------------------
    call get_rlm_waterbal( tile(:,:)%soil%phy )

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
    call writeout_nc_forcing()
    call writeout_nc_gpp()
    call writeout_nc_waterbal()

    ! xxx insignificant
    c_uptake = 0.0

  end function biosphere_annual

end module md_biosphere
