module md_biosphere

  use md_params_core
  use md_classdefs
  use md_plant, only: plant_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, psoilphystype, waterbal, getsolar, initdaily_waterbal, initglobal_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_gpp, only: outtype_pmodel, getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, getout_daily_gpp, getout_annual_gpp, writeout_ascii_gpp
  use md_vegdynamics, only: vegdynamics
  use md_tile, only: tile_type, initglobal_tile

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  type( tile_type ) , dimension(nlu,maxgrid)     :: tile
  type( plant_type ), dimension(npft,maxgrid)    :: plant ! npft counts over PFTs in all land units (tiles)
  type( solartype )                              :: solar
  type( psoilphystype ), dimension(nlu,maxgrid)  :: psoilphys
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
    use md_interface
  
    ! return variable
    real :: c_uptake   ! annual net global C uptake by biosphere (gC/yr)

    ! local variables
    integer :: dm, moy, jpngr, day

    ! xxx verbose
    logical, parameter :: verbose = .false.
    real            :: cbal1, cbal2
    type( orgpool ) :: orgtmp1, orgtmp2, orgbal1, orgbal2
    real :: eps = 9.999e-11

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (interface%steering%init) then

      !----------------------------------------------------------------
      ! GET MODEL PARAMETERS
      ! read model parameters that may be varied for optimisation
      !----------------------------------------------------------------
      call getpar_modl_plant()
      call getpar_modl_waterbal()
      call getpar_modl_gpp()

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      call initglobal_tile( tile(:,:) )
      call initglobal_plant( plant(:,:) )
      call initglobal_waterbal( psoilphys(:,:) )

      !----------------------------------------------------------------
      ! Open input/output files
      !----------------------------------------------------------------
      call initio_waterbal()
      call initio_gpp()
      call initio_plant()

    endif 

    !----------------------------------------------------------------
    ! Initialise output variables for this year
    !----------------------------------------------------------------
    call initoutput_waterbal()
    call initoutput_gpp()
    call initoutput_plant()

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
      solar = getsolar( &
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
      out_pmodel(:,:) = getlue( &
        interface%pco2, & 
        interface%climate(jpngr)%dtemp(:), & 
        interface%climate(jpngr)%dvpd(:), & 
        interface%grid(jpngr)%elv & 
        )

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
          call initdaily_plant()
          call initdaily_waterbal()
          call initdaily_gpp()

          !----------------------------------------------------------------
          ! get soil moisture, and runoff
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling waterbal() ... '
          call waterbal( &
            psoilphys(:,jpngr), &
            day, & 
            interface%grid(jpngr)%lat, & 
            interface%grid(jpngr)%elv, & 
            interface%climate(jpngr)%dprec(day), & 
            interface%climate(jpngr)%dtemp(day), & 
            interface%climate(jpngr)%dfsun(day)  &
            )
          if (verbose) write(0,*) '... done'

          !----------------------------------------------------------------
          ! update canopy and tile variables and simulate daily 
          ! establishment / sprouting
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling vegdynamics() ... '
          call vegdynamics( tile(:,jpngr), plant(:,jpngr), solar, out_pmodel(:,:), interface%mfapar_field(moy,jpngr) )
          if (verbose) write(0,*) '... done'

          !/////////////////////////////////////////////////////////////////
          ! calculate GPP
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling gpp() ... '
          ! print*,'nind', tile(1,1)%nind
          ! print*,'acrown', plant(1,1)%acrown
          ! print*,'in biosphere: fapar', plant(1,1)%fapar_ind
          ! print*,'in biosphere: acrown', plant(1,1)%acrown
          ! print*,'in biosphere: nind', tile(1,1)%nind
          call gpp( &
            out_pmodel(:,moy), solar, plant(:,jpngr), tile(:,jpngr), day, moy, &
            interface%climate(jpngr)%dtemp(day) &
            )
          if (verbose) write(0,*) '... done'

          !----------------------------------------------------------------
          ! collect from daily updated state variables for annual variables
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling getout_daily() ... '
          call getout_daily_waterbal( jpngr, moy, day, solar, psoilphys(:,jpngr) )
          call getout_daily_gpp( out_pmodel(:,moy), jpngr, day )
          call getout_daily_plant( plant(:,jpngr), jpngr, moy, day )
          if (verbose) write(0,*) '... done'

        end do

      end do

      !----------------------------------------------------------------
      ! collect annual output
      !----------------------------------------------------------------
      call getout_annual_plant( plant(:,jpngr), jpngr )
      call getout_annual_gpp( jpngr )

      !----------------------------------------------------------------
      ! Write to output
      !----------------------------------------------------------------
      call writeout_ascii_waterbal()
      call writeout_ascii_gpp()
      call writeout_ascii_plant()

    end do

    ! xxx insignificant
    c_uptake = 0.0

  end function biosphere_annual

end module md_biosphere