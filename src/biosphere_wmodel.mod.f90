module md_biosphere

  use md_params_core
  use md_params_siml
  use md_params_site
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: solartype, waterbal, getsolar, initdaily_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_tile, only: tile_type, initglobal_tile

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  type( tile_type ) , dimension(nlu,maxgrid)     :: tile
  type( solartype )                              :: solar

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
    integer :: dm, moy, jpngr, day

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
      call getpar_modl_waterbal()

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      call initglobal_tile( tile(:,:) )

      !----------------------------------------------------------------
      ! Open input/output files
      !----------------------------------------------------------------
      call initio_waterbal()

    endif 

    !----------------------------------------------------------------
    ! Initialise output variables for this year
    !----------------------------------------------------------------
    call initoutput_waterbal()

    !----------------------------------------------------------------
    ! LOOP THROUGH GRIDCELLS
    !----------------------------------------------------------------
    gridcellloop: do jpngr=1,maxgrid

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
        interface%climate(jpngr)%dfsun(:) & 
        )
      if (verbose) write(0,*) '... done'

      !----------------------------------------------------------------
      ! LOOP THROUGH MONTHS
      !----------------------------------------------------------------
      day=0
      monthloop: do moy=1,nmonth

        !----------------------------------------------------------------
        ! LOOP THROUGH DAYS
        !----------------------------------------------------------------
        dayloop: do dm=1,ndaymonth(moy)
          day=day+1

          if (verbose) write(0,*) '----------------------'
          if (verbose) write(0,*) 'YEAR, DAY ', interface%steering%year, day
          if (verbose) write(0,*) '----------------------'

          !----------------------------------------------------------------
          ! initialise daily updated variables 
          !----------------------------------------------------------------
          call initdaily_waterbal()

          !----------------------------------------------------------------
          ! get soil moisture, and runoff
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling waterbal() ... '
          call waterbal( &
            tile(:,jpngr)%soil%phy, &
            day, & 
            interface%grid(jpngr)%lat, & 
            interface%grid(jpngr)%elv, & 
            interface%climate(jpngr)%dprec(day), & 
            interface%climate(jpngr)%dtemp(day), & 
            interface%climate(jpngr)%dfsun(day), &
            interface%climate(jpngr)%dnetrad(day)&
            )
          if (verbose) write(0,*) '... done'

          !----------------------------------------------------------------
          ! collect from daily updated state variables for annual variables
          !----------------------------------------------------------------
          if (verbose) write(0,*) 'calling getout_daily() ... '
          call getout_daily_waterbal( jpngr, moy, day, solar, tile(:,jpngr)%soil%phy )
          if (verbose) write(0,*) '... done'

        end do dayloop

      end do monthloop

      !----------------------------------------------------------------
      ! collect annual output
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      ! Write to output
      !----------------------------------------------------------------
      call writeout_ascii_waterbal()

    end do gridcellloop

    ! xxx insignificant
    c_uptake = 0.0

  end function biosphere_annual

end module md_biosphere
