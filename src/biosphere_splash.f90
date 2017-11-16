subroutine biosphere( &
  year, lon, lat, elv &
  , params_soil_field, lu_area, pco2 &
  , dtemp_field, dprec_field &
  , dfsun_field, dvpd_field, dndep_field &
  , c_uptake &
  , mfapar_field &
  ) 

  !////////////////////////////////////////////////////////////////
  ! Subroutine BIOSPHERE calculates net ecosystem exchange (nee)
  ! in response to environmental boundary conditions (atmospheric 
  ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
  ! LPJ, also formulated as subroutine.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core
  use md_params_siml
  use md_params_site
  use md_soiltemp, only: soiltemp, ((interface%steering%init))output_soiltemp, ((interface%steering%init))io_soiltemp, getout_daily_soiltemp, writeout_ascii_soiltemp
  use md_waterbal, only: waterbal, getsolar_alldays, ((interface%steering%init))daily_waterbal, ((interface%steering%init))global_waterbal, ((interface%steering%init))io_waterbal, getout_daily_waterbal, ((interface%steering%init))output_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_params_soil, only: paramtype_soil

  implicit none

  ! arguments
  integer, intent(in)                           :: year       ! simulation year
  real, intent(in), dimension(maxgrid)          :: lon        ! longitude vector/field (degrees E)              
  real, intent(in), dimension(maxgrid)          :: lat        ! latitude vector/field (degrees N)             
  real, intent(in), dimension(maxgrid)          :: elv        ! elevation (altitude) vector/field (m above sea level)                  
  type(paramtype_soil), intent(in), dimension(maxgrid) :: params_soil_field
  real, dimension(3)                            :: lu_area    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
  real, intent(in)                              :: pco2
  real, intent(in), dimension(ndayyear,maxgrid) :: dtemp_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dprec_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dfsun_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dvpd_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dndep_field
  real, intent(out)                             :: c_uptake   ! annual net global C uptake by biosphere

  ! optional arguments
  real, intent(in), dimension(nmonth,maxgrid)   :: mfapar_field

  ! local variables
  integer :: dm, moy, jpngr, day

  ! ! XXX PMODEL_TEST
  ! write(0,*) 'WARNING: FAPAR = 1.00 USED IN PMODEL'

  !----------------------------------------------------------------
  ! INITIALISATIONS
  !----------------------------------------------------------------
  if (((interface%steering%init))) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    ! print*,'getting model parameters'
    call getpar_modl_waterbal()
    ! call getpar_modl_soil()

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    call ((interface%steering%init))global_waterbal()

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    call ((interface%steering%init))io_waterbal()
    call ((interface%steering%init))io_soiltemp()

  endif 

  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  call ((interface%steering%init))output_waterbal()
  call ((interface%steering%init))output_soiltemp()

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
    call getsolar_alldays( lat(jpngr), elv(jpngr), dfsun_field(:,jpngr) )

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
        ! ((interface%steering%init))ialise daily updated variables 
        !----------------------------------------------------------------
        call ((interface%steering%init))daily_waterbal()

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        ! write(0,*) 'calling waterbal() ... '
        call waterbal( jpngr, day, lat(jpngr), elv(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr) )
        ! write(0,*) '... done'

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        ! write(0,*) 'calling soiltemp() ... '
        call soiltemp( jpngr, moy, day, dtemp_field(:,jpngr), params_soil_field(jpngr) )
        ! write(0,*) '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        ! write(0,*) 'calling getout_daily_waterbal() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        ! write(0,*) '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    call writeout_ascii_waterbal( year )
    call writeout_ascii_soiltemp( year )

  end do

end subroutine biosphere

