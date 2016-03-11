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
  use _params_core
  use _params_siml
  use _params_site
  use _params_modl
  use _vars_core, only: initannual, initdaily, initglobal, initpft
  use _soiltemp, only: soiltemp, initoutput_soiltemp, initio_soiltemp, getout_daily_soiltemp, writeout_ascii_soiltemp
  use _params_soil, only: paramtype_soil
  use _waterbal, only: waterbal, getsolar_alldays, initdaily_waterbal, initglobal_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use _gpp, only: getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, getout_daily_gpp, writeout_ascii_gpp

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
  if (init) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    ! print*,'getting model parameters'
    call getpar_modl()
    call getpar_modl_waterbal()
    call getpar_modl_gpp()

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    call initglobal()
    call initglobal_waterbal()

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    call initio_waterbal()
    call initio_soiltemp()
    call initio_gpp()

  endif 

  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  call initoutput_waterbal()
  call initoutput_soiltemp()
  call initoutput_gpp()

  !----------------------------------------------------------------
  ! LOOP THROUGH GRIDCELLS
  !----------------------------------------------------------------
  do jpngr=1,maxgrid

    !----------------------------------------------------------------
    ! initialise annually updated variables
    !----------------------------------------------------------------
    call initannual()

    !----------------------------------------------------------------
    ! Get monthly light use efficiency, and Rd per unit of light absorbed
    ! Photosynthetic parameters acclimate at monthly time scale
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a monthly loop within 'getlue'!
    !----------------------------------------------------------------
    call getlue( jpngr, pco2, dtemp_field(:,jpngr), dvpd_field(:,jpngr), elv(jpngr) )

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
        ! initialise daily updated variables 
        !----------------------------------------------------------------
        call initdaily()
        call initdaily_waterbal()
        call initdaily_gpp()

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
        ! calculate GPP
        !----------------------------------------------------------------
        ! write(0,*) 'calling gpp() ... '
        call gpp( jpngr, day, moy, mfapar_field(moy,jpngr) )
        ! call gpp( jpngr, day, moy, 1.00 )
        ! write(0,*) '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        ! write(0,*) 'calling getout_daily_waterbal() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        call getout_daily_gpp( jpngr, moy, day )
        ! write(0,*) '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    call writeout_ascii_waterbal( year, spinup )
    call writeout_ascii_soiltemp( year, spinup )
    call writeout_ascii_gpp( year, spinup )

  end do

end subroutine biosphere

