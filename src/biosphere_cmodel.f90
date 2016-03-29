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
  use md_plant, only: initdaily_plant, initglobal_plant, initpft, getout_daily_plant, getpar_modl_plant, initoutput_plant, writeout_ascii_plant, initio_plant, getout_annual_plant
  use md_soiltemp, only: soiltemp, initoutput_soiltemp, initio_soiltemp, getout_daily_soiltemp, writeout_ascii_soiltemp
  use md_params_soil, only: paramtype_soil
  use md_waterbal, only: waterbal, getsolar_alldays, initdaily_waterbal, initglobal_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, writeout_ascii_waterbal
  use md_phenology, only: gettempphenology, getpar_modl_phenology
  use md_gpp, only: getpar_modl_gpp, initio_gpp, initoutput_gpp, initdaily_gpp, getlue, gpp, getout_daily_gpp, writeout_ascii_gpp
  use md_npp, only: npp
  use md_turnover, only: turnover
  use md_allocation, only: allocation_daily
  use md_vegdynamics, only: vegdynamics
  
  ! xxx debug
  use md_plant, only: ispresent, canopy

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
  ! print*, 'WARNING: FAPAR = 1.00 USED IN PMODEL'

  !----------------------------------------------------------------
  ! INITIALISATIONS
  !----------------------------------------------------------------
  if (init) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    ! print*,'getting model parameters ...'
    call getpar_modl_plant()
    call getpar_modl_waterbal()
    call getpar_modl_gpp()
    call getpar_modl_phenology()
    ! call getpar_modl_soil()
    ! print*,'... done'

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    ! print*,'initialising variables ...'
    call initglobal_plant()
    call initglobal_waterbal()
    ! call initglobal_soil()
    ! print*,'... done'

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    ! print*,'initialising IO ...'
    call initio_waterbal()
    call initio_soiltemp()
    call initio_gpp()
    call initio_plant()
    ! call initio_soil()
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
  ! call initoutput_soil()
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
    call getsolar_alldays( lat(jpngr), elv(jpngr), dfsun_field(:,jpngr) )

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
    call gettempphenology( jpngr, dtemp_field(:,jpngr) )

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

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        ! print*, 'calling waterbal() ... '
        call waterbal( jpngr, day, lat(jpngr), elv(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr) )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        ! print*, 'calling soiltemp() ... '
        call soiltemp( jpngr, moy, day, dtemp_field(:,jpngr), params_soil_field(jpngr) )
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
        call gpp( jpngr, day, moy, dtemp_field(day,jpngr), mfapar_field(moy,jpngr) )
        ! call gpp( jpngr, day, moy, 1.00 )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        ! print*, 'calling npp() ... '
        call npp( jpngr, dtemp_field(day,jpngr), day )
        ! print*, '... done'

        ! amount of NPP added to reproduction (xxx ignore at this point)
        !call reproduction( jpngr )

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        ! print*, 'calling turnover() ... '
        call turnover( jpngr, day )
        ! print*, '... done'

        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        ! print*, 'calling allocation() ... '
        call allocation_daily( jpngr, day, moy, dm )
        ! print*, '... done'

        ! !----------------------------------------------------------------
        ! ! litter and soil decomposition and N mineralisation
        ! !----------------------------------------------------------------
        ! ! print*, 'calling littersom() ... '
        ! call littersom( jpngr, day )
        ! ! print*, '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        ! print*, 'calling getout_daily_*() ... '
        call getout_daily_waterbal( jpngr, moy, day )
        call getout_daily_soiltemp( jpngr, moy, day )
        call getout_daily_gpp( jpngr, moy, day )
        call getout_daily_plant( jpngr, moy, day )
        ! call getout_daily_soil( jpngr, moy, day )
        ! print*, '... done'

      end do

    end do

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    call writeout_ascii_waterbal( year )
    call writeout_ascii_soiltemp( year )
    call writeout_ascii_gpp( year )
    call writeout_ascii_plant( year )
    ! call writeout_ascii_soil( year )

  end do

end subroutine biosphere

