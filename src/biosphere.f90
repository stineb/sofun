subroutine biosphere( &
  year, lon, lat, elv &
  , params_soil_field, lu_area, pco2 &
  , dtemp_field, dprec_field &
  , dfsun_field, dvpd_field, dndep_field &
  , c_uptake &
  , fapar_field &
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
  use md_params_modl
  use md_vars_core, only: ((interface%steering%init))annual, ((interface%steering%init))daily, ((interface%steering%init))global, ((interface%steering%init))pft, solar
  use md_outvars, only: ((interface%steering%init))io, getout_annual, getout_daily, ((interface%steering%init))output
  use md_soiltemp, only: soiltemp, ((interface%steering%init))output_soiltemp, ((interface%steering%init))io_soiltemp, getout_daily_soiltemp, writeout_ascii_soiltemp
  use md_phenology, only: gettempphenology, getpar_phenology  !, phenology
  use md_vegdynamics, only: estab_daily
  use md_gpp, only: gpp, getlue, ((interface%steering%init))daily_gpp, writeout_ascii_gpp, getout_annual_gpp, ((interface%steering%init))io_gpp, ((interface%steering%init))output_gpp
  use md_npp, only: npp, getpar_modl_npp
  use md_allocation, only: allocation_daily, ((interface%steering%init))io_allocation, ((interface%steering%init))output_allocation, getout_daily_allocation, writeout_ascii_allocation
  use md_turnover, only: turnover
  use md_waterbal, only: waterbal, getsolar_alldays, outdcpa, ((interface%steering%init))io_waterbal, getout_daily_waterbal, ((interface%steering%init))output_waterbal, getpar_modl_waterbal
  use md_nuptake, only: ((interface%steering%init))daily_nuptake, ((interface%steering%init))io_nuptake, getout_daily_nuptake, ((interface%steering%init))output_nuptake, getpar_modl_nuptake
  use md_ntransform, only: ntransform, ((interface%steering%init))_global_ntransform, ((interface%steering%init))io_ntransform, ((interface%steering%init))output_ntransform, getout_daily_ntransform, writeout_ascii_ntransform, getpar_modl_ntransform
  use md_littersom, only: littersom, ((interface%steering%init))io_littersom, ((interface%steering%init))output_littersom, getout_annual_littersom, getpar_modl_littersom !, getout_daily_littersom
  use md_params_soil, only: paramtype_soil

  use md_waterbal, only: writeout_ascii_waterbal

  ! xxx try
  ! use md_waterbal, only: getpar_modl_waterbal
  ! use md_nuptake, only: getpar_modl_nuptake
  ! use md_npp, only: getpar_modl_npp
  ! use md_phenology, only: getpar_phenology
  ! use md_ntransform, only: getpar_modl_ntransform
  ! use md_littersom, only: getpar_modl_littersom
  ! use md_sofunutils, only: getparreal


  ! use md_nuptake, only: writeout_ascii_nuptake
  ! use md_littersom, only: writeout_ascii_littersom

  ! xxx debug
  ! use md_vars_core
  ! use md_vegdynamics
  ! use md_vars_core

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
  real, intent(in), dimension(nmonth)           :: fapar_field

  ! local variables
  integer :: dm, mo, jpngr, day


  ! XXX uncomment this again XXX
  ! This is necessary because of keyword-specification of arguments to 'gpp'
  ! interface
  !   subroutine gpp( jpngr, ppfd, co2, tc, vpd, elv, dgpp_data )
                    
  !     ! optional set of arguments if 'dgpp' is to be predicted
  !     integer, intent(in), optional :: jpngr    ! gridcell number
  !     real, intent(in), optional    :: ppfd     ! photon flux density (xxx)
  !     real, intent(in), optional    :: co2      ! atmospheric CO2 (ppm)
  !     real, intent(in), optional    :: tc       ! air temperature (deg C)
  !     real, intent(in), optional    :: vpd      ! vapour pressure deficit (Pa)
  !     real, intent(in), optional    :: elv      ! elevation above sea level (m)

  !     ! optional set of arguments if 'dgpp' is to be prescribed
  !     real, dimension(npft), intent(in), optional :: dgpp_data
    
  !   end subroutine gpp
  ! end interface


  ! !----------------------------------------------------------------
  ! ! Declare interface here for subroutines with optional arguments
  ! !----------------------------------------------------------------
  ! interface
  !   subroutine gpp( jpngr, doy, moy, dgpp_data )
  !     integer, intent(in), optional               :: jpngr     ! gridcell number
  !     integer, intent(in), optional               :: doy       ! day of year and month of year
  !     integer, intent(in), optional               :: moy       ! month of year and month of year
  !     real, dimension(npft), intent(in), optional :: dgpp_data ! prescribed daily GPP for this gridcell
  !   end subroutine gpp
  ! end interface

  ! if (present(fapar).and..not.fapar(1)==dummy) write(0,*) 'süpi: fapar ', fapar
  ! stop


  !----------------------------------------------------------------
  ! INITIALISATIONS
  !----------------------------------------------------------------
  if (((interface%steering%init))) then

    !----------------------------------------------------------------
    ! GET MODEL PARAMETERS
    ! read model parameters that may be varied for optimisation
    !----------------------------------------------------------------
    ! print*,'getting model parameters'
    call getpar_modl()
    call getpar_modl_waterbal()
    call getpar_modl_npp()
    call getpar_modl_nuptake()
    call getpar_modl_ntransform()
    call getpar_modl_littersom()
    call getpar_phenology()

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    call ((interface%steering%init))global()
    call ((interface%steering%init))_global_ntransform()

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    call ((interface%steering%init))io()
    call ((interface%steering%init))io_waterbal()
    call ((interface%steering%init))io_soiltemp()
    call ((interface%steering%init))io_gpp()
    call ((interface%steering%init))io_nuptake()
    call ((interface%steering%init))io_ntransform()
    call ((interface%steering%init))io_littersom()
    call ((interface%steering%init))io_allocation()

  endif 


  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  call ((interface%steering%init))output()
  call ((interface%steering%init))output_waterbal()
  call ((interface%steering%init))output_soiltemp()
  call ((interface%steering%init))output_gpp()
  call ((interface%steering%init))output_nuptake()
  call ((interface%steering%init))output_ntransform()
  call ((interface%steering%init))output_littersom()
  call ((interface%steering%init))output_allocation()

  !----------------------------------------------------------------
  ! LOOP THROUGH GRIDCELLS
  !----------------------------------------------------------------
  do jpngr=1,maxgrid

    !----------------------------------------------------------------
    ! ((interface%steering%init))ialise annually updated variables
    !----------------------------------------------------------------
    call ((interface%steering%init))annual

    !call snow

    !----------------------------------------------------------------
    ! get temperature-driven phenology (drought-driven phenology is 
    ! calculated after waterbalance)
    !----------------------------------------------------------------
    call gettempphenology( jpngr, dtemp_field(:,jpngr) )

    !call climate20
    !call bioclim -> SR contained in 'establishment'
    !call conversion( lu_area )

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

    ! !----------------------------------------------------------------
    ! ! Get leaf traits based on photosynthetic parameters (Vcmax25) 
    ! ! for this year given this year's climate (temp., VPD) and CO2.
    ! !----------------------------------------------------------------
    ! call gettraits( jpngr )

    ! !----------------------------------------------------------------
    ! ! sapling addition on bare ground, tree geometry update
    ! !----------------------------------------------------------------
    ! call establishment( jpngr, fpc_grid_data(:,jpngr) )

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    day=0
    do mo=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
      do dm=1,ndaymonth(mo)
        day=day+1

        write(0,*) '-----------'
        write(0,*) 'DAY ',day

        !----------------------------------------------------------------
        ! ((interface%steering%init))ialise daily updated variables 
        !----------------------------------------------------------------
        call ((interface%steering%init))daily()
        call ((interface%steering%init))daily_nuptake()
        call ((interface%steering%init))daily_gpp()

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !call waterbal( day, lon(jpngr), lat(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr), elv(jpngr) )
        ! xxx try
        !----------------------------------------------------------------
        write(0,*) 'calling waterbal() ... '
        call waterbal( jpngr, day, lat(jpngr), elv(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr) )
        write(0,*) '... done'

        !----------------------------------------------------------------
        ! simulate daily establishment / sprouting
        !----------------------------------------------------------------
        write(0,*) 'calling estab_daily() ... '
        call estab_daily( jpngr, day ) 
        write(0,*) '... done'
        ! XXX this is no longer executed like this. Instead, provide a certain C and N free of charge in allocation XX

        !----------------------------------------------------------------
        ! calculate GPP
        !----------------------------------------------------------------
        write(0,*) 'calling gpp() ... '
        call gpp( jpngr, day, mo, fapar_field(mo) )
        write(0,*) '... done'

        ! xxx        dgpp creates FPE      xxx

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        write(0,*) 'calling soiltemp() ... '
        write(0,*) 'with arguments ... '
        write(0,*) 'dtemp_field(:,jpngr) =', dtemp_field(:,jpngr)
        write(0,*) 'outdcpa(:,:,jpngr)   =', outdcpa(:,:,jpngr)
        call soiltemp( jpngr, mo, day, dtemp_field(:,jpngr) )
        write(0,*) '... done'

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        write(0,*) 'calling npp() ... '
        call npp( jpngr, dtemp_field(day,jpngr), day )
        write(0,*) '... done'

        ! write(0,*) 'A pninorg(lu,jpngr)%n14',pninorg
        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        write(0,*) 'calling allocation() ... '
        call allocation_daily( jpngr, day, mo, dm )
        ! write(0,*) 'B pninorg(lu,jpngr)%n14',pninorg
        write(0,*) '... done'

        ! amount of NPP added to reproduction (xxx ignore at this point)
        !call reproduction( jpngr )

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        write(0,*) 'calling turnover() ... '
        call turnover( jpngr, day )
        ! write(0,*) 'C pninorg(lu,jpngr)%n14',pninorg
        write(0,*) '... done'

        ! !----------------------------------------------------------------
        ! ! light competition
        ! !----------------------------------------------------------------
        ! call light( jpngr )

        !----------------------------------------------------------------
        ! litter and soil decomposition and N mineralisation
        !----------------------------------------------------------------
        write(0,*) 'calling littersom() ... '
        call littersom( jpngr, day )
        write(0,*) '... done'

        ! write(0,*) 'D pninorg(lu,jpngr)%n14',pninorg
        !----------------------------------------------------------------
        ! inorganic soil N dynamics
        !----------------------------------------------------------------
        write(0,*) 'calling ntransform() ... '
        ! ! xxx try
        ! dndep_field(day,jpngr) = dndep_field(day,jpngr) * max( 1.0, 100.0 - real( realyear ) )
        call ntransform( dm, mo, jpngr, dndep_field(day,jpngr), sum(dprec_field(:,jpngr)) )
        ! write(0,*) 'E pninorg(lu,jpngr)%n14',pninorg
        write(0,*) '... done'

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        write(0,*) 'calling getout_daily() ... '
        call getout_daily( jpngr, mo, day )
        write(0,*) 'calling getout_daily_waterbal() ... '
        call getout_daily_waterbal( jpngr, mo, day )
        write(0,*) 'calling getout_daily_nuptake() ... '
        call getout_daily_nuptake( jpngr, mo, day )
        write(0,*) 'calling getout_daily_ntransform() ... '
        ! call getout_daily_littersom( jpngr, mo, day )
        call getout_daily_ntransform( jpngr, mo, day )
        write(0,*) 'calling getout_daily_allocation() ... '
        call getout_daily_allocation( jpngr, mo, day )
        write(0,*) '... done'
        call getout_daily_soiltemp( jpngr, mo, day )

        ! stop 'end of day'
        
      end do

    end do

    ! !----------------------------------------------------------------
    ! ! allocate biomass increment
    ! ! precompiler flag used here because option requires allocation to be called daily
    ! !----------------------------------------------------------------
    ! call allocation_annual( jpngr )

    ! reduce part of NPP for reproduction
    !call reproduction

    ! light competition by limiting to maximum FPC
    ! call light( jpngr )
    ! print*,'fpc_grid', fpc_grid(1,1)
    !call mortality
    !call fire

    !----------------------------------------------------------------
    ! collect annually updated output variables
    !----------------------------------------------------------------
    call getout_annual( jpngr )
    call getout_annual_gpp( jpngr )
    call getout_annual_littersom( jpngr )

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    call writeout_ascii( year, dtemp_field(:,jpngr) )
    call writeout_ascii_gpp( year, spinup )
    call writeout_ascii_waterbal( year, spinup )
    call writeout_ascii_allocation( year, spinup )
    call writeout_ascii_ntransform( year, spinup )
    call writeout_ascii_soiltemp( year, spinup )
    
    ! XXX Mysterious problem: when called from here, model does weird things (check e.g. gpp and soil C)
    ! call writeout_ascii_nuptake( year, spinup )
    ! call writeout_ascii_littersom( year, spinup )

  end do

  ! try xxx
  c_uptake = 1.e15 ! in gC


  ! write output
  !call write_output

end subroutine biosphere

