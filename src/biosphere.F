subroutine biosphere( year, realyear, lon, lat, elv, lu_area, pco2 &
  , dtemp_field, dprec_field, dfsun_field, dvpd_field, dndep_field &
  , c_uptake &
  , dgpp_data &
  , fpc_grid_data &
  )
  !////////////////////////////////////////////////////////////////
  !  Subroutine BIOSPHERE calculates net ecosystem exchange (nee)
  !  in response to environmental boundary conditions (atmospheric 
  !  CO2, temperature, Nitrogen deposition. This SR "replaces" 
  !  LPJ, also formulated as subroutine.
  !----------------------------------------------------------------
  use _params_core
  use _params_siml
  use _params_site
  use _params_modl
  use _vars_core, only: initannual, initdaily, initglobal, initpft
  use _outvars, only: initio, getout_annual, getout_daily, initoutput
  use _soiltemp, only: soiltemp
  use _phenology, only: summerphenology
  use _vegdynamics, only: estab_daily
  use _gpp, only: gpp, getlue, gettraits, initdaily_gpp
  use _npp, only: npp
  use _allocation, only: allocation_daily
  use _turnover, only: turnover
  use _waterbal, only: waterbal, getrad, outdcpa, initdaily_waterbal, initio_waterbal, getout_daily_waterbal, initoutput_waterbal, getout_monthly_waterbal
  use _nuptake, only: initdaily_nuptake, initio_nuptake, getout_daily_nuptake, initoutput_nuptake
  use _ntransform, only: ntransform, init_global_ntransform, initio_ntransform, initoutput_ntransform
  use _littersom, only: littersom, initio_littersom, getout_daily_littersom, initoutput_littersom, getout_annual_littersom

  use _waterbal, only: writeout_ascii_waterbal

  ! xxx debug
  ! use _vars_core
  ! use _vegdynamics
  ! use _vars_core

  implicit none

  ! arguments
  integer, intent(in)                           :: year       ! simulation year
  integer, intent(in)                           :: realyear   ! year AD (=first simulation year )
  real, intent(in), dimension(maxgrid)          :: lon        ! longitude vector/field (degrees E)              
  real, intent(in), dimension(maxgrid)          :: lat        ! latitude vector/field (degrees N)             
  real, intent(in), dimension(maxgrid)          :: elv        ! elevation (altitude) vector/field (m above sea level)                  
  real, dimension(3)                            :: lu_area    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
  real, intent(in)                              :: pco2
  real, intent(in), dimension(ndayyear,maxgrid) :: dtemp_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dprec_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dfsun_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dvpd_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dndep_field
  real, intent(out)                             :: c_uptake   ! annual net global C uptake by biosphere

  ! optional arguments
  real, intent(in), dimension(ndayyear,npft), optional :: dgpp_data       ! for runs with prescribed GPP
  real, intent(in), dimension(nlu,maxgrid), optional   :: fpc_grid_data   ! for runs with prescribed GPP

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

  !----------------------------------------------------------------
  ! INITIALISATIONS
  !----------------------------------------------------------------
  if (init) then

    ! xxx check: if (present(fpc_grid_data)) .and. _dynveg then abort

    !----------------------------------------------------------------
    ! Initialise pool variables and/or read from restart file (not implemented)
    !----------------------------------------------------------------
    call initglobal()
    call init_global_ntransform()

    !----------------------------------------------------------------
    ! Open input/output files
    !----------------------------------------------------------------
    call initio()
    call initio_waterbal()
    call initio_nuptake()
    call initio_ntransform() ! nothing yet xxx
    call initio_littersom()

  endif 


  !----------------------------------------------------------------
  ! Initialise output variables for this year
  !----------------------------------------------------------------
  call initoutput()
  call initoutput_waterbal()
  call initoutput_nuptake()
  call initoutput_ntransform()  ! nothing yet
  call initoutput_littersom()

  !----------------------------------------------------------------
  ! LOOP THROUGH GRIDCELLS
  !----------------------------------------------------------------
  do jpngr=1,maxgrid

    !----------------------------------------------------------------
    ! initialise annually updated variables
    !----------------------------------------------------------------
    call initannual

    !call snow

    !----------------------------------------------------------------
    ! get temperature-driven phenology (drought-driven phenology is calculated after waterbalance)
    !----------------------------------------------------------------
    call summerphenology(jpngr, dtemp_field(:,jpngr))

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
    ! Get radiation, evaporation, and soil moisture based on daily precipitation,
    ! temperature, sunshine fraction, and elevation (previously known!)
    ! This is not compatible with a daily biosphere-climate coupling. I.e., 
    ! there is a daily loop within 'getrad'!
    !----------------------------------------------------------------
    call getrad( lat(jpngr), elv(jpngr), dfsun_field(:,jpngr) )

    !----------------------------------------------------------------
    ! Get leaf traits based on photosynthetic parameters (Vcmax25) 
    ! for this year given this year's climate (temp., VPD) and CO2.
    !----------------------------------------------------------------
    call gettraits( jpngr )

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

        ! write(0,*) '-----------'
        ! write(0,*) 'DAY ',day

        !----------------------------------------------------------------
        ! initialise daily updated variables 
        !----------------------------------------------------------------
        call initdaily()
        call initdaily_nuptake()
        call initdaily_waterbal()
        call initdaily_gpp()

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !call waterbal( day, lon(jpngr), lat(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr), elv(jpngr) )
        ! xxx try
        !----------------------------------------------------------------
        call waterbal( jpngr, day, lat(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr), elv(jpngr) )

        !----------------------------------------------------------------
        ! simulate daily establishment / sprouting
        !----------------------------------------------------------------
        call estab_daily( jpngr, day, fpc_grid_data(:,jpngr) )

        !----------------------------------------------------------------
        ! calculate GPP
        !----------------------------------------------------------------
        if ( present(dgpp_data) ) then
          call gpp( jpngr, day, mo, dgpp_data(day,:) )
        else 
          call gpp( jpngr, day, mo )
        end if   

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        call soiltemp( jpngr, mo, day, dtemp_field(:,jpngr), outdcpa(:,:,jpngr) )

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        call npp( jpngr, dtemp_field(day,jpngr) )

        ! write(0,*) 'A pninorg(lu,jpngr)%n14',pninorg
        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        call allocation_daily( jpngr, day, mo )
        ! write(0,*) 'B pninorg(lu,jpngr)%n14',pninorg


        ! amount of NPP added to reproduction (xxx ignore at this point)
        !call reproduction( jpngr )

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        call turnover( jpngr )
        ! write(0,*) 'C pninorg(lu,jpngr)%n14',pninorg

        !----------------------------------------------------------------
        ! litter and soil decomposition and N mineralisation
        !----------------------------------------------------------------
        call littersom( jpngr, day )

        ! write(0,*) 'D pninorg(lu,jpngr)%n14',pninorg
        !----------------------------------------------------------------
        ! inorganic soil N dynamics
        !----------------------------------------------------------------
        call ntransform( dm, mo, jpngr, sum(dprec_field(:,jpngr)) )
        ! write(0,*) 'E pninorg(lu,jpngr)%n14',pninorg

        !----------------------------------------------------------------
        ! collect from daily updated state variables for annual variables
        !----------------------------------------------------------------
        call getout_daily( jpngr, mo, day )
        call getout_daily_waterbal( jpngr, mo, day )
        call getout_daily_nuptake( jpngr, mo, day )
        call getout_daily_littersom( jpngr, mo, day )

        ! stop 'end of day'
        
      end do

      !----------------------------------------------------------------
      ! MODULE-SPECIFIC VARIABLES
      !----------------------------------------------------------------
      call getout_monthly_waterbal( jpngr, mo )

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
    call getout_annual_littersom( jpngr )

    !----------------------------------------------------------------
    ! Write to output
    !----------------------------------------------------------------
    call writeout_ascii( year, dtemp_field(:,jpngr) )
    call writeout_ascii_waterbal( year, spinup )

  end do

  ! try xxx
  c_uptake = 1.e15 ! in gC


  ! write output
  !call write_output

end subroutine biosphere

