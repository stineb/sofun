subroutine biosphere( year, realyear, lon, lat, elv, lu_area, pco2 &
  , dtemp_field, dprec_field, dfsun_field, dndep_field, c_uptake &
  , dgpp_data &
  , fpc_grid_data &
  )
  !////////////////////////////////////////////////////////////////
  !  Subroutine BIOSPHERE calculates net ecosystem exchange (nee)
  !  in response to environmental boundary conditions (atmospheric 
  !  CO2, temperature, Nitrogen deposition. This SR "replaces" 
  !  LPJ, also formulated as subroutine.
  !----------------------------------------------------------------
  use params_core
  use params_siml
  use params_site
  use params_modl
  use waterbal_stash, only: waterbal
  use fixveg
  use soiltemp_sitch, only: soiltemp
  use waterbal_stash, only: outdcpa ! used for soiltemp xxx make this global?
  use phenology_lpx, only: summerphenology
  use ntransform_xuri, only: ntransform

  ! xxx debug
  use fluxes
  use treegeometry
  use pools

  implicit none

  ! arguments
  integer, intent(in)                   :: year       ! simulation year
  integer, intent(in)                   :: realyear   ! year AD (=first simulation year )
  real, intent(in), dimension(maxgrid)  :: lon        ! longitude vector/field (degrees E)              
  real, intent(in), dimension(maxgrid)  :: lat        ! latitude vector/field (degrees N)             
  real, intent(in), dimension(maxgrid)  :: elv        ! elevation (altitude) vector/field (m above sea level)                  
  real, intent(in)                      :: pco2
  real, intent(in), dimension(ndayyear,maxgrid) :: dtemp_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dprec_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dfsun_field
  real, intent(in), dimension(ndayyear,maxgrid) :: dndep_field
  real, dimension(3)                    :: lu_area    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
  real, intent(out)                     :: c_uptake   ! annual net global C uptake by biosphere

  ! optional arguments
  real, intent(in), dimension(ndayyear,npft), optional :: dgpp_data       ! for runs with prescribed GPP
  real, intent(in), dimension(nlu,maxgrid), optional   :: fpc_grid_data   ! for runs with prescribed GPP

  ! local variables
  integer :: dm, mo, jpngr, day


  ! This is necessary because of keyword-specification of arguments to 'gpp'
  interface
    subroutine gpp( jpngr, pco2, dgpp_data )
      
      ! optional set of arguments if 'dgpp' is to be predicted
      integer, intent(in), optional :: jpngr                  ! gridcell number
      real, intent(in), optional    :: pco2                    ! atmospheric CO2 concentration [ppm]

      ! optional set of arguments if 'dgpp' is to be prescribed
      real, dimension(npft), intent(in), optional :: dgpp_data
    
    end subroutine gpp
  end interface


  if (init) then

    ! xxx check: if (present(fpc_grid_data)) .and. dynveg then abort

    ! First simulation year
    ! xxx debug
    !print*,'A'
    !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
    !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
    !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

    ! Initialise pool variables and/or read from restart file (not implemented)
    call initglobal

    !print*,'B'
    !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
    !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
    !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

    ! Open input/output files
    call initio(runname)

  endif 

  !print*,'C'
  !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
  !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
  !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

  ! Initialise output variables for this year
  call initoutput

  !print*,'D'
  !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
  !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
  !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

  ! LOOP THROUGH GRIDCELLS
  do jpngr=1,maxgrid

    ! initialise annually updated variables
    call initannual

    !call snow

    ! get temperature-driven phenology (drought-driven phenology is calculated after waterbalance)
    call summerphenology(jpngr, dtemp_field(:,jpngr))

    !call climate20
    !call bioclim -> SR contained in 'establishment'
    !call conversion( lu_area )

    ! sapling addition on bare ground, tree geometry update
    call establishment( jpngr, fpc_grid_data(:,jpngr) )


    ! LOOP THROUGH MONTHS
    day=0
    do mo=1,nmonth

      ! LOOP THROUGH DAYS
      do dm=1,ndaymonth(mo)
        day=day+1

        ! initialise daily updated variables 
        call initdaily

        !print*,'E'
        !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
        !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
        !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12
        !print*,'nh4 ',nh4(1,1)%n14 

        ! get soil moisture, and runoff
        !call waterbal( day, lon(jpngr), lat(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr), elv(jpngr) )
        ! xxx try
        call waterbal( jpngr, day, lon(jpngr), lat(jpngr), dprec_field(day,jpngr), dtemp_field(day,jpngr), dfsun_field(day,jpngr), elv(jpngr) )

        ! get GPP from input
        if (present(dgpp_data)) then
          call gpp( dgpp_data=dgpp_data(day,:) )
        else 
          call gpp( jpngr=jpngr, pco2=pco2 )
        end if   

        ! use SR 'sitchtemp'
        call soiltemp( jpngr, mo, day, dtemp_field(:,jpngr), outdcpa(:,:,jpngr) )

        ! substract autotrophic respiration to get NPP, remainder is added to
        ! labile pool (plabl)
        ! xxx try: second argument is dtemp_soil
        call npp( jpngr, dtemp_field(day,jpngr) )

        !print*,'F'
        !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
        !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
        !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 


        !#if _allocation_pipe
        ! allocate biomass increment
        ! precompiler flag used here because option requires allocation to be called daily
        !call allocation_daily( jpngr )
        !#endif

        ! amount of NPP added to reproduction (xxx ignore at this point)
        !call reproduction( jpngr )

        ! leaf, sapwood, and fine-root turnover
        call turnover( jpngr )

        ! litter and soil turnover
        call littersom( jpngr, day )
        !print*,'F'
        !print*,'nh4 ',nh4(1,1)%n14 


        ! inorganic soil N chemistry
        call ntransform( dm, mo, jpngr, sum(dprec_field(:,jpngr)) )

        !print*,'G'
        !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
        !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
        !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 
        !print*,'nh4 ',nh4(1,1)%n14 

        ! collect from daily updated state variables for annual variables
        call getout_daily( jpngr, mo, day )

        !print*,'H'
        !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
        !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
        !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

      end do

    end do

    !print*,'I'
    !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
    !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
    !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 


    ! allocate biomass increment
    ! precompiler flag used here because option requires allocation to be called daily
    call allocation_annual( jpngr )

    !print*,'K'
    !print*,'plabl(1,1)%c%c12 ',plabl(1,1)%c%c12 
    !print*,'pleaf(1,1)%c%c12 ',pleaf(1,1)%c%c12
    !print*,'proot(1,1)%c%c12 ',proot(1,1)%c%c12 

    ! reduce part of NPP for reproduction
    !call reproduction

    ! light competition by limiting to maximum FPC
    ! call light( jpngr )
    !print*,'fpc_grid', fpc_grid(1,1)
    !call mortality
    !call fire

    ! collect annually updated output variables
    call getout_annual( jpngr )

    ! Write to output
    call writeout_ascii( year, dtemp_field(:,jpngr) )

  end do

  ! try xxx
  c_uptake = 1.e15 ! in gC


  ! write output
  !call write_output

end subroutine biosphere

