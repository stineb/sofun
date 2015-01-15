subroutine biosphere( year, realyear, pco2, dtemp, dprec, dndep, landuse, c_uptake )
  !////////////////////////////////////////////////////////////////
  !  Subroutine BIOSPHERE calculates net ecosystem exchange (nee)
  !  in response to environmental boundary conditions (atmospheric 
  !  CO2, temperature, Nitrogen deposition. This SR "replaces" 
  !  LPJ, also formulated as subroutine.
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use params_core
  use params_siml
  use params_site
  use params_modl
  use waterbal_stash
  
  ! xxx debug
  use fluxes
  use treegeometry
  use pools

  implicit none

  ! arguments
  integer, intent(in)                   :: year       ! simulation year
  integer, intent(in)                   :: realyear   ! year AD (=first simulation year )
  real, intent(in), dimension(ndayyear) :: dtemp
  real, intent(in), dimension(ndayyear) :: dprec
  real, intent(in), dimension(ndayyear) :: dndep
  real, intent(in)                      :: pco2
  real, dimension(3)                    :: landuse    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
  real, intent(out)                     :: c_uptake   ! annual net global C uptake by biosphere

  ! local variables
  integer :: dm, mo, jpngr, day


  if (init) then
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
    !call summerphenology (drough-driven phenology is calculated after waterbalance)
    !call climate20
    !call bioclim -> SR contained in 'establishment'
    !call conversion( lu_area )

    ! sapling addition on bare ground, tree geometry update
    call establishment( jpngr )


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
        call waterbal( lon, lat, dprec, dtemp, (1.0-dccov), elv )

        ! get GPP from input
        call gpp( jpngr, pco2, year, realyear, mo, dm, day )

        ! calculate water balance
        !call waterbalance

        ! use SR 'sitchtemp'
        !call soiltemp

        ! substract autotrophic respiration to get NPP, remainder is added to
        ! labile pool (plabl)
        ! xxx try: second argument is dtemp_soil
        call npp( jpngr, dtemp(day), 1.0 )

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
        call littersom( jpngr )
        !print*,'F'
        !print*,'nh4 ',nh4(1,1)%n14 


        ! inorganic soil N chemistry
        call ntransform( dm, mo, jpngr, dprec )

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

#if _fix_veg == 0
    ! light competition by limiting to maximum FPC
    call light( jpngr )
    print*,'fpc_grid', fpc_grid(1,1)
    !call mortality
    !call fire
#endif

    ! collect annually updated output variables
    call getout_annual( jpngr )

    ! Write to output
    call writeout_ascii( year )

  end do

  ! try xxx
  c_uptake = 1.e15 ! in gC

  ! write output
  !call write_output

end subroutine biosphere

