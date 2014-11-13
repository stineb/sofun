      program main
      !////////////////////////////////////////////////////////////////
      !  Main program for site scale simulations, here used for 
      !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
      !  Nitrogen)
      !----------------------------------------------------------------
#include "sofun_module_control.inc"
      use params_core
      use params_siml
      use params_site
      use params_modl
      use forcing
      use input_site
      use gridvars
      use statevars

      implicit none

      integer :: year
      integer :: realyear

      real, dimension(ndayyear) :: dtemp
      real, dimension(ndayyear) :: dprec
      real, dimension(ndayyear) :: dndep
      real                      :: pco2
      real, dimension(3)        :: landuse      ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'biosphere'
      real                      :: c_uptake     ! annual net global C uptake by biosphere


      !----------------------------------------------------------------
      ! READ RUNNAME FROM STANDARD INPUT
      !----------------------------------------------------------------
      read (*,*) runname
      ! make sure runname length is smaller/equal than maxlen_runname
      if (len_trim(runname)>=maxlen_runname) then
        stop'runname too long'
      endif

      ! write simulation name to standard output (screen)
      write(0,*) '------------SOFUN : '//trim(runname)//'-------------'
!       write (*,200)
!200    format ('SIMULATION NAME: '$)
!      write (*,*) runname


      !----------------------------------------------------------------
      ! GET SIMULATION PARAMETERS FROM FILE <runname>.sofun.parameter
      ! SR getpar_siml is defined in params_siml.mod.F
      !----------------------------------------------------------------
      call getpar_siml(trim(runname))

      !----------------------------------------------------------------
      ! GET SITE PARAMETERS
      ! site location (lon,lat), soil type, vegetation type
      ! SR getpar_site is defined in params_site.mod.F. 
      ! 'sitename' is global variable
      !----------------------------------------------------------------
      call getpar_site


      ! LOOP THROUGH YEARS
      write(0,*) '------------START OF SIMULATION-------------'
      do year=1,runyears

        !----------------------------------------------------------------
        ! Define simulations "steering" variables (realyear, etc.)
        !----------------------------------------------------------------
        call getsteering(year,realyear)

        if (year == spinupyears+1 ) then
          write(0,*) '-----------TRANSIENT SIMULATION-------------'
        endif

        !----------------------------------------------------------------
        ! Get external (environmental) forcing
        !----------------------------------------------------------------
        call getco2( trim(runname), realyear, pco2 )
        call getclimate( trim(runname), realyear, dtemp, dprec )
        call getndep( trim(runname), realyear, dprec, dndep )
        call getlanduse( trim(runname), realyear, landuse )

        !----------------------------------------------------------------
        ! Call SR biosphere at an annual time step but with vectors 
        ! containing data for each day of this year.
        !----------------------------------------------------------------
        call biosphere( year, realyear, pco2, dtemp, dprec, dndep, landuse, c_uptake ) 

        write(0,100) 'year AD, pco2',realyear,pco2

100   format(A,I6,F8.2,F8.2)

      enddo

      write(0,*) '--------------END OF SIMULATION---------------'

      contains
      !////////////////////////////////////////////////////////////////
      !  Subroutine BIOSPHERE calculates net ecosystem exchange (nee)
      !  in response to environmental boundary conditions (atmospheric 
      !  CO2, temperature, Nitrogen deposition. This SR "replaces" 
      !  LPJ, also formulated as subroutine.
      !----------------------------------------------------------------

      subroutine biosphere( year, realyear, pco2, dtemp, dprec, dndep, landuse, c_uptake )

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
         
        ! Initialise pool variables and/or read from restart file (not implemented)
        call initglobal

        ! Open input/output files
        call initio(runname)

        ! Get soil parameters (if not defined in <sitename>.parameter)
        !call getsoilpar

        !----------------------------------------------------------------
        ! GET MODEL PARAMETERS
        ! read model parameters that may varied for optimisation
        !----------------------------------------------------------------
        call getpar_modl

        ! Obtain land unit dependent parameters, define decomposition rates
        !call luparameters

        !----------------------------------------------------------------
        ! GET MODEL INPUT DATA
        ! Prescribe variables instead of calculating these parts online
        ! This reads daily values for each day of the transient simulation
        !----------------------------------------------------------------
        call getinput_site

      endif

      ! Initialise output variables
      call initoutput


      ! LOOP THROUGH GRIDCELLS
      do jpngr=1,maxgrid

        !call snow
        !call summerphenology (drough-driven phenology is calculated after waterbalance)
        !call climate20
        !call bioclim -> SR contained in 'establishment'
        !call conversion( lu_area )

        ! LOOP THROUGH MONTHS
        day=0
        do mo=1,nmonth

          ! LOOP THROUGH DAYS
          do dm=1,ndaymonth(mo)
            day=day+1

            ! initialise daily updated variables 
            call initdaily

            ! get GPP from input
            call gpp( pco2, year, realyear, mo, dm, day )

            ! calculate water balance
            !call waterbalance

            ! use SR 'sitchtemp'
            !call soiltemp

            ! substract autotrophic respiration to get NPP, remainder is added to
            ! labile pool (plabl)
            ! xxx try: second argument is dtemp_soil
            call npp( jpngr, dtemp(day), 1.0 )

#if _allocation_pipe
            ! allocate biomass increment
            ! precompiler flag used here because option requires allocation to be called daily
            call allocation_daily( jpngr )
#endif

            ! amount of NPP added to reproduction (xxx ignore at this point)
            !call reproduction( jpngr )

            ! leaf, sapwood, and fine-root turnover
            call turnover( jpngr )

            ! litter and soil turnover
            call littersom( jpngr )

            ! inorganic soil N chemistry
            call ntransform( dm, mo, jpngr, dprec )

            ! collect from daily updated state variables for annual variables
            call getout_daily( jpngr )

            ! write to daily output
            !call writeout_daily( year, day ) 

          end do

        end do

#if _allocation_pipe==0
        ! allocate biomass increment
        ! precompiler flag used here because option requires allocation to be called daily
        call allocation_annual( jpngr )
#endif

        ! reduce part of NPP for reproduction
        !call reproduction

        !call light
        !call mortality
        !call fire
        call establishment(jpngr)

        ! collect output field variables
        !call getout_fields

      end do

      ! try xxx
      c_uptake = 1.e15 ! in gC

      ! write output
      !call write_output

      end subroutine biosphere

      end program
