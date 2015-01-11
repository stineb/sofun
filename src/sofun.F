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
        use fluxes

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

        enddo

        write(0,*) '--------------END OF SIMULATION---------------'

100     format(A,I6,F8.2,F8.2)

      end program
