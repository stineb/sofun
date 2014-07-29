       program main
       !////////////////////////////////////////////////////////////////
       !  Main program for site scale simulations, here used for 
       !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
       !  Nitrogen)
       !----------------------------------------------------------------
       use params_core
       use params_siml
       use params_site
       use params_modl
       use forcing
       use gridvars
       use dailyvars

       implicit none
       character(len=maxlen_runname) :: runname

       integer :: year
       integer :: realyear

       real, dimension(ndayyear) :: dtemp
       real, dimension(ndayyear) :: dprec
       real, dimension(ndayyear) :: dndep
       real                      :: pco2
       real, dimension(3)        :: landuse       ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'biosphere'
       real                      :: c_uptake      ! annual net global C uptake by biosphere


       !----------------------------------------------------------------
       ! READ RUNNAME FROM STANDARD INPUT
       !----------------------------------------------------------------
       read (*,*) runname
       ! make sure runname length is smaller/equal than maxlen_runname
       if (len_trim(runname)>=maxlen_runname) then
         print*,'abort'
       endif
       ! write simulation name to standard output (screen)
       write(0,*) '-------------SOFUN : '//trim(runname)//'--------------'
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
       write(0,*) '-------------START OF SIMULATION--------------'
       do year=1,runyears

         !----------------------------------------------------------------
         ! Define simulations "steering" variables (realyear, etc.)
         !----------------------------------------------------------------
         call getsteering(year,realyear)

         if (year == spinupyears+1 ) then
           write(0,*) '------------TRANSIENT SIMULATION--------------'
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

100      format(A,I6,F8.2,F8.2)

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
       integer, intent(in)                   :: year
       integer, intent(in)                   :: realyear   ! used here because GPP is input from obs.
       real, intent(in), dimension(ndayyear) :: dtemp
       real, intent(in), dimension(ndayyear) :: dprec
       real, intent(in), dimension(ndayyear) :: dndep
       real, intent(in)                      :: pco2
       real, dimension(3)                    :: landuse    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
       real, intent(out)                     :: c_uptake   ! annual net global C uptake by biosphere

       ! local variables
       integer :: dm, mo, jpngr

       if (init) then
         
         ! Initialise variables and/or read from restart file (not implemented)
         call initcell

         ! Open input/output files
         !call initio

         ! Get soil parameters (if not defined in <sitename>.parameter)
         !call getsoilpar

         !----------------------------------------------------------------
         ! GET MODEL PARAMETERS
         ! read model parameters that may varied for optimisation
         !----------------------------------------------------------------
         call getpar_modl

         ! Obtain land unit dependent parameters, define decomposition rates
         !call luparameters

       endif

       !call conversion( lu_area )

       ! initialise daily variables
       call initdaily

       ! LOOP THROUGH GRIDCELLS
       do jpngr=1,maxgrid

         ! LOOP THROUGH MONTHS
         do mo=1,nmonth

           ! LOOP THROUGH DAYS
           do dm=1,ndaymonth(mo)


             ! get GPP from input
             call gpp( dummy, realyear, mo, dm )

             ! calculate water balance
             !call waterbalance

             ! use SR 'sitchtemp'
             !call heatdiff

             ! calculate autotrophic respiration and exudation
             ! exudation is calculated by SOFUN model
             ! call npp

             ! reduce part of NPP for reproduction
             !call reproduction

             ! plant biomass turnover
             !call turnover

             ! litter and soil turnover
             call littersom( dm, mo, jpngr )

             ! inorganic soil N chemistry
             !call ntransform

             ! collect daily output variables from daily updated state variables
             !call getoutdaily

           end do

         end do

       end do

       ! allocate biomass increment 
       !call allocation

       ! collect annual output variables
       !call outannual

       ! write output
       !call write_output

       end subroutine biosphere

       end program main
