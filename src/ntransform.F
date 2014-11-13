      subroutine ntransform( dm, mo, jpngr, dprec )
      !////////////////////////////////////////////////////////////////
      !  Litter and SOM decomposition and nitrogen mineralisation.
      !  1st order decay of litter and SOM pools, governed by temperature
      !  and soil moisture following LPJ (Sitch et al., 2003) and 
      !  Xu-Ri & Prentice (XXX).
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_core
      use classdefs
      !use params_siml, only: init
      use params_modl
      use rates
      use pools, only: nh4, no3, no2, no_w, no_d, n2o_w, n2o_d, n2_w
      use statevars, only: dn2o, dnloss, ddenitr, dnitr, dvol, dleach, dn2o, dn2, dno, doc

      implicit none

      ! ARGUMENTS
      integer, intent(in) :: mo                              ! month
      integer, intent(in) :: dm                              ! day of the current month
      integer, intent(in) :: jpngr                           ! grid cell number
      integer, dimension(ndayyear), intent(in) :: dprec      ! daily precipitation [mm/d]

      ! LOCAL VARIABLES
      integer :: lu                                          ! gridcell unit counter variable
      integer, save :: day

      real, save :: ph_soil
      real, save :: nh3max

      real :: dnmax                                          ! labile carbon availability modifier
      real :: ftemp_ninorg                                   ! temperature modifier

      real :: no3_inc, n2o_inc, no_inc, no2_inc, n2_inc      ! temporary variables
      real :: tmp                                            ! temporary variable
      
      real :: nh4_w, no3_w, no2_w                            ! anaerobic pools
      real :: nh4_d, no3_d, no2_d                            ! aerobic pools
      real :: doc_w, no, n2o, n2                             ! anaerobic pools

      real :: doc_d                                          ! aerobic pools

      ! xxx try
      real :: soilwater                                      ! auxiliary variable
      real :: soiltemp                                       ! auxiliary variable
      real, dimension(ndayyear,nlu) :: dfleach

      
      if ( dm==1 .and. mo==1 ) then
        !///////////////////////////////////////////////////////////////////////
        ! ANNUAL INITIALIZATION 
        !-----------------------------------------------------------------------
        day = 0

        ! Calculate soil PH using empirical relationship with annual precip
        ! Eq.5, Tab.5, XP08 (ntransform.cpp:65) (c++:aprec in mm/yr; F: mm/yr)
        !------------------------------------------------------------------
        ph_soil = 3810./(762.+sum(dprec(:)))+3.8

        ! Deprotonation of NH4 to NH3 depends on soil pH
        !------------------------------------------------------------------
        if (ph_soil>6.) then
          nh3max = 1.
        else
          nh3max = 0.00001
        endif
        
      endif
      
        
      !/////////////////////////////////////////////////////////////////////////
      ! DO EVERY DAY (in daily, as well as in monthly mode!)
      !-------------------------------------------------------------------------
      day=day+1


      !///////////////////////////////////////////////////////////////////////
      ! DAILY INITIALIZATION 
      !-----------------------------------------------------------------------
      dleach(:) = 0.
      dvol(:)   = 0.
      dn2o(:)   = 0.
      dn2(:)    = 0.
      dno(:)    = 0.


      ! LOOP OVER GRIDCELL LAND UNITS
      do lu=1,nlu

                
        ! Daily updated soil moisture and soil temperature are required for
        ! ntransform.
        ! xxx try
        soiltemp = 10.0
        soilwater = 0.5 ! XXX beni: this is wrong: dw1 is only plant available water. should be water-filled pore space = ( (porosity - ice) - (total fluid water volume) ) / dz
        dfleach(day,lu) = 0.05

        ! must rather be wtot_up which includes water below permanent wilting point (see waterbalance.F).
        !------------------------------------------------------------------
        ! reference temperature: 25°C
        ftemp_ninorg = min( 1., ftemp( soiltemp, "lloyd_and_taylor", ref_temp=25. ) )   


        !///////////////////////////////////////////////////////////////////////
        ! AMMONIUM VOLATILIZATION (ntransform.cpp:41)
        !-----------------------------------------------------------------------
        ! use mw1 for monthly timestep and wpool for daily, because this is updated daily
        ! XXX nh3max is not considered in the equations presented in the paper! XXX
        dvol(lu) = ( nh3max*(soilwater*(1.0-soilwater)) * ftemp_ninorg*ftemp_ninorg*exp(2.*(ph_soil-10.)) ) * nh4(lu,jpngr)%n14
        nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - dvol(lu)

        
        !///////////////////////////////////////////////////////////////////////
        ! NITRATE LEACHING
        !-----------------------------------------------------------------------
        ! Reduce NO3 by fraction dleach(lu)
        !------------------------------------------------------------------
        dleach(lu) = dfleach(day,lu) * no3(lu,jpngr)%n14
        no3(lu,jpngr)%n14 = no3(lu,jpngr)%n14 - dleach(lu)

        !///////////////////////////////////////////////////////////////////////
        ! SUBSTRATE PARTITIONING (ntransform.cpp:95)
        !------------------------------------------------------------------
        ! Nitrification (aerobic) and denitrification (anaerobic) can occur
        ! simulataneously in different microsites. Substrate is thus parti-
        ! tioned according to the water content.
        
        ! wet (anaerobic) fraction
        !------------------------------------------------------------------
        nh4_w = soilwater/3.3 * nh4(lu,jpngr)%n14
        no3_w = soilwater/3.3 * no3(lu,jpngr)%n14
        no2_w = soilwater/3.3 * no2(lu,jpngr)%n14

        doc_w = doc(lu)*soilwater/3.3
        
        ! dry (aerobic) fraction
        !------------------------------------------------------------------
        nh4_d = (1.0-soilwater/3.3) * nh4(lu,jpngr)%n14
        no3_d = (1.0-soilwater/3.3) * no3(lu,jpngr)%n14
        no2_d = (1.0-soilwater/3.3) * no2(lu,jpngr)%n14

        doc_d = doc(lu)*(1.0-soilwater/3.3)

        
        !///////////////////////////////////////////////////////////////////////
        ! NITRIFICATION in aerobic microsites (ntransform.cpp:123)
        !------------------------------------------------------------------
        ftemp_ninorg = max( min( (((70.-soiltemp)/(70.-38.))**12.) * exp(12.*(soiltemp-38.)/(70.-38.)), 1.), 0.)

        
        ! gross nitrification rate (Eq.1, Tab.8, XP08)
        !------------------------------------------------------------------
        no3_inc = MAXNITR*ftemp_ninorg*nh4_d
        dnitr(lu) = no3_inc
        nh4_d = nh4_d-no3_inc   

        
        ! NO from nitrification (Eq.3, Tab.8, XP08)
        !------------------------------------------------------------------
        no_inc = NON*no3_inc
        no3_inc = no3_inc-no_inc

        no_d(lu,jpngr) = no_d(lu,jpngr) + no_inc

        
        ! N2O from nitrification (Eq.4, Tab.8, XP08)
        !------------------------------------------------------------------
        n2o_inc = N2ON*no3_inc
        no3_inc = no3_inc-n2o_inc

        n2o_d(lu,jpngr) = n2o_d(lu,jpngr) + n2o_inc
        no3_d = no3_d + no3_inc
        
        
        !///////////////////////////////////////////////////////////////////////
        ! DENITRIFICATION (ntransform.cpp:177) in anaerobic microsites
        !------------------------------------------------------------------
        ! reference temperature: 22°C
        ftemp_ninorg = ftemp( soiltemp, "lloyd_and_taylor", ref_temp=22. )

        
        ! Effect of labile carbon availability on denitrification (Eq.2, Tab.9, XP08)
        ! doc is last year's doc because it is only available at the end of the month
        ! while this SR is calculated daily, even when _dailymode==0.
        !------------------------------------------------------------------
        dnmax = DOCMAX*doc_w/(KDOC+doc_w)                     ! dnmax < 1 for all doc_w 

        
        ! Denitrification ratio, NO3->NO2 (Eq.3, Tab.9, XP08)
        !------------------------------------------------------------------
        no2_inc = min(dnmax*ftemp_ninorg*no3_w/(KN+no3_w)*1000.,no3_w)
        if (no2_inc>no3_w) stop 'no2_inc > no3_w'

        ddenitr(lu) = no2_inc
        no3_w = no3_w - no2_inc
        no2_w = no2_w + no2_inc

        
        ! Transformation NO2->N2 (Eq.4., Tab.9, XP08)
        !------------------------------------------------------------------
        n2_inc = min(dnmax*ftemp_ninorg*no2_w/(KN+no2_w)*1000.,no2_w)
        if (n2_inc>no2_w) stop 'n2_inc > no2_w'

        no2_w = no2_w-n2_inc

        
        ! N2O from denitrification (Eq.6, Tab.9, XP08)
        !------------------------------------------------------------------
        ! n2o_inc = 0.018d0*ftemp_ninorg*(1.01d0-0.21d0*soilwater)*n2_inc  !Colin says 0.018 was used here. Code I got had 0.015
        ! Factor reduced from 1.8% to 1.2% to get ~6.5 TgN/yr N2O emissions
        n2o_inc = DNITR2N2O*ftemp_ninorg*(1.01-0.21*soilwater)*n2_inc
        n2_inc = n2_inc-n2o_inc

        n2o_w(lu,jpngr) = n2o_w(lu,jpngr) + n2o_inc


        ! NO from denitrification (Eq.5, Tab.9, XP08)
        !------------------------------------------------------------------
        no_inc = 0.0001*ftemp_ninorg*(1.01-0.21*soilwater)*n2_inc
        n2_inc = n2_inc-no_inc

        no_w(lu,jpngr) = no_w(lu,jpngr) + no_inc

        
        ! N2 from denitrification
        !------------------------------------------------------------------
        n2_w(lu,jpngr) = n2_w(lu,jpngr) + n2_inc


        !///////////////////////////////////////////////////////////////////////
        ! UPDATE POOLS (ntransform.cpp:389)
        !------------------------------------------------------------------
        ! nh4, no3 and no2 (sum of sub-pools in wet and dry microsites) are
        ! defined as global variables. They are split into sub-pools at the
        ! beginning of ntransform (see "substrate partitioning"), processed
        ! through "denitrification" and "nitrification" and added up again
        ! here. In contrast, each sub-pool for wet and dry conditions of no,
        ! n2o and n2 is defined as a global (common) variable, while the sum
        ! of the sup-pools (no, n2o and n2) is defined locally and only used
        ! for the diffusion/emission (see below).
        !------------------------------------------------------------------
        nh4(lu,jpngr)%n14 = nh4_w + nh4_d
        no3(lu,jpngr)%n14 = no3_w + no3_d
        no2(lu,jpngr)%n14 = no2_w + no2_d

        no  = no_w(lu,jpngr) + no_d(lu,jpngr)
        n2o = n2o_w(lu,jpngr) + n2o_d(lu,jpngr)
        n2  = n2_w(lu,jpngr)

                
        !///////////////////////////////////////////////////////////////////////
        ! Diffusion of NO, N2O and N2 from the soil (ntransform.cpp:281)
        !------------------------------------------------------------------
        ! reference temperature: 25°C
        ftemp_ninorg = min( 1., ftemp( soiltemp, "lloyd_and_taylor", ref_temp=25. ))


        ! Nitrification fluxes
        !------------------------------------------------------------------
        tmp = ftemp_ninorg*(1.0-soilwater)*no_d(lu,jpngr)
        no_d(lu,jpngr) = no_d(lu,jpngr)-tmp
        
        tmp = ftemp_ninorg*(1.0-soilwater)*n2o_d(lu,jpngr)
        n2o_d(lu,jpngr) = n2o_d(lu,jpngr)-tmp

        
        ! Denitrification fluxes
        !------------------------------------------------------------------
        tmp = ftemp_ninorg*(1.0-soilwater)*no_w(lu,jpngr)
        no_w(lu,jpngr) = no_w(lu,jpngr)-tmp
        
        tmp = ftemp_ninorg*(1.0-soilwater)*n2o_w(lu,jpngr)
        n2o_w(lu,jpngr) = n2o_w(lu,jpngr)-tmp
                   

        ! Total fluxes (XXX should be equal to sum of nitrification and
        ! denitrification fluxes. XXX)
        !------------------------------------------------------------------
        dno(lu) = ftemp_ninorg*(1.0-soilwater)*no
        dn2o(lu)    = ftemp_ninorg*(1.0-soilwater)*n2o
        dn2(lu) = ftemp_ninorg*(1.0-soilwater)*n2

        n2_w(lu,jpngr) = n2_w(lu,jpngr)-dn2(lu)


        ! Record for output 
        !------------------------------------------------------------------
        dnloss(lu)  = dn2(lu) + dno(lu) + dn2o(lu) + dleach(lu) + dvol(lu)
        ddenitr(lu) = ddenitr(lu)

        
      enddo                                                 ! lu

      end subroutine ntransform
           
