      subroutine npp(jpngr,dtemp,dtemp_soil)
      !/////////////////////////////////////////////////////////////////////////
      ! NET PRIMARY PRODUCTIVITY
      ! Calculate maintenance and growth respiration and substract this from GPP 
      ! to get NPP before additional root respiration for nutrient uptake (see 
      ! SR nuptake).
      !-------------------------------------------------------------------------

      use classdefs
      use params_core, only: npft, nlu
      use params_modl, only: rresp_growth, lu_category
      use rates
      use dailyvars, only: dgpp, dnpp, dra
      use pools, only: pleaf, proot
      use annvars, only: nind, ispresent

      implicit none

      ! ARGUMENTS
      integer, intent(in)              :: jpngr
      real, intent(in)                 :: dtemp      ! air temperature at this day
      real, dimension(nlu), intent(in) :: dtemp_soil ! soil temperature temperature at this day

      ! LOCAL VARIABLES
      integer :: pft
      integer :: lu

      real :: ftemp_air
      real :: ftemp_soil

      real :: dresp_leaf = 0.0        ! daily leaf mainenance respiration (=lresp_d)
      real :: dresp_root = 0.0        ! daily root mainenance respiration (=rresp_d)
      real :: dresp_maint = 0.0      ! daily total maintenance respiration
      real :: dresp_growth = 0.0      ! daily growth respiration
      
      !-------------------------------------------------------------------------
      ! PFT LOOP
      !-------------------------------------------------------------------------
      do pft=1,npft
        if (ispresent(pft,jpngr)) then
          
          lu=lu_category(pft)
          
          ! reference temperature: 10°C
          ftemp_air  = ftemp( dtemp, "lloyd_and_taylor" )   ! xxx leave dtemp really with day-dimension?
          ftemp_soil = ftemp( dtemp_soil(lu), "lloyd_and_taylor" )

          
          !/////////////////////////////////////////////////////////////////////////
          ! MAINTENANCE RESPIRATION
          ! use function 'resp_main'
          !-------------------------------------------------------------------------
          dresp_leaf = resp_maint( pleaf(pft,jpngr)%n%n14*nind(pft,jpngr), pft, ftemp=ftemp_air )
          dresp_root = resp_maint( proot(pft,jpngr)%n%n14*nind(pft,jpngr), pft, ftemp=ftemp_soil )
          ! done here: sapwood respiration for trees
          
          ! daily total maintenance repsiration
          dresp_maint = dresp_leaf + dresp_root

          !/////////////////////////////////////////////////////////////////////////
          ! GROWTH RESPIRATION
          ! XXX address this: here growth respiration is calculated before exudation
          ! is substracted from (GPP-Rmaint). This does not make sense as not all of
          ! (GPP-Rmaint) will later be available for growth. This raises the question
          ! whether N-uptake should be calculated here (after maintentance resp.).
          ! Further, this raises the question whether allocation should be calculated
          ! daily in order to calculate growth respiration on only the part that is 
          ! allocated to actual growth and not on the part that is diverted to 'plabl'.
          !-------------------------------------------------------------------------
          dresp_growth = max( (dgpp(pft)%c12-dresp_maint)*rresp_growth(pft), 0.0 )

          !/////////////////////////////////////////////////////////////////////////
          ! DAILY NPP 
          ! NPP is the sum of C available for growth and for N uptake 
          !-------------------------------------------------------------------------
          dnpp(pft) = carbon( dgpp(pft)%c12 - dresp_maint - dresp_growth )
          dra(pft)  = carbon( dresp_maint + dresp_growth )
    
        endif                                               !(present(pft,jpngr))
      enddo                                                 !pft

      contains

      function resp_maint( nmass, pft, ftemp )
      !*************************************************************************
      ! MAINTENANCE RESPIRATION
      ! Calculate tissue maintenance respiration values today [Eqn (7)]
      ! incorporating daily phenology. This is using the new leaf respiration method
      ! (spahni 2009-12-11: set leaf_resp_rs==0 for old respiration method)
      !-------------------------------------------------------------------------
      use params_modl, only: rresp_maint, kresp

      implicit none
#include "sofun_module_control.inc"
      real :: resp_maint                    ! return value: maintenance respiration [gC/m2]
      real, intent(in) :: nmass             ! N mass per unit area [gN/m2]
      integer, intent(in) :: pft            ! PFT number
      real, intent(in), optional :: ftemp   ! temperature modifier

#if _resp_maint_makela
      resp_maint = rresp_maint(pft) * nmass
#else
      resp_maint = rresp_maint(pft) * kresp * nmass * ftemp ! xxx try: *dphen(day,pft)
#endif

      ! alternative method in LPX: directly from photosynthesis ('dresp' in LPX-Bern)

      end function resp_maint

      end subroutine npp


