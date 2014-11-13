      subroutine npp( jpngr, dtemp, dtemp_soil )
      !/////////////////////////////////////////////////////////////////////////
      ! NET PRIMARY PRODUCTIVITY
      ! Calculate maintenance and growth respiration and substract this from GPP 
      ! to get NPP before additional root respiration for nutrient uptake (see 
      ! SR nuptake). NPP is defined so as to include C allocated to growth as 
      ! well as to exudation. Thus, exudation is not part of autotrophic respir-
      ! ation, but diverted to the quickly decaying exudates pool. Exudates decay
      ! is calculated in SR 'littersom' and is kept track of as soil respiration 
      ! ('rsoil'). This implies that growth respiration is "paid" also on exu-
      ! dates. 
      !-------------------------------------------------------------------------
      use classdefs
      use params_core, only: npft, nlu
      use params_modl, only: rresp_growth, lu_category, cton_pro
      use rates
      use statevars, only: dgpp, dnpp, drauto, drleaf, drroot, drhet, drgrow, dcex, dnup, nind, ispresent
      use pools, only: pleaf, proot, pexud, plabl

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
      real :: dresp_maint = 0.0    ! daily total maintenance respiration
      
      !-------------------------------------------------------------------------
      ! PFT LOOP
      !-------------------------------------------------------------------------
      do pft=1,npft
        if (ispresent(pft,jpngr)) then

          ! xxx consider this: order of PFT within this loop is relevant for N uptake
          ! maybe permute order? maybe divide N availability up for each PFT?
          
          lu=lu_category(pft)
          
          ! reference temperature: 10°C
          ftemp_air  = ftemp( dtemp, "lloyd_and_taylor" )   ! xxx leave dtemp really with day-dimension?
          ftemp_soil = ftemp( dtemp_soil(lu), "lloyd_and_taylor" )

          
          !/////////////////////////////////////////////////////////////////////////
          ! MAINTENANCE RESPIRATION
          ! use function 'resp_main'
          !-------------------------------------------------------------------------
          ! argument 'ftemp' not used when _resp_maint_makela activated
          ! fine roots should have a higher repsiration coefficient than other tissues (Franklin et al., 2007).
          !drleaf(pft) = carbon( resp_maint( pleaf(pft,jpngr)%n%n14*nind(pft,jpngr), pft, ftemp=ftemp_air ) )
          !drroot(pft) = carbon( resp_maint( proot(pft,jpngr)%n%n14*nind(pft,jpngr), pft, ftemp=ftemp_soil ) )
          drleaf(pft) = resp_maint( pleaf(pft,jpngr)%n%n14 * nind(pft,jpngr), pft, ftemp=ftemp_air )
          drroot(pft) = resp_maint( proot(pft,jpngr)%n%n14 * nind(pft,jpngr), pft, ftemp=ftemp_soil )
          ! done here: sapwood respiration for trees
          
          ! daily total maintenance repsiration
          !dresp_maint = cplus( drleaf(pft), drroot(pft) )
          dresp_maint = drleaf(pft) + drroot(pft)


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
          !drgrow(pft) = cfrac( rresp_growth(pft), cminus( dgpp(pft), dresp_maint ) )
          drgrow(pft) = rresp_growth(pft) * ( dgpp(pft) - dresp_maint )


          !/////////////////////////////////////////////////////////////////////////
          ! DAILY NPP 
          ! NPP is the sum of C available for growth and for N uptake 
          ! This is where isotopic signatures are introduced because only 'dbminc'
          ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
          ! trophic respiration is immediate, it makes thus no sense to calculate 
          ! full isotopic effects of gross exchange fluxes.
          ! XXX NPP may get negative here!
          !-------------------------------------------------------------------------
          drauto(pft) = dresp_maint + drgrow(pft)
          dnpp(pft)   = carbon( max(dgpp(pft) - drauto(pft), 0.) ) 

          !if (dnpp(pft)%c12<0.) then
            !print*,'gpp, ra ',dgpp(pft), drauto(pft)
            !print*,'negative NPP'
          !else
            !print*,'positive NPP'
          !endif 

          !/////////////////////////////////////////////////////////////////////////
          ! EXUDATION FOR N UPTAKE
          ! This calculates exudation 'dcex', N uptake 'dnup', ...
          ! Labile C exuded for N uptake in interaction with mycorrhiza.
          ! Calculate otpimal C expenditure for N uptake (FUN approach).
          ! PFT loop has to be closed above to get total NPP over all PFTs in each 
          ! LU. 
          !-------------------------------------------------------------------------
          !call nuptake( jpngr, pft )

          ! xxx try
          dcex(pft) = carbon(0.0)
          dnup(pft) = nitrogen( dnpp(pft)%c12/cton_pro(pft) )

          ! Add exuded C to exudates pool (fast decay)
          call ccp( dcex(pft), pexud(pft,jpngr) )

          !/////////////////////////////////////////////////////////////////////////
          ! TO LABILE POOL
          ! NPP available for growth first enters the labile pool ('plabl ').
          ! XXX Allocation is called here without "paying"  growth respir.?
          !-------------------------------------------------------------------------
          !print*,'dnpp, dcex, dnup ', dnpp(pft), dcex(pft), dnup(pft)
          call orgcp( orgpool( cminus( dnpp(pft), dcex(pft) ), dnup(pft) ), plabl(pft,jpngr) )

          !print*,'plabl ',plabl(pft,jpngr)%c%c12

        endif

      end do

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




