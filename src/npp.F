      subroutine npp(jpngr,dtemp,dtemp_soil)

      use classdefs
      use params_core, only: npft, nlu
      use params_modl, only: respcoeff, kresp, lu_category
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
      real :: dresp_maintn = 0.0      ! daily total maintenance respiration
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
          ! Calculate tissue maintenance respiration values today [Eqn (7)]
          ! incorporating daily phenology. This is using the new leaf respiration method
          ! (spahni 2009-12-11: set leaf_resp_rs==0 for old respiration method)
          !-------------------------------------------------------------------------
          
          ! done here: alternative method in LPX: directly from photosynthesis ('dresp' in LPX-Bern)

          ! daily leaf and root maintenance respiration
          dresp_leaf = respcoeff(pft)*kresp*pleaf(pft,jpngr)%c%c12*nind(pft,jpngr)*ntoc( pleaf(pft,jpngr) )*ftemp_air   ! xxx try: *dphen(day,pft)

          dresp_root = respcoeff(pft)*kresp*proot(pft,jpngr)%c%c12*nind(pft,jpngr)*ntoc( proot(pft,jpngr) )*ftemp_soil  ! xxx try: *dphen(d,pft) !(canexch.cpp:976)

          ! done here: sapwood respiration for trees
          
          ! daily total maintenance repsiration
          dresp_maintn = dresp_leaf + dresp_root


          !/////////////////////////////////////////////////////////////////////////
          ! GROWTH RESPIRATION
          ! 25% of what's left after maintenance respiration
          !-------------------------------------------------------------------------
          ! xxx put '0.25' to parameters!
          dresp_growth = max( (dgpp(pft)%c12-dresp_maintn)*0.25, 0.0 )


          !/////////////////////////////////////////////////////////////////////////
          ! DAILY NPP 
          ! NPP is the sum of C available for growth and for N uptake 
          !-------------------------------------------------------------------------
          dnpp(pft) = carbon( dgpp(pft)%c12 - dresp_maintn - dresp_growth )
          dra(pft)  = carbon( dresp_maintn + dresp_growth )
    
        endif                                               !(present(pft,jpngr))
      enddo                                                 !pft


      end subroutine npp
      
      
      ! xxxxxxxxxxxxxxxxxxx


      ! xxx try
!      do pft=1,npft

!        call ccp( cfrac( 0.5, dgpp(pft) ), dnpp(pft) )
!        call ccp( cfrac( 0.3, dgpp(pft) ), dra(pft) )
!        call ccp( cfrac( 0.2, dgpp(pft) ), pexud(pft,jpngr) )

!        ! accolate remaining NPP to labile pool before allocation
!        call ccp( dnpp(pft), plabl(pft,jpngr)%c )

!        ! xxx try: take up as much N as needed for NPP
!        dnup(pft) = nitrogen( dnpp(pft)%c12/cton_pro(pft) )

!        call ncp( dnup(pft), plabl(pft,jpngr)%n )

!        !print*,'dgpp, dnpp ',dgpp,dnpp

!      end do  

!      end subroutine npp