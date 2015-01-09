#if _fix_veg
      subroutine establishment( jpngr )
      !//////////////////////////////////////////////////////////////////////////////
      ! GET FIXED PRESCRIBED VEGETATION COVER
      ! Subroutine defines 'fpc_grid' and updates 'nind' following geo-
      ! metry constraints.
      ! b.stocker@imperial.ac.uk
      !-----------------------------------------------------------------------------
      use params_site, only: lpft1, lpft2, lpft3, lpft4, lpft5, lpft6,&
       lpft7, lpft8, lpft9, in_fpc_grid_1, in_fpc_grid_2, in_fpc_grid_3, &
       in_fpc_grid_4, in_fpc_grid_5, in_fpc_grid_6, in_fpc_grid_7, &
       in_fpc_grid_8, in_fpc_grid_9

      implicit none

      ! ARGUMENTS
      integer, intent(in)       :: jpngr

      ! LOCAL VARIABLES
      logical, dimension(npft) :: estab

      !-----------------------------------------------------------------------------
      ! define which PFTs to establish
      !-----------------------------------------------------------------------------
      estab(:) = .false.
      if (lpft1) estab(1) = .true.
      if (lpft2) estab(2) = .true.
      if (lpft3) estab(3) = .true.
      if (lpft4) estab(4) = .true.
      if (lpft5) estab(5) = .true.
      if (lpft6) estab(6) = .true.
      if (lpft7) estab(7) = .true.
      if (lpft8) estab(8) = .true.
      if (lpft9) estab(9) = .true.

      !-----------------------------------------------------------------------------
      ! Get fractional plant cover from site input file
      !-----------------------------------------------------------------------------
      if (lpft1) fpc_grid(1,jpngr) = in_fpc_grid_1
      if (lpft2) fpc_grid(2,jpngr) = in_fpc_grid_2
      if (lpft3) fpc_grid(3,jpngr) = in_fpc_grid_3
      if (lpft4) fpc_grid(4,jpngr) = in_fpc_grid_4
      if (lpft5) fpc_grid(5,jpngr) = in_fpc_grid_5
      if (lpft6) fpc_grid(6,jpngr) = in_fpc_grid_6
      if (lpft7) fpc_grid(7,jpngr) = in_fpc_grid_7
      if (lpft8) fpc_grid(8,jpngr) = in_fpc_grid_8
      if (lpft9) fpc_grid(9,jpngr) = in_fpc_grid_9

      end subroutine establishment

#else
      subroutine establishment( jpngr )
      !////////////////////////////////////////////////////////////////
      ! DEFINE DYNAMIC VEGETATION COVER
      ! Subroutine defines 'fpc_grid' and updates 'nind' following geo-
      ! metry constraints.
      ! b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl
      use params_site, only: lpft1, lpft2, lpft3, lpft4, lpft5, lpft6,&
       lpft7, lpft8, lpft9      
      use statevars, only: aestab, ispresent, nind, crownarea
      use pools
      use gridvars, only: lu_area
      use classdefs, only: ntoc

      implicit none

      ! ARGUMENTS
      integer,intent(in)       :: jpngr

      ! LOCAL VARIABLES
      integer                  :: pft
      integer                  :: lu
      real, dimension(nlu)     :: fpc_total
      real, dimension(nlu)     :: fpc_tree_total   !total grid FPC for tree PFTs
      real, dimension(nlu)     :: fpc_grass_total  !total grid FPC for tree PFTs
      integer, dimension(nlu)  :: ntrees
      integer, dimension(nlu)  :: ngrass
      real, dimension(nlu)     :: bare
      type(orgpool)            :: sm_ind_temp

      logical, dimension(npft) :: estab

      real                     :: estab_rate       ! sapling establishment rate over area available for establishment (indiv/m2)
      real, dimension(nlu)     :: nind_add         ! grid-level establishment rate (indiv/m2) (=estab_grid)
      real                     :: nind_old         ! previous (year's) nind

      ! was done here: Check if 'survive' is false or nind<nind_min. If so, kill PFT.

      ! was done here: if aprec(jpngr)<aprec_min_estab; then estab=.false.
      ! minimum precipitation for establishment to bioclim? (100 mm/yr used in LPX)

      !-----------------------------------------------------------------------------
      ! DYNAMIC ESTABLISHMENT ON BARE GROUND
      !-----------------------------------------------------------------------------
      ! define which PFTs to establish
      estab(:) = .false.
      if (lpft1) estab(1) = .true.
      if (lpft2) estab(2) = .true.
      if (lpft3) estab(3) = .true.
      if (lpft4) estab(4) = .true.
      if (lpft5) estab(5) = .true.
      if (lpft6) estab(6) = .true.
      if (lpft7) estab(7) = .true.
      if (lpft8) estab(8) = .true.
      if (lpft9) estab(9) = .true.
      
      !-----------------------------------------------------------------------------
      ! Get total area of trees and grassses (FPC = fractional plant cover)
      !------------------------------------------------------------------------------
      fpc_tree_total(:)=0.0      
      fpc_total(:)=0.0
      fpc_grass_total(:)=0.0

      fpc_tree_total = sum( fpc_grid(:,jpngr), mask=tree )


      do pft=1,npft
        lu=lu_category(pft)
        if(ispresent(pft,jpngr)) then
          if (tree(pft)) then
            fpc_tree_total(lu) = fpc_tree_total(lu)+fpc_grid(pft,jpngr)
          elseif (grass(pft)) then ! grass
            fpc_grass_total(lu) = fpc_grass_total(lu)+fpc_grid(pft,jpngr)
          endif
          fpc_total(lu) = fpc_total(lu)+fpc_grid(pft,jpngr)
        endif
      enddo                     ! pft

      !------------------------------------------------------------------------------
      ! Set 'estab' to false if FPC exceeds 1
      ! xxx this allows FPC>1 when trees and grasses coexist.
      !------------------------------------------------------------------------------
      do pft=1,npft
        lu=lu_category(pft)
        if(lu_area(lu,jpngr)<1.0e-12) estab(pft)=.false.
        if(tree(pft)) then
          if (fpc_tree_total(lu)>=1.0) estab(pft)=.false.
        elseif (grass(pft)) then 
          if (fpc_grass_total(lu)>=1.0) estab(pft)=.false.
        endif 
      enddo


      !------------------------------------------------------------------------------
      ! introduce new PFTs if any
      !------------------------------------------------------------------------------
      do pft=1,npft
        if (.not.ispresent(pft,jpngr).and.estab(pft)) then !note aprec condition is included in estab
          ispresent(pft,jpngr) = .true.
          ! initialise pools with 0. initialisation with sapling pool size is done below.
          call initpft(pft,jpngr)
        endif
      enddo


      !------------------------------------------------------------------------------
      ! count number of establishing PFTs
      !------------------------------------------------------------------------------
      ntrees(:)=0
      ngrass(:)=0
      do pft=1,npft
        if(ispresent(pft,jpngr)) then
          if(estab(pft)) then
            lu=lu_category(pft)
            if(tree(pft)) then  !trees
              ntrees(lu)=ntrees(lu)+1 
            elseif (grass(pft)) then !grasses
              ngrass(lu)=ngrass(lu)+1         
            endif
          endif
        endif
      enddo
      
      !------------------------------------------------------------------------------
      ! calculate establishment rate, for grasses ~ bare ground
      !------------------------------------------------------------------------------
      do lu=1,nlu

        ! Calculate grid-level establishment rate per woody PFT
        ! Space available for woody PFT establishment is proportion of grid cell
        ! not currently occupied by woody PFTs.
        ! xxx parameters 0.24 (=estab_max) and 5.0 (no explicit name in LPX) are
        ! adopted from LPX. Maybe add to modelparameters? Really necessary?
        if (ntrees(lu)>0.0) then
          estab_rate = max(0.0,0.24*(1.0-exp(5.0* &
              (fpc_tree_total(lu)-1.0))))/real(ntrees(lu))
          nind_add(lu) = estab_rate*(1.0-fpc_tree_total(lu))
        else
          nind_add(lu) = 0.0
        endif

        ! Grasses and mosses can establish in non-vegetated areas
        if (ngrass(lu)>0.0) then
          bare(lu) = (1.0-fpc_total(lu))/(real(ngrass(lu)))
        endif

      enddo


      !------------------------------------------------------------------------------
      ! initialise with sapling pool size
      !------------------------------------------------------------------------------
      do pft=1,npft
        lu=lu_category(pft)

        if (ispresent(pft,jpngr).and.estab(pft)) then

          if (tree(pft)) then
            ! //////////////////////////////////////////////////////////////////////
            ! TREE ESTABLISHMENT
            ! ----------------------------------------------------------------------                        
            ! Add new saplings to current population
            
            nind_old=nind(pft,jpngr)
            nind(pft,jpngr)=nind_old+nind_add(lu)
            

            pleaf(pft,jpngr) = distributed_to_nind( &
                pleaf(pft,jpngr), nind_old          &
              , lm_sapl(pft),     nind_add(lu)      &
              , nind(pft,jpngr)                     &
              )

            pwood(pft,jpngr) = distributed_to_nind( &
                pwood(pft,jpngr), nind_old          &
              , hm_sapl(pft),     nind_add(lu)      &
              , nind(pft,jpngr)                     &
              )

            proot(pft,jpngr) = distributed_to_nind( &
                proot(pft,jpngr), nind_old          &
              , rm_sapl(pft),     nind_add(lu)      &
              , nind(pft,jpngr)                     &
              )

            sm_ind_temp = distributed_to_nind( &
                psapw(pft,jpngr), nind_old          &
              , sm_sapl(pft),     nind_add(lu)      &
              , nind(pft,jpngr)                     &
              )

            ! record implicit C and N fixation
            call orgcp( orgfrac(nind_add(lu), lm_sapl(pft)), aestab(lu) )
            call orgcp( orgfrac(nind_add(lu), sm_sapl(pft)), aestab(lu) )
            call orgcp( orgfrac(nind_add(lu), hm_sapl(pft)), aestab(lu) )
            call orgcp( orgfrac(nind_add(lu), rm_sapl(pft)), aestab(lu) )

           
            ! Calculate height, diameter and crown area for new average
            ! individual such that the basic allometric relationships (A-C below)
            ! are satisfied.
            call update_tree_geometry( pft, jpngr, sm_ind_temp, ntoc(psapw(pft,jpngr)) )


          else if (grass(pft)) then
            !------------------------------------------------------------------------------
            ! GRASS ESTABLISHMENT
            !------------------------------------------------------------------------------
            call orgcp( orgfrac(bare(lu),lm_sapl(pft)), pleaf(pft,jpngr) )            
            call orgcp( orgfrac(bare(lu),rm_sapl(pft)), proot(pft,jpngr) )

            ! record implicit C and N fixation
            call orgcp( orgfrac(bare(lu),lm_sapl(pft)), aestab(lu) )        
            call orgcp( orgfrac(bare(lu),rm_sapl(pft)), aestab(lu) )        

            ! if PFT is not yet ispresent:
            ! add C13 signature for saplings: C4: co2(2)-3.6,    C3: co2(2)-17.8
            ! add C14 signature for saplings: C4: co2(2)-3.6*2., C3: co2(2)-17.8*2.

            ! if PFT is ispresent: 
            ! add C13 and C14 signature from annual NPP

          endif
        endif
      enddo                     !pft

      !------------------------------------------------------------------------------
      ! update LAI and FPC
      !------------------------------------------------------------------------------
      do pft=1,npft
        if (ispresent(pft,jpngr)) then
          call update_fpc( pft, jpngr )
        endif
      enddo  

      contains

      subroutine update_tree_geometry( pft, jpngr, sm_ind_temp, afntoc_sm )
        !////////////////////////////////////////////////////////////////
        ! Subroutine updating tree geometry.
        ! Requires global variables:
        !   - psapw, pwood
        ! Requires global model parameters:
        !   - wooddens, allom, crownarea_max, reinickerp, latosa, sla
        ! Updates global variables (geometry-related):
        !   - height, crownarea, psapw
        !----------------------------------------------------------------
        use pools, only: pleaf, psapw, pwood
        use params_modl, only: wooddens, allom1, allom2, allom3, &
          crownarea_max, reinickerp, latosa, sla
        use treegeometry, only: height, crownarea

        ! arguments
        integer, intent(in)          :: pft, jpngr
        type(orgpool), intent(in)    :: sm_ind_temp
        real, intent(in)             :: afntoc_sm

        ! local variables
        real    :: stemdiam

        ! (A) (leaf area) = latosa * (sapwood xs area)
        ! (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
        ! (B) (leaf mass) = lmtorm * (root mass)
        ! (C) height = allom2 * (stem diameter)**allom3
        ! (source?)
        ! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
        ! crownarea_max)

        ! From (A),
        ! (1) sap_xsa = lm_ind * sla / latosa
        ! (2) wooddens = (sm_ind + hm_ind) / stemvolume
        ! (3) stemvolume = stem_xsa * height
        ! From (1), (2) & (3),
        ! (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
        ! (5) stem_xsa = pi * (stemdiam**2) / 4
        ! From (5),
        ! (6) stemdiam = ( 4 * stem_xsa / pi )**0.5
        ! From (4) & (6),
        ! (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
        ! pi )**0.5
        ! From (C) & (7),
        ! (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
        ! ( allom2 * stemdiam**allom3 ) / pi )**0.5
        ! From (8),
        ! (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / pi /
        ! allom2 )**( 1 / (2 + allom3) )

        stemdiam = (4.0*(sm_ind_temp%c%c12+pwood(pft,jpngr)%c%c12)/ &
          wooddens(pft)/pi/allom2(pft))**(1.0/(2.0+allom3)) !Eqn 9

        height(pft,jpngr) = allom2(pft)*stemdiam**allom3 !Eqn C

        crownarea(pft,jpngr) = min(crownarea_max, allom1(pft)*stemdiam**reinickerp) !Eqn D

        ! note kuno 2005-01-18: allometry is independent of sapwood/heartwood ratio,
        ! only sum of the two is used (eq. 9)
        
        ! Recalculate sapwood mass, transferring excess sapwood to heartwood
        ! compartment, if necessary to satisfy Eqn A
        psapw(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12*height(pft,jpngr)*wooddens(pft)*sla(pft)/latosa(pft)
        
        ! xxx really do this? N budget may not be conserved!
        !psapw(pft,jpngr)%n%n14 = psapw(pft,jpngr)%n%n14*afntoc_sm

        ! change sg 21.10.2001: prevent negative allocation to heartwood
        call orgcp( orgminus( sm_ind_temp, psapw(pft,jpngr)), pwood(pft,jpngr) )

      end subroutine update_tree_geometry


      function distributed_to_nind( pool1, nind1, pool2, nind2, nind_new )
        !////////////////////////////////////////////////////////////////
        !  Generic function to return variable of type 'orgpool' and size
        !  of a fraction 'frac' of source pool ('from')
        !----------------------------------------------------------------
        implicit none

       ! arguments
        real, intent(in)           :: nind1, nind2, nind_new
        type (orgpool), intent(in) :: pool1, pool2
        
        ! function return variable
        type (orgpool) :: distributed_to_nind

        distributed_to_nind%c%c12 = (pool1%c%c12*nind1 + pool2%c%c12*nind2) &
          / nind_new

        distributed_to_nind%n%n14 = (pool1%n%n14*nind1 + pool2%n%n14*nind2) &
          / nind_new

      end function distributed_to_nind


      end subroutine establishment
#endif