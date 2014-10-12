      subroutine establishment(jpngr)
      !////////////////////////////////////////////////////////////////
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl
      use params_site, only: lpft1
      use annvars
      use pools
      use gridvars, only: lu_area

      implicit none

      ! ARGUMENTS
      integer,intent(in)       :: jpngr

      ! LOCAL VARIABLES
      integer                  :: pft
      integer                  :: lu
      real, dimension(nlu)     :: fpc_total
      integer, dimension(nlu)  :: npft_estab
      integer, dimension(nlu)  :: ngrass
      real, dimension(nlu)     :: bare

      logical, dimension(npft) :: estab

      ! was done here: Check if 'survive' is false or nind<nind_min. If so, kill PFT.

      ! was done here: if aprec(jpngr)<aprec_min_estab; then estab=.false.
      ! minimum precipitation for establishment to bioclim? (100 mm/yr used in LPX)

      ! define which PFTs to establish
      estab(:) = .false.
      do pft=1,npft
        if (pft==1) then
          if (lpft1) estab(pft) = .true.
        endif
      end do


      !------------------------------------------------------------------------------
      ! check if there is enough space for establishment
      !------------------------------------------------------------------------------
      fpc_total(:)=0.d0      
      do pft=1,npft
        lu=lu_category(pft)
        if(ispresent(pft,jpngr)) then
          fpc_total(lu)=fpc_total(lu)+fpc_grid(pft,jpngr)
        endif
      enddo                     ! pft

      do pft=1,npft
        lu=lu_category(pft)
        if(lu_area(lu,jpngr)<1.0d-12) estab(pft)=.false.
        if(fpc_total(lu)>1.0) estab(pft)=.false.
      enddo


      !------------------------------------------------------------------------------
      ! introduce new PFTs if any
      !------------------------------------------------------------------------------
      do pft=1,npft
        if (.not.ispresent(pft,jpngr).and.estab(pft)) then !note aprec condition is included in estab

          ispresent(pft,jpngr) = .true.

          ! initialise pools with 0. initialisation with sapling pool size is done below.
          call initpft(pft,jpngr)

          if (tree(pft)) then
            nind(pft,jpngr)=0.0
            crownarea(pft,jpngr)=0.0
          else  !grass or moss              
            nind(pft,jpngr)=1.0
            crownarea(pft,jpngr)=1.0
          endif

        endif
      enddo


      !------------------------------------------------------------------------------
      ! count number of establishing PFTs
      !------------------------------------------------------------------------------
      npft_estab(:)=0
      ngrass(:)=0
      do pft=1,npft
        if(ispresent(pft,jpngr)) then
          if(estab(pft)) then
            lu=lu_category(pft)
            if(tree(pft)) then  !trees
              npft_estab(lu)=npft_estab(lu)+1 
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
        ! belongs here: if npft_estab > 0, calculate estab_rate and estab_grid 
        if (ngrass(lu)>0) then ! .or. nmoss>0
          ! Grasses and mosses can establish in non-vegetated areas
          bare(lu)=(1.0-fpc_total(lu))/(real(ngrass(lu)))
        endif
      enddo


      !------------------------------------------------------------------------------
      ! initialise with sapling pool size
      !------------------------------------------------------------------------------
      do pft=1,npft
        lu=lu_category(pft)
        if(ispresent(pft,jpngr).and.estab(pft)) then

          if  (grass(pft)) then
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


      end subroutine establishment
