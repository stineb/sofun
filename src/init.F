      subroutine initglobal
      !////////////////////////////////////////////////////////////////
      !  Initialisation of all pools on all gridcells at the beginning
      !  of the simulation.
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use pools
      use annvars

      ! initialise all pools with zero
      pleaf(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))
      proot(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))
      plabl(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))

      pexud(:,:)    = carbon(0.0)

      plitt_af(:,:) = orgpool(carbon(0.0),nitrogen(0.0))
      plitt_as(:,:) = orgpool(carbon(0.0),nitrogen(0.0))
      plitt_bg(:,:) = orgpool(carbon(0.0),nitrogen(0.0))

      psoil_sl(:,:) = orgpool(carbon(0.0),nitrogen(0.0))
      psoil_fs(:,:) = orgpool(carbon(0.0),nitrogen(0.0))

      nh4(:,:)      = nitrogen(0.0)
      no3(:,:)      = nitrogen(0.0)
      no2(:,:)      = nitrogen(0.0)
    
      no_w(:,:)     = 0.0
      no_d(:,:)     = 0.0
      n2o_w(:,:)    = 0.0
      n2o_d(:,:)    = 0.0
      n2_w(:,:)     = 0.0

      ! initialise other properties
      fpc_grid(:,:)  = 0.0
      ispresent(:,:) = .false.

      end subroutine initglobal


      subroutine initpft(pft,jpngr)
      !////////////////////////////////////////////////////////////////
      !  Initialisation of specified PFT on specified gridcell
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use pools
      use annvars

      implicit none
      integer, intent(in) :: pft
      integer, intent(in) :: jpngr

      ! initialise all pools with zero
      pleaf(pft,jpngr) = orgpool(carbon(0.),nitrogen(0.))
      proot(pft,jpngr) = orgpool(carbon(0.),nitrogen(0.))
      plabl(pft,jpngr) = orgpool(carbon(0.),nitrogen(0.))

      ! initialise other properties
      fpc_grid(pft,jpngr)  = 0.0
      nind(pft,jpngr)      = 0.0
      lai_ind(pft,jpngr)   = 0.0
      crownarea(pft,jpngr) = 0.0

      !if (.not.grass(pft)) then
      !sm_ind(pft,jpngr,1)=orgpool(carbon(0.),nitrogen(0.))
      !hm_ind(pft,jpngr,1)=orgpool(carbon(0.),nitrogen(0.))
      ! end if
    
      ! xxx put this somewhere else
      !if(.not.tree(pft)) crownarea(pft,jpngr)=1.0d0
      !leafon(pft,jpngr)=.true.
      !leafondays(pft,jpngr)=0.0d0
      !leafoffdays(pft,jpngr)=0.0d0

      end subroutine initpft


      subroutine initdaily
      !////////////////////////////////////////////////////////////////
      !  Initialises all daily variables with zero.
      !----------------------------------------------------------------
      use dailyvars
      use pools

      implicit none

      dnpp(:)    = carbon(0.0)
      drsoil(:)  = carbon(0.0)
      drhet(:)   = carbon(0.0)
      dnep(:)    = carbon(0.0)
      dcex(:)    = carbon(0.0)
      dreco(:)   = carbon(0.0)
      dgpp(:)    = 0.0
      drauto(:)  = 0.0
      drleaf(:)  = 0.0
      drroot(:)  = 0.0
      drgrow(:)  = 0.0

      dnup(:)    = nitrogen(0.0)
      dnfix_free(:) = 0.0
      dnup_pas(:) = 0.0
      dnup_act(:) = 0.0
      dnup_fix(:) = 0.0
      dnup_ret(:) = 0.0

      ! xxx try
      ! hold soil carbon and litter pools at a constant size
      ! here example gridcell 9.5 E, 47.5 N
      psoil_fs = orgpool( carbon(4000.), nitrogen(150.) )
      psoil_sl = orgpool( carbon(10000.), nitrogen(550.) )

      plitt_af = orgpool( carbon(500.), nitrogen(10.))
      plitt_as = orgpool( carbon(1500.), nitrogen(30.))
      plitt_bg = orgpool( carbon(400.), nitrogen(8.))


      end subroutine initdaily


      subroutine initannual
      !////////////////////////////////////////////////////////////////
      !  Initialises all daily variables with zero.
      !----------------------------------------------------------------
      use annvars

      implicit none

      aestab(:) = orgpool(carbon(0.0),nitrogen(0.0))
      anpp(:)   = carbon(0.0)
      agpp(:)   = 0.0

      end subroutine initannual
