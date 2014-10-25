      subroutine getout_annual
      !////////////////////////////////////////////////////////////////
      !  SR called daily to sum up output variables.
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl, only: npft
      use dailyvars
      use annvars

      ! LOCAL VARIABLES
      integer :: pft

      do pft=1,npft

        call ccp( dgpp(pft), agpp(pft) )
        call ccp( dnpp(pft), anpp(pft) )

        call ncp( dnup(pft), anup(pft) )

      end do

      end subroutine getout_annual