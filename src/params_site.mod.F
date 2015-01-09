      module params_site
            !////////////////////////////////////////////////////////////////
            !  Module contains site-specific parameters read by getpar_site
            !----------------------------------------------------------------
            use params_siml
            implicit none
            real :: lon_site, lat_site, alt_site
            logical :: lpft1,lpft2,lpft3,lpft4,lpft5,lpft6,lpft7,lpft8,lpft9
            real    :: in_fpc_grid_1,in_fpc_grid_2,in_fpc_grid_3,in_fpc_grid_4,&
                  in_fpc_grid_5,in_fpc_grid_6,in_fpc_grid_7,in_fpc_grid_8,in_fpc_grid_9

            contains

            subroutine getpar_site
                  !////////////////////////////////////////////////////////////////
                  !  SR reads site-specific parameters from <sitename>.parameter
                  !----------------------------------------------------------------
                  implicit none

                  write(0,*) 'reading site parameter file ', sitename//".parameter ..."
                  call getreal(sitename//".parameter",'longitude',lon_site)
                  call getreal(sitename//".parameter",'latitude',lat_site)
                  call getreal(sitename//".parameter",'altitude',alt_site)

                  ! Get prescribed PFT selection 
                  ! xxx try: make npft a definable variable depending on this selection 
                  call getlogical(sitename//".parameter",'lpft1',lpft1)
                  call getlogical(sitename//".parameter",'lpft2',lpft2)
                  call getlogical(sitename//".parameter",'lpft3',lpft3)
                  call getlogical(sitename//".parameter",'lpft4',lpft4)
                  call getlogical(sitename//".parameter",'lpft5',lpft5)
                  call getlogical(sitename//".parameter",'lpft6',lpft6)
                  call getlogical(sitename//".parameter",'lpft7',lpft7)
                  call getlogical(sitename//".parameter",'lpft8',lpft8)
                  call getlogical(sitename//".parameter",'lpft9',lpft9)

                  ! Get prescribed fractional plant cover (FPC) for each PFT
                  call getreal(sitename//".parameter",'in_fpc_grid_1',in_fpc_grid_1)
                  call getreal(sitename//".parameter",'in_fpc_grid_2',in_fpc_grid_2)
                  call getreal(sitename//".parameter",'in_fpc_grid_3',in_fpc_grid_3)
                  call getreal(sitename//".parameter",'in_fpc_grid_4',in_fpc_grid_4)
                  call getreal(sitename//".parameter",'in_fpc_grid_5',in_fpc_grid_5)
                  call getreal(sitename//".parameter",'in_fpc_grid_6',in_fpc_grid_6)
                  call getreal(sitename//".parameter",'in_fpc_grid_7',in_fpc_grid_7)
                  call getreal(sitename//".parameter",'in_fpc_grid_8',in_fpc_grid_8)
                  call getreal(sitename//".parameter",'in_fpc_grid_9',in_fpc_grid_9)

            end subroutine getpar_site

      end module params_site
