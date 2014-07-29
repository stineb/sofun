       module pools
       !////////////////////////////////////////////////////////////////
       !  Module contains all pool variables. These are handled as state
       !  variables of the model and contain a dummy spatial dimension
       !  (to be extended). Daily output for pool variables are done by
       !  copying value of pool to a daily variable. Always use 'p' as
       !  a prefix of variable name.
       !----------------------------------------------------------------
       use params_core
       use classdefs

       implicit none

       type (orgpool), dimension(npft,maxgrid)  :: pleaf    ! (=lm_ind)
       type (orgpool), dimension(npft,maxgrid)  :: proot    ! (=rm_ind)
       type (orgpool), dimension(npft,maxgrid)  :: plabl    ! labile pool (=bm_inc but contains also N) 

       type (carbon),  dimension(npft,maxgrid)  :: pexud

       type (orgpool), dimension(npft,maxgrid)  :: plitt_af
       type (orgpool), dimension(npft,maxgrid)  :: plitt_as
       type (orgpool), dimension(npft,maxgrid)  :: plitt_bg

       type (orgpool), dimension(nlu,maxgrid)   :: psoil_sl
       type (orgpool), dimension(nlu,maxgrid)   :: psoil_fs

       type (nitrogen), dimension(nlu,maxgrid)  :: nh4
       type (nitrogen), dimension(nlu,maxgrid)  :: no3
       type (nitrogen), dimension(nlu,maxgrid)  :: no2

       real, dimension(nlu,maxgrid)             :: no_w,  no_d
       real, dimension(nlu,maxgrid)             :: n2o_w, n2o_d
       real, dimension(nlu,maxgrid)             :: n2_w


       end module pools

