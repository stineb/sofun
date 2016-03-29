subroutine writeout_ascii( year, dtemp )
  !/////////////////////////////////////////////////////////////////////////



  ! XXXX DISCONTINUED XXXX




  
  ! Write daily ASCII output
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !-------------------------------------------------------------------------
  use md_params_siml, only: spinup, daily_out_startyr, &
    daily_out_endyr, outyear
  use md_params_core, only: ndayyear
  use md_params_modl, only: lu_category
  use md_outvars
  use md_nuptake, only: writeout_ascii_nuptake
  use md_littersom, only: writeout_ascii_littersom
  use md_params_siml, only: & 
      loutdgpp       &
    , loutdtransp    &
    , loutdnpp       &
    , loutdnup       &
    , loutdCleaf     &
    , loutdCroot     &
    , loutdClabl     &
    , loutdNlabl     &
    , loutdClitt     &
    , loutdNlitt     &
    , loutdCsoil     &
    , loutdNsoil     &
    , loutdlai       &
    , loutdninorg

  ! arguments
  integer, intent(in) :: year       ! simulation year
  real, dimension(ndayyear), intent(in) :: dtemp  ! xxx un-nice to have it passed as argument

  ! local variables
  real :: itime
  integer :: day, moy, jpngr

  ! xxx implement this: sum over gridcells? single output per gridcell?
  if (maxgrid>1) stop 'writeout_ascii: think of something ...'
  jpngr = 1


  !-------------------------------------------------------------------------
  ! Collect variables to output variables
  !-------------------------------------------------------------------------
  !do lu=1,nlu
  !  do pft=1,npft
  !    if (lu==lu_category(pft)) then
  !    end if
  !  end do
  !end do
  if (nlu>1) stop 'Output only for one LU category implemented.'


  !-------------------------------------------------------------------------
  ! DAILY OUTPUT
  ! Write daily value, summed over all PFTs / LUs
  ! xxx implement taking sum over PFTs (and gridcells) in this land use category
  !-------------------------------------------------------------------------
  if ( .not. spinup .and. outyear>=daily_out_startyr .and. outyear<=daily_out_endyr ) then

    ! Write daily output only during transient simulation
    do day=1,ndayyear

      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(outyear) + real(day-1)/real(ndayyear)
      
      if (loutdnpp      ) write(102,999) itime, sum(outdnpp(:,day,jpngr))
      if (loutdCleaf    ) write(103,999) itime, sum(outdCleaf(:,day,jpngr))
      if (loutdnup      ) write(104,999) itime, sum(outdnup(:,day,jpngr))
      if (loutdninorg   ) write(107,999) itime, sum(outdninorg(:,day,jpngr))
      if (loutdCroot    ) write(111,999) itime, sum(outdCroot(:,day,jpngr))
      if (loutdClabl    ) write(112,999) itime, sum(outdClabl(:,day,jpngr))
      if (loutdClitt    ) write(113,999) itime, sum(outdClitt(:,day,jpngr))
      if (loutdtransp   ) write(114,999) itime, sum(outdtransp(:,day,jpngr))
      if (loutdNlabl    ) write(115,999) itime, sum(outdNlabl(:,day,jpngr))
      if (loutdCsoil    ) write(118,999) itime, sum(outdCsoil(:,day,jpngr))
      if (loutdNlitt    ) write(119,999) itime, sum(outdNlitt(:,day,jpngr))
      if (loutdNsoil    ) write(120,999) itime, sum(outdNsoil(:,day,jpngr))
      if (loutdlai      ) write(121,999) itime, sum(outdlai(:,day,jpngr))
      write(110,999) itime, dtemp(day)
        
    end do
  end if

  !-------------------------------------------------------------------------
  ! ANNUAL OUTPUT
  ! Write annual value, summed over all PFTs / LUs
  ! xxx implement taking sum over PFTs (and gridcells) in this land use category
  !-------------------------------------------------------------------------
  itime = real(outyear)

  write(307,999) itime, sum(outaCveg2lit(:,jpngr))
  write(308,999) itime, sum(outaNveg2lit(:,jpngr))
  write(311,999) itime, sum(outanpp(:,jpngr))
  write(312,999) itime, sum(outaCveg(:,jpngr))
  write(316,999) itime, sum(outaninorg(:,jpngr))
  write(317,999) itime, sum(outanup(:,jpngr))
  write(318,999) itime, sum(outalai(:,jpngr))
  write(319,999) itime, sum(outanarea_mb(:,jpngr))
  write(320,999) itime, sum(outanarea_cw(:,jpngr))
  write(321,999) itime, sum(outacton_lm(:,jpngr))
  write(322,999) itime, sum(outalma(:,jpngr))
  write(323,999) itime, sum(outavcmax(:,jpngr))

  !-------------------------------------------------------------------------
  ! MODULE-SPECIFIC OUTPUT
  ! XXX WEIRD: can't remove these calls from here and I don't know why XXX
  !-------------------------------------------------------------------------
  call writeout_ascii_nuptake( year, spinup )
  call writeout_ascii_littersom( year, spinup )

  return

999     format (F20.8,F20.8)

end subroutine writeout_ascii
