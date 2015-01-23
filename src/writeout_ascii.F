subroutine writeout_ascii( year )
  !/////////////////////////////////////////////////////////////////////////
  ! Write daily ASCII output
  !-------------------------------------------------------------------------
  use params_siml, only: firstyeartrend, spinupyears, spinup
  use params_core, only: ndayyear
  use params_modl, only: lu_category
  use fluxes
  use pools
  use outvars
  use waterbal_stash, only: writeout_ascii_waterbal
  use nuptake_fun, only: writeout_ascii_nuptake

  ! Arguments
  integer, intent(in) :: year       ! simulation year

  ! Local variables
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
  if (.not.spinup) then

    print*,'year',year

    ! Write daily output only during transient simulation
    do day=1,ndayyear

      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

      write(101,999) itime, sum(outdgpp(:,day,jpngr))
      write(102,999) itime, sum(outdnpp(:,day,jpngr))
      write(103,999) itime, sum(outdCleaf(:,day,jpngr))
      write(104,999) itime, sum(outdnup(:,day,jpngr))
      write(105,999) itime, sum(outdcex(:,day,jpngr))
      write(106,999) itime, sum(outdnetmin(:,day,jpngr))
      write(107,999) itime, sum(outdninorg(:,day,jpngr))
      write(108,999) itime, sum(outdnfixfree(:,day,jpngr))

    end do
  end if


  !-------------------------------------------------------------------------
  ! MONTHLY OUTPUT
  ! Write monthly value
  !-------------------------------------------------------------------------
  !        do moy=1,nmonth        !

  !          itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(moy-1)/real(nmonth)        !

  !        end do


  !-------------------------------------------------------------------------
  ! ANNUAL OUTPUT
  ! Write annual value, summed over all PFTs / LUs
  ! xxx implement taking sum over PFTs (and gridcells) in this land use category
  !-------------------------------------------------------------------------
  itime = real(year) + real(firstyeartrend) - real(spinupyears)
  
  write(301,999) itime, sum(outaClitt(:,jpngr))
  write(302,999) itime, sum(outaCsoil(:,jpngr))
  write(303,999) itime, sum(outaCalloc(:,jpngr))
  write(304,999) itime, sum(outanreq(:,jpngr))
  write(305,999) itime, sum(outaClit2soil(:,jpngr))
  write(306,999) itime, sum(outaNlit2soil(:,jpngr))
  write(307,999) itime, sum(outaCveg2lit(:,jpngr))
  write(308,999) itime, sum(outaNveg2lit(:,jpngr))
  write(309,999) itime, sum(outaNalloc(:,jpngr))


  !-------------------------------------------------------------------------
  ! MODULE-SPECIFIC OUTPUT
  !-------------------------------------------------------------------------
  call writeout_ascii_waterbal( year, spinup )
  call writeout_ascii_nuptake( year, spinup )

  return

999     format (F20.8,F20.8)

end subroutine writeout_ascii
