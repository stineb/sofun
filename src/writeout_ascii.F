subroutine writeout_ascii( year, dtemp )
  !/////////////////////////////////////////////////////////////////////////
  ! Write daily ASCII output
  !-------------------------------------------------------------------------
  use _params_siml, only: firstyeartrend, spinupyears, spinup
  use _params_core, only: ndayyear
  use _params_modl, only: lu_category
  use _outvars
  use _waterbal, only: writeout_ascii_waterbal
  use _nuptake, only: writeout_ascii_nuptake
  use _littersom, only: writeout_ascii_littersom

  ! Arguments
  integer, intent(in) :: year       ! simulation year
  real, dimension(ndayyear), intent(in) :: dtemp  ! xxx un-nice to have it passed as argument

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

    ! Write daily output only during transient simulation
    do day=1,ndayyear

      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

      write(101,999) itime, sum(outdgpp(:,day,jpngr))
      write(102,999) itime, sum(outdnpp(:,day,jpngr))
      write(103,999) itime, sum(outdCleaf(:,day,jpngr))
      write(104,999) itime, sum(outdnup(:,day,jpngr))
      write(107,999) itime, sum(outdninorg(:,day,jpngr))
      write(109,999) itime, sum(outdtemp_soil(:,day,jpngr))
      write(110,999) itime, dtemp(day)
      
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

  write(303,999) itime, sum(outaCalloc(:,jpngr))
  write(307,999) itime, sum(outaCveg2lit(:,jpngr))
  write(308,999) itime, sum(outaNveg2lit(:,jpngr))
  write(309,999) itime, sum(outaNalloc(:,jpngr))
  write(310,999) itime, sum(outagpp(:,jpngr))
  write(311,999) itime, sum(outanpp(:,jpngr))
  write(312,999) itime, sum(outaCveg(:,jpngr))
  write(316,999) itime, sum(outaninorg(:,jpngr))

  !-------------------------------------------------------------------------
  ! MODULE-SPECIFIC OUTPUT
  !-------------------------------------------------------------------------
  call writeout_ascii_waterbal( year, spinup )
  call writeout_ascii_nuptake( year, spinup )
  call writeout_ascii_littersom( year, spinup )

  return

999     format (F20.8,F20.8)

end subroutine writeout_ascii
