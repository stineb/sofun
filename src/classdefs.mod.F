module _classdefs
  !////////////////////////////////////////////////////////////////
  !  Module contains Fortran 90 derived-type declarations to define
  !  objects
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  implicit none

  ! Minimum precision
  real, parameter :: epsilon = 1.0e-5 

  ! Carbon, so far contains only c12 (to be extended for c13)
  type carbon
   real :: c12
  end type carbon

  type nitrogen
   real :: n14
  end type nitrogen

  ! Organic pools, contain carbon (c12) and nitrogen (n14)
  type orgpool
   type(carbon)   :: c
   type(nitrogen) :: n
  end type orgpool

  !! Plants, contain leaves and roots
  !type plantclass
   !type(orgpool) :: lm ! leafmass
   !type(orgpool) :: rm ! rootmass
  !end type plantclass

  ! Soil, contains a pool with fast and slow turnover
  !type soilclass
   !type(orgpool) :: fs ! fast
   !type(orgpool) :: sl ! slow
  !end type soilclass

  ! Litter, contains ...
  !type litterclass
   !type(orgpool) :: af ! above-ground fast
   !type(orgpool) :: as ! above-ground slow
   !type(orgpool) :: bg ! below-ground
  !end type litterclass

contains
!=========================LOW-LEVEL================================


!--------------------------ORGANIC---------------------------------

  subroutine orgcp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" organic mass to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------

    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: to
    real, optional, intent(in) :: scale
    !real, optional, intent(in) :: d13C

    if (present(scale )) then
      call ccp( amount%c,to%c,scale )
      call ncp( amount%n,to%n,scale )
    else
      call ccp( amount%c,to%c)
      call ncp( amount%n,to%n)
    end if

    !         if (present(d13C)) then
    !           to%c%c12 = amount%c%c12 + to%c%c12
    !           to%n%n14 = amount%n%n14 + to%n%n14
    !         else
    !           to%c%c12 = amount%c%c12 + to%c%c12
    !           to%n%n14 = amount%n%n14 + to%n%n14
    !         end if

  end subroutine orgcp


  subroutine orgcpRec( amount, to, outc, outn, d13C)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" organic mass to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source.
    !  Additionally records amount copied (adding to outc and outn).
    !----------------------------------------------------------------    
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(inout) :: outn
    real, optional, intent(in) :: d13C

    outc = outc + amount%c%c12
    outn = outn + amount%n%n14

    call ccp( amount%c,to%c)
    call ncp( amount%n,to%n)

  end subroutine orgcpRec


  subroutine orgsub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" organic mass ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from

#if _check_sanity
    if ( amount%c%c12>from%c%c12+epsilon) then
      stop 'in ORGSUB: attempting to remove C amount > from-pool'
    else if ( amount%n%n14>from%n%n14+epsilon) then
      stop 'in ORGSUB: attempting to remove N amount > from-pool'
    else if (from%c%c12<0.0) then
      stop 'in ORGSUB: C in from-pool negative'
    else if (from%n%n14<0.0) then
      stop 'in ORGSUB: N in from-pool negative'
    endif
#endif


    call csub( amount%c,from%c)
    call nsub( amount%n,from%n)

  end subroutine orgsub


  subroutine orgmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" organic mass ('amount') from source pool 
    !  ('from') to destination pool ('to'). Substracts amount moved 
    !  ('amount') from source. 'orgmv' is the combination of 'orgcp' 
    !  and 'orgsub'
    !----------------------------------------------------------------
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from
    type(orgpool), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if (present(scale )) then
      call orgcp(orgfrac(scale,amount),to) 
      call orgsub( amount, from )       
    else
      call orgcp( amount, to)
      call orgsub( amount, from )
    endif  

  end subroutine orgmv


  subroutine orgmvRec( amount, from, to, outc, outn, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" organic mass ('amount') from source pool 
    !  ('from') to destination pool ('to') and to additionally record
    !  amount moved. Substracts amount moved ('amount') from source. 
    ! 'orgmvRec' is the combination of 'orgcp' and 'orgsub'.
    !----------------------------------------------------------------    
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from
    type(orgpool), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(inout) :: outn
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if (present(scale )) then
      outc = outc + amount%c%c12 * scale
      outn = outn + amount%n%n14 * scale
      call orgcp(orgfrac(scale,amount),to) 
      call orgsub( amount, from )
    else
      outc = outc + amount%c%c12
      outn = outn + amount%n%n14
      call orgcp( amount, to)
      call orgsub( amount, from )
    endif  

  end subroutine orgmvRec


  subroutine orginit( pool )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(orgpool), intent(inout) :: pool

    call cinit(pool%c)
    call ninit(pool%n)

  end subroutine orginit


  !--------------------------CARBON----------------------------------

  subroutine cmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'cmv' is the combination 
    !  of 'ccp' and 'csub'. 
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from
    type(carbon), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if (present(scale )) then
      call ccp(cfrac(scale,amount),to)
      call csub( amount, from )
    else
      call ccp( amount, to)
      call csub( amount, from )
    endif

  end subroutine cmv


  subroutine cmvRec( amount, from, to, outc, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'cmvRec' is the combination 
    !  of 'ccp' and 'csub'. Additionally adds addition to outc/outn.
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from
    type(carbon), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')


    if (present(scale )) then
      outc = outc + amount%c12 * scale
      call ccp(cfrac(scale,amount),to)
      call csub( amount, from )
    else
      outc = outc + amount%c12
      call ccp( amount, to)
      call csub( amount, from )
    endif

  end subroutine cmvRec


  subroutine ccp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" carbon to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: to
    real, optional, intent(in) :: scale
    !real, optional, intent(in) :: d13C

    if (present(scale )) then
      to%c12 = to%c12 + amount%c12 * scale
    else
      to%c12 = to%c12 + amount%c12 
    end if

    ! if (present(d13C)) then
    !   to%c%c12 = amount%c%c12 + to%c%c12
    !   to%n%n14 = amount%n%n14 + to%n%n14
    ! else
    !   to%c%c12 = amount%c%c12 + to%c%c12
    !   to%n%n14 = amount%n%n14 + to%n%n14
    ! end if

  end subroutine ccp


  subroutine ccpRec( amount, to, outc, d13C)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" carbon to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: to
    real, intent(inout) :: outc
    real, optional, intent(in) :: d13C

    to%c12 = amount%c12 + to%c12
    outc = outc + amount%c12

    ! if (present(d13C)) then
    !   to%c%c12 = amount%c%c12 + to%c%c12
    !   to%n%n14 = amount%n%n14 + to%n%n14
    ! else
    !   to%c%c12 = amount%c%c12 + to%c%c12
    !   to%n%n14 = amount%n%n14 + to%n%n14
    ! end if

  end subroutine ccpRec


  subroutine csub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" organic mass ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from

#if _check_sanity
    if ( amount%c12 > from%c12+epsilon) then
      write(0,*) 'amount', amount%c12
      write(0,*) 'from  ', from%c12
      write(0,*) 'in CSUB: attempting to remove amount > from-pool'
      stop
    endif
#endif
    from%c12 = from%c12 - amount%c12
     
  end subroutine csub


  subroutine cinit(pool)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(carbon), intent(inout) :: pool

    pool%c12 = 0.

  end subroutine cinit


  !--------------------------NITROGEN--------------------------------

  subroutine nmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'nmv' is the combination 
    !  of 'ccp' and 'csub'. 
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from
    type(nitrogen), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if (present(scale )) then
      call ncp(nfrac(scale,amount),to)
      call nsub( amount, from )
    else
      call ncp( amount, to)
      call nsub( amount, from )
    endif

  end subroutine nmv


  subroutine nmvRec( amount, from, to, outn, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'nmvRec' is the combination 
    !  of 'ccp' and 'csub'. Additionally adds addition to outc/outn.
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from
    type(nitrogen), intent(inout) :: to
    real, intent(inout) :: outn
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if (present(scale )) then
      outn = outn + amount%n14 * scale
      call ncp(nfrac(scale,amount),to)
      call nsub( amount, from )
    else
      outn = outn + amount%n14
      call ncp( amount, to)
      call nsub( amount, from )
    endif

  end subroutine nmvRec


  subroutine ncp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" nitrogen to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: to
    real, intent(in), optional :: scale

    if (present(scale )) then
      to%n14 = to%n14 + amount%n14 * scale
    else
      to%n14 = to%n14 + amount%n14 
    end if

  end subroutine ncp


  subroutine ncpRec( amount, to, outn )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" nitrogen to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: to
    real, intent(inout) :: outn

    to%n14 = amount%n14 + to%n14
    outn = outn + amount%n14

  end subroutine ncpRec


  subroutine nsub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" nitrogen ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from

#if _check_sanity
    if ( amount%n14>from%n14+epsilon) then
      stop 'in NSUB: attempting to remove amount > from-pool'
    endif
#endif
    from%n14 = from%n14 - amount%n14

    return

  end subroutine nsub


  subroutine ninit(pool)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(nitrogen), intent(inout) :: pool

    pool%n14 = 0.0

  end subroutine ninit


  !--------------------------FUNCTIONS--------------------------------

  function orgfrac( frac, from )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable of type 'orgpool' and size
    !  of a fraction 'frac' of source pool ('from')
    !----------------------------------------------------------------
 
    ! arguments
    real, intent(in)          :: frac
    type(orgpool), intent(in) :: from

    ! function return variable
    type(orgpool), intent(out) :: orgfrac

    orgfrac%c%c12 = frac * from%c%c12
    orgfrac%n%n14 = frac * from%n%n14

  end function orgfrac


  function cfrac( frac, from )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable of type 'carbon' and size
    !  of a fraction 'frac' of source pool ('from')
    !----------------------------------------------------------------
 
    ! arguments
    real, intent(in) :: frac
    type(carbon), intent(in) :: from

    ! function return variable
    type(carbon), intent(out) :: cfrac

    cfrac%c12 = frac * from%c12

  end function cfrac


  function nfrac( frac, from )
   !////////////////////////////////////////////////////////////////
   !  Generic function to return variable of type 'nitrogen' and size
   !  of a fraction 'frac' of source pool ('from')
   !----------------------------------------------------------------

   ! arguments
   real, intent(in) :: frac
   type(nitrogen), intent(in) :: from
   
   ! function return variable
   type(nitrogen), intent(out) :: nfrac

   nfrac%n14 = frac * from%n14

  end function nfrac


  function orgplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'orgpool'. Sum is of type 'orgpool' as well.
    !----------------------------------------------------------------

    ! arguments
    type(orgpool), intent(in) :: pool1
    type(orgpool), intent(in) :: pool2
    type(orgpool), intent(in),optional :: pool3
    type(orgpool), intent(in),optional :: pool4
    type(orgpool), intent(in),optional :: pool5
    type(orgpool), intent(in),optional :: pool6
    type(orgpool), intent(in),optional :: pool7
    type(orgpool), intent(in),optional :: pool8
    type(orgpool), intent(in),optional :: pool9
    type(orgpool), intent(in),optional :: pool10

    ! function return variable
    type(orgpool), intent(out) :: orgplus

    orgplus%c = cplus(pool1%c,pool2%c)
    orgplus%n = nplus(pool1%n,pool2%n)

    if (present(pool3)) then
      orgplus%c = cplus(orgplus%c,pool3%c)
      orgplus%n = nplus(orgplus%n,pool3%n)
      if (present(pool4)) then
        orgplus%c = cplus(orgplus%c,pool4%c)
        orgplus%n = nplus(orgplus%n,pool4%n)
        if (present(pool5)) then
          orgplus%c = cplus(orgplus%c,pool5%c)
          orgplus%n = nplus(orgplus%n,pool5%n)
          if (present(pool6)) then
            orgplus%c = cplus(orgplus%c,pool6%c)
            orgplus%n = nplus(orgplus%n,pool6%n)
            if (present(pool7)) then
              orgplus%c = cplus(orgplus%c,pool7%c)
              orgplus%n = nplus(orgplus%n,pool7%n)
              if (present(pool8)) then
                orgplus%c = cplus(orgplus%c,pool8%c)
                orgplus%n = nplus(orgplus%n,pool8%n)
                if (present(pool9)) then
                  orgplus%c = cplus(orgplus%c,pool9%c)
                  orgplus%n = nplus(orgplus%n,pool9%n)
                  if (present(pool10)) then
                    orgplus%c = cplus(orgplus%c,pool10%c)
                    orgplus%n = nplus(orgplus%n,pool10%n)
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if


  end function orgplus


  function cplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------

    ! arguments
    type(carbon), intent(in) :: pool1
    type(carbon), intent(in) :: pool2
    type(carbon), intent(in), optional :: pool3
    type(carbon), intent(in), optional :: pool4
    type(carbon), intent(in), optional :: pool5
    type(carbon), intent(in), optional :: pool6
    type(carbon), intent(in), optional :: pool7
    type(carbon), intent(in), optional :: pool8
    type(carbon), intent(in), optional :: pool9
    type(carbon), intent(in), optional :: pool10

    ! function return variable
    type(carbon), intent(out) :: cplus

    cplus%c12 = pool1%c12 + pool2%c12

    if (present(pool3)) then
      cplus%c12 = cplus%c12 + pool3%c12
      if (present(pool4)) then
        cplus%c12 = cplus%c12 + pool4%c12
        if (present(pool5)) then
          cplus%c12 = cplus%c12 + pool5%c12
          if (present(pool6)) then
            cplus%c12 = cplus%c12 + pool6%c12
            if (present(pool7)) then
              cplus%c12 = cplus%c12 + pool7%c12
              if (present(pool8)) then
                cplus%c12 = cplus%c12 + pool8%c12
                if (present(pool9)) then
                  cplus%c12 = cplus%c12 + pool9%c12
                  if (present(pool10)) then
                    cplus%c12 = cplus%c12 + pool10%c12
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if

  end function cplus


  function nplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'nitrogen'. Sum is of type 'nitrogen' as well.
    !----------------------------------------------------------------

    ! arguments
    type(nitrogen), intent(in) :: pool1
    type(nitrogen), intent(in) :: pool2
    type(nitrogen), intent(in), optional :: pool3
    type(nitrogen), intent(in), optional :: pool4
    type(nitrogen), intent(in), optional :: pool5
    type(nitrogen), intent(in), optional :: pool6
    type(nitrogen), intent(in), optional :: pool7
    type(nitrogen), intent(in), optional :: pool8
    type(nitrogen), intent(in), optional :: pool9
    type(nitrogen), intent(in), optional :: pool10

    ! function return variable
    type(nitrogen), intent(out) :: nplus

    nplus%n14 = pool1%n14 + pool2%n14

    if (present(pool3)) then
      nplus%n14 = nplus%n14 + pool3%n14
      if (present(pool4)) then
        nplus%n14 = nplus%n14 + pool4%n14
        if (present(pool5)) then
          nplus%n14 = nplus%n14 + pool5%n14
          if (present(pool6)) then
            nplus%n14 = nplus%n14 + pool6%n14
            if (present(pool7)) then
              nplus%n14 = nplus%n14 + pool7%n14
              if (present(pool8)) then
                nplus%n14 = nplus%n14 + pool8%n14
                if (present(pool9)) then
                  nplus%n14 = nplus%n14 + pool9%n14
                  if (present(pool10)) then
                    nplus%n14 = nplus%n14 + pool10%n14
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if

  end function nplus


  function orgminus( pool1, pool2 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'orgpool'. Sum is of type 'orgpool' as well.
    !----------------------------------------------------------------

    ! arguments
    type(orgpool), intent(in) :: pool1
    type(orgpool), intent(in) :: pool2

    ! function return variable
    type(orgpool), intent(out) :: orgminus

    orgminus%c = cminus( pool1%c, pool2%c )
    orgminus%n = nminus( pool1%n, pool2%n )

  end function orgminus


  function cminus( pool1, pool2 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------
    ! arguments
    type(carbon), intent(in) :: pool1
    type(carbon), intent(in) :: pool2

    ! function return variable
    type(carbon), intent(out) :: cminus

    cminus%c12 = pool1%c12 - pool2%c12

  end function cminus


  function nminus( pool1, pool2 )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------
    ! arguments
    type(nitrogen), intent(in) :: pool1
    type(nitrogen), intent(in) :: pool2

    ! function return variable
    type(nitrogen), intent(out) :: nminus

    nminus%n14 = pool1%n14 - pool2%n14

  end function nminus


  function cton( pool, default )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return the C:N ratio of an organic pool.
    !----------------------------------------------------------------
    ! arguments
    type(orgpool), intent(in) :: pool
    real, intent(in), optional :: default

    ! function return variable
    real, intent(out) :: cton

    if (present(default)) then
      if (pool%n%n14==0.0) then
        cton = default
      else
        cton = pool%c%c12 / pool%n%n14
      end if
    else
#if _check_sanity
    if (pool%n%n14==0.) then
      stop 'in CTON: N is zero'
    endif
    if (pool%n%n14<0.0 .or. pool%c%c12<0.0) then
      stop 'in CTON: C and/or N is negative'
    endif
#endif
    cton = pool%c%c12 / pool%n%n14
    end if

  end function cton


  function ntoc( pool, default )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return the N:C ratio of an organic pool.
    !  This is equal to the inverse of the 'cton' function.
    !----------------------------------------------------------------
    ! arguments
    type(orgpool), intent(in) :: pool
    real, intent(in), optional :: default

    ! function return variable
    real :: ntoc

    if (present(default)) then
      if (pool%c%c12==0.0) then
        ntoc = default
      else
        ntoc = pool%n%n14 / pool%c%c12
      end if
    else
#if _check_sanity
      if (pool%c%c12==0.) then
        stop 'in NTOC: C is zero'
      endif
      if (pool%n%n14<0.0 .or. pool%c%c12<0.0) then
        stop 'in NTOC: C and/or N is negative'
      endif
#endif
      ntoc = pool%n%n14 / pool%c%c12
    end if

  end function ntoc

end module _classdefs
