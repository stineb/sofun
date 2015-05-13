module _findroot_fzeroin
  !************************************************************************
  !*                                                                      *
  !* Compute a zero of a continuous real valued function with the         *
  !* ------------------------------------------------------------         *
  !* Zeroin method.                                                       *
  !* --------------                                                       *
  !*                                                                      *
  !* exported funktion:                                                   *
  !*   - zeroin():  Zeroin- method for real functions                     *
  !*                                                                      *
  !* Programming language: ANSI C                                         *
  !* Compiler:             Borland C++ 2.0                                *
  !* Computer:             IBM PS/2 70 mit 80387                          *
  !* Author:               Siegmar Neuner                                 *
  !* Modifications:        Juergen Dietel, Rechenzentrum, RWTH Aachen     *
  !* Source:               FORTRAN source code                            *
  !* Date:                 11. 27. 1992                                   *
  !* -------------------------------------------------------------------- *
  !* Ref.: "Numerical Algorithms with C By G. Engeln-Mueller and F. Uhlig,*
  !*        Springer-Verlag, 1996" [BIBLI 11].                            *
  !*                                                                      *
  !*                                 F90 version by J-P Moreau, Paris.    *
  !*                                        (www.jpmoreau.fr)             *
  !************************************************************************
  ! use _allocation, only: eval_imbalance

  parameter(FOUR=4.0,HALF=0.5d0,ONE=1.0,THREE=3.0,TWO=2.0)
  parameter(ZERO=0.0, XMACHEPS = 1.2e-16)   !machine smallest real number

  type outtype_zeroin
    real :: root
    real :: froot
    integer :: niter
    integer :: error
  end type outtype_zeroin

contains

  function zeroin( func, myabserr, myrelerr, fmax, x1, x2, protnam ) result( out_zeroin )
    !************************************************************************
    !* Given a real valued function f on an interval [a,b] with             *
    !* f(a) * f(b) < 0, we compute a root of f in [a,b].                    *
    !* The Zeroin method combines bisection and secant methods with inverse *
    !* quadratic interpolation.                                             *
    !*                                                                      *
    !* Input parameters:                                                    *
    !* =================                                                    *
    !* abserr\  Error bounds with   abserr >= 0  and  relerr >= 0. Their    *
    !* relerr/  sum must exceed zero. For break-off we use the test         *
    !*              |xm|  <=  0.5 * (abserr + |b| * relerr),                *
    !*          where xm denotes half the length of the final inclusion     *
    !*          interval. For relerr = 0 we test only for absolute error,   *
    !*          while for  abserr = 0, we test the relative error only.     *
    !*          abserr and relerr are used as put in only when both exceed  *
    !*          four times the machine constant. Otherwise we set them to   *
    !*          this value.                                                 *
    !* fmax     upper limit of calls of function                            *
    !* protnam  Name of a file used for intermediate results. If the pointer*
    !*          is set to zero, we do not use this file.                    *
    !* a,b      end points of the interval, that includes a root            *
    !*                                                                      *
    !* Output parameters:                                                   *
    !* =================                                                    *
    !* abserr\  actual error bounds used                                    *
    !* relerr/                                                              *
    !* b        approximate root                                            *
    !* fb       value of f at the root b                                    *
    !* fanz     number of calls of  func()                                   *
    !*                                                                      *
    !* Error code rc:                                                       *
    !* =============                                                        *
    !* = -2: abserr or relerr is negative, or both are zero, or fmax < 1.   *
    !* = -1: The necessary assumption  f(a) * f(b) < 0  is violated.        *
    !* =  0: Desired accuracy has been reached :                            *
    !*           |xm|  <=  0.5 * (abserr + |b| * relerr).                   *
    !* =  1: b on output is a root with  func(b) = 0.                        *
    !* =  2: Either a or b is a root at input.                              *
    !* =  3: After fmax calls of func the algorithm has not found a root     *
    !* =  4: Error when opening intermediate result file                    *
    !*                                                                      *
    !************************************************************************
    ! arguments
    real, external :: func
    real, intent(in)    :: myabserr, myrelerr
    integer, intent(in) :: fmax
    real, intent(in)    :: x1
    real, intent(in)    :: x2
    character(len=*), intent(in), optional :: protnam

    ! function return variable
    type(outtype_zeroin), intent(out) :: out_zeroin

    ! local variables
    real :: fa         ! Function value func(a)
    real :: fc         ! Function value func(c)                    
    real :: eps        ! minimal value for error bounds abserr and relerr
    real :: tol1       ! auxiliary variable for mixed error bound 
    real :: xm         ! half of the current interval length      
    real :: c          ! value inside [a,b]                       
    real :: d          ! Distance to the nearest approximate root 
    real :: e          ! previous value of d                      
    real :: p
    real :: q
    real :: r
    real :: s
    real :: tmp                ! auxiliary variable to check inclusion    
                               !   f(a) * f(b) < 0                        
    integer fehler             ! error code of this function              

    real    :: abserr, relerr
    real    :: a
    real    :: b
    
    real    :: fb
    integer :: fanz
    integer :: rc

    ! ------------------ initialize variables -----------------------------
    a = x1
    b = x2

    abserr = myabserr
    relerr = myrelerr

    c=ZERO; d=ZERO; e=ZERO; rc=0; fehler=0
    eps = FOUR * XMACHEPS

    fa = func(a)
    fb = func(b)

    fanz = 2

    ! ----------- check   f(a) * f(b) < 0  --------------------------------
    tmp = fa * fb
    if (tmp > ZERO) then
      rc = -1
    else if (tmp == ZERO) then
      rc = 2
    end if
    if (rc.ne.0)  goto 999

    ! ----------- check usability of given error bounds -------------------
    if (abserr < ZERO.or.relerr < ZERO.or.abserr + relerr <= ZERO.or.fmax < 1) then
      rc = -2
      goto 999
    end if

    if (relerr==ZERO.and.abserr < eps) then
      abserr = eps
    else if (abserr==ZERO.and.relerr < eps) then
      relerr = eps
    else
      if (abserr < eps)  abserr = eps
      if (relerr < eps)  relerr = eps
    end if

    if (present(protnam)) open(unit=2,file='tzeroin.log',status='unknown')  
    
    fc=fb
    do while (2>1)                                        ! start iteration

      if (fb * (fc / abs(fc)) > ZERO) then ! no inclusion of a root
        c  = a                              ! between b and c ?
        fc = fa                             ! alter c so that b and c
        e  = b - a                          ! include the root of f
        d  = e
      end if

      if (abs(fc) < abs(fb)) then         ! If fc has the smaller modulus
        a   = b                             ! interchange interval end
        b   = c                             ! points
        c   = a
        fa  = fb
        fb  = fc
        fc  = fa
      end if

      if (present(protnam)) then
        write(2,100)  a, b, c
        write(2,110)  fa, fb, fc
      end if

      tol1 = HALF * (abserr + relerr * abs(b));
      xm   = HALF * (c - b);

      if (present(protnam)) write(2,120)  xm, tol1

      if (abs(xm) <= tol1) then       ! reached desired accuracy ?
        fehler = 0
        goto 999                         ! normal exit
      end if

      if (fb == ZERO) then             ! Is the best approximate root b
                                       ! a precise root of f ?
        fehler = 1
        goto 999                         ! other normal exit, b is a root
      end if

      r = ZERO

      if (abs(e) < tol1.or.abs(fa) <= abs(fb)) then
        e = xm;
        d = e;
        if (present(protnam)) write(2,*) ' Bisection'
      else
        if (a.ne.c) then       ! if a is not equal to c
          q = fa / fc          ! With a, b and c we have 3 points for
          r = fb / fc          ! an inverse quadratic interpolation
          s = fb / fa
          p = s * (TWO * xm * q * (q - r) - (b - a) * (r - ONE))
          q = (q - ONE) * (r - ONE) * (s - ONE)
        else                   ! If a equals  c :
          s = fb / fa          ! Use the secant method or linear
          p = TWO * xm * s     ! interpolation
          q = ONE - s
        end if

        if (p > ZERO) then     ! alter the sign of  p/q for the
          q = -q               ! subsequent division
        else
          p = abs(p)
        end if 

        if (TWO * p  >=  THREE * xm * q - abs(tol1 * q).or.p >= abs(HALF * e * q)) then
          e = xm; d = e
          if (present(protnam)) write(2,*) ' Bisection'
        else
          e = d        ! Compute the quotient p/q for both iterations
          d = p / q    ! which will be used to modify b
          if (present(protnam)) then
            if (r == ZERO) then
              write(2,*) ' Secant method'
            else
              write(2,*) ' Inverse quadratic interpolation'
            end if
          end if
        end if
      end if

      a  = b           ! store the best approximation b and its
      fa = fb          ! function value fb in a and fa

      if (abs(d) > tol1) then                    ! d large enough?
        b = b + d                                 ! add d to b
        if (present(protnam)) write(2,*) ' Difference d from new b: ', d
      else                                        ! d too small?
        b = b + Sign(tol1, xm)                    ! add tol1 to b
        if (present(protnam)) then
          if (xm < ZERO) then
            write(2,*) ' Subtract error bound: d = ', -tol1
          else
            write(2,*) ' Add error bound: d = ', tol1
          end if
        end if
      end if

      fb = func( b )
      fanz=fanz+1                        !      up evaluation counter

      if (present(protnam)) write(2,*) ' b = ',b,' fb= ', fb,' Number of functional evaluations = ', fanz
      if (fanz > fmax) then              ! too many function evaluations?
        fehler = 3
        goto 50
      end if
    end do                                                ! end iteration

    50 close(unit=2)

    rc = fehler

    ! construct function return type
    999 out_zeroin%root  = b
    out_zeroin%froot = fb
    out_zeroin%niter = fanz
    out_zeroin%error = rc

    return

    100 format(' a  = ', F20.14, ' b  = ', F20.14, ' c  = ', F20.14)
    110 format(' fa = ', F20.14, ' fb = ', F20.14, ' fc = ', F20.14)
    120 format(' xm = ', F20.14, ' tol1 = ', F20.14) 

  end function zeroin

end module _findroot_fzeroin
