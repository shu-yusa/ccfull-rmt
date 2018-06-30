module special_fct
  real(8), parameter, private :: PI=3.141592653589793d0

  public :: gammaln, arg_gamma, jL, PL, YL0, dYL0
  interface gammaln
    module procedure gammaln_d, gammaln_c
  end interface

  contains
!******************************************************************!
  function gammaln_d(x) result(g)
    implicit none
    real(8), intent(in) :: x
    real(8), dimension(6) :: coef = (/76.18009172947146d0, &
                -86.50532032941677d0, 24.01409824083091d0, &
                -1.231739572450155d0, 0.1208650973866179d-2, &
                -0.5395239384953d-5/)
    real(8), parameter :: stp = 2.5066282746310005d0
    real(8) :: tmp, g

    if (x < 0) stop 'negative argument in gammaln'
    tmp = x + 5.5d0
    tmp = (x + 0.5d0) * log(tmp) - tmp
    g = tmp + log(stp/x*(1.000000000190015d0 &
       + coef(1)/(x+1.0d0) + coef(2)/(x+2.0d0) + coef(3)/(x+3.0d0) &
       + coef(4)/(x+4.0d0) + coef(5)/(x+5.0d0) + coef(6)/(x+6.0d0)))

  end function
!******************************************************************!
  function gammaln_c(z) result(g)
    implicit none
    complex(8), intent(in) :: z
    complex(8) :: tmp, g
    real(8), dimension(6) :: coef = (/76.18009172947146d0, &
                -86.50532032941677d0, 24.01409824083091d0, &
                -1.231739572450155d0, 0.1208650973866179d-2, &
                -0.5395239384953d-5/)
    real(8), parameter :: stp = 2.5066282746310005d0

    if (dble(z) < 0) stop 'negative argument in gammaln'
    tmp = z + 5.5d0
    tmp = (z + 0.5d0) * log(tmp) - tmp
    g = tmp + log(stp/z*(1.000000000190015d0 &
       + coef(1)/(z+1.0d0) + coef(2)/(z+2.0d0) + coef(3)/(z+3.0d0) &
       + coef(4)/(z+4.0d0) + coef(5)/(z+5.0d0) + coef(6)/(z+6.0d0)))

  end function
!******************************************************************!
  function arg_gamma(z) result(t)
    implicit none
    real(8) :: t
    complex(8), intent(in) :: z
    complex(8) :: gam

    gam = exp(gammaln_c(z))
    t = atan(aimag(gam)/dble(gam))
!   t = aimag(gammaln_c(z))
  
  end function
!**********************************************************************!
  elemental function  jL(L, x)  result(f)
!----------------------------------------------------------------------
!   spherical bessel function
!----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: L
    integer :: i
    real(8), intent(in) :: x
    real(8) :: f, j0, j1

    if (x == 0.0d0) then
      if (L == 0) then
        f = 1.0d0
      else
        f = 0.0d0
      end if
      return
    end if

    if (L == 0) then
      f = sin(x) / x
      return
    else if (L == 1) then
      f = (sin(x) - x * cos(x)) / (x * x)
      return
    end if

    j0 = sin(x) / x
    j1 = (sin(x) - x * cos(x)) / (x * x)
    do i=2, L
      f = dble(2 * i - 1) * j1 / x - j0
      j0 = j1
      j1 = f
    end do
  end function
!----------------------------------------------------------------------
  elemental function PL(L, x) result(P)
!----------------------------------------------------------------------
!   Legendre polynomial
!----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: L
    integer :: i
    real(8), intent(in) :: x
    real(8) :: P0, P1, P

    if (L == 0) then
      P = 1.0d0
      return
    else if (L == 1) then
      P = x
      return
    end if

    P0 = 1.0d0
    P1 = x
    do i=2, L
      P = (dble(2 * i - 1) * x * P1 - dble(i - 1) * P0) / dble(i)
      P0 = P1
      P1 = P
    end do
  end function
!----------------------------------------------------------------------
  elemental function YL0(L, t) result(f)
    implicit none
    integer, intent(in) :: L
    real(8), intent(in) :: t
    real(8) :: f

    f = sqrt(dble(2*L+1)/(4.0d0*PI)) * PL(L, cos(t))
  end function
!----------------------------------------------------------------------
  elemental function dYL0(L, t) result(dY)
    implicit none
    integer, intent(in) :: L
    real(8), intent(in) :: t
    real(8) :: dY, cost

    if (t == 0.0d0 .or. t == PI) then
      dY = 0.0d0
      return
    end if
    
    cost = cos(t)
    dY = dble(L) / sin(t) * sqrt(dble(2*L+1)/(4.0d0*PI)) &
       * (cost * PL(L,cost) - PL(L-1,cost))
  end function
!----------------------------------------------------------------------
end module

