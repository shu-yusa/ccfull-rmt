      module gaus_int
      implicit none
      public :: int_herm, int_lagu, int_legen
      contains
!**********************************************************************!
      subroutine int_herm(n, m, w, fp, fn, S)
      implicit none
      integer, intent(in) :: n, m
      integer :: j
      real(8), intent(in) :: w(m), fp(m), fn(m)
      real(8), intent(out) :: S

      if (mod(n,2) == 1) then
        S = w(m) * fp(m)
      else
        S = 2.0d0 * w(m) * fp(m)
      end if

      do j=1, m-1
        S = S + w(j) * (fp(j) + fn(j))
      end do

      return
      end subroutine
!**********************************************************************!
      subroutine int_lagu(n, w, f, S)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in), dimension(n) :: w, f
      real(8), intent(out) :: S

      S = sum(w(1:n) * f(1:n))

      return
      end subroutine
!**********************************************************************!
      subroutine int_legen(n, m, a, b, w, fp, fn, S)
      implicit none
      integer, intent(in) :: n, m
      integer :: j
      real(8), intent(in) :: a, b, w(m), fp(m), fn(m)
      real(8), intent(out) :: S

      if (mod(n,2) == 1) then
        S = w(m) * fp(m)
      else
        S = w(m) * (fp(m) + fn(m))
      end if

      do j=1, m-1
        S = S + w(j) * (fp(j) + fn(j))
      end do

      S = S * 0.5d0 * (b - a)

      return
      end subroutine
!**********************************************************************!
      end module
