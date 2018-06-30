!======================================================================!
      module gammaln_fct
      implicit none
      public :: gammaln, arg_gamma
      interface gammaln
        module procedure gammaln_d, gammaln_c
      end interface

      contains
!**********************************************************************!
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
!**********************************************************************!
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
!**********************************************************************!
      function arg_gamma(z) result(t)
      implicit none
      real(8) :: t
      complex(8), intent(in) :: z
      complex(8) :: gam
      
!     gam = exp(gammaln_c(z))
!     t = atan(aimag(gam)/dble(gam))
      t = aimag(gammaln_c(z))
      
      end function
!**********************************************************************!
      end module
