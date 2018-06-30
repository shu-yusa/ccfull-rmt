      module Random_num
      public :: gaus_ran_num, unifm_ran_num
      contains
!======================================================================!
      function gaus_ran_num(idum) result(y)
!----------------------------------------------------------------------!
!     This subroutine generates a standard normally distributed random !
!     number.                                                          !
!     idum : Seed                                                      !
!     y    : Random number.                                            !
!----------------------------------------------------------------------!
      implicit none
      logical, save :: gaus_store=.false.
      integer, intent(inout) :: idum
      real(8) :: y
      real(8), parameter :: PI=3.141592653589793d0
      real(8), save :: ystore
      real(8) :: x1, x2, w

      if (gaus_store) then
        y = ystore
        gaus_store = .false.
      else
        x1 = unifm_ran_num(idum)
        x2 = unifm_ran_num(idum)
        w = sqrt(-2.0d0 * log(x1))
        y = w * cos(2.0d0 * PI * x2)
        ystore = w * sin(2.0d0 * PI * x2)
        gaus_store = .true.
      end if

      return
      end function
!======================================================================!
      function unifm_ran_num(idum) result(r)
!----------------------------------------------------------------------!
!     "Minimal random number generator of Park and Miller combined     !
!     with a Marsaglia shift sequence. Returns a uniform random        !
!     deviate between 0.0 and 1.0 (exclusive of the endpoint values).  !
!     This fully portable, scalar generator has the "traditional"      !
!     (not Fortran 90) calling squence with a random deviate as the    !
!     returned function value: call with idum a negative integer to    !
!     initialize; thereafter, do not alter idum except to reinitialize.!
!     The period of this generator is about 3.1 X 10^18.               !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(inout) :: idum
      integer, parameter :: IA=16807,IM=2147483647,IQ=12773,IR=2836
      integer, save :: ix=-1, iy=-1, k
      real(8) :: r
      real(8), save :: am

      if (idum <= 0 .or. iy < 0) then
        am = nearest(1.0d0, -1.0d0) / IM
        iy = ior(ieor(888889999, abs(idum)), 1)
        ix = ieor(777755555, abs(idum))
        idum = abs(idum) + 1
      end if
      ix = ieor(ix, ishft(ix,13))
      ix = ieor(ix, ishft(ix,-17))
      ix = ieor(ix, ishft(ix,5))
      k = iy / IQ
      iy = IA * (iy - k * IQ) - IR * k
      if (iy < 0) iy = iy + IM
      r = am * ior(iand(IM, ieor(ix,iy)), 1)

      end function
!======================================================================!
      end module
