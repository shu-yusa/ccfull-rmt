module level_density
  use global_constant
  use print_array
  use angular_coup

  type, public :: level_dens
    integer :: Amass
    real(8) :: a
    real(8) :: eta
    real(8) :: Epair

    contains

    procedure :: level_dens_
    procedure :: sigma
    procedure :: rho
    procedure :: rho_J
    procedure :: rho2
    procedure :: rho_J2
  end type

  contains

  subroutine level_dens_(this, Amass, a, eta, Epair)
    implicit none
    class(level_dens), intent(out) :: this
    real(8), intent(in) :: a, eta, Epair
    integer, intent(in) :: Amass

    this%a = a
    this%eta = eta
    this%Epair = Epair
    this%Amass = Amass

  end subroutine

  function sigma(this, E) result(s)
    implicit none
    class(level_dens), intent(in) :: this
    real(8), intent(in) :: E
    real(8) :: s

    s = sqrt(0.0888d0 * this%Amass**(2.0d0/3.0d0) * sqrt(this%a*(E-this%Epair)))

  end function

  function rho(this, E) result(r)
    implicit none
    class(level_dens), intent(in) :: this
    real(8), intent(in) :: E
    real(8) :: r, U

    U = E - this%Epair
    r = this%eta * exp(2.0d0*sqrt(this%a*U))    &
      / (12.0d0*sqrt(2.0d0)*this%a**(0.25d0)*U**(1.25d0)*sigma(this,E))

  end function

  function rho2(this, E) result(r)
    implicit none
    class(level_dens), intent(in) :: this
    real(8), intent(in) :: E
    real(8) :: r, E0
! 92Zr
    real(8), parameter :: l = 0.659459d0
    real(8), parameter :: b = 540.705d0
    real(8), parameter :: c = -737.154d0
    real(8), parameter :: d = 366.72d0
    real(8), parameter :: h = -87.0023d0
    real(8), parameter :: f = 10.0678d0
    real(8), parameter :: g = -0.458924d0

! 90Zr
!   real(8), parameter :: b = 182.509d0
!   real(8), parameter :: c = -286.525d0
!   real(8), parameter :: d = 119.671d0
!   real(8), parameter :: h = -22.649d0
!   real(8), parameter :: f = 2.05739d0
!   real(8), parameter :: g = -0.0727757d0

! 208Pb
!   real(8), parameter :: b = 6968.88d0
!   real(8), parameter :: c = -2612.28d0
!   real(8), parameter :: d = 497.522d0
!   real(8), parameter :: h = -49.5888d0
!   real(8), parameter :: f = 2.34719d0
!   real(8), parameter :: g = -0.0363248d0

     if (E < 3.0d0) then
       E0 = 3.0d0
       r = ((((6.0d0*g*E0 + 5.0d0*f)*E0 + 4.0d0*h)*E0 + 3.0d0*d)*E0 + 2.0d0*c)*E0+b
     else if (E > 4.0d0) then
       E0 = 4.0d0
       r = ((((6.0d0*g*E0 + 5.0d0*f)*E0 + 4.0d0*h)*E0 + 3.0d0*d)*E0 + 2.0d0*c)*E0+b
     else
       r = ((((6.0d0*g*E + 5.0d0*f)*E + 4.0d0*h)*E + 3.0d0*d)*E + 2.0d0*c)*E+b
     end if
!   write(6,*) r

  end function

  function rho_J(this, E, J) result(r)
    implicit none
    class(level_dens), intent(in) :: this
    integer, intent(in) :: J
    real(8), intent(in) :: E
    real(8) :: r, U

    U = E - this%Epair
    if (E < 2.5d0) then
      U = 2.5d0
    else
      U = E
    end if
    r = rho(this,U) * dble(2*J+1)*exp(-(dble(J)+0.5d0)**2/(2.0d0*sigma(this,U)**2)) &
      / (2.0d0 * sigma(this,U)**2)
  end function

  function rho_J2(this, E, J) result(r)
    implicit none
    class(level_dens), intent(in) :: this
    integer, intent(in) :: J
    real(8), intent(in)  :: E
    real(8) :: r, U

    U = E - this%Epair
    if (E < 2.5d0) then
      U = 2.5d0
    else
      U = E
    end if
    r = rho2(this, E) * dble(2*J+1)*exp(-(dble(J)+0.5d0)**2/(2.0d0*sigma(this,U)**2)) &
      / (2.0d0 * sigma(this,U)**2)
  end function

end module
