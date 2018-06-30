module potentials
  use input_data

  type, public ::  potential
    type(inp), pointer :: ip
    real(8) :: rb, rp
    real(8) :: Vmin, Vmax
    contains
    procedure :: potential_
    procedure :: find_rmin
  end type

  contains

!----------------------------------------------------------------------
  subroutine potential_(this, ip)
    implicit none
    class(potential), intent(out) :: this
    type(inp), intent(in), target :: ip

    this%ip => ip
  end subroutine
!----------------------------------------------------------------------
  elemental function Vc(ip, r) result(V)
    use global_constant, only : e2
    implicit none
    real(8), intent(in) :: r
    real(8) :: V
    type(inp), intent(in) :: ip

    if (r > ip%Rc) then
      V = ip%Zp * ip%Zt * e2 / r
    else
      V = 0.5d0*ip%Zp*ip%Zt*e2 * (3.0d0*ip%Rc*ip%Rc - r*r) / ip%Rc**3
    end if

    return
  end function
!----------------------------------------------------------------------!
  elemental function Vn(ip, r) result(V)
    implicit none
    real(8), intent(in) :: r
    real(8) :: V
    type(inp), intent(in) :: ip

    if (r >= ip%rcut) then
      V = 0.0d0
      return
    end if
    V = - ip%V0 / (1.0d0 + exp((r - ip%Rn) / ip%a))

    return
  end function
!----------------------------------------------------------------------!
  elemental function Wn(ip, r) result(W)
    implicit none
    real(8), intent(in) :: r
    real(8) :: W
    type(inp), intent(in) :: ip

    if (r >= ip%rcut) then
      W = 0.0d0
      return
    end if
    W = - ip%W0 / (1.0d0 + exp((r - ip%Rw) / ip%aw))

    return
  end function
!----------------------------------------------------------------------!
  elemental function Vcnt(ip, J, r) result(V)
    use global_constant, only : hbar
    implicit none
    integer, intent(in) :: J
    real(8), intent(in) :: r
    real(8) :: V
    type(inp), intent(in) :: ip

    V = 0.5d0 * hbar * hbar / ip%rmass * dble(J * (J + 1)) / (r * r)

    return
  end function
!----------------------------------------------------------------------
  function find_rmin(this, J) result(pocket)
    implicit none
    class(potential), intent(inout) :: this
    logical :: pocket
    integer, intent(in) :: J
    integer :: i
    real(8), parameter :: dr = 0.1d0
    real(8), parameter :: r0 = 3.0d0, rmax = 20.0d0
    real(8), parameter :: epsr = 1.0d-15, epsa = 1.0d-200
    real(8) :: V1, V2, r1, r2, rh

    r1 = r0
    r2 = r1 + dr
    do while(r1 < this%ip%rmax)
      if (sign(1.0d0,dV_dr(J,r1)) * sign(1.0d0,dV_dr(J,r2)) > 0.0d0) then
        r1 = r2
        r2 = r1 + dr
      else
        exit
      end if
    end do
    if (r1 >= this%ip%rmax) then
      this%rp = 0.0d0
      this%rb = 0.0d0
      pocket = .false.
    end if


    rh = 0.5d0 * (r1 + r2)
    V2 = 0.0d0
    i = 0
    pocket = .true.
    do
      V1 = dV_dr(J,rh)
      if (abs((V1 - V2)) < abs(V1) * epsr + epsa) exit
      if (V1 < 0.0d0) then
        r1 = rh
      else
        r2 = rh
      end if
      rh = 0.5d0 * (r1 + r2)
      V2 = V1
      if (i == 50 .or. rh > rmax) then
        pocket = .false.
        return
      end if
      i = i + 1
    end do
    this%rp = rh
    this%Vmin = Vn(this%ip, this%rp) + Vc(this%ip, this%rp)

    r1 = rh + 1.0d-14
    r2 = r1 + dr
    do 
      if (dV_dr(J,r1) * dV_dr(J,r2) >= 0.0d0) then
        r1 = r2
        r2 = r1 + dr
      else
        exit
      end if
    end do

    rh = 0.5d0 * (r1 + r2)
    V2 = 0.0d0
    i = 0
    do
      V1 = dV_dr(J,rh)
      if (abs((V1 - V2)) < abs(V1) * epsr + epsa) exit
      if (V1 > 0.0d0) then
        r1 = rh
      else
        r2 = rh
      end if
      rh = 0.5d0 * (r1 + r2)
      V2 = V1
      if (i == 60) then
        exit
      end if
      i = i + 1
    end do
    this%rb = rh
    this%Vmax = Vn(this%ip, this%rb) + Vc(this%ip, this%rb)

    contains
!----------------------------------------------------------------------
    elemental function dV_dr(J, r) result(dV)
      use global_constant, only : e2, hbar
      implicit none
      integer, intent(in) :: J
      real(8), intent(in) :: r
      real(8) :: dV, f

      f = exp((r - this%ip%Rn) / this%ip%a)
      dV = this%ip%V0*f / (this%ip%a*(1.0d0+f)**2) - this%ip%Zp*this%ip%Zt*e2 / (r*r)  &
           - dble(J * (J + 1)) * hbar * hbar / (this%ip%rmass*r**3)
    end function
  end function
end module

