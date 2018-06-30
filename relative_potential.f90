module relative_potential
  use potentials
  use input_data
  private

  type, public, extends(potential) :: rel_pot
    complex(8), allocatable, dimension(:,:) :: Vrel
    contains
    procedure :: rel_pot_
    procedure :: get_Vrel
    procedure :: make_Vrel
    procedure :: destruct_rel_pot
  end type
  contains

!----------------------------------------------------------------------
  subroutine rel_pot_(this, ip)
    implicit none
    class(rel_pot), intent(out) :: this
    type(inp), intent(in), target :: ip

    call this%potential%potential_(ip)
    allocate(this%Vrel(0:ip%Jmax,ip%rgrid+2))
  end subroutine
!----------------------------------------------------------------------
  elemental subroutine destruct_rel_pot(this)
    implicit none
    class(rel_pot), intent(inout) :: this

    if (allocated(this%Vrel)) deallocate(this%Vrel)
  end subroutine
!----------------------------------------------------------------------
  elemental function get_Vrel(this, J, ir) result(V)
    implicit none
    class(rel_pot), intent(in) :: this
    integer, intent(in) :: J, ir
    complex(8) :: V

    V = this%Vrel(J,ir)

  end function
!----------------------------------------------------------------------
  elemental subroutine make_Vrel(this)
    implicit none
    class(rel_pot), intent(inout) :: this
    integer :: i, J, ist
    real(8) :: r
    complex(8), parameter :: ai=(0.0d0,1.0d0)

    if (this%ip%rmin == 0.0d0) then
      ist = 2
      this%Vrel(:,1) = 0.0d0
    else
      ist = 1
    end if
    do i=ist, this%ip%rgrid+2
      r = this%ip%rmin + dble(i-1) * this%ip%dr
      this%Vrel(:,i) = (Vn(this%ip,r) + ai * Wn(this%ip,r)) + Vc(this%ip,r)
    end do
    do J=0, this%ip%Jmax
      do i=ist, this%ip%rgrid+2
        r = this%ip%rmin + dble(i-1) * this%ip%dr
        this%Vrel(J,i) = this%Vrel(J,i) + Vcnt(this%ip,J,r)
      end do
    end do

    return
  end subroutine
end module

