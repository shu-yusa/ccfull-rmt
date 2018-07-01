!======================================================================!
module fitting
  private
  integer, parameter :: MaxIter = 5000

  type, abstract, public :: fit
    integer :: npara
    integer :: ndata
    integer :: iter=0
    real(8), allocatable :: pms(:)
    real(8), allocatable :: exdata(:,:)
    contains
    procedure(th_cal_proto), deferred, pass :: th_cal
    procedure, private :: fit_1
    procedure, private :: fit_2
    procedure :: read_data => read_data1
    procedure :: set_npara
    procedure :: set_exdata
    procedure :: set_pms
    procedure :: get_npara
    procedure :: get_ndata
    procedure :: get_pms
    procedure :: OROF_method
    procedure, private :: kai2
    procedure, private :: kai22
    generic :: kai => kai2, kai22
    generic :: fit_ => fit_1, fit_2
  end type

  abstract interface
    function th_cal_proto(this, pms) result(th_val)
      import :: fit
      class(fit), intent(inout) :: this
      real(8), intent(in) :: pms(:)
      real(8) :: th_val(this%ndata)
    end function
  end interface

contains

  subroutine fit_1(this, npara, ndata)
    implicit none
    class(fit), intent(out) :: this
    integer, intent(in) :: npara, ndata

    this%npara = npara
    this%ndata = ndata

    allocate(this%pms(npara))
    allocate(this%exdata(3,ndata))

    this%pms = 0.0d0
    this%exdata = 0.0d0
  end subroutine
!**********************************************************************
  subroutine fit_2(this, pms, exdata_x, exdata_y, exdata_z)
    implicit none
    class(fit), intent(out) :: this
    real(8), intent(in), dimension(:) :: pms,exdata_x,exdata_y,exdata_z

    this%npara = size(pms)
    this%ndata = size(exdata_x)  

    if (.not. allocated(this%pms)) allocate(this%pms(this%npara))
    if (.not. allocated(this%exdata)) allocate(this%exdata(3,this%ndata))

    this%pms = pms
    this%exdata(1,:) = exdata_x
    this%exdata(2,:) = exdata_y
    this%exdata(3,:) = exdata_z
  end subroutine
!**********************************************************************
  subroutine read_data1(this, unitn, Fname)
    implicit none
    class(fit), intent(inout), target :: this
    integer, intent(in) :: unitn
    integer :: l, ios
    character(len=*), intent(in) :: Fname

    open(unitn,file=trim(Fname))
      l = 0
      do
        read(unitn,*,iostat=ios)
        if (ios < 0) then
          exit
        else
          l = l + 1
        end if
      end do

      this%ndata = l       ! add g.s. channel

      if (.not. allocated(this%exdata)) then
        allocate(this%exdata(3,this%ndata))
      else if (size(this%exdata,dim=2) /= this%ndata) then
        write(6,*) 'read_data: warning: inconsistent data number'
        deallocate(this%exdata)
        allocate(this%exdata(3,this%ndata))
      end if

      rewind(unitn)

      do l=1, this%ndata
        read(unitn,*) this%exdata(1,l),this%exdata(2,l),this%exdata(3,l)
      end do
    close(unitn)

  end subroutine
!**********************************************************************
  function get_npara(this) result(npara)
    implicit none
    class(fit), intent(in) :: this
    integer :: npara

    npara = this%npara

  end function
!**********************************************************************
  function get_ndata(this) result(ndata)
    implicit none
    class(fit), intent(in) :: this
    integer :: ndata

    ndata = this%ndata

  end function
!**********************************************************************
  subroutine set_npara(this, n)
    use, intrinsic :: iso_fortran_env
    implicit none
    class(fit), intent(out) :: this
    integer, intent(in) :: n

    if (.not. allocated(this%pms)) then
      this%npara = n
      allocate(this%pms(n))
    else
      write(output_unit,*) 'npara already set'
    end if

    end subroutine
!**********************************************************************
  subroutine set_exdata(this, exdata_x, exdata_y, exdata_z)
    implicit none
    class(fit), intent(out) :: this
    real(8), intent(in), dimension(:) :: exdata_x, exdata_y, exdata_z

    if (.not. allocated(this%exdata)) then
      allocate(this%exdata(3,size(exdata_x)))
    end if
    this%exdata(1,:) = exdata_x
    this%exdata(2,:) = exdata_y
    this%exdata(3,:) = exdata_z
  end subroutine
!**********************************************************************
  subroutine set_pms(this, pms)
    implicit none
    class(fit), intent(inout) :: this
    real(8), intent(in) :: pms(:)

    if (.not. allocated(this%pms)) then
      allocate(this%pms(size(pms)))
      this%npara = size(pms)
    else if (size(this%pms,1) /= size(pms,1)) then
      stop "set_params: wrong dimension"
    end if

    this%pms = pms
  end subroutine
!**********************************************************************
  function get_pms(this) result(pms)
    implicit none
    class(fit), intent(in) :: this
    real(8) :: pms(this%npara)

    pms = this%pms

  end function
!**********************************************************************
  function kai2(this) result(kai)
    implicit none
    class(fit), intent(inout) :: this
    integer :: i
    real(8) :: kai, th_val(this%ndata)

    th_val = this%th_cal(this%pms)

    kai = 0.0d0
    do i=1, this%ndata
      kai = kai + ((th_val(i) - this%exdata(2,i)) / this%exdata(3,i)) ** 2
    end do

    kai = kai / dble(this%ndata)

  end function
!**********************************************************************
  function kai22(this, pms) result(kai)
    implicit none
    class(fit), intent(inout) :: this
    integer :: i
    real(8), intent(in) :: pms(this%npara)
    real(8) :: kai, th_val(this%ndata)

    th_val = this%th_cal(pms)

    kai = 0.0d0
    do i=1, this%ndata
      kai = kai + ((th_val(i) - this%exdata(2,i)) / this%exdata(3,i)) ** 2
    end do

    kai = kai / dble(this%ndata)
  end function
!**********************************************************************
  subroutine OROF_method(this, epsr)
    use, intrinsic :: iso_fortran_env
    implicit none
    class(fit), intent(inout) :: this
    integer :: i, j, iter
    real(8), intent(in) :: epsr
    real(8) :: kai0, th_val0(this%ndata), kai2
    character(len=11), parameter :: FM='(1x,a,i5)'

    do i=1, MaxIter
      iter = i
      th_val0 = this%th_cal(this%pms)
      kai0 = 0.0d0
      do j=1, this%ndata
        kai0 = kai0 + ((th_val0(j) - this%exdata(2,j)) / this%exdata(3,j)) ** 2
      end do
      kai0 = kai0 / dble(this%ndata)
      call renew_pms(this, th_val0, kai0)
      kai2 = this%kai2()
      write(output_unit,*) "kai square = ", kai2
      if ((abs(kai2 - kai0)) < epsr * kai0) then
        write(output_unit,FM) 'Iteration :',i
        exit
      end if
    end do
    if (i==MaxIter) print *,  'fitting failed'
    write(output_unit,*) kai2

  contains
!**********************************************************************
    subroutine renew_pms(this, th_val0, kai0)
!     use mkl95_lapack
      use lin_coup
      implicit none
      class(fit), intent(inout) :: this
      real(8), intent(in) :: th_val0(:), kai0
      real(8), target :: dp(this%npara)
      real(8) :: ddF(this%npara,this%npara)
      real(8) :: phi, z
      real(8), pointer :: p_dp(:,:)

      call Deriv(this, th_val0, dp, ddF)
      dp = - dp
      p_dp(1:this%npara,1:1) => dp
!     call sysv(ddF, p_dp, "L")
      call lin_sy(ddF, p_dp)
      phi = maxval(dp / this%pms)
      if (phi < 0.1d0) then
         z = 1.0d0
      else if (phi < 0.3d0) then
         z = 0.5d0
      else if (phi < 0.5d0) then
         z = 0.2d0
      else if (phi < 1.0d0) then
         z = 0.1d0
      else if (phi >= 1.0d0) then
         z = 0.05d0 / phi
      end if
      call Parab_Approx(this, kai0, dp, z)

    end subroutine
!**********************************************************************
    subroutine Deriv(this, th_val0, dF, ddF)
      implicit none
      class(fit), intent(inout) :: this
      integer :: i, j, k
      real(8), intent(out) :: dF(this%npara), ddF(this%npara,this%npara)
      real(8), dimension(this%ndata) :: kai, th_val, th_val0
      real(8), parameter :: del_p = 0.01d0
      real(8) :: pms(this%npara), S, rf(this%npara,this%ndata)

      pms = this%pms
      do i=1, this%npara
        pms(i) = (1.0d0 + del_p) * this%pms(i)
        th_val = this%th_cal(pms)
        do k=1, this%ndata
          kai(k) = (th_val0(k) - this%exdata(2,k)) / this%exdata(3,k) ** 2
          rf(i,k) = th_val(k) - th_val0(k)
        end do
        dF(i) = 2.0d2 * dot_product(kai,rf(i,:)) / this%pms(i)
        do j=1, i
          S = 0.0d0
          do k=1, this%ndata
            S = S + rf(i,k) * rf(j, k) / this%exdata(3,k) ** 2
          end do
          ddF(i,j) = 2.0d4 * S / (this%pms(j) * this%pms(i))
        end do
        pms(i) = this%pms(i)
      end do

    end subroutine
!**********************************************************************
    subroutine Parab_Approx(this, y0, dp, z)
      implicit none
      class(fit), intent(inout) :: this
      real(8), value :: z
      real(8), intent(in) :: dp(this%npara), y0
      real(8), dimension(this%npara) :: pa, pb
      real(8) :: ya, yb, z0, za, zb, zm, ym

      pa = this%pms + z * dp
      ya = this%kai22(pa)
      if (ya <= y0) then
        za = z
        zb = 2.0d0 * z
        pb = this%pms + zb * dp   
        yb = this%kai22(pb) 
      else
        za = 0.5d0 * z
        zb = z
        pb = pa
        yb = ya
        pa = this%pms + za * dp
        ya = this%kai22(pa)
      end if
      if ((yb - y0) * za + (y0 - ya) * zb > 0.0d0) then
        zm = 0.5d0*za * (3.0d0*y0 - 4.0d0*ya + yb) / (y0 - 2.0d0*ya+yb)
        if (zm <= - zb) then
          z = - zb
        else if (zm <= zb) then
          z = zm
        else if (zm > zb) then
          z = 2.0d0 * zb
        end if
      else
        if (yb < y0) then
          z = 3.0d0 * za
        else if (yb >= y0) then
          z = - za
        end if
      end if
      this%pms = this%pms + z * dp

      return
    end subroutine
  end subroutine

end module

