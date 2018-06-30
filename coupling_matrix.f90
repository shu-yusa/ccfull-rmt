module coupling_matrix
  use input_data, only : inp
  use potentials
  use print_array
  use global_constant, only : PI
  use level_density
  use random_num

  integer, parameter :: idum0=-123456789

  type, public, extends(potential) :: coup_mat
    type(level_dens) :: ld
    real(8), allocatable, dimension(:) :: e_n, rhoa
    real(8), allocatable, dimension(:,:) :: Vcp_linear
    real(8), allocatable, dimension(:,:,:) :: wn
    real(8), allocatable :: w_l(:)
    integer, allocatable, dimension(:) :: lambda
    integer :: lambda_min, lambda_max, idum
    integer :: ist, ied
    real(8) :: rho0
    real(8) :: alevel
    real(8) :: alpha
    real(8) :: sigma
    real(8) :: Delta
    real(8) :: w_2
    real(8) :: sqrt_dr
    contains
    procedure :: coup_mat_
    procedure :: set_rmt_params
    procedure :: make_Vcp
    procedure :: destruct_coup_mat
    procedure :: white_noise
    procedure :: record_strength_dist
    procedure :: getVcp
    procedure :: noise_regenerate
    procedure :: h
  end type

  contains
!----------------------------------------------------------------------
  subroutine coup_mat_(this, ip, Fname, Fname2, nprocs, myrank, idum)
    implicit none
    class(coup_mat), intent(out) :: this
    type(inp), target, intent(in) :: ip
    integer, intent(in) :: idum, nprocs, myrank
    integer :: n, i, Amass, step, rem
    integer :: Fnum = 7
    integer :: lambda(ip%Nch_t_noncoll)
    real(8), dimension(ip%Nch_t_noncoll) :: beta, eps_t
    real(8) :: eps
    real(8) :: e, a, eta, Epair, w_l
    character(len=*), intent(in) :: Fname, Fname2

    call this%potential%potential_(ip)
    n = ip%Nch * (ip%Nch + 1) / 2
    allocate(this%Vcp_linear(n,ip%rgrid+2))
    allocate(this%e_n(ip%Nch))
    allocate(this%lambda(max(4,ip%Nch_t)))
    n = ip%Nch_p * (ip%Nch_p + 1) / 2
    allocate(this%wn(ip%Nch_t_coll+1:ip%Nch_t,n,ip%rgrid+2))

    this%idum = idum

    Amass = 208
    a = 10.33d0
    eta = 1.15d0
    Epair = 2.226d0
    this%sqrt_dr = sqrt(this%ip%dr)

    call this%ld%level_dens_(Amass, a, eta, Epair)

    open(Fnum, file=Fname, action="read", status="old")
    read(Fnum,*) this%Delta, this%sigma, this%alpha
    read(Fnum,*) this%lambda_min, this%lambda_max
    read(Fnum,*) w_l
    close(Fnum)
    allocate(this%w_l(this%lambda_min:this%lambda_max))
    this%w_l = w_l



! noncoll levels input
    call read_noncoll_data(this, Fname2, eps_t, lambda, beta) 
! spin
    this%lambda(1) = 0
    this%lambda(2) = this%ip%lambda_t
    if (this%ip%Nphonon_t == 2) then
      this%lambda(3) = 2
      if (this%ip%Nphonon_t2 == 1) then
        this%lambda(4) = this%ip%lambda_t2
      end if
    else if (this%ip%Nphonon_t2 == 1) then
      this%lambda(3) = this%ip%lambda_t2
    end if
    this%lambda(this%ip%Nch_t_coll+1:) = lambda

! excitation energy
    if (this%ip%coup_t == - 1 .and. this%ip%coup_p == - 1) then
      this%e_n = 0.0d0
    else if (this%ip%coup_t == - 1) then
      do i=1, this%ip%Nch
        this%e_n(i) = eps(this,i)
      end do
    else if (this%ip%coup_p == - 1) then
      do i=1, this%ip%Nch_t_coll
        this%e_n(i) = eps(this,i)
      end do
      this%e_n(this%ip%Nch_t_coll+1:) = eps_t
    else
      call eps_pro_tar(this, eps_t, this%e_n)
    end if


! process assinged
    if (nprocs == 1) then
      this%ist = 1
      this%ied = this%ip%rgrid+2
    else
!     step = nint(this%ip%rmax / this%ip%dr)
      step = nint(this%ip%rcut / this%ip%dr)
      rem = mod(step, nprocs-1)
      if (rem == 0 .or. myrank < rem) then
        this%ist = ((step-1)/ (nprocs-1)+1) * myrank + 1
        this%ied = ((step-1)/ (nprocs-1)+1) * (myrank + 1)
      else
        this%ist = ((step-1)/ (nprocs-1)) * myrank + 1 + rem 
        this%ied = ((step-1)/ (nprocs-1)) * (myrank + 1) + rem
      end if
      if (myrank == nprocs-1) then
        this%ist = step+1
        this%ied = this%ip%rgrid+2
      end if
    end if

    contains
      function count_lines(Fname) result(n)
      implicit none
      character(len=*), intent(in) :: Fname
      integer :: n, ios

      open(4723, file=trim(Fname), action="read")

      n = 0
      do
        read(4723,*, iostat=ios)
        if (ios < 0) exit
        n = n + 1
      end do


      close(4723)

      end function
  end subroutine
!----------------------------------------------------------------------
  subroutine noise_regenerate(this, idum)
    implicit none
    class(coup_mat), intent(inout) :: this
    integer, intent(in), optional :: idum
    integer :: n

    if (present(idum)) this%idum = idum
    if (.not. allocated(this%wn)) then
      n = this%ip%Nch_p * (this%ip%Nch_p + 1) / 2
      allocate(this%wn(this%ip%Nch_t_coll+1:this%ip%Nch_t,n,this%ip%rgrid+2))
    end if
    call white_noise(this)

  end subroutine
!----------------------------------------------------------------------
  elemental subroutine set_rmt_params(this, w_l, Delta, sigma, alpha)
    implicit none
    class(coup_mat), intent(out) :: this
    real(8), intent(in) :: Delta, sigma, alpha, w_l

    this%w_l = w_l
    this%Delta = Delta
    this%sigma = sigma
    this%alpha = alpha
  end subroutine
!----------------------------------------------------------------------
  elemental subroutine destruct_coup_mat(this)
    implicit none
    class(coup_mat), intent(inout) :: this

    if (allocated(this%Vcp_linear)) deallocate(this%Vcp_linear)
  end subroutine
!----------------------------------------------------------------------
  function getVcp(this, ir) result(V)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, nst, nlen
    real(8) :: V(this%ip%Nch,this%ip%Nch)

    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      V(n:this%ip%Nch,n) = this%Vcp_linear(nst:nst-1+nlen,ir)
      V(n,2:this%ip%Nch) = V(2:this%ip%Nch,n)
      nst = nst + nlen
    end do

    return
  end function
!----------------------------------------------------------------------
  subroutine make_Vcp(this)
    implicit none
    class(coup_mat), intent(inout) :: this
!   integer, intent(in) :: nprocs, myrank
    integer :: i, n, m, ist, ied, step, rem
    integer :: idum = -123456789
    real(8) :: r
    real(8), allocatable :: Vnm(:,:)
    character(len=40), parameter :: FM='(1x,10f8.3)'

!   call record_strength_dist(this, "strength_dist/strength.dat")
!   stop

!   open(7,file="Vnm_coll.dat")
!   open(8,file="Vnm.dat")
!   open(9,file="Vnm2.dat")
!   allocate(Vnm(this%ip%Nch, this%ip%Nch))

! coupling matrix element
    do i=this%ist, this%ied
!   do i=1, this%ip%rgrid+2
      r = this%ip%rmin + dble(i-1) *  this%ip%dr
!     write(6,*) "r1 = ",r
      if (this%ip%coup_t == - 1 .and. this%ip%coup_p == - 1) then
        this%Vcp_linear(:,i) = 0.0d0
      else if (this%ip%coup_t == - 1 .or. this%ip%coup_p == - 1) then
        call Vcoup(this, r, i)
      else
        call Vcoup_pro_tar(this, r, i)
      end if
!     Vnm = this%getVcp(i)
!     write(7,FM) r, Vnm(2,1),Vnm(3,1)
!     write(8,FM) r, Vnm(4,1),Vnm(5,1),Vnm(6,1),Vnm(7,1),Vnm(8,1),Vnm(9,1),Vnm(10,1),Vnm(11,1)
!     write(9,FM) r, Vnm(12,1),Vnm(13,1),Vnm(14,1),Vnm(15,1),Vnm(16,1),Vnm(17,1),Vnm(18,1),Vnm(19,1)
!     call print_mat(Vnm(:8,:8))
!     if (r > this%ip%rcut) exit
    end do
!   close(7)
!   close(8)
!   close(9)
!   stop
    if (allocated(this%wn)) deallocate(this%wn)


    return
  end subroutine
!----------------------------------------------------------------------!
  subroutine read_noncoll_data(this, Fname, eps_t, lambda, beta)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(out) :: lambda(this%ip%Nch_t_noncoll)
    integer :: i, ios
    integer, parameter :: Nfile = 7
    character(len=*), intent(in) :: Fname
    real(8), dimension(this%ip%Nch_t_noncoll), intent(out) :: eps_t, beta

    open(Nfile,file=Fname,status="old",action="read",iostat=ios)
      do i=1, this%ip%Nch_t_noncoll
!       read(Nfile,*)  lambda(i), eps_t(i), beta(i)
        read(Nfile,*)  beta(i), eps_t(i), lambda(i)
      end do
      eps_t = eps_t / 1000.0d0
    close(Nfile)

  end subroutine
!----------------------------------------------------------------------
  subroutine Vcoup(this, r, ir)
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: n, m, nst, nlen
    integer, intent(in) :: ir
    real(8), intent(in) :: r
    real(8), allocatable, dimension(:,:) :: Vn_cp, Vc_cp
    real(8) :: V

    allocate(Vn_cp(this%ip%Nch,this%ip%Nch))
    allocate(Vc_cp(this%ip%Nch,this%ip%Nch))


    if (this%ip%coup == 0) then
      call Vn_vib(this, ir, r, Vn_cp)
      call Vc_vib(this, r, Vc_cp)
    else if (this%ip%coup == 1) then
      call Vn_rot(this, r, Vn_cp)
      call Vc_rot(this, r, Vc_cp)
    end if
    do n=2, this%ip%Nch_t
      do m=max(this%ip%Nch_t_coll+1,n), this%ip%Nch_t
        Vn_cp(m,n) = 0.0d0
      end do
    end do

    V = Vn(this%ip,r)
    nst = 1
    if (this%ip%coup_t /= -1) then
      do n=1, this%ip%Nch_t_coll
        nlen = this%ip%Nch - n + 1
        this%Vcp_linear(nst,ir) = (Vn_cp(n,n) - V) + Vc_cp(n,n)
        nst = nst + nlen
      end do
    else if (this%ip%coup_p /= -1) then
      do n=1, this%ip%Nch
        nlen = this%ip%Nch - n + 1
        this%Vcp_linear(nst,ir) = (Vn_cp(n,n) - V) + Vc_cp(n,n)
        nst = nst + nlen
      end do
    end if

    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      this%Vcp_linear(nst+1:nst-1+nlen,ir) = Vn_cp(n+1:this%ip%Nch,n) + Vc_cp(n+1:this%ip%Nch,n)
      nst = nst + nlen
    end do

    deallocate(Vn_cp, Vc_cp)
    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vn_rot(this, r, Vpot)
    use mkl95_lapack
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: i, n, m, I1, I2
    real(8), intent(in) :: r
    real(8), intent(out) :: Vpot(this%ip%Nch,this%ip%Nch)
    real(8), parameter :: PI = 3.141592653589793d0
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = sqrt(9.0d0/(4.0d0*PI))
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), dimension(this%ip%Nch) :: Vn
    real(8) :: wg1, wg2, w1, w2, c, f

    if (.not. allocated(O)) then
      allocate(O(this%ip%Nch,this%ip%Nch), Oa(this%ip%Nch))
      w1 = CONST1 * this%ip%beta2 * this%ip%Rm
      w2 = CONST2 * this%ip%beta4 * this%ip%Rm

      if (this%ip%coup_p == 1) then
        do m=1, this%ip%Nrot_p+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%Nrot_p+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            O(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
          end do
        end do

!  2nd mode in projectile
        c = this%ip%betaL_p2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch
          do n=this%ip%Nrot_p+2, this%ip%Nch
            if (m == 1 .and. n == this%ip%Nrot_p+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nrot_p+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nrot_p))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

      else if (this%ip%coup_t == 1) then
        do m=1, this%ip%nrot_t+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%nrot_t+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            O(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
          end do
        end do

!  2nd mode in target
        c = this%ip%betaL_t2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch
          do n=this%ip%nrot_t+2, this%ip%Nch
            if (m == 1 .and. n == this%ip%nrot_t+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%nrot_t+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%nrot_t))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do
      end if
      call syev(O, Oa, 'V','L')
    end if

!   Vpot = 0.0d0
!   f = exp((r - this%ip%Rn)/this%ip%a)
!   Vpot = -O * this%ip%V0/this%ip%a * f / (1.0d0+f)**2
!   forall (n=1:this%ip%Nch)
!     Vpot(n,n) = -this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn)/this%ip%a))
!   end forall


    Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))

    Vpot = 0.0d0
    do m=1, this%ip%Nch
      do n=m, this%ip%Nch
        do i=1, this%ip%Nch
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do

!   call print_mat(Vpot)

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vc_rot(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I1, I2
    real(8), intent(in) :: r
    real(8), intent(out) :: V(this%ip%Nch,this%ip%Nch)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV  = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, x2, wg1, wg2, w1, w2, ZZe2

    ZZe2 = this%ip%Zp * this%ip%Zt * e2
    w1 = 0.6d0 * ZZe2 * CONST1 * this%ip%beta2 * (1.0d0 + CONST2 * this%ip%beta2)
    w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4 + CONST3 * this%ip%beta2 * this%ip%beta2)
    if (r > this%ip%Rm) then
      x = this%ip%Rm / r
      x2 = x * x / r
    else
      x = r / this%ip%Rm
      x2 = x * x / this%ip%Rm
    end if

    if (this%ip%coup_p == 1) then
      do m=1, this%ip%Nrot_p+1
        I2 = 2 * (m - 1)
        do n=m, this%ip%Nrot_p+1
          I1 = 2 * (n - 1)
          call wig3j(I1, 2, I2, 0, 0, 0, wg1)
          call wig3j(I1, 4, I2, 0, 0, 0, wg2)
          V(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                        * sqrt(dble((2*I1+1) * (2*I2+1)))
        end do
      end do

!   The second mode of phonon excitation
      w1 = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w1 * (this%ip%Rm / r) ** this%ip%lambda_p2 / r
      else
        x = w1 * (r / this%ip%Rm) ** this%ip%lambda_p2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch
        do n=this%ip%nrot_p+2, this%ip%Nch
          if (m == 1 .and. n == this%ip%nrot_p+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%nrot_p+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%nrot_p))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

    else if (this%ip%coup_t == 1) then
      do m=1, this%ip%nrot_t+1
        I2 = 2 * (m - 1)
        do n=m, this%ip%nrot_t+1
          I1 = 2 * (n - 1)
          call wig3j(I1, 2, I2, 0, 0, 0, wg1)
          call wig3j(I1, 4, I2, 0, 0, 0, wg2)
          V(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                        * sqrt(dble((2*I1+1) * (2*I2+1)))
        end do
      end do

!   The second mode of phonon excitation
      w1 = this%ip%betaL_t2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w1 * (this%ip%Rm / r) ** this%ip%lambda_t2 / r
      else
        x = w1 * (r / this%ip%Rm) ** this%ip%lambda_t2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch
        do n=this%ip%nrot_t+2, this%ip%Nch
          if (m == 1 .and. n == this%ip%nrot_t+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%nrot_t+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%nrot_t))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do
    end if

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine  Vn_vib(this, ir, r, Vpot)
    use mkl95_lapack
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, i, n2, idum
    real(8), intent(in) :: r
    real(8), intent(out) :: Vpot(:,:)
!   real(8), intent(in) :: beta(this%ip%Nch_t_noncoll)
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), allocatable, save :: O_sub(:,:), Oa_sub(:)
    real(8), dimension(this%ip%Nch) :: Vn
    real(8), parameter :: PI = 3.141592653589793d0
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: c

    if (.not. allocated(O)) then

      c = this%ip%betaL * this%ip%Rm * SQRT4PI_INV 

      if (this%ip%coup_p == 0) then
        allocate(O(this%ip%Nch,this%ip%Nch))
        allocate(Oa(this%ip%Nch))
        do m=1, this%ip%Nphonon_p+1
          do n=m, this%ip%Nphonon_p+1
            if (n == m + 1) then
              O(n,m) = c * sqrt(dble(m))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

!  2nd mode in projectile
        c = this%ip%betaL_p2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch_p
          do n=max(m,this%ip%Nphonon_p+2), this%ip%Nch_p
            if (m == 1 .and. n == this%ip%Nphonon_p+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nphonon_p+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nphonon_p))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

        call syev(O, Oa, 'V','L')
      else if (this%ip%coup_t == 0) then
        allocate(O(this%ip%Nch_t_coll,this%ip%Nch_t_coll))
        allocate(Oa(this%ip%Nch_t_coll))
        do m=1, this%ip%Nphonon_t+1
          do n=m, this%ip%Nphonon_t+1
            if (n == m + 1) then
              O(n,m) = c * sqrt(dble(m))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

! 2nd mode in target
        c = this%ip%betaL_t2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch_t_coll
          do n=this%ip%Nphonon_t+2, this%ip%Nch_t_coll
            if (m == 1 .and. n == this%ip%Nphonon_t+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nphonon_t+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nphonon_t))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

        call syev(O, Oa, 'V','L')

      end if

    end if
    
    if (this%ip%coup_p == 0) then
      n2 = this%ip%Nch
      Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))
    else if (this%ip%coup_t == 0) then
      n2 = this%ip%Nch_t_coll
      Vn(1:n2) = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa(1:n2)) / this%ip%a))
      call Vn_noncoll(this, ir, r, Vpot)
    end if

    Vpot(1:n2,1:n2) = 0.0d0
    do m=1, n2
      do n=m, n2
        do i=1, n2
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vc_vib(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m
    real(8), intent(in) :: r
    real(8), intent(out) :: V(:,:)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, w, ZZe2
    real(8) :: wa(this%ip%Nch_t_noncoll)

    ZZe2 = this%ip%Zp * this%ip%Zt * e2
    w = this%ip%betaLC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda+1) * ZZe2
    if (r > this%ip%Rm) then
      x = w * (this%ip%Rm / r) ** this%ip%lambda / r
    else
      x = w * (r / this%ip%Rm) ** this%ip%lambda / this%ip%Rm
    end if

    if (this%ip%coup_p == 0) then
      do m=1, this%ip%Nphonon_p+1
        do n=m, this%ip%Nphonon_p+1
          if (n == m + 1) then
            V(n,m) = x * sqrt(dble(m))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

! 2nd mode in projectile
      w = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w * (this%ip%Rm / r) ** this%ip%lambda_p2 / r
      else
        x = w * (r / this%ip%Rm) ** this%ip%lambda_p2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch_p
        do n=max(m,this%ip%Nphonon_p+2), this%ip%Nch_p
          if (m == 1 .and. n == this%ip%Nphonon_p+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%Nphonon_p+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%Nphonon_p))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

    else if (this%ip%coup_t == 0) then
      do m=1, this%ip%Nphonon_t+1
        do n=m, this%ip%Nphonon_t+1
          if (n == m + 1) then
            V(n,m) = x * sqrt(dble(m))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

! 2nd mode in target
      w = this%ip%betaL_t2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w * (this%ip%Rm / r) ** this%ip%lambda_t2 / r
      else
        x = w * (r / this%ip%Rm) ** this%ip%lambda_t2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch_t_coll
        do n=this%ip%Nphonon_t+2, this%ip%Nch_t_coll
          if (m == 1 .and. n == this%ip%Nphonon_t+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%Nphonon_t+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%Nphonon_t))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

! non-collective
!     wa = beta * SQRT4PI_INV * 3.0d0/dble(2*lambda + 1) * ZZe2
!     if (r > this%ip%Rt) then
!       V(this%ip%Nch_t_coll+1:,1) = wa * (this%ip%Rt / r) ** lambda / r
!     else
!       V(this%ip%Nch_t_coll+1:,1) = wa * (r / this%ip%Rt) ** lambda / this%ip%Rt
!     end if

!     do n=2, this%ip%Nch_t_coll
!       V(this%ip%Nch_t_coll+1:,n) = 0.0d0
!     end do

      do n=1, this%ip%Nch_t
        V(max(n,this%ip%Nch_t_coll+1):,n) = 0.0d0
      end do

    end if

    return
  end subroutine
!----------------------------------------------------------------------
  function eps(this, n) result(e)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: n
    integer :: I
    real(8) :: e

    if (this%ip%coup_p == 0) then
      if (n < this%ip%Nphonon_p+2) then
        e = dble(n-1) * this%ip%omega
      else
        e = dble(n-this%ip%Nphonon_p-1) * this%ip%omega_p2
      end if
    else if (this%ip%coup_t == 0) then
      if (n < this%ip%Nphonon_t+2) then
        e = dble(n-1) * this%ip%omega
      else
        e = dble(n-this%ip%Nphonon_t-1) * this%ip%omega_t2
      end if
    else if (this%ip%coup_p == 1) then
      if (n < this%ip%nrot_p+2) then 
        I = 2 * (n - 1)
        e = dble(I * (I + 1)) / 6.0d0 * this%ip%E2
      else
        e = dble(n - this%ip%nrot_p - 2) * this%ip%omega_p2
      end if
    else if (this%ip%coup_t == 1) then
      if (n < this%ip%nrot_t+2) then 
        I = 2 * (n - 1)
        e = dble(I * (I + 1)) / 6.0d0 * this%ip%E2
      else
        e = dble(n - this%ip%nrot_t - 2) * this%ip%omega_t2
      end if
    end if


    return
  end function
!----------------------------------------------------------------------!
  subroutine Vcoup_pro_tar(this, r, ir)
    use potentials, only : Vn
    implicit none
    type(coup_mat), intent(inout) :: this
    integer, intent(in) :: ir
    integer :: n, m, nst, nlen
    real(8), intent(in) :: r
    real(8), allocatable, dimension(:,:) :: Vn_cp, Vc_cp
    real(8) :: V

    allocate(Vn_cp(this%ip%Nch,this%ip%Nch), Vc_cp(this%ip%Nch,this%ip%Nch))
    call Vn_pro_tar(this, ir, r, Vn_cp)
    call Vc_pro_tar(this, r, Vc_cp)
    V = Vn(this%ip,r)
    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      this%Vcp_linear(nst,ir) = (Vn_cp(n,n) - V)+ Vc_cp(n,n)
      this%Vcp_linear(nst+1:nst-1+nlen,ir) = Vn_cp(n+1:this%ip%Nch,n) + Vc_cp(n+1:this%ip%Nch,n)
      nst = nst + nlen
    end do

    deallocate(Vn_cp, Vc_cp)

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vn_pro_tar(this, ir, r, Vpot)
    use mkl95_lapack
    use global_constant, only : PI
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: n, m, i, I1, I2, n1st, idum, st, ed
    integer, intent(in) :: ir
    real(8), intent(in) :: r
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), intent(out) :: Vpot(this%ip%Nch,this%ip%Nch)
    real(8), dimension(this%ip%Nch_t, this%ip%Nch_t) :: Q
    real(8), dimension(this%ip%Nch) :: Vn
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = sqrt(9.0d0/(4.0d0*PI))
    real(8) :: ct, wg1, wg2, w1, w2

    if (r > this%ip%rcut) then
      Vpot = 0.0d0
      return
    end if

    if (.not. allocated(O)) then
      allocate(O(this%ip%Nch,this%ip%Nch))
      allocate(Oa(this%ip%Nch))

      do n=1, this%ip%Nch_t
        Q(n:this%ip%Nch_t,n) = 0.0d0
      end do

! 1st mode
      select case(this%ip%coup_t)
        case(0)
          ct = this%ip%betaL_t * this%ip%Rt * SQRT4PI_INV
          do m=1, this%ip%Nphonon_t+1
            do n=m, this%ip%Nphonon_t+1
              if (n == m + 1) then
                Q(n,m) = ct * sqrt(dble(m))
              else
                Q(n,m) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_t
        case(1)
          w1 = CONST1 * this%ip%beta2t * this%ip%Rt
          w2 = CONST2 * this%ip%beta4t * this%ip%Rt
          do m=1, this%ip%Nrot_t+1
            I2 = 2 * (m - 1)
            do n=m, this%ip%Nrot_t+1
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              Q(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                          * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
            end do
          end do
          n1st = this%ip%Nrot_t
      end select

! 2nd mode
      ct = this%ip%betaL_t2 * this%ip%Rt * SQRT4PI_INV 
      do m=1, this%ip%Nch_t_coll
        do n=max(m,n1st+2), this%ip%Nch_t_coll
          if (m == 1 .and. n == n1st+2) then
            Q(n,m) = ct
          else if (n == m + 1 .and. m > n1st+1) then
            Q(n,m) = ct * sqrt(dble(m-n1st))
          else
            Q(n,m) = 0.0d0
          end if
        end do
      end do

! projectile excitations
      call sub_Vn(Q, O)

      call syev(O, Oa, 'V','L')
    end if

    Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))

    Vpot = 0.0d0
    do m=1, this%ip%Nch
      do n=m, this%ip%Nch
        do i=1, this%ip%Nch
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do

    call Vn_noncoll_mutual(this, ir, r, Vpot)

    contains
!----------------------------------------------------------------------
    subroutine sub_Vn(Q, O)
      implicit none
      integer :: i, I1, I2, n, m, nt, mt, n1st
      real(8), intent(in) :: Q(:,:)
      real(8), intent(out) :: O(:,:)
      real(8) :: cp, d, w1, w2, wg1, wg2, Nch_t2

      select case(this%ip%coup_p)
! 1st mode in projectile
        case(0)
          cp = this%ip%betaL_p * this%ip%Rp * SQRT4PI_INV
          do m=1, this%ip%Nphonon_p+1
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nphonon_p+1
              nt = (n - 1) * this%ip%Nch_t
              if (n == m + 1) then
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
                do i=1, this%ip%Nch_t
!               do i=1, this%ip%Nch_t_coll
                  O(nt+i,mt+i) = cp * sqrt(dble(m))
                end do
              else if (n == m) then
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Q
              else
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_p
        case(1)
          w1 = CONST1 * this%ip%beta2p * this%ip%Rp
          w2 = CONST2 * this%ip%beta4p * this%ip%Rp
          do m=1, this%ip%Nrot_p+1
            I2 = 2 * (m - 1)
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nrot_p+1
              nt = (n - 1) * this%ip%Nch_t
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              if (n == m) then
                O(nt+1:nt+this%ip%Nch_t, mt+1:mt+this%ip%Nch_t) = Q
              else
                O(nt+1:nt+this%ip%Nch_t, mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
              d = (w1 * wg1 * wg1 + w2 * wg2 * wg2)           &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
              do i=1, this%ip%Nch_t
!             do i=1, this%ip%Nch_t_coll
                O(nt+i,mt+i) = O(nt+i,mt+i) + d
              end do
            end do
          end do
          n1st = this%ip%Nrot_p
      end select

! 2nd mode in projectile
      if (this%ip%Nphonon_p2 > 0) then
        cp = this%ip%betaL_p2 * this%ip%Rp * SQRT4PI_INV
        do m=1, this%ip%Nch_p
          mt = (m - 1) * this%ip%Nch_t
          do n=max(m,n1st+2), this%ip%Nch_p
            nt = (n - 1) * this%ip%Nch_t
            if (m == 1 .and. n == n1st+2) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
  !           do i=1, this%ip%Nch_t_coll
                O(nt+i,mt+i) = cp
              end do
            else if (n == m + 1 .and. m > n1st+1) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
  !           do i=1, this%ip%Nch_t_coll
                O(nt+i,mt+i) = cp * sqrt(dble(m-n1st))
              end do
            else if (n == m) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Q
            else
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
            end if
          end do
        end do
      end if

      return
    end subroutine
  end subroutine
!----------------------------------------------------------------------!
  subroutine Vc_pro_tar(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I1, I2, n1st
    real(8), intent(in) :: r
    real(8), intent(out) :: V(:,:)
    real(8), dimension(this%ip%Nch_t,this%ip%Nch_t) :: Vt
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, x2, w, wg1, wg2, w1, w2, ZZe2
!   real(8) :: wa(this%ip%Nch_t_noncoll)

    ZZe2 = this%ip%Zt * this%ip%Zp * e2

! non-collective
!   wa = 3.0d0 / dble(2 * lambda + 1) * beta * SQRT4PI_INV * ZZe2
!   if (r > this%ip%Rt) then
!     Vt(this%ip%Nch_t_coll+1:,1) = wa * (this%ip%Rt / r) ** lambda / r
!   else
!     Vt(this%ip%Nch_t_coll+1:,1) = wa * (r / this%ip%Rt) ** lambda / this%ip%Rt
!   end if

! clear non-coll part
    do n=1, this%ip%Nch_t
      Vt(max(this%ip%Nch_t_coll+1,n):this%ip%Nch_t,n) = 0.0d0
    end do

! 1st-mode
    select case (this%ip%coup_t)
      case(0)
        w = this%ip%betaL_tC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t+1) * ZZe2
        if (r > this%ip%Rt) then
          x = w * (this%ip%Rt / r) ** this%ip%lambda_t / r
        else
          x = w * (r / this%ip%Rt) ** this%ip%lambda_t / this%ip%Rt
        end if

        do m=1, this%ip%Nphonon_t+1
          do n=m, this%ip%Nphonon_t+1
            if (n == m + 1) then
              Vt(n,m) = x * sqrt(dble(m))
            else
              Vt(n,m) = 0.0d0
            end if
          end do
        end do
        n1st = this%ip%Nphonon_t
      case(1)
        if (r > this%ip%Rt) then
          x = this%ip%Rt / r
          x2 = x * x / r
        else
          x = r / this%ip%Rt
          x2 = x * x / this%ip%Rt
        end if
        w1 = 0.6d0 * ZZe2 * CONST1* this%ip%beta2t * (1.0d0 + CONST2*this%ip%beta2t)
        w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4t + CONST3*this%ip%beta2t*this%ip%beta2t)
        do m=1, this%ip%Nrot_t+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%Nrot_t+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            Vt(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                           * sqrt(dble((2*I1+1) * (2*I2+1)))
          end do
        end do
        n1st = this%ip%Nrot_t
    end select

! 2nd-mode in target
    w = this%ip%betaL_t2C * SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
    if (r > this%ip%Rt) then
      x = w * (this%ip%Rt / r) ** this%ip%lambda_t2 / r
    else
      x = w * (r / this%ip%Rt) ** this%ip%lambda_t2 / this%ip%Rt
    end if

    do m=1, this%ip%Nch_t_coll
      do n=max(m,n1st+2), this%ip%Nch_t_coll
        if (m == 1 .and. n == n1st+2) then
          Vt(n,m) = x
        else if (n == m + 1 .and. m > n1st+1) then
          Vt(n,m) = x * sqrt(dble(m-n1st))
        else
          Vt(n,m) = 0.0d0
        end if
      end do
    end do

    call sub_Vc(r, Vt, V)

    contains
!----------------------------------------------------------------------!
    subroutine sub_Vc(r, Vt, V)
      use global_constant, only : e2
      implicit none
      integer :: i, I1, I2, m, n, nt, mt, n1st
      real(8), intent(in) :: r, Vt(:,:)
      real(8), intent(out) :: V(:,:)
      real(8) :: x, x2, w, w1, w2, wg1, wg2, d, ZZe2

      ZZe2 = this%ip%Zt * this%ip%Zp * e2
      select case(this%ip%coup_p)
        case(0)
          w = this%ip%betaL_pC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p+1) * ZZe2
          if (r > this%ip%Rp) then
            x = w * (this%ip%Rp / r) ** this%ip%lambda_p / r
          else
            x = w * (r / this%ip%Rp) ** this%ip%lambda_p / this%ip%Rp
          end if
          do m=1, this%ip%Nphonon_p+1
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nphonon_p+1
              nt = (n - 1) * this%ip%Nch_t
              if (n == m + 1) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
                do i=1, this%ip%Nch_t
                  V(nt+i, mt+i) = x * sqrt(dble(m))
                end do
              else if (n == m) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
              else
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_p
        case(1)
          if (r > this%ip%Rp) then
            x = this%ip%Rp / r
            x2 = x * x / r
          else
            x = r / this%ip%Rp
            x2 = x * x / this%ip%Rp
          end if
          w1 = 0.6d0 * ZZe2 * CONST1 * this%ip%beta2p *(1.0d0 + CONST2*this%ip%beta2p)
          w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4p + CONST3*this%ip%beta2p*this%ip%beta2p)
          do m=1, this%ip%Nrot_p+1
            mt = (m - 1) * this%ip%Nch_t
            I2 = 2 * (m - 1)
            do n=m, this%ip%Nrot_p+1
              nt = (n - 1) * this%ip%Nch_t
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              if (n == m) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
              else
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
              d = (w1 * wg1 * wg1 + w2 * (wg2*x) ** 2) * x2       &
                     * sqrt(dble((2*I1+1) * (2*I2+1)))
              do i=1, this%ip%Nch_t
                V(nt+i,mt+i) = V(nt+i,mt+i) + d
              end do
            end do
          end do
          n1st = this%ip%Nrot_p
      end select

! 2nd mode in projectile
      if (this%ip%Nphonon_p2 > 0) then
        w = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
        if (r > this%ip%Rp) then
          x = w * (this%ip%Rp / r) ** this%ip%lambda_p2 / r
        else
          x = w * (r / this%ip%Rp) ** this%ip%lambda_p2 / this%ip%Rp
        end if
        do m=1, this%ip%Nch_p
          mt = (m - 1) * this%ip%Nch_t
          do n=max(m,n1st+2), this%ip%Nch_p
            nt = (n - 1) * this%ip%Nch_t
            if (m == 1 .and. n == n1st+2) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                V(nt+i, mt+i) = x
              end do
            else if (n == m + 1 .and. m > n1st+1) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                V(nt+i, mt+i) = x * sqrt(dble(m-n1st))
              end do
            else if (n == m) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
            else
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
            end if
          end do
        end do
      end if

    end subroutine
  end subroutine
!----------------------------------------------------------------------!
  subroutine eps_pro_tar(this, eps_t, eps)
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I, n1st
    real(8), intent(out) :: eps(this%ip%Nch)
    real(8), intent(in) :: eps_t(this%ip%Nch_t_noncoll)
    real(8), dimension(this%ip%Nch_t) :: eps_t_all
    real(8) :: e(4)

    select case(this%ip%coup_t)
      case(0)
        forall (n=1:this%ip%Nphonon_t+1) eps_t_all(n) = dble(n-1)*this%ip%omega_t
        n1st = this%ip%Nphonon_t
      case(1)
        do n=1, this%ip%Nrot_t + 1
          I = 2 * (n - 1)
          eps_t_all(n) = dble(I*(I+1)) / 6.0d0 * this%ip%E2t
        end do
        n1st = this%ip%Nrot_t
    end select

!$  2nd phonon
    forall (n=n1st+2:this%ip%Nch_t) 
      eps_t_all(n) = dble(n-n1st-1) * this%ip%omega_t2
    end forall

!$  noncollective
    eps_t_all(this%ip%Nch_t_coll+1:) = eps_t
!       e(1) = 0.0d0
!       e(2) = 1.6337d0
!       e(3) = 4.2477d0
!       e(4) = 8.7776d0
      
    select case(this%ip%coup_p)
      case(0)
        do n=1, this%ip%Nphonon_p+1
          m = (n - 1) * this%ip%Nch_t
          eps(m+1:m+this%ip%Nch_t) = dble(n - 1) * this%ip%omega_p + eps_t_all
        end do
        n1st = this%ip%Nphonon_p
      case(1)
        do n=1, this%ip%Nrot_p+1
          m = (n - 1) * this%ip%Nch_t
          I = 2 * (n - 1)
          eps(m+1:m+this%ip%Nch_t) = dble(I*(I+1))/6.0d0*this%ip%E2p + eps_t_all
!         eps(m+1:m+this%ip%Nch_t) = e(n) + eps_t_all
        end do
        n1st = this%ip%Nrot_p
    end select
! 2nd mode in projectile
    do n=n1st+2, this%ip%Nch_p
      m = (n - 1) * this%ip%Nch_t
      eps(m+1:m+this%ip%Nch_t) = dble(n-n1st-1)*this%ip%omega_p2 + eps_t_all
    end do

    return
  end subroutine
!----------------------------------------------------------------------!
  function rho1(this, e) result(r)
    implicit none
    class(coup_mat), intent(in) :: this
    real(8), intent(in) :: e
    real(8) :: r

    r = this%rho0 * exp(2.0d0 * sqrt(0.125d0 * this%alevel * e))

    r = 100.0d0
    if (e > 6.0d0) then
      r = 200.0d0
    else if (e < 4.0d0) then
      r = 100.0d0
    end if


  end function
!----------------------------------------------------------------------!
! function rho2(this, n) result(r)
!   implicit none
!   class(coup_mat), intent(in) :: this
!   integer, intent(in) :: n
!   real(8) :: r

!   r = this%rhoa(n)

! end function
!----------------------------------------------------------------------!
  function g(this, m, n, r1, r2) result(f)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: m, n
    integer :: L, I1, I2
    real(8), intent(in) :: r1, r2
    real(8) :: f, s, wg, e1, e2

    e1 = this%e_n(m)
    e2 = this%e_n(n)
    I1 = this%lambda(m)
    I2 = this%lambda(n)

    s = 0.0d0
    do L=2, this%lambda_max
      call wig3j(I1,L,I2,0,0,0,wg)
      s = s + wg*wg*this%w_l(L)
    end do
    f = sqrt(s * sqrt(dble((2*I1+1)*(2*I2+1)) / (this%ld%rho_J(e1,I1)*this%ld%rho_J(e2,I2))))  &
      * exp(- 0.25d0*((e1-e2)/this%Delta)**2)         &
      * (2.0d0/(PI*this%sigma*this%sigma)) ** 0.25d0  &
      * exp(- ((r1-r2)/this%sigma)**2)                &
      * exp(- 0.5d0*(r1/this%alpha)**2)

  end function   
!----------------------------------------------------------------------!
  function h(this, r) result(f)
    implicit none
    class(coup_mat), intent(in) :: this
    real(8), intent(in) :: r
    real(8) :: f, g

    g = exp((r - this%ip%Rn)/this%ip%a)
!   g = exp((r - this%alpha)/this%ip%a)
    f = g / (1.0d0 + g)**2

  end function
!----------------------------------------------------------------------!
  function g2(this, m, r1, r2) result(f)
!----------------------------------------------------------------------!
!   g2 is nearly zero for large r1.                                    !
!----------------------------------------------------------------------!
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: m
    integer :: L, I1, I2
    real(8), intent(in) :: r1, r2
    real(8) :: f, s, wg, e1

    e1 = this%e_n(m)
    I1 = this%lambda(m)

    s = 0.0d0
    do L=0, this%lambda_max
      call wig3j(I1,L,0,0,0,0,wg)
      s = s + wg*wg*this%w_l(L)
    end do
    f = sqrt(s * sqrt(dble(2*I1+1) / this%ld%rho2(e1)))    &
      * exp(- 0.25d0*(e1/this%Delta)**2)                 &
      * (2.0d0/(PI*this%sigma*this%sigma)) ** 0.25d0          &
      * exp(- ((r1-r2)/this%sigma)**2)                        &
!     * exp(- 0.5d0*(r1/this%alpha)**2)
      * h(this, r1)

  end function   
!----------------------------------------------------------------------!
  subroutine white_noise(this)
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: ir, n, m, n_max
!   integer :: idum
    real(8) :: hoge

    hoge = 0.1d0
    n_max = this%ip%Nch_p * (this%ip%Nch_p+1) / 2
    do ir=1, this%ip%rgrid+2
!   do ir=1, this%ip%rgrid_cut+2
      do n=1, n_max
!     do n=1, this%ip%Nch_t
        do m=this%ip%Nch_t_coll+1, this%ip%Nch_t
!       do m=n, this%ip%Nch_t
!           hoge = gaus_ran_num(idum)
!           hoge = gaus_ran_num(idum)
            this%wn(m,n,ir) = gaus_ran_num(this%idum)
        end do
      end do
    end do
  end subroutine
!----------------------------------------------------------------------!
  function get_wn(this, m, n, ir) result(w)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: m,n,ir
    real(8) :: w

    if (ir > 2*this%ip%rgrid_cut) then
      w = 0.0d0
    else
      w = this%wn(m,n,ir)
    end if

    return
  end function
!----------------------------------------------------------------------!
  subroutine Vn_noncoll(this, ir, r1, Vn)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, i, j, lambda
    integer :: L, I1, I2
    real(8), intent(in) :: r1
    real(8), intent(inout) :: Vn(:,:)
    real(8) :: r2, wg, cof, b(this%ip%Nch_t), str(this%ip%Nch_t)
    real(8) :: e1, e2, s, f
    character(len=40), parameter :: FM='(1x,10f8.3)'


    do n=1, this%ip%Nch_t
      do m=max(this%ip%Nch_t_coll+1,n), this%ip%Nch_t
        Vn(m,n) = 0.0d0
      end do
    end do

!   do n=1, this%ip%Nch_t
!     do m=max(this%ip%Nch_t_coll+1,n), this%ip%Nch_t
    do n=1, 1
      do m=max(this%ip%Nch_t_coll+1,n), this%ip%Nch_t
        do j=1, this%ip%rgrid+2
          r2 = this%ip%rmin + dble(j-1) * this%ip%dr
!         Vn(m,n) = Vn(m,n) + g(this,m,n,r1,r2)*this%wn(m,n,j)*this%sqrt_dr
          Vn(m,n) = Vn(m,n) + g2(this,m,r1,r2)*get_wn(this,m,n,j)*this%sqrt_dr
        end do
      end do
    end do

  end subroutine
!----------------------------------------------------------------------!
  subroutine Vn_noncoll_mutual(this, ir, r1, Vn)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, i, j, lambda, nb, mb, m2, k, p
    integer :: L, I1, I2
    real(8), intent(in) :: r1
    real(8), intent(inout) :: Vn(:,:)
    real(8) :: r2, wg, cof, b(this%ip%Nch_t), str(this%ip%Nch_t)
    real(8) :: e1, e2, s, f
    character(len=40), parameter :: FM='(1x,10f8.3)'

    p = 0
    do n=1, this%ip%Nch_p
      nb = (n - 1) * this%ip%Nch_t + 1
      do m=n, this%ip%Nch_p
        mb = (m - 1)*this%ip%Nch_t
        p = p + 1
        do k=this%ip%Nch_t_coll+1, this%ip%Nch_t
          do j=1, this%ip%rgrid+2
!         do j=1, this%ip%rgrid_cut+2
            r2 = this%ip%rmin + dble(j-1) * this%ip%dr
            Vn(mb+k,nb) = Vn(mb+k,nb) + g2(this,k,r1,r2)*get_wn(this,k,p,j)*this%sqrt_dr
          end do
        end do
      end do
    end do
!   stop

  end subroutine
!----------------------------------------------------------------------!
  subroutine record_strength_dist(this, Fname)
    implicit none
    class(coup_mat), intent(in) :: this
    character(len=*), intent(in) :: Fname
    integer :: n, L, I1
    real(8) :: e1, s, f, wg, g, eff_str, w
    character(len=50), parameter :: FM="(1x,3f8.3)"

    open(7,file=Fname)
    w = 200.0d0
    do n=this%ip%Nch_t_coll+1, this%ip%Nch_t
      e1 = this%e_n(n)
      I1 = this%lambda(n)
      s = 0.0d0
      do L=2, this%lambda_max
        call wig3j(0,L,I1,0,0,0,wg)
!       s = s + wg*wg*this%w_l(L)
        s = s + wg*wg
      end do
      g = w * s * sqrt(dble(2*I1+1) / (this%ld%rho2(e1))) &
            * exp(-0.5d0*(e1/this%Delta)**2) 
      write(7,*) e1, sqrt(g)
    end do
    close(7)
    write(6,*) sqrt(g)

  end subroutine
end module

