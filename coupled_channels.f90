module coupled_channels
  use input_data, only : inp
  use relative_potential
  use coupling_matrix
  private

  type, public :: cc_scat
    type(inp), pointer, private :: ip
    type(rel_pot), pointer :: vpot_rel
    type(coup_mat), pointer :: vcoup_mat
    integer :: num_stab_pt
    integer :: ibar
    integer, allocatable, dimension(:) :: stab_pt
    real(8), private, allocatable, dimension(:,:) :: sigC
    real(8), private, allocatable, dimension(:,:) :: IE
    real(8), private :: w
    complex(8), private, allocatable, dimension(:,:) :: Cnm, Dnm
    complex(8), private, allocatable, dimension(:,:) :: kai0, kai2
    complex(8), private, allocatable, dimension(:,:) :: S
    complex(8), private, allocatable, dimension(:,:) :: Sl
    contains
    procedure :: cc_scat_
    procedure :: cc_scattering
    procedure, private:: Numerov
    procedure, private:: match
    procedure, private:: Penet
    procedure :: scat_amp_Coul
    procedure :: scat_amp_nucl
    procedure :: Ruthefrd_sig
    procedure :: destruct_cc_scat
!   final :: destruct
  end type

  contains
!----------------------------------------------------------------------
  subroutine cc_scat_(this, ip, vpot_rel, vcoup_mat)
    use global_constant
    implicit none
    class(cc_scat), intent(out) :: this
    type(inp), intent(in), target :: ip
    type(rel_pot), target, intent(in) :: vpot_rel
    type(coup_mat), target, intent(in) :: vcoup_mat
    integer :: i

    this%ip => ip
    this%vpot_rel => vpot_rel
    this%vcoup_mat => vcoup_mat

    allocate(this%Cnm(ip%Nch,ip%Nch),this%kai0(ip%Nch,ip%Nch))
    allocate(this%Dnm(ip%Nch,ip%Nch),this%kai2(ip%Nch,ip%Nch))
    allocate(this%S(ip%Nch,ip%Nch))
    allocate(this%Sl(ip%Nch,0:ip%Jmax))
    allocate(this%sigC(0:ip%Jmax,ip%Nch))
    allocate(this%IE(ip%Nch,ip%Nch))

    this%IE = 0.0d0
    forall (i=1:ip%Nch) this%IE(i,i) = 1.0d0

    this%w = 0.5d0 * hbar * hbar * PI / this%ip%rmass

    this%ibar = nint((this%vcoup_mat%rb - this%ip%rmin) / this%ip%dr)
    this%num_stab_pt = ip%num_stab_pt
    if (this%num_stab_pt <= 1) then
      this%num_stab_pt = 1
      allocate(this%stab_pt(1))
      this%stab_pt(1) = this%ibar
    else if (this%num_stab_pt > ip%rgrid_cut) then
      this%num_stab_pt = ip%rgrid_cut
      allocate(this%stab_pt(this%num_stab_pt))
      do i=1, this%num_stab_pt
        this%stab_pt(i) = (i*this%ip%rgrid_cut) / this%num_stab_pt
      end do
    else
      allocate(this%stab_pt(this%num_stab_pt))
      do i=1, this%num_stab_pt
        this%stab_pt(i) = (i*this%ip%rgrid_cut) / this%num_stab_pt
      end do
    end if


  end subroutine
!----------------------------------------------------------------------
  elemental subroutine destruct_cc_scat(this)
    implicit none
    class(cc_scat), intent(inout) :: this

    if (allocated(this%Cnm)) deallocate(this%Cnm)
    if (allocated(this%Dnm)) deallocate(this%Dnm)
    if (allocated(this%kai0)) deallocate(this%kai0)
    if (allocated(this%kai2)) deallocate(this%kai2)
    if (allocated(this%S)) deallocate(this%S)
    if (allocated(this%Sl)) deallocate(this%Sl)
    if (allocated(this%sigC)) deallocate(this%sigC)

  end subroutine
!----------------------------------------------------------------------
  subroutine cc_scattering(this, E, spin, sig_fus, sig_iel_n, sig_iel)
    implicit none
    class(cc_scat), intent(inout) :: this
    integer :: J
    real(8), intent(in) :: E
    real(8), intent(out) :: sig_fus, spin, sig_iel
    real(8), dimension(:), intent(out) :: sig_iel_n
    real(8) :: Rn(this%ip%Nch), P

    sig_fus = 0.0d0
    spin = 0.0d0
    sig_iel_n = 0.0d0
    do J=0, this%ip%Jmax
      call Numerov(this, J, E)
      call match(this, J, E)
      call Penet(this, J, E, Rn, P)
      sig_fus = sig_fus + dble(2 * J + 1) * P
      spin = spin + dble(J) * dble(2 * J + 1) * P
      sig_iel_n = sig_iel_n + dble(2 * J + 1) * Rn
    end do
    if (sig_fus /= 0.0d0) then
      spin = spin / sig_fus
    else
      spin = 0.0d0
    end if
    sig_fus = sig_fus * this%w / E * 10.0d0
    sig_iel_n = sig_iel_n * this%w / E * 10.0d0
    sig_iel = sum(sig_iel_n(2:this%ip%Nch))

  end subroutine
!----------------------------------------------------------------------
  pure subroutine Numerov(this, J, E)
!   use mkl95_lapack
    use lin_coup
    implicit none
    class(cc_scat), intent(inout) :: this
    integer, intent(in) :: J
    integer :: i, k, L, istab, ibar
    real(8), intent(in) :: E!, rb
    real(8), parameter :: w3 = 1.732050807568877d0  ! sqrt(3.0d0)
    real(8) :: w1, w2, norm(this%ip%Nch)
    complex(8), allocatable, dimension(:,:) :: psi1, A0, A1, T, U
    
    w1 = this%ip%dr * this%ip%dr / 12.0d0
    w2 = this%ip%dr * this%ip%dr / sqrt(12.0d0)

    allocate(A0(this%ip%Nch,this%ip%Nch))
    allocate(A1(this%ip%Nch,this%ip%Nch))
    allocate(T(this%ip%Nch,this%ip%Nch))
    allocate(psi1(this%ip%Nch,this%ip%Nch))
    allocate(U(this%ip%Nch,this%ip%Nch))

    call Anm(2, J, E, A1)
    this%kai0 = 0.0d0
    psi1 = 0.01d0 * this%IE 
    psi1 = matmul(this%IE-w1*A1, psi1)
    istab = 1
    do i=3, this%ip%rgrid+1
      T = w2 * A1 + w3 * this%IE
      do k=1, this%ip%Nch
        U(k,k) = sum(T(:,k) * T(:,k)) - 1.0d0
        do L=k+1, this%ip%Nch
          U(L,k) = sum(T(:,L) * T(:,k))
          U(k,L) = U(L,k)
        end do
      end do
      this%kai2 = matmul(U, psi1) - this%kai0
      if (i == this%stab_pt(istab)) then
        call stabilize(A1, psi1, this%kai2, this%S)
        istab = min(istab + 1,10)
      end if
      norm = maxval(abs(this%kai2), dim=1)
      do k=1, this%ip%Nch
        this%kai0(:,k) = psi1(:,k) / norm(k)
        psi1(:,k) = this%kai2(:,k) / norm(k)
      end do
      if (i == this%ip%rgrid+1) A0 = A1
      call Anm(i, J, E, A1)
    end do

    T = w2 * A1 + w3 * this%IE
    do k=1, this%ip%Nch
      U(k,k) = sum(T(:,k) * T(:,k)) - 1.0d0
      do L=k+1, this%ip%Nch
        U(L,k) = sum(T(:,L) * T(:,k))
        U(k,L) = U(L,k)
      end do
    end do
    this%kai2 = matmul(U, psi1) - this%kai0

    call Anm(this%ip%rgrid+2, J, E, U)
    A1 = this%IE - w1 * U
    T = this%IE
!   call sysv(A1, T)  ! LAPACK
    call lin_sy(A1, T)
    this%kai2 = matmul(T, this%kai2)
    A1 = this%IE - w1 * A0
    T = this%IE
!   call sysv(A1, T)  ! LAPACK
    call lin_sy(A1, T)
    this%kai0 = matmul(T, this%kai0)

    deallocate(A0, A1, T, U, psi1)

    contains
!----------------------------------------------------------------------!
    pure subroutine Anm(ir, J, E, A)
      use global_constant, only : hbar
      implicit none
      integer, intent(in) :: J, ir
      integer :: n
      real(8), intent(in) :: E
      complex(8), intent(out) :: A(:,:)
      
      A = this%vcoup_mat%getVcp(ir)
      forall (n=1:this%ip%Nch)
        A(n,n) = A(n,n) + this%vpot_rel%get_Vrel(J,ir) + this%vcoup_mat%e_n(n) - E
      end forall
      A = 2.0d0 * this%ip%rmass / (hbar * hbar) * A

      return
    end subroutine
!----------------------------------------------------------------------!
    pure subroutine stabilize(A1, u1, u2, S)
!     use mkl95_lapack
      use lin_coup
      implicit none
      real(8) :: w
      complex(8), intent(in) :: A1(:,:)
      complex(8), intent(inout), dimension(:,:) :: u1, u2
      complex(8), intent(out) :: S(:,:)
      complex(8), allocatable, dimension(:,:) :: Sinv, T

      w = this%ip%dr * this%ip%dr / 12.0d0
      allocate(T(this%ip%Nch,this%ip%Nch))
      allocate(Sinv(this%ip%Nch,this%ip%Nch))

      S = this%IE - w1 * A1
      T = this%IE
!     call gesv(S, T)  ! LAPACK
      call lin_ge(S, T)
      S = matmul(T, u1)
      Sinv = this%IE
!     call gesv(S, Sinv)  ! LAPACK
      call lin_ge(S, Sinv)
      u1 = matmul(u1, Sinv)
      u2 = matmul(u2, Sinv)

      deallocate(T, Sinv)

      return
    end subroutine
  end subroutine
!----------------------------------------------------------------------!
  subroutine match(this, J, E)
    use global_constant, only : hbar, e2
    use Coulomb, only : dfcoul
    implicit none
    class(cc_scat), intent(inout) :: this
    integer, intent(in) :: J
    integer :: n, iexp(0:this%ip%Jmax)
    real(8), intent(in) :: E
    real(8), dimension(0:this%ip%Jmax) :: fc, fpc, gc, gpc, sig0
    real(8) :: eta, rho, E_C, w1, w2
    complex(8), dimension(this%ip%Nch) :: Hp0, Hp2, Hn0, Hn2
    complex(8), parameter :: ai = (0.0d0,1.0d0)
    complex(8) :: denom

    w1 = this%ip%Zp * this%ip%Zt * e2 * sqrt(0.5d0 * this%ip%rmass) / hbar
    w2 = sqrt(2.0d0 * this%ip%rmass) / hbar

    do n=1, this%ip%Nch
      E_C = E - this%vcoup_mat%e_n(n)
      eta = w1 / sqrt(E_C)
      rho = w2 * sqrt(E_C) * (this%ip%rmax - this%ip%dr)
      call dfcoul(eta, rho, fc, fpc, gc, gpc, sig0, J, iexp)
      Hp0(n) = gc(J) + ai * fc(J)
      Hn0(n) = conjg(Hp0(n))

      rho = w2 * sqrt(E_C) * (this%ip%rmax + this%ip%dr)
      call dfcoul(eta, rho, fc, fpc, gc, gpc, this%sigC(:,n), J, iexp)
      Hp2(n) = gc(J) + ai * fc(J)
      Hn2(n) = conjg(Hp2(n))

      denom = Hp0(n) * Hn2(n) - Hp2(n) * Hn0(n)
      this%Cnm(n,:) = (Hp0(n) * this%kai2(n,:) - Hp2(n) * this%kai0(n,:)) / denom
      this%Dnm(n,:) = (Hn2(n) * this%kai0(n,:) - Hn0(n) * this%kai2(n,:)) / denom
    end do
    this%Cnm = matmul(this%Cnm,this%S)
    this%Dnm = matmul(this%Dnm,this%S)

    return
  end subroutine
!----------------------------------------------------------------------!
  pure subroutine Penet(this, J, E, Rn, P)
    use global_constant, only : hbar
!   use mkl95_lapack
    use lin_coup
    implicit none
    integer :: n, i
    integer, intent(in) :: J
    class(cc_scat), intent(inout) :: this
    real(8), intent(in) :: E
    real(8), intent(out) :: Rn(this%ip%Nch), P
    real(8) :: R
    complex(8), allocatable :: Cinv(:,:)
    complex(8) :: ref_cof

    allocate(Cinv(this%ip%Nch,this%ip%Nch))
    Cinv = this%IE
!   call gesv(this%Cnm, Cinv)  ! LAPACK
    call lin_ge(this%Cnm, Cinv)
    do n=1, this%ip%Nch
      ref_cof = sum(this%Dnm(n,:) * Cinv(:,1))
      this%Sl(n,J) = - (max(E - this%vcoup_mat%e_n(n),0.0d0) / E) ** 0.25d0 * ref_cof
    end do
    Rn = abs(this%Sl(:,J)) ** 2
    R = sum(Rn)
    P = max(0.0d0, 1.0d0 - R)
    deallocate(Cinv)

    return
  end subroutine
!----------------------------------------------------------------------!
  pure subroutine scat_amp_Coul(this, E, angl, fc)
    use global_constant, only : e2, hbar
    implicit none
    class(cc_scat), intent(inout) :: this
    integer :: n
    real(8), intent(in) :: E, angl(this%ip%nth)
    real(8) :: eta, k0, sin_2
    complex(8), intent(out) :: fc(this%ip%nth)
    complex(8), parameter :: ai = (0.0d0,1.0d0)

    k0 = sqrt(2.0d0 * this%ip%rmass * E) / hbar
    eta = this%ip%Zp * this%ip%Zt * e2 * sqrt(0.5d0 * this%ip%rmass / E) / hbar
    do n=1, this%ip%nth
      sin_2 = sin(0.5d0 * angl(n)) ** 2
      fc(n) = - eta / (2.0d0 * k0 * sin_2)                &
              * exp(- ai * (eta * log(sin_2) - 2.0d0 * this%sigC(0,1)))
    end do

    return
    end subroutine
!----------------------------------------------------------------------!
  pure subroutine scat_amp_nucl(this, E, angl, fN)
    use global_constant, only : hbar
    use special_fct
    implicit none
    class(cc_scat), intent(inout) :: this
    integer :: i, n, J
    real(8), intent(in) :: E, angl(this%ip%nth)
    real(8) :: k0
    complex(8), intent(out) :: fN(this%ip%Nch,this%ip%nth)
    complex(8), parameter :: ai = (0.0d0,1.0d0)

    k0 = sqrt(2.0d0 * this%ip%rmass * E) / hbar
    fN = 0.0d0
    do i=1, this%ip%nth
      do J=0, this%ip%Jmax
        fN(1,i) = fN(1,i) + dble(2 * J + 1) * (this%Sl(1,J) - 1.0d0)   &
                  * exp(2.0d0 * ai * this%sigC(J,1))                &
                  * PL(J,cos(angl(i)))
      end do
      do n=2, this%ip%Nch
        do J=0, this%ip%Jmax
          fN(n,i) = fN(n,i) + dble(2 * J + 1) * this%Sl(n,J)           &
                    * exp(ai * (this%sigC(J,1)+this%sigC(J,n)))  &
                    * PL(J,cos(angl(i)))
        end do
      end do
    end do
    fN = fN / (2.0d0 * ai * k0)
  end subroutine
!----------------------------------------------------------------------!
  pure subroutine Ruthefrd_sig(this, E, nth, angl, dsig_R)
    use global_constant, only : e2, hbar
    implicit none
    class(cc_scat), intent(inout) :: this
    integer, intent(in) :: nth
    real(8), intent(in) :: E, angl(nth)
    real(8), intent(out) :: dsig_R(nth)
    real(8) :: eta, k

    k = sqrt(2.0d0 * this%ip%rmass * E) / hbar
    eta = this%ip%Zp * this%ip%Zt * e2 * sqrt(0.5d0 * this%ip%rmass / E) / hbar
    dsig_R = (0.5d0 * eta / (k * sin(0.5d0*angl) ** 2)) ** 2

    return
  end subroutine
end module

