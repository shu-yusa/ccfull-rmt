      module lin_coup
      public :: lin_ge, lin_sy, lin_he
      private :: LU_Decomp_real_ge, Par_Pivot_real_ge,        &
                 LU_Decomp_cmplx_ge ,Par_Pivot_cmplx_ge,      &
                 LDLT_Decomp_real_sy, Pivvoting_real_sy,      &
                 LDLT_Decomp_cmplx_sy, Pivvoting_cmplx_sy,    &
                 LDLT_Decomp_he, Pivvoting_he

      interface lin_ge
        module procedure lin_real_ge, lin_cmplx_ge
      end interface

      interface lin_sy
         module procedure lin_real_sy, lin_cmplx_sy
      end interface

      contains
!**********************************************************************!
      pure subroutine lin_real_ge(a, b)
      implicit none
      integer :: N, M
      integer :: i, k
      integer, allocatable :: mm(:)
      real(8), intent(inout), dimension(:,:) :: a, b
      real(8), allocatable, dimension(:) :: x

      N = size(a,1)
      M = size(b,2)
      allocate(x(N),mm(N))
      call LU_Decomp_real_ge(N, a, mm)
      do i=1, M
        do k=1, N
          x(k) = (b(mm(k),i) - dot_product(a(k,1:k-1),x(1:k-1))) / a(k,k)
        end do
        do k=N, 1, -1
          b(k,i) = x(k) - dot_product(a(k,k+1:N), b(k+1:N,i))
        end do
      end do
      deallocate(x,mm)

      return
      end subroutine
!**********************************************************************!
      pure subroutine  LU_Decomp_real_ge(N, a, mm)
      implicit none
      integer :: i, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      real(8), intent(inout) :: a(N,N)
      real(8) :: p

      forall (i=1:N) mm(i) = i
      do k=1, N-1
        call Par_Pivot_real_ge(N, a, k, mm, p)
        a(k, k+1:N) = a(k, k+1:N) / p
        do i=k+1, N
          a(i, k+1:N) = a(i, k+1:N) - a(i, k) * a(k, k+1:N)
        end do
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine  Par_Pivot_real_ge(N, a, k, mm, p)
      implicit none
      integer, intent(in) :: k, N
      integer, intent(inout) :: mm(N)
      integer :: iw, j, L
      real(8), intent(out) :: p
      real(8), intent(inout) :: a(N,N)
      real(8) :: w(N)

      L = k
      p = a(k,k)

      do j=k+1, N
        if (abs(p) < abs(a(j,k))) then
          L = j
          p = a(j,k)
        end if
      end do

      if (L /= k) then
        w(:)   = a(k,:)
        a(k,:) = a(L,:)
        a(L,:) = w(:)
        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!**********************************************************************!
      pure subroutine lin_cmplx_ge(a, b)
      implicit none
      integer :: N, M
      integer :: i, j, k
      integer, allocatable :: mm(:)
      complex(8), intent(inout), dimension(:,:) :: a, b
      complex(8), allocatable, dimension(:) :: x

      N = size(a,1)
      M = size(b,2)
      allocate(mm(N),x(N))
      call LU_Decomp_cmplx_ge(N, a, mm)
      do i=1, M
        do k=1, N
          x(k) = b(mm(k),i) 
          do j=1, k-1
            x(k) = x(k) - a(k,j) * x(j)
          end do
          x(k) = x(k) / a(k,k)
        end do
        do k=N, 1, -1
          b(k,i) = x(k) - sum(a(k,k+1:N) * b(k+1:N,i))
        end do
      end do
      deallocate(mm,x)

      return
      end subroutine
!**********************************************************************!
      pure subroutine  LU_Decomp_cmplx_ge(N, a, mm)
      implicit none
      integer :: i, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      complex(8), intent(inout) :: a(N,N)
      complex(8) :: p

      do i=1, N
        mm(i) = i
      end do
      
      do k=1, N-1
        call Par_Pivot_cmplx_ge(N, a, k, mm, p)
        a(k, k+1:N) = a(k, k+1:N) / p
        do i=k+1, N
          a(i, k+1:N) = a(i, k+1:N) - a(i, k) * a(k, k+1:N)
        end do
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine  Par_Pivot_cmplx_ge(N, a, k, mm, p)
      implicit none
      integer, intent(in) :: k, N
      integer, intent(inout) :: mm(N)
      integer :: iw, j, L
      complex(8), intent(out) :: p
      complex(8), intent(inout) :: a(N, N)
      complex(8) :: w(N)

      L = k
      p = a(k, k)

      do j=k+1, N
        if (abs(p) < abs(a(j, k))) then
          L = j
          p = a(j, k)
        end if
      end do

      if (L /= k) then
        w(:)    = a(k, :)
        a(k, :) = a(L, :)
        a(L, :) = w(:)
        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!**********************************************************************!
      subroutine lin_real_sy(a, b)
      implicit none
      integer :: N, M
      integer :: i, k
      integer, allocatable :: mm(:)
      real(8), intent(inout), dimension(:,:) :: a, b
      real(8), allocatable, dimension(:) :: x, D

      N = size(a,1)
      M = size(b,2)
      allocate(x(N),D(N),mm(N))
      call LDLT_Decomp_real_sy(N, a, D, mm)
      do i=1, M
        do k=1, N
          x(k) = b(mm(k),i) - dot_product(a(k,1:k-1), x(1:k-1))
        end do
        do k=N, 1, -1
          x(k) = x(k) / D(k) - dot_product(a(k+1:N,k), x(k+1:N))
        end do
        forall (k=1:N) b(mm(k),i) = x(k)
      end do
      deallocate(x,D,mm)

      return
      end subroutine
!**********************************************************************!
      subroutine LDLT_Decomp_real_sy(N, L, D, mm)
      implicit none
      integer :: i, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      real(8), intent(inout) :: L(N,N)
      real(8), intent(out) :: D(N)
      real(8) :: u(N)
      
      forall (i=1:N) mm(i) = i

      do k=1, N
        call  Pivotting_real_sy(N, k, L, mm, D(k))
        do i=1, k-1
          u(i) = L(k,i) - dot_product(u(1:i-1), L(i,1:i-1))
          L(k,i) = u(i) / D(i)
        end do
        D(k) = D(k) - dot_product(u(1:k-1), L(k,1:k-1))
      end do

      return
      end subroutine
!**********************************************************************!
      subroutine Pivotting_real_sy(N, k, a, mm, p)
      implicit none
      integer, intent(in) :: N
      integer :: j, L, iw
      integer, intent(in) :: k
      integer, intent(inout) :: mm(N)
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: p
      real(8) :: w(N), w2(N), w3

      L = k
      p = a(k,k)

      do j=k+1, N
        if (abs(p) < abs(a(j,j))) then
          L = j
          p = a(j,j)
        end if
      end do

!     L > k
      if (L /= k) then
!       w(:)   = a(k,:)
!       a(k,:) = a(L,:)
!       a(L,:) = w(:)
!           
!       w(:)   = a(:,k)
!       a(:,k) = a(:,L)
!       a(:,L) = w(:)

        w(1:k) = a(k,1:k)
        w(k+1:N) = a(k+1:N,k)
        w2(1:L) = a(L,1:L)
        w2(L+1:N) = a(L+1:N,L)

        a(L,:) = w
        a(k,:) = w2

        a(:,k) = w2
        w3 = a(k,k)
        a(k,k) = a(L,k)
        a(L,k) = w3

        a(:,L) = w 
        w3 = a(k,L)
        a(k,L) = a(L,L)
        a(L,L) = w3

        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!**********************************************************************!
      pure subroutine lin_cmplx_sy(a, b)
      implicit none
      integer :: N, M
      integer :: i, j, k
      integer, allocatable :: mm(:)
      complex(8), intent(inout), dimension(:,:) :: a, b
      complex(8), allocatable, dimension(:) :: x, D

      N = size(a,1)
      M = size(b,2)
      allocate(mm(N),x(N), D(N))
      call LDLT_Decomp_cmplx_sy(N, a, D, mm)
      do i=1, M
        do k=1, N
          x(k) = b(mm(k),i)
          do j=1, k-1
            x(k) = x(k) - a(k,j) * x(j)
          end do
        end do
        do k=N, 1, -1
          x(k) = x(k) / D(k)
          do j=k+1, N
            x(k) = x(k) - a(j,k) * x(j)
          end do
        end do
        forall(k=1:N) b(mm(k),i) = x(k)
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine LDLT_Decomp_cmplx_sy(N, L, D, mm)
      implicit none
      integer :: i, j, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      complex(8), intent(out) :: L(N,N), D(N)
      complex(8) :: u(N)
      
      forall (i=1:N) mm(i) = i

      do k=1, N
        call  Pivotting_cmplx_sy(N, k, L, mm, D(k))
        do i=1, k-1
          u(i) = L(k,i)
          do j=1, i-1
            u(i) = u(i) - u(j) * L(i,j)
          end do
          L(k,i) = u(i) / D(i)
        end do
        do i=1, k-1
          D(k) = D(k) - u(i) * L(k,i)
        end do
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine Pivotting_cmplx_sy(N, k, a, mm, p)
      implicit none
      integer, intent(in) :: N
      integer :: j, L, iw
      integer, intent(in) :: k
      integer, intent(out) :: mm(N)
      complex(8), intent(inout) :: a(N,N)
      complex(8), intent(out) :: p
      complex(8) :: w(N), w2(N), w3

      L = k
      p = a(k,k)

      do j=k+1, N
        if (abs(p) < abs(a(j,j))) then
          L = j
          p = a(j,j)
        end if
      end do

      if (L /= k) then
!       w(:)   = a(k,:)
!       a(k,:) = a(L,:)
!       a(L,:) = w(:)
!           
!       w(:)   = a(:,k)
!       a(:,k) = a(:,L)
!       a(:,L) = w(:)

        w(1:k) = a(k,1:k)
        w(k+1:N) = a(k+1:N,k)
        w2(1:L) = a(L,1:L)
        w2(L+1:N) = a(L+1:N,L)

        a(L,:) = w
        a(k,:) = w2

        a(:,k) = w2
        w3 = a(k,k)
        a(k,k) = a(L,k)
        a(L,k) = w3

        a(:,L) = w 
        w3 = a(k,L)
        a(k,L) = a(L,L)
        a(L,L) = w3

        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!**********************************************************************!
      subroutine lin_he(a, b)
      implicit none
      integer :: N, M
      integer :: i, j, k
      integer, allocatable :: mm(:)
      complex(8), intent(inout), dimension(:,:) :: a, b
      complex(8), allocatable, dimension(:) :: x, D

      N = size(a,1)
      M = size(b,2)
      allocate(mm(N),x(N), D(N))
      call LDLT_Decomp_he(N, a, D, mm)
      do i=1, M
        do k=1, N
          x(k) = b(mm(k),i)
          do j=1, k-1
            x(k) = x(k) - a(k,j) * x(j)
          end do
        end do
        do k=N, 1, -1
          x(k) = x(k) / D(k) - dot_product((a(k+1:N,k)),x(k+1:N))
        end do
        forall(k=1:N) b(mm(k),i) = x(k)
      end do

      return
      end subroutine
!**********************************************************************!
      subroutine LDLT_Decomp_he(N, L, D, mm)
      implicit none
      integer :: i, j, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      complex(8), intent(out) :: L(N,N), D(N)
      complex(8) :: u(N)
      
      forall (i=1:N) mm(i) = i

      do k=1, N
        call  Pivotting_he(N, k, L, mm, D(k))
        do i=1, k-1
          u(i) = conjg(L(k,i))
          do j=1, i-1
            u(i) = u(i) - u(j) * L(i,j)
          end do
          L(k,i) = conjg(u(i)) / D(i)
        end do
        do i=1, k-1
          D(k) = D(k) - u(i) * L(k,i)
        end do
      end do

      return
      end subroutine
!**********************************************************************!
      subroutine Pivotting_he(N, k, a, mm, p)
      implicit none
      integer, intent(in) :: N
      integer :: j, L, iw
      integer, intent(in) :: k
      integer, intent(out) :: mm(N)
      complex(8), intent(inout) :: a(N,N)
      complex(8), intent(out) :: p
      complex(8) :: w(N), w2(N), w3

      L = k
      p = a(k,k)

      do j=k+1, N
        if (abs(p) < abs(a(j,j))) then
          L = j
          p = a(j,j)
        end if
      end do

      if (L /= k) then
!       w(:)   = a(k,:)
!       a(k,:) = a(L,:)
!       a(L,:) = w(:)
!           
!       w(:)   = a(:,k)
!       a(:,k) = a(:,L)
!       a(:,L) = w(:)

        w(1:k) = a(k,1:k)
        w(k+1:N) = conjg(a(k+1:N,k))
        w2(1:L) = a(L,1:L)
        w2(L+1:N) = conjg(a(L+1:N,L))

        a(L,:) = w
        a(k,:) = w2

        a(:,k) = conjg(w2)
        w3 = a(k,k)
        a(k,k) = a(L,k)
        a(L,k) = w3

        a(:,L) = conjg(w)
        w3 = a(k,L)
        a(k,L) = a(L,L)
        a(L,L) = w3

        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!**********************************************************************!
      end module

