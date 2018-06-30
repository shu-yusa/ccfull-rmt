      module print_array
      interface print_mat
        module procedure print_real_mat, print_cmplx_mat, &
                         print_logical_mat
      end interface

      interface print_vec
        module procedure print_real_vec, print_cmplx_vec, &
                         print_logical_vec
      end interface
      
      contains
!**********************************************************************!
      subroutine print_real_mat(a, p)
      implicit none
      integer, optional, intent(in) :: p
      integer :: i, j, m, L, L2
      real(8), intent(in) :: a(:,:)
      character(len=30) :: FM, c, d, e, f
      character(len=20), allocatable :: g(:)

      L = size(a,1)
      L2 = size(a,2)

      if (min(L,L2) <= 0) then
        print *, 'print_mat : wrong dimension of array'
        return
      else if (max(L,L2) > 30) then
        print *, 'print_mat : too large dimension for printing'
        return
      end if


      write(c,*) L2
      m = max(int(log10(maxval(abs(a*(1.0d0+1.0d-14))))) + 1,1)

      if (present(p)) then
        if (p > 0) then
          write(e,*) p
          write(d,*) p + 8
          e = 'es'//trim(adjustl(d))//'.'//trim(adjustl(e))
          FM = '(1x,a,'//trim(adjustl(c))//trim(e)//',a)'
        else
          print *, 'print_real_mat : wrong input value for format.'
          return
        end if
      else
        if (m >= 5) then
          FM = '(1x,a,'//trim(adjustl(c))//'es12.4,a)'
        else
          write(d,*) m + 7
          e = 'f'//trim(adjustl(d))//'.4'
          FM = '(1x,a,'//trim(adjustl(c))//trim(e)//',a)'
        end if
      end if

      allocate(g(L))
      do i=1, L
        write(g(i),*) i
        g(i) = adjustl(g(i))
      end do
      m = len_trim(g(L))

      do i=1, L
        write(6,FM) adjustr(g(i)(1:m))//' |', (a(i,j), j=1, L2), ' |'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      subroutine  print_cmplx_mat(a, p)
      implicit none
      integer, optional, intent(in) :: p
      integer :: i, j, m, L, L2
      complex(8), intent(in) :: a(:,:)
      character(len=30), parameter :: FM1='(1x,a)'
      character(len=40) :: FM2, c, d, f
      character(len=20), allocatable :: g(:)

      L = size(a,1)
      L2 = size(a,2)

      if (min(L,L2) <= 0) then
        print *, 'print_mat : wrong dimension of array'
        return
      else if (max(L,L2) > 15) then
        print *, 'print_mat : too large dimension for printing'
        return
      end if
      m = max(max(int(log10(maxval(abs(dble(a))))),int(log10(maxval(abs(aimag(a)))))) + 1,1)

      if (present(p)) then
        if (p > 0) then
          write(c,*) p
          write(d,*) p + 7
          c = trim(adjustl(d))//'.'//trim(adjustl(c))
          FM2 = '(a,es'//trim(c)//',",",es'//trim(c)//',a)' 
        else
          print *, 'print_mat : wrong input value for format.'
          return
        end if
      else
        if (m >= 5) then
          FM2 = '(a,es11.4,",",es11.4,a)'
        else
          write(d,*) m + 6
          c = 'f'//trim(adjustl(d))//'.4'
          FM2 = '(a,'//trim(c)//',",",'//trim(c)//',a)'
        end if
      end if

      allocate(g(L))
      do i=1, L
        write(g(i),*) i
        g(i) = adjustl(g(i))
      end do
      m = len_trim(g(L))

      do i=1, L
        write(6,FM1,advance='no') adjustr(g(i)(1:m))//' | '
        do j=1, L2
          write(6,FM2,advance='no') '(',a(i,j), ') '
        end do
        write(6,'(a)') '|'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      subroutine print_logical_mat(a)
      implicit none
      logical, intent(in) :: a(:,:)
      integer :: i, j, L, L2
      character(len=20), allocatable :: c(:)
      character(len=20) :: FM

      L = size(a,1)
      L2 = size(a,2)

      if (min(L,L2) <= 0) then
        print *, 'print_mat : wrong dimension of array'
        return
      else if (max(L,L2) > 80) then
        print *, 'print_mat : too large dimension for printing'
        return
      end if

      allocate(c(max(L,L2)))
      do i=1, max(L,L2)
        write(c(i),*) i
        c(i) = adjustl(c(i))
      end do

      FM = '(5x,'//c(L2)(1:2)//'a)'

      write(6,FM) (' '//adjustr(c(i)(1:2)),i=1,L2)

      FM = '(1x,a,'//c(L2)(1:2)//'l3,a)'

      do i=1, L
        write(6,FM) adjustr(c(i)(1:2))//' |', (a(i,j), j=1, L2), ' |'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      subroutine print_real_vec(v, p)
      implicit none
      integer, optional, intent(in) :: p
      integer :: i, m, L
      real(8), intent(in) :: v(:)
      character(len=30) :: FM, c, d, f
      character(len=20), allocatable :: g(:)

      L = size(v)

      if (L <= 0) then
        print *, 'print_vec : wrong dimension of array'
        return
      end if
      m = max(int(log10(maxval(abs(v)))) + 1,1)
!     m = max(int(log10(maxval(abs(v*(1.0d0+1.0d-14))))) + 1,1)

      if (present(p)) then
        if (p > 0) then
          write(c,*) p
          write(d,*) p + 7
          c = trim(adjustl(d))//'.'//trim(adjustl(c))
          FM = '(1x,a,es'//trim(c)//',a)' 
        else
          print *, 'print_real_vec : wrong input value for format.'
          return
        end if
      else
        if (m > 6) then
          FM = '(1x,a,es12.5,a)' 
        else
          write(d,*) m + 7
          d = 'f'//trim(adjustl(d))//'.5'
          FM = '(1x,a,'//trim(d)//',a)'
        end if
      end if

      allocate(g(L))
      do i=1, L
        write(g(i),*) i
        g(i) = adjustl(g(i))
      end do
      m = len_trim(g(L))

      do i=1, L
        write(6,FM) adjustr(g(i)(1:m))//' | ',dble(v(i)),' |'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      subroutine print_cmplx_vec(v, p)
      implicit none
      integer, optional, intent(in) :: p
      integer :: i, m, L
      complex(8), intent(in) :: v(:)
      character(len=30) :: FM, c, d, f
      character(len=20), allocatable :: g(:)

      L = size(v)

      if (L <= 0) then
        print *, 'print_vec : wrong dimension of array'
        return
      end if
      m = max(max(int(log10(maxval(abs(dble(v))))),int(log10(maxval(abs(aimag(v)))))) + 1,1)

      if (present(p)) then
        if (p > 0) then
          write(c,*) p
          write(d,*) p + 7
          c = trim(adjustl(d))//'.'//trim(adjustl(c))
          FM = '(1x,a,es'//trim(c)//',a,es'//trim(c)//',a)' 
        else
          print *, 'print_vec : wrong input value for format.'
          return
        end if
      else
        if (m > 5) then
          FM = '(1x,a,es12.5,a,es12.5,a)'
        else
          write(d,*) m + 6
          c = 'f'//trim(adjustl(d))//'.4'
          FM = '(1x,a,'//trim(c)//',a,'//trim(c)//',a)'
        end if
      end if

      allocate(g(L))
      do i=1, L
        write(g(i),*) i
        g(i) = adjustl(g(i))
      end do
      m = len_trim(g(L))

      do i=1, L
        write(6,FM) adjustr(g(i)(1:m))//' | (',dble(v(i)),',',aimag(v(i)), ') |'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      subroutine print_logical_vec(v)
      implicit none
      logical, intent(in) :: v(:)
      integer :: i, j, L, k
      character(len=20), allocatable :: c(:)
      character(len=20) :: FM

      L = size(v)

      if (L <= 0) then
        print *, 'print_vec : wrong dimension of array'
        return
      end if

      allocate(c(L))
      do i=1, L
        write(c(i),*) i
        c(i) = adjustl(c(i))
      end do
      k = len_trim(c(L))

      FM = '(1x,a,l3,a)'

      do i=1, L
        write(6,FM) adjustr(c(i)(1:k))//' |', v(i), '  |'
      end do
      write(6,*)

      return
      end subroutine
!**********************************************************************!
      end module

