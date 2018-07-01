!======================================================================!
      module diff
!----------------------------------------------------------------------!
!     This is a module subroutine which computes a derivative of a     !
!     given function f. To do this, we compute a polynomial of order n !
!     which coincides with f at n+1 points x1, x2, ..., x(n+1). The    !
!     polynomial can be rearranged as a polynomial of (x - x0) where   !
!     x0 is a point at which we want to compute the derivative. By     !
!     comparing the polynomial and the n-th term of the Talyor         !
!     expansion of the function f, we obtains the formula to compute   !
!     the derivative of f.                                             !
!                                                                      !
!     ** USAGE **                                                      !
!     use diff, only : derive                                          !
!     ** imput parameters **                                           !
!     N ... The number of data points. integer, intent(in) .           !
!     Nd ... The number of points at where we compute the derivative.  !
!            integer, intent(in).                                      !
!     Lmax ... Maximum order of the derivative. integer, intent(in).   !
!     m ... The number of points neighboring to the x0.                !
!           integer, intent(in).                                       !
!     x(1:N) ... position of data points. real(8), intent(in)          !
!     f(1:N) ... The value of the function at data points x(1:N).      !
!                real(8), intent(in)                                   !
!     d(1:Nd) ... The points at where we compute the derivative.       !
!                 real(8), intent(in)                                  !
!     df(Nd,0:Lmax) ... The derivative of the function f of order L at !
!                       the points d. real(8), intent(out).            !
!----------------------------------------------------------------------!
      contains
!**********************************************************************!
      subroutine derive(N, Nd, Lmax, m, x, f, d, df)
      implicit none
      integer, intent(in) :: N, Nd, Lmax, m
      integer, parameter :: sentinel = 1.0d2
      integer :: i, j, k, L
      real(8), intent(in) :: x(N), f(N), d(Nd)
      real(8), intent(out) :: df(Nd,0:Lmax)
      real(8), allocatable :: dvdf(:,:)
      real(8) :: x0, xk, Zk, fact, dd, b(-1:Lmax)

      allocate(dvdf(0:N,0:m))
      dvdf(1:N,0) = f(1:N)
      do L=1, m
        do k=1, N-L
          dvdf(k,L) = (dvdf(k+1,L-1) - dvdf(k,L-1)) / (x(k+L) - x(k))
        end do
      end do
      do L=1, m
        dvdf(0,L) = sentinel * max(abs(dvdf(1,L)), 1.0d0)
        dvdf(N-L+1,L) = sentinel * max(abs(dvdf(N-L,L)), 1.0d0)
      end do

      df = 0.0d0
      do i=1, Nd
        x0 = d(i)
        if (x0 < x(1) .or. x0 > x(N)) stop 'impossible to derive'
        do j=1, N
          if (x(j) > x0) then
            exit
          end if
        end do
        j = j - 1  
        b(-1) = 0.0d0
        b(0) = 1.0d0
        b(2:Lmax) = 0.0d0
        df(i,0) = f(j)
        xk = x(j)
        do k=1, m
          Zk = x0 - xk
          if (abs(dvdf(j-1,k)) < abs(dvdf(j,k))) then
            j = j - 1
            xk = x(j)
          else
            xk = x(j+k)
          end if
          dd = dvdf(j,k)          
!         if (abs(dd) > 1.0d2) stop 'Large divided difference'
          if (abs(dd) > 1.0d2) write(6,*) 'Large divided difference'
          do L=min(k,Lmax), 0, -1
            b(L) = b(L-1) + b(L) * Zk
            df(i,L) = df(i,L) + b(L) * dd
          end do
        end do
        fact = 1.0d0
        do L=2, Lmax
          fact = dble(L) * fact
          df(i,L) = fact * df(i,L)
        end do
      end do
      deallocate(dvdf)

      return
      end subroutine
!**********************************************************************!
      end module
!======================================================================!
