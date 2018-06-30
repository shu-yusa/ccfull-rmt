!======================================================================!
      module angular_coup
      public :: Wig3j, CGcof, W_Racah
      private :: fact, nPm, log_fact, log_nPm, Del_R, log_Del_R
      contains
!**********************************************************************!
      pure subroutine Wig3j(j1, j2, j3, m1, m2, m3, wg)
      implicit none
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      real(kind=8), intent(out) :: wg
      real(kind=8) :: CG

      call CGcof(j1, j2, j3, m1, m2, -m3, CG)
      wg = (- 1.0d0) ** (j1 - j2 - m3) * CG / sqrt(dble(2 * j3 + 1))

      return
      end subroutine
!**********************************************************************!
      pure subroutine CGcof(j1, j2, j3, m1, m2, m3, CG)
      implicit none
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      integer :: k, ka, kb, kc, kd
      real(kind=8), intent(out) :: CG
      real(kind=8) :: CG2, fa, fb, fc, fd

!     if (j1 < 0 .or. j2 < 0 .or. j3 < 0) then
!       stop 'negative angular momentum'
!     else if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) then
!       stop 'wrong angular momentum'
!     end if
! special cases
      if (m3 /= m1+m2 .or. j3 < abs(j1-j2) .or. j3 > j1+j2) then
        CG = 0.0d0
        return
      end if
      if (m1==0 .and. m2==0 .and. m3==0) then
        if (mod(j1 + j2 + j3, 2) /= 0) then
          CG = 0.0d0
        else
          kd = j1 + j2 + j3 + 1
          ka = max(j1, j2, j3)
          kc = min(j1, j2, j3)
          kb = kd - ka - kc - 1
          fa = fact(kc+kb-ka)
          fb = fa * nPm(kc+ka-kb, 2*(ka-kb))
          fc = fact((-ka+kb+kc)/2)
          fd = fc * nPm((ka+kc-kb)/2, ka-kb)
          CG = (- 1) ** ((j1+j2-j3)/2) * sqrt(dble(j3 + j3 + 1)   & 
     &         * fa * fb / nPm(kd, 2*kc+1))                       &
     &         * nPm((kd-1)/2, kc) / (fc * fd)
        end if
        return
      end if
      if (j1 == j2 .and. m1 == m2 .and. mod(j1+j1+j3, 2) /= 0) then
        CG = 0.0d0
        return
      end if
      if ((j2 == 0 .and. m2 == 0) .or. (j1 == 0 .and. m1 == 0)) then
        if (j1 + j2 == j3 .and. m1 + m2 == m3) then
          CG = 1.0d0
        else
          CG = 0.0d0
        end if
        return
      end if

      CG2 = 0.0d0
      do k=0, min(j1+j2-j3,j1-m1,j2+m2)
        ka = j3 - j2 + m1 + k
        if (ka < 0) cycle
        kb = j3 - j1 - m2 + k
        if (kb < 0) cycle
        fa = fact(ka) * fact(kb) * fact(j1+j2-j3-k)         &
     &                * fact(j1-m1-k) * fact(j2+m2-k)
        CG2 = CG2 + (-1.0d0) ** k / (fact(k) * fa)
      end do

      kd = j1 + j2 + j3 + 1
      ka = max(j1, j2, j3)
      kc = min(j1, j2, j3)
      kb = kd - ka - kc - 1
      fa = fact(kc+kb-ka)
      fb = fa * nPm(kc+ka-kb, 2*(ka-kb))
      CG = CG2 * sqrt(dble(j3 + j3 + 1) * fa * fb / nPm(kd, 2*kc+1))

      select case(m1)
        case(0)      ; CG = CG * fact(j1) 
        case default ; CG = CG * sqrt(fact(j1-m1) * fact(j1+m1))
      end select
      select case(m2)
        case(0)      ; CG = CG * fact(j2) 
        case default ; CG = CG * sqrt(fact(j2-m2) * fact(j2+m2))
      end select
      select case(m3)
        case(0)      ; CG = CG * fact(j3) 
        case default ; CG = CG * sqrt(fact(j3-m3) * fact(j3+m3))
      end select

      return
      end subroutine
!**********************************************************************!
      pure function fact(n) result(g)
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(kind=8) :: g

!     if (n < 0) stop 'negative in factorial'
      g = 1.0d0
      do i=2, n
        g = g * dble(i)
      end do

      return
      end function
!**********************************************************************!
      pure function nPm(n, m) result(f)
      implicit none
      integer, intent(in) :: n, m
      integer :: i
      real(kind=8) :: f

!     if (n < 0 .or. m < 0) stop 'negative in permutation'
      if (m == 0) then
        f = 1.0d0
        return
      end if
      f = n - m + 1
      do i=n-m+2, n
        f = f * dble(i)
      end do

      return
      end function
!**********************************************************************!
    function log_fact(n) result(f)
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(8) :: f

      f = 0.0d0
      do i=2, n
        f = f + log(dble(i))
      end do

    end function

    function log_nPm(n,m) result(f)
      implicit none
      integer, intent(in) :: n,m
      integer :: i
      real(8) :: f

      if (m == 0) then
        f = 0.0d0
        return
      end if
      f = log(dble(n-m+1))
      do i=n-m+2, n
        f = f + log(dble(i))
      end do

    end function

    function Del_R(a,b,c) result(d)
      implicit none
      integer, intent(in) :: a, b, c
      real(8) :: d

      d = exp(log_Del_R(a,b,c))
    end function

    function log_Del_R(a,b,c) result(d)
      implicit none
      integer, intent(in) :: a, b, c
      integer :: ka, kb, kc
      real(8) :: d

      ka = max(a, b, c)
      kc = min(a, b, c)
      kb = a + b + c - ka - kc
      d = 0.5d0*(log_fact(ka-kb+kc)+log_fact(-ka+kb+kc)-log_nPm(a+b+c+1,2*kc+1))
    end function

    function W_racah(a,b,c,d,e,f) result(W)
      implicit none
      integer, intent(in) :: a, b, c, d, e, f
      integer :: k, kmax, kmin
      real(8) :: W, s, logD

!     if (e == 0) then
!       if (a /= b .or. c /= d) then
!         W = 0.0d0
!       else
!         W = (-1.0d0)**(f-b-d)/sqrt(dble((2*b+1)*(2*d+1)))
!       end if
!       return
!     else if (f == 0) then
!       if (a /= c .or. b /= d) then
!         W = 0.0d0
!       else
!         W = (-1.0d0)**(e-c-d)/sqrt(dble((2*c+1)*(2*d+1)))
!       end if
!       return
!     else if (a == 0) then
!       if (e /= b .or. c /= f) then
!         W = 0.0d0
!       else
!         W = (-1.0d0)**(e-b)/sqrt(dble((2*b+1)*(2*f+1)))
!       end if
!       return
!     else if (b == 0) then
!       if (a /= e .or. f /= d) then
!         W = 0.0d0
!       else
!         W = (-1.0d0)**(f-b-d)/sqrt(dble((2*e+1)*(2*d+1)))
!       end if
!       return
!     end if
      kmin = max(a+b+e,c+d+e,a+c+f,b+d+f)
      kmax = min(a+b+c+d,a+d+e+f,b+c+e+f)

      logD = log_Del_R(a,b,e)+ log_Del_R(c,d,e)+ log_Del_R(a,c,f)+ log_Del_R(b,d,f)
      W = 0.0d0
      do k=kmin, kmax
        s = logD + log_fact(k+1)   &
           - log_fact(k-a-b-e)-log_fact(k-c-d-e)  &
           - log_fact(k-a-c-f)-log_fact(k-b-d-f)  &
           - log_fact(a+b+c+d-k)-log_fact(a+d+e+f-k)-log_fact(b+c+e+f-k)
        W = W + (-1.0d0)**(k+a+b+c+d)*exp(s)
      end do

    end function
      end module
