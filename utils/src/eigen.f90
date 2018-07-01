      module eigen
      public :: tred2, eigval_u, eigval_d, triql, houshld_ql, jacobi
      private :: set_limit, num_ch_sign, bisec, pythag
!----------------------------------------------------------------------!
!     This program is a module subroutine which seeks eigenvalues      !
!     (and eigenvectors if you want) by Householder-QL method.         !
!     ** USAGE **                                                      !
!     use Eigen, only : Houshld_QL                                     !
!     -- Input parameters --                                           !
!     N  : Dimension of the matrix                                     !
!     A  : Matrix                                                      !
!     d  : Eigen values                                                !
!     x  : Eigen vectors                                               !
!----------------------------------------------------------------------!
      contains
!**********************************************************************!
      subroutine houshld_ql(a, d, c)
      implicit none
      integer :: N
      real(8), intent(inout) :: a(:,:)
      real(8), intent(out) :: d(:)
      real(8), allocatable :: sub_diag(:)
      character(len=1), intent(in) :: c

      N = size(a,1)
      allocate(sub_diag(N))
      if (c == 'V' .or. c == 'v') then
        call tred2(a, d, sub_diag)
        call triql(d, sub_diag, a)
      else if (c == 'N' .or. c == 'n') then
        call tred2(a, d, sub_diag, .true.)
        call triql(d, sub_diag)
      else
        stop "houshld_ql : third augument must be 'V' or 'N' !"
      end if
      deallocate(sub_diag)

      return
      end subroutine
!**********************************************************************!
      subroutine tred2(a, d, e, novectors)
      implicit none
      real(8), intent(inout) :: a(:,:)
      real(8), intent(out), dimension(:) :: d, e
      logical, optional, intent(in) :: novectors
      integer :: i, j, l, n
      real(8) :: f, g, h, hh, scale
      real(8), dimension(size(a,1)) :: gg
      logical, save :: yesvec = .true.

      n = size(a,1)
      if (present(novectors)) yesvec = .not. novectors
      do i=n, 2, -1
        l = i - 1
        h = 0.0d0
        if (l > 1) then
          scale = sum(abs(a(i,1:l)))
          if (scale == 0.0d0) then
            e(i) = a(i,l)
          else
            a(i,1:l) = a(i,1:l) / scale
            h = sum(a(i,1:l) ** 2)
            f = a(i,l)
            g = - sign(sqrt(h),f)
            e(i) = scale * g
            h = h - f * g
            a(i,l) = f - g
            if (yesvec) a(1:l,i) = a(i,1:l) / h
            do j=1, l
              e(j) = (dot_product(a(j,1:j),a(i,1:j))        &
                   + dot_product(a(j+1:l,j),a(i,j+1:l))) / h
            end do
            f = dot_product(e(1:l),a(i,1:l))
            hh = f / (h + h)
            e(1:l) = e(1:l) - hh * a(i,1:l)
            do j=1,l
              a(j,1:j) = a(j,1:j) - a(i,j) * e(1:j) - e(j) * a(i,1:j)
            end do
          end if
        else
          e(i) = a(i,l)
        end if
        d(i) = h
      end do

      if (yesvec) d(1) = 0.0d0
      e(1) = 0.0d0
      do i=1, n
        if (yesvec) then
          l = i - 1
          if (d(i) /= 0.0d0) then
            gg(1:l) = matmul(a(i,1:l),a(1:l,1:l))
            a(1:l,1:l) = a(1:l,1:l) - spread(a(1:l,i),2,l) * spread(gg(1:l),1,l)
          end if
          d(i) = a(i,i)
          a(i,i) = 1.0d0
          a(i,1:l) = 0.0d0
          a(1:l,i) = 0.0d0
        else
          d(i) = a(i,i)
        end if
      end do

      end subroutine
!**********************************************************************!
      pure subroutine eigval_d(N, Nreq, epsr, Diag, Sub_Diag, Eig)
      implicit none
      integer, intent(in) :: N, Nreq
      integer :: Nsl, Nsr, Nl, Nr, Nm, m
      real(8), intent(in) :: Diag(N), Sub_Diag(N), epsr
      real(8), intent(out) :: Eig(N)
      real(8) :: a, b, Xsr, x, Xsl, G
      real(8) :: Xr, Xl, Sub_Diag2(N), Xm, NCS_a, Xl0

      call  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      m = 0
      Xsl = a
      Xsr = b
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, a, Nsl, G)
      if (G == 0.0d0) then
        Eig(N) = a
        m = m + 1
      end if
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, b, Nsr, G)
      if (G == 0.0d0) then
        Eig(1) = b
        m = m + 1
      end if
      NCS_a = Nsl
       
      if (m == Nreq) return
      do
        Xl = Xsl ; Xr = Xsr
        Nl = Nsl ; Nr = Nsr
        if (Nsl - Nsr == 1) then
          Xl0 = Xl
          call  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
          m = m + 1
          EIG(m) = x
          Xr = Xl0
          Xl = a
          Nr = Nl
          Nl = NCS_a
          if (m == Nreq) return
        end if
         
        Xm = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, Xm, Nm, G)

        if (Nl > Nm) then
          Xsl = Xl ; Xsr = Xm
          Nsl = Nl ; Nsr = Nm
        end if
        if (Nm > Nr) then
          Xsl = Xm ; Xsr = Xr
          Nsl = Nm ; Nsr = Nr
        end if
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine eigval_u(N, Nreq, epsr, Diag, Sub_Diag, Eig)
      implicit none
      integer, intent(in) :: N, Nreq
      integer :: Nsl, Nsr, Nl, Nr, Nm, m
      real(8), intent(in)  :: Diag(N), Sub_Diag(N), epsr
      real(8), intent(out) :: Eig(N)
      real(8) :: a, b, Xsr, x, Xsl, G
      real(8) :: Xr, Xl, Sub_Diag2(N), Xm, NCS_b, Xr0

      call  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      m = 0
      Xsl = a
      Xsr = b
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, a, Nsl, G)
      if (G == 0.0d0) then
        Eig(1) = a
        m = m + 1
      end if
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, b, Nsr, G)
      if (G == 0.0d0) then
        Eig(N) = b
        m = m + 1
      end if
      NCS_b = Nsr
       
      if (m == Nreq) return
      do
        Xl = Xsl ; Xr = Xsr
        Nl = Nsl ; Nr = Nsr
        if (Nsl - Nsr == 1) then
          Xr0 = Xr
          call  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
          m = m + 1
          EIG(m) = x
          Xr = b
          Xl = Xr0
          Nl = Nr
          Nr = NCS_b
          if (m == Nreq) return
        end if
        Xm = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, Xm, Nm, G)
        if (Nm > Nr) then
          Xsl = Xm ; Xsr = Xr
          Nsl = Nm ; Nsr = Nr
        end if
        if (Nl > Nm) then
          Xsl = Xl ; Xsr = Xm
          Nsl = Nl ; Nsr = Nm
        end if
      end do

      return
      end subroutine
!**********************************************************************!
      pure subroutine bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
      implicit none
      integer, intent(in) :: N, Nl, Nr
      integer :: Nm
      real(8), intent(inout) :: Xl, Xr
      real(8), intent(in) :: Diag(N), Sub_Diag2(N), epsr
      real(8), intent(out) :: x
      real(8), parameter :: epsa=1.0d-300
      real(8) :: G
        
      do
        x = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, x, Nm, G)
        if (Nm == Nr) then
          Xr = x
        else if (Nm == Nl) then
          Xl = x
        end if
        if (Xr - Xl < epsa + epsr * (abs(Xr) + abs(Xl))) exit
      end do
      
      return
      end subroutine
!**********************************************************************!
      pure subroutine set_limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real(8), intent(in) :: Diag(N), Sub_Diag(N)
      real(8), intent(out) :: Sub_Diag2(N), a, b

      a = Diag(1) - abs(Sub_Diag(1))
      b = Diag(1) + abs(Sub_Diag(1))
      
      do i=2, N-1
         a = Min(a, Diag(i) - abs(Sub_Diag(i-1)) - abs(Sub_Diag(i)))
         b = Max(b, Diag(i) + abs(Sub_Diag(i-1)) + abs(Sub_Diag(i)))
         Sub_Diag2(i-1) = Sub_Diag(i-1) * Sub_Diag(i-1)
      end do

      a = Min(a, Diag(N) - abs(Sub_Diag(N-1)))
      b = Max(b, Diag(N) + abs(Sub_Diag(N-1)))
      Sub_Diag2(N-1) = Sub_Diag(N-1) * Sub_Diag(N-1)

      return
      end subroutine
!**********************************************************************!
      pure subroutine num_ch_sign(N, Diag, Sub_Diag2, x, NCS, G)
      implicit none
      integer, intent(in) :: N
      integer, intent(inout) :: NCS
      integer :: i
      real(8), intent(in) :: Diag(N), Sub_Diag2(N), x
      real(8), intent(out) :: G
      real(8), parameter :: eps = 1.0d-10
      real(8) :: G0
      
      G0 = x - Diag(1)

      if (G0 < 0.0d0) then
        NCS = 1      
      else 
        NCS = 0
      end if

      do i=2, N
        if (G0 == 0.0d0) G0 = eps
        G = x - Diag(i) - Sub_Diag2(i-1) / G0
        if (G < 0.0d0) NCS = NCS + 1
        G0 = G
      end do

      return
      end subroutine
!**********************************************************************!
      subroutine triql(d, e, z)
      implicit none
      logical :: vec
      integer :: n
      integer :: i, iter, l, m
      real(8), intent(inout), dimension(:) :: d, e
      real(8), optional, intent(inout) :: z(:,:)
      real(8) :: b, c, dd, f, g, p, r, s, ff(size(d))

      n = size(d)
      vec = present(z)
      e = eoshift(e, 1)
      do l=1, n
        iter = 0
 it:    do
          do m=l, n-1
            dd = abs(d(m)) + abs(d(m+1))
            if (abs(e(m)) + dd == dd) exit
          end do
          if (m == l) exit it
          if (iter == 30) stop 'too many iterations in triql'
          iter = iter + 1
          g = (d(l+1) - d(l)) / (2.0d0 * e(l))      ! shift of origin
          r = pythag(g,1.0d0)                       ! absolutly smaller
          g = d(m) - d(l) + e(l) / (g + sign(r,g))  ! eigenvalue
          s = 1.0d0
          c = 1.0d0
          p = 0.0d0
          do i=m-1, l, -1
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r == 0.0d0) then
              d(i+1) = d(i+1) - p
              e(m) = 0.0d0
              cycle it
            end if
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
            if (vec) then
              ff = z(:,i+1)
              z(:,i+1) = s * z(:,i) + c * ff
              z(:,i) = c * z(:,i) - s * ff
            end if
          end do
          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0d0
        end do it
      end do

      return
      end subroutine
!**********************************************************************!
      pure function pythag(a, b) result(c)
      implicit none
      real(8), intent(in) :: a, b
      real(8) :: c, absa, absb

      absa = abs(a)
      absb = abs(b)
      if (absa > absb) then
        c = absa * sqrt(1.0d0 + (absb / absa) ** 2)
      else if (absb /= 0.0d0) then
        c = absb * sqrt(1.0d0 + (absa / absb) ** 2)
      else
        c = 0.0d0
      end if
        
      return
      end function
!**********************************************************************!
      subroutine jacobi(N, a, epsr, x)
      implicit none
      integer :: i, j, p, q
      integer, intent(in)  :: N
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: x(N,N)
      real(8), intent(in) :: epsr
      real(8) :: amax, Dmax
      real(8) :: alpha, c, s, t, tau, h
      real(8) :: apj, ajp, ajq, aqj, xpj, xjp, xqj, xjq
      
      x = 0.0d0
      do i=1, N
        x(i,i) = 1.0d0
      end do

      do 
        p = 1
        q = 2
        amax = abs(a(q,p))
        Dmax = abs(a(N,N))
        do j=1, N-1
          do i=j+1, N
            if (abs(a(i,j)) > amax) then
              q = i
              p = j
              amax = abs(a(i,j))
            end if
          end do
          Dmax = max(Dmax, a(j,j))
        end do
        if (amax < epsr * Dmax)  exit

        alpha = (a(q,q) - a(p,p)) / (2.0d0 * a(q,p))
        t = sign(1.0d0,alpha) / (abs(alpha) + sqrt(1.0d0+alpha*alpha))
        c = 1.0d0 / sqrt(1.0d0 + t * t)
        s = t * c
        tau = s / (1.0d0 + c)
        h = t * a(q,p)
    
        a(p,p) = a(p,p) - h
        a(q,q) = a(q,q) + h
        a(q,p) = 0.0d0
 
        do j=1, p-1
          ajp = a(p,j)
          ajq = a(q,j)
          a(p,j) = ajp - s * (ajq + tau * ajp)
          a(q,j) = ajq + s * (ajp - tau * ajq)
        end do
        do j=p+1, q-1
          apj = a(j,p)
          ajq = a(q,j)
          a(j,p) = apj - s * (ajq + tau * apj)
          a(q,j) = ajq + s * (apj - tau * ajq)
        end do
        do j=q+1, n
          apj = a(j,p)
          aqj = a(j,q)
          a(j,p) = apj - s * (aqj + tau * apj)
          a(j,q) = aqj + s * (apj - tau * aqj)
        end do
        do j=1, n
          xjp = x(j,p)
          xjq = x(j,q)
          x(j,p) = xjp - s * (xjq + tau * xjp)
          x(j,q) = xjq + s * (xjp - tau * xjq)
        end do
      end do

      return 
      end subroutine
!**********************************************************************!
      end module
!======================================================================!
