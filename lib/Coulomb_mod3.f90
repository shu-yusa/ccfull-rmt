!======================================================================!
      module Coulomb2
      public :: dfcoul
      private :: jflgam, yfclen, yfasym, yfireg, yfrica, dfcz0, jfdelg
!======================================================================!
!     * no implicit definition
!     * all internal functions are used as generic name
!     * modified constants(attached 'd0')
!======================================================================!
!----------------------------------------------------------------------!
! Subroutine for the Coulomb wave function                             !
!                                                                      !
!     eta  --- Sommerfeld parameter. real(8), intent(in)               !
!     rho  --- Independent variable. real(8), intent(in)               !
!     fcw(0:L) --- Regular Coulomb wave function of angular momentum L.!
!                  real(8), intent(out)                                !
!     fpcw(0:L) --- Derivative of regular Coulomb wave function.       !
!                  real(8), intent(out)                                !
!     gcw(0:L)  --- Irregular Coulomb wave function.                   !
!                   real(8), intent(out)                               !
!     gpcw(0:L) --- Derivative of irregular Coulomb wave function.     !
!                   real(8), intent(out)                               !
!     sigmad(0:L) --- Coulomb phase shift. real(8), intent(out)        !
!     L --- Angular momentum. We can compute above quantities up to    !
!           this angular momentum. integer, intent(in)                 !
!     iexp(0:L) --- ???. integer, intent(out)                          !
!----------------------------------------------------------------------!
      contains
!**********************************************************************!
      subroutine jflgam(xd, yd, par, pai, nbchif)
      implicit none
      integer :: i, k, kr
      integer, intent(out) :: nbchif
      real(8), intent(in) :: xd, yd
      real(8), intent(out) :: par, pai
      real(8), parameter :: hlo2pi=0.918938533204672d0
      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: pis2=1.570796326794897d0
      real(8), parameter :: pis4=0.785398163397448d0
      real(8), parameter :: alo2pi=1.837877066409345d0
      real(8), parameter :: rac2=0.3465735902799726d0
      real(8), parameter :: depi=6.283185307179586d0
      real(8), parameter :: alopi=1.1447298858494002d0
      real(8), parameter :: supint=2147483647.0d0
      real(8) :: test(7)=(/2.9152d+7, 2.2958d+3, 1.4124d+2, 3.9522d+1, &
     &                     19.6611d0, 12.791d0,-10.d0/)
      real(8) :: c(6)=(/8.333333333333333d-2, -2.777777777777777d-3,   &
     &                  7.936507936507937d-4, -5.952380952380952d-4,   &
     &                  8.417508417508418d-4, -1.917526917526918d-3/)
      real(8) :: x, y, u, v, tra, tra1, cosi, sini, cos2i, sin2i
      real(8) :: zmod, xx
      
      nbchif = 15
      x = abs(xd)
      xx = x
      if (yd) 1,2,1
    1 y = abs(yd)
      kr = 1
      i = mod(10.99d0-x,supint)
!     translation
      if (i) 3,3,4
    4 tra = i
      x = x + tra
!     logarithme(x+iy) (x,y positifs)
    3 if (x - y) 5,6,7
    5 tra1 = x / y
      if (tra1) 8,8,9
    8 u = log(y)
      v = pis2
      go to 10
    6 u = rac2 + log(x)
      v = pis4
      go to 10
    9 tra = y * sqrt(1.d0 + tra1 * tra1)
      tra1 = y / x
   11 u = log(tra)
      v = atan(tra1)
   10 go to (12,19,23),kr
    7 tra1 = y / x
      tra=x*sqrt(1.d0+tra1*tra1)
      go to 11
   12 kr=2
!     developpement asymptotique ( x superieur a 10 )
      tra=x-0.5d0
      pai=v*tra+y*(u-1.d0)
      par=-x+hlo2pi+u*tra-v*y
      zmod=x+y
      if(zmod-test(1))13,13,14
   13 tra=x*x+y*y
      cosi=x/tra
      sini=y/tra
      sin2i=(sini*cosi)+(sini*cosi)
      cos2i=(cosi+sini)*(cosi-sini)
      k=1
      go to 15
   16 tra=cosi*cos2i-sini*sin2i
      sini=sini*cos2i+cosi*sin2i
      cosi=tra
   15 par=par+c(k)*cosi
      pai=pai-c(k)*sini
      k=k+1
      if(zmod-test(k))16,16,14
!     translation inverse
   17 i=i-1
      x=i
      x=xx+x
      go to 3
   19 par=par-u
      pai=pai-v
   14 if(i-1)18,60,17
   60 if(xd)17,61,17
!     controle du quadrant
   18 if(xd)20,61,21
   61 tra=pi*y
      if(tra-1.d-2)300,300,301
  300 par= tra*(2.d0+tra*(-2.d0+tra*(1.333333333333333d0+tra*(&
     &-0.6666666666666666d0+tra*(0.2666666666666666d0+tra*(   &
     &-0.08888888888888888d0+tra*0.02539682539682540d0))))))
      tra1=-log(y)-log(par)
      go to 302
  301 par=1.d0-exp(-tra-tra)
      tra1=-log(y*par)
  302 par=0.5d0*(alo2pi-tra+tra1)
      pai=pai-pis2
   21 if(yd)28,100,100
!     x+iy change en -x-iy
   20 tra=pi*y
      par=alo2pi-u-par-tra
      pai=pi-v-pai
      tra=exp(-tra-tra)
      x=pi*mod(x,2.d0)
      sini=(1.d0-tra)*cos(x)
      cosi=(1.d0+tra)*sin(x)
      kr=3
      x=abs(cosi)
      y=abs(sini)
      go to 3
   23 if(cosi)24,25,25
   24 v=pi-sign(v,sini)
      go to 26
   25 if(sini)27,26,26
   27 v=-v
   26 par=par-u
      pai=pai-v
      if(yd)100,100,28
   28 pai=-pai
!     argument dans -pi,pi
  100 tra=abs(pai/depi)
      if(tra-1.d+15)203,204,204
  204 nbchif=0
      pai=0.d0
      go to 29
  203 if(tra-1.d0)205,205,206
  206 nbchif=log10(tra)
      nbchif=14-nbchif
      tra=mod(tra,supint)
      pai=mod(tra,1.d0)*sign(depi,pai)
  205 if(abs(pai)-pi)29,29,207
  207 pai=pai-sign(depi,pai)
      go to 29
!     jflgam reel
    2 pai=0.d0
      if(xd)31,32,33
!     conditions d existence
   32 write (6,1000)
 1000 format (21h jflgam(0) est infini)
      go to 50
   31 if(x-4503599627370496.d0)34,35,35
   35 write (6,1001)
 1001 format (30h argument de jflgam trop grand)
      go to 50
   34 y=mod(x,supint)
      if(y)400,36,400
  400 if(y-0.99d0)33,33,405
  405 tra=int(y+0.1d0)
      if(abs(y-tra)-5.d-15)36,36,33
   36 write (6,1002)
 1002 format (28h jflgam (-entier) est infini)
   50 par=1.d+74
      nbchif=0
      go to 29
!     translation
   33 i=mod(10.99d0-x,supint)
      if(i)37,37,38
   38 tra=i
      x=x+tra
!     developpement asymptotique
   37 y=log(x)
      par=-x+hlo2pi+y*(x-0.5d0)
      if(x-test(1))39,39,43
   39 cosi=1.d0/x
      cos2i=cosi*cosi
      k=1
      go to 41
   42 cosi=cosi*cos2i
   41 par=par+c(k)*cosi
      k=k+1
      if(x-test(k))42,42,40
!     translation inverse
   44 x=x-1.d0
   48 y=log(x)
      par=par-y
      i=i-1
   40 if(i-1)43,49,44
   49 x=abs(xd)
      go to 48
!     x negatif
   43 if(xd)45,29,29
   45 par=alopi-par-y
      y=pi*mod(x,2.d0)
      y=-sin(y)
      if(y)46,36,47
   46 y=-y
      pai=pi
   47 par=par-log(y)
      entry jflgv1
   29 return
      end subroutine
!====================================================================!
      subroutine yfclen(eta, ro, u, up, v, vp, sigma0, idiv, nn)
      implicit none
      integer, intent(in) :: nn
      integer, intent(out) :: idiv
      integer :: indg, nb, n, m, k, j
      real(8), intent(in) :: eta, ro, sigma0
      real(8), intent(out) :: u, up, v, vp
      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: xa=0.577215664901533d0 
      real(8) :: ro2, etap, pieta, z, pieta2, par, pai, r, x, xx
      real(8) :: z1, u0, v0, v1, xn, u2, v2, pp, w, p, u1, xn1, wp
      real(8) :: eta2, t0, t1, tj, tm, tl, tk
      complex(8) :: fa, ak, am, al, bj, bn, bm, bl, bk, f, g, gp
      complex(8) :: c1, c2, c3, c4, c5, c6, d, hk, axpo
!
      if (nn == 1) then
!
!          serie origine
!
        ro2 = ro * ro
        etap = eta + eta
        pieta = pi * eta
        z = 138.15510557964276d0
        idiv = 0
        if (abs(pieta) > z) then
          indg = int(pieta / z)
          idiv = 60 * indg
          if (eta < 0) idiv = 0
          pieta = pieta - z * dble(indg)
        end if
        pieta2 = 0.5d0 * pieta
        p = exp(pieta2) * sqrt(sinh(pieta) / pieta)
        call jfdelg(1.0d0, eta, par, pai, nb)
        z1 = etap * (xa + xa + log(2.0d0) - 1.0d0 + par)
        u0 = 0.0d0
        u1 = ro
        v0 = 1.0d0
        v1 = z1 * ro
        u = u0 + u1
        v = v0 + v1
        up = 1.0d0
        vp = z1
        xn = 2.0d0
        do n=2, 10000
          xn1 = xn * (xn - 1.0d0)
          u2 = (etap * ro * u1 - ro2 * u0) / xn1
          u = u + u2
          v2 = (etap*ro*v1 - ro2*v0 - etap*(xn+xn-1.0d0)*u2) / xn1
          v = v + v2
          up = up + xn * u2 / ro
          vp = vp + xn * v2 / ro
          if (abs(u2/u) <= 1.0d-14 .and. abs(v2/v) <= 1.0d-14) exit
          u0 = u1
          u1 = u2
          v0 = v1
          v1 = v2
          xn = xn + 1.0d0
        end do
        pp = v + etap * u * log(ro)
        w = u / p
        wp = up / p
        v = p * pp
        vp = p * (vp + etap * (up * log(ro) + u / ro))
        u = w
        up = wp
        return
      end if
!
      eta2 = eta * eta
      fa = cmplx(1.0d0, eta, kind=8)
      m = 0.25 * eta + 4.0d1
!
!          polynomes de tchebichev jusqu'au rang m
!
      k = m + 2
      x = (eta + eta) / ro
      xx = x + x - 1.0d0
      t0 = 1.0d0
      t1 = xx
      xx = xx + xx
      do j=2, k
        tj = xx * t1 - t0
        t0 = t1
        t1 = tj
      end do
      tm = t1
      tl = t0
!
!          initialisation
!
      am = (0.0d0, 0.0d0)
      al = (1.0d-40, 1.0d-40)
      bn = (0.0d0, 0.0d0)
      bm = (1.0d-40, 1.0d-40)
      bl = (0.0d0, 0.0d0)
      bk = cmplx(4.0d0 * dble(m+1), 0.0d0, kind=8) * al + bm
      f = (0.0d0, 0.0d0)
      g = (0.0d0, 0.0d0)
      gp = (0.0d0, 0.0d0)
!
!          recurrence descendante
!
      k = m
      do 
        r = k
        tk = xx * tl - tm
        tm = tl
        tl = tk
        hk = cmplx(tk, 0.0d0, kind=8)
        c1 = cmplx(r*(r+1.0d0)-eta2, eta*(r+r+1.0d0), kind=8)
        c2 = (4.0d0, 0.0d0) * cmplx(r+1.0d0, 0.0d0, kind=8)    &
     &      * cmplx(-r-1.0d0, eta*3.0d0, kind=8)
        c3 = fa * cmplx(-r-r-4.0d0, eta, kind=8)
        c4 = cmplx((7.0d0*r+5.0d0)*0.25d0, 0.0d0, kind=8)
        c5 = cmplx(r+r+2.0d0, 0.0d0, kind=8)
        c6 = cmplx((r+3.0d0)*0.25d0, 0.0d0, kind=8)
        ak = (c2 * al + c3 * am - c4 * bl - c5 * bm - c6 * bn) / c1
        j = k / 2
        j = k - j - j
        if (j == 0) then
          f = f + ak
        else
          f = f - ak
        end if
        z = abs(ak)
        g = g + hk * ak
        gp = gp + hk * bk
!
!          f=a0/2-a1+a2-a3+a4-a5+...
!
!          congruence modulo 10**60
!
        if (z >= 1.0d60) then
          d = (1.0d60, 0.0d0)
          ak = ak / d
          al = al / d
          am = am / d
          bk = bk / d
          bl = bl / d
          bm = bm / d
          bn = bn / d
          f = f / d
          g = g / d
          gp = gp / d
        end if
        if (k <= 0) exit
        d = cmplx(4.0d0*r,0.0d0, kind=8)
        bj = d * ak + bl
        am = al
        al = ak
        bn = bm
        bm = bl
        bl = bk
        bk = bj
        k = k - 1
      end do
!
!          normalisation et calcul de z(ro)
!
      d = (0.5d0, 0.0d0) * ak
      f = f - d
      g = g - d
      gp = gp -(0.5d0, 0.0d0) * bk
      d = cmplx(0.0d0, -eta*log(2.0d0)+sigma0, kind=8)
      axpo = exp(d)
      f = f / axpo
      g = g / f
      gp = gp / f
!
!          calcul de f et g
!
      d = cmplx(0.0d0, ro-eta*log(ro), kind=8)
      axpo = exp(d)
      d = g * axpo
      gp = axpo * (cmplx(0.0d0, 1.0d0-eta/ro, kind=8) * g   &
     &           - cmplx(x/ro, 0.0d0, kind=8) * gp)
      v = d
      d = (0.0d0, -1.0d0) * d
      u = d
      vp = gp
      gp = (0.0d0, -1.0d0) * gp
      up = gp
      idiv = 0

      return
      end subroutine
!====================================================================!
      subroutine yfasym(eta, rau, fo, fpo, go, gpo, sigo, iexp)
      implicit none
      integer, intent(out) :: iexp
      integer :: n, ntruc
      real(8), intent(in) :: eta, rau
      real(8), intent(out) :: fo, fpo, go, gpo, sigo
      real(8) :: trb, rau2, etac, tra, ps, pt, gs, gt, sf, sg, spg
      real(8) :: spf, denom, an, bn, ps1, gs1, pt1, gt1, test, tetao

      iexp = 0
      trb = 0.0d0
      rau2 = rau + rau
      etac = eta * eta
      call jflgam(1.0d0, eta, tra, sigo, ntruc)
      n = 0
      ps = 1.0d0
      gs = 0.0d0
      pt = 0.0d0
      gt = 1.0d0 - eta / rau
      sf = ps
      sg = gs
      spf = pt
      spg = gt
      do
        denom = dble(n + 1) * rau2
        an = dble(n + n + 1) * eta / denom
        bn = (etac - dble(n * (n + 1))) / denom
        ps1 = an * ps - bn * pt
        gs1 = an * gs - bn * gt - ps1 / rau
        pt1 = an * pt + bn * ps
        gt1 = an * gt + bn * gs - pt1 / rau
        sf = sf + ps1
        sg = sg + gs1
        spf = spf + pt1
        spg = spg + gt1
        n = n + 1
        if (n == 17) then
          tra = ps * ps + pt * pt
          trb = ps1 * ps1 + pt1 * pt1
          test = tra - trb
          if (test <= 0.0d0) exit
        else if (n > 17) then
          trb = ps1 * ps1 + pt1 * pt1
          test = tra - trb
          if (test <= 0.0d0) exit
        end if
        ps = ps1
        gs = gs1
        pt = pt1
        gt = gt1
        tra = trb
      end do
      tetao = rau - eta * log(rau2) + sigo
      tra = sin(tetao)
      trb = cos(tetao)
      go = sf * trb - spf * tra
      gpo = sg * trb - spg * tra
      fo = spf * trb + sf * tra
      fpo = spg * trb + sg * tra

      return
      end subroutine
!====================================================================!
      subroutine dfcoul(eta, ro, f, fp, g, gp, sigma, l, iexp)
      implicit none
      integer, intent(in) :: l
      integer, intent(out) :: iexp(l+1)
      real(8), intent(in) :: eta, ro
      real(8), intent(out), dimension(l+1) ::f, fp, g, gp, sigma
      real(8), parameter :: depi=6.283185307179586d0
      real(8) :: etac, f0, fp0, g0, gp0, z, roinf, zig, zag, sig, xm
      real(8) :: ftest, fptest, fi, fpi, fl, fpl, fact, factp, a, b, s
      integer :: i, j, l1, ind, linf, lin, ipi, lmax, indice, i1, i2

      etac = eta * eta
      l1 = l + 1

      call dfcz0(eta, ro, f0, fp0, g0, gp0, s, i)
      
      f(1) = f0
      fp(1) = fp0
      g(1) = g0
      gp(1) = gp0
      iexp(1) = i
      sigma(1) = s

      if (l <= 0) return   ! angular momentum must be greater than 0.
      linf = 0
      ind = 0
      if (eta > 0.0d0 .and. ro < (eta + eta)) then
        lin = 1
      else
        z = eta + sqrt(etac + dble(l * (l + 1)))
        if (ro < z) then
          do 
            roinf = eta + sqrt(etac + dble(linf * (linf + 1)))
            if (ro < roinf) then
              ind = 1
              lin = linf + 1
              exit
            end if
            if (linf < l) then
              linf = linf + 1
            else
              lin = linf + 1
              exit
            end if
          end do
        end if

        xm = 1.0d0
        if (ind == 0) lin = l1

        do j=2, lin
          zig = sqrt(etac + xm * xm) / xm
          zag = eta / xm + xm / ro
          f(j) = (zag * f(j-1) - fp(j-1)) / zig
          fp(j) = zig * f(j-1) - zag * f(j)
          g(j) = (zag * g(j-1) - gp(j-1)) / zig
          gp(j) = zig * g(j-1) - zag * g(j)
          iexp(j) = i
          sig = sigma(j-1) + atan(eta / (j-1))
          ipi = int(sig / depi)
          sig = sig - dble(ipi) * depi
          if (sig < 0.0d0 .and. sig < - 0.5d0 * depi) sig = sig + depi
          if (sig > 0.0d0 .and. sig >   0.5d0 * depi) sig = sig - depi
          sigma(j) = sig
          xm = xm + 1.0d0
        end do
        if (ind == 0) return
      end if

      ftest = f(lin)
      fptest = fp(lin)
      lmax = linf + 25 + int(5.0d0 * abs(eta))
      if (lmax < l) lmax = l
      fi = 1.0d0
      fpi = 1.0d0

      do
        xm = dble(lmax + 1)
        zig = sqrt(etac + xm * xm) / xm
        zag = eta / xm + xm / ro
        fl = (zag * fi + fpi) / zig
        fpl = zag * fl - zig * fi
        if (abs(fl) >= 1.0d15) then
          fl = fl * 1.0d-15
          fpl = fpl * 1.0d-15
        else if (abs(fpl) >= 1.0d15) then
          fl = fl * 1.0d-15
          fpl = fpl * 1.0d-15
        end if
        fi = fl
        fpi = fpl
        if (lmax - l <= 0) then
          f(lmax+1) = fl
          fp(lmax+1) = fpl
          if (lmax <= linf) exit
        end if
        lmax = lmax - 1
      end do

      fact = ftest / f(lin)
      factp = fptest / fp(lin)
      indice = i / 60
      xm = linf

      do j=lin, l1
        f(j) = f(j) * fact
        fp(j) = fp(j) * factp
        if (j == 1) then
          xm = xm + 1
          cycle
        end if
        zig = sqrt(etac + xm * xm) / xm
        zag = eta / xm + xm / ro
        g(j) = (zag * g(j-1) - gp(j-1)) / zig
        gp(j) = zig * g(j-1) - zag * g(j)
        if (abs(g(j)) >= 1.0d60) then
          g(j) = g(j) / 1.0d60      !divide -> multiple  better? 
          gp(j) = gp(j) / 1.0d60
          indice = indice + 1
        else if (abs(gp(j)) >= 1.0d60) then
          g(j) = g(j) / 1.0d60
          gp(j) = gp(j) / 1.0d60
          indice = indice + 1
        end if
        iexp(j) = indice * 60
        a = fp(j) * g(j)
        b = - f(j) * gp(j)
        if (a < 1.0d0) then
          i1 = int(log10(a))
          i2 = int(log10(b))
        else
          i1 = int(log10(a)) + 1
          i2 = int(log10(b)) + 1
        end if
        f(j) = f(j) * 10.0d0 ** (- i2)
        fp(j) = fp(j) * 10.0d0 ** (- i1)
        sig = sigma(j-1) + atan(eta / (j-1))
        ipi = sig / depi
        sig = sig - ipi * depi
        if (sig < 0.0d0 .and. sig < - 0.5d0 * depi) then 
          sig = sig + depi
        else if (sig > 0.0d0 .and. sig > 0.5d0 * depi) then
          sig = sig - depi
        end if
        sigma(j) = sig
        xm = xm + 1.0d0
      end do

      return
      end subroutine
!====================================================================!
      subroutine yfireg(eta,ro,g0,gp0)
      implicit none
      integer :: iexp, n, nb
      real(8), intent(in) :: eta, ro
      real(8), intent(out) :: g0, gp0
      real(8) :: rau0, f0, fp0, sigma0, x, x2, x3, unr, etr0, u0, u1
      real(8) :: u2, s, v1, v2, t, xn, xn1, u3, v3, pi, ga, eta2, ro2
      real(8) :: etap, pieta, pieta2, b, par, pai, c1, v0, u, v, up
      real(8) :: vp, gp

      if (eta <= 0.0d0) then
        if(ro <= 0.5d0 * eta + 9.0d0)goto 200
        goto  300
      else if (eta <= 3.0d0) then
        if(ro <= 2.25d0 + 7.35d0 * (3.0d0 - eta))goto 200
        goto  300
      else if (eta <= 1.0d1) then
        if(ro <= 1.2d0 + 1.5d-1 * (1.0d1 - eta))goto 200
        goto  300
      else if (eta <= 18.0d0) then
        if(ro <= 0.6d0 + 0.75d-1 * (18.0d0 - eta))goto 200
        goto  300
      else if (eta <= 22.0d0) then
        if(ro <= 0.4d0 + 0.5d-1 * (22.0d0 - eta))goto 200
        goto  300
      end if
      if(ro <= 0.3d0 + (3.0d1 - eta) / 8.0d1)goto 200
!   serie de taylor depart rau0

  300 continue
      rau0 = 1.666666666666667d0 * abs(eta) + 7.5d0
      call yfasym(eta, rau0, f0, fp0, g0, gp0, sigma0, iexp)
      x = rau0 - ro
      x2 = x * x
      x3 = x * x2
      unr = 1.0d0 / rau0
      etr0 = 1.0d0 - 2.0d0 * eta * unr
      u0 = g0
      u1 = - x * gp0
      u2 = - 0.5d0 * etr0 * x2 * u0
      s = u0 + u1 + u2
      v1 = u1 / x
      v2 = 2.0d0 * u2 / x
      t = v1 + v2
      xn = 3.0d0
      do n=3, 10000
!       n=n
        xn1 = xn - 1.0d0
        xn1 = xn * xn1
        u3 = x*u2*unr*(1.0d0-2.d0/xn) - etr0*u1*x2/xn1 + x3*u0*unr/xn1
        s = s + u3
        v3 = xn * u3 / x
        t = t + v3
        if (abs(u3/s) <= 1.0d-11 .and. abs(v3/t) <= 1.0d-11) exit
        u0 = u1
        u1 = u2
        u2 = u3
        xn = xn + 1.0d0
      end do
      g0 = s
      gp0 = - t
      return

!   serie  origine
  200 continue
      pi = 3.141592653589793d0
      ga = 0.577215664901533d0
      eta2 = eta * eta
      ro2 = ro * ro
      etap = eta + eta
      pieta = pi * eta
      pieta2 = 0.5d0 * pieta
      b = exp(pieta2) * sqrt(sinh(pieta) / pieta)
      call jfdelg(1.0d0, eta, par, pai, nb)
      c1 = etap * (ga + ga + log(2.0d0) - 1.0d0 + par)
      u0 = 0.0d0
      u1 = ro
      v0 = 1.0d0
      v1 = c1 * ro
      u = u0 + u1
      v = v0 + v1
      up = 1.0d0
      vp = c1
      xn = 2.0d0
      do n=2, 10000
        xn1 = xn * (xn - 1.0d0)
        u2 = (etap * ro * u1 - ro2 * u0) / xn1
        u = u + u2
        v2 = (etap*ro*v1 - ro2*v0 - etap*(xn+xn-1.0d0)*u2) / xn1
        v = v + v2
        up = up + xn * u2 / ro
        vp = vp + xn * v2 / ro
        if (abs(u2/u) <= 1.0d-14 .and. abs(v2/v) <= 1.0d-14) exit
        u0 = u1
        u1 = u2
        v0 = v1
        v1 = v2
        xn = xn + 1.0d0
      end do
      gp = v + etap * u * log(ro)
      g0 = b * gp
      gp0 = b * (vp + etap * (up * log(ro) + u / ro))

      return
      end subroutine
!====================================================================!
      subroutine yfrica(eta, ro, fo, fpo, go, gpo, sigma0, idiv)
      implicit none
      integer, intent(out) :: idiv
      integer :: n, ind, jnd, ig, nn, indice, indg
      real(8), intent(in) :: eta
      real(8) :: ro
      real(8), intent(out) :: fo, fpo, go, gpo, sigma0
      real(8), parameter :: g61=0.1159057617187498d-1
      real(8), parameter :: g62=0.3863525390624998d-1
      real(8), parameter :: g63=0.4660034179687498d-1
      real(8), parameter :: g64=0.4858398437499998d-1
      real(8), parameter :: g65=0.1156514485677080d1
      real(8), parameter :: g66=0.5687475585937496d1
      real(8), parameter :: g67=0.1323888288225445d2
      real(8), parameter :: g68=0.1713083224826384d2
      real(8), parameter :: g69=0.1269003295898436d2
      real(8), parameter :: g610=0.5055236816406248d1
      real(8), parameter :: g611=0.8425394694010415d0
      real(8), parameter :: g81=0.1851092066083633d-01
      real(8), parameter :: g82=0.8638429641723630d-01
      real(8), parameter :: g83=0.1564757823944092d0
      real(8), parameter :: g84=0.1430139541625977d0
      real(8), parameter :: g85=0.1924622058868408d0
      real(8), parameter :: g86=0.8500803152720129d1
      real(8), parameter :: g87=0.7265429720878595d2
      real(8), parameter :: g88=0.3057942376817972d3
      real(8), parameter :: g89=0.7699689544836672d3
      real(8), parameter :: g810=0.1254157054424285d4
      real(8), parameter :: g811=0.1361719536066055d4
      real(8), parameter :: g812=0.9831831171035763d3
      real(8), parameter :: g813=0.4547869927883148d3
      real(8), parameter :: g814=0.1222640538215636d3
      real(8), parameter :: g815=0.1455524450256709d2
      real(8), parameter :: gp61=0.2897644042968748d-01
      real(8), parameter :: gp62=0.2318115234375000d0
      real(8), parameter :: gp63=0.8056640625000000d0
      real(8), parameter :: gp64=0.1601562499999998d1
      real(8), parameter :: gp65=0.3046875000000000d0
      real(8), parameter :: gp66=0.5624999999999998d1
      real(8), parameter :: gp81=0.6478822231292720d-01
      real(8), parameter :: gp82=0.6910743713378906d0
      real(8), parameter :: gp83=0.3322952270507811d1
      real(8), parameter :: gp84=0.9483032226562498d1
      real(8), parameter :: gp85=0.1769653320312499d2
      real(8), parameter :: gp86=0.3478710937499998d2
      real(8), parameter :: gp87=0.5020312499999999d2
      real(8), parameter :: gp88=0.7874999999999999d2
      real(8) :: tra, rau2, rauc, etac, eta2, etaro, etaro2, pieta, ro1
      real(8) :: rau0, x, u, x2, ru, rx, tre, trb, fi, tr1, tr2, tr3, s
      real(8) :: tr4, tr5, tr6, tr7, tr8, psip, xxx, psi, fip, eta0, et
      real(8) :: etad, et0, et1, et2, et3, et4, et5, x3, unr, etr0, u0
      real(8) :: u1, u2, u3, v1, v2, v3, t, xn, xn1, ho, hpo, trd, trc
      real(8) :: rau
      real(8) :: q(5)=(/0.4959570165d-1, 0.8888888889d-2,              &
     &              0.2455199181d-2, 0.9108958061d-3, 0.2534684115d-3/)
      real(8) :: qp(5)=(/0.1728260369d0,0.3174603174d-3,               &
     &              0.3581214850d-2, 0.3117824680d-3, 0.9073966427d-3/)
  
      rau = ro
      call jflgam(1.0d0, eta, tra, sigma0, ind)
      rau2 = rau + rau
      rauc = rau * rau
      etac = eta * eta
      eta2 = eta + eta
      etaro = eta * rau
      etaro2 = etaro + etaro
      pieta = 3.141592653589793d0 * eta
      ind = 0
      jnd = 0
      ig = 0

      if (eta <= 0.0d0) then
        if (etaro > - 14.0625d0) then 
          nn = 1
          call yfclen(eta, rau, fo, fpo, go, gpo, sigma0, idiv, nn)
          return
        else
          indice = 1
          idiv = 0
        end if
        goto 2                                    ! **********
      end if

      if (abs(rau - eta2) <= 1.0d-9) goto 18      ! **********
      if (rau == eta2) then
        goto 18                                   ! **********
      else if (rau > eta2) then
        if (rau >= eta2 + 2.0d1 * (eta ** 0.25d0)) then
          indice = 0
          idiv = 0
          goto 2                                  ! **********
        end if    
        nn = 0
      end if

      if (etaro-12.d0 < 0.0d0) then
        nn = 1
        call yfclen(eta,rau,fo,fpo,go,gpo,sigma0,idiv,nn)
        return
      end if
      tra = eta2 - 6.75d0 * (eta ** 0.4d0)
      if (rau > tra) then
        ind = 1
        jnd = 1
        ro1 = rau
        rau = tra
        rau0 = tra
      end if
!!              riccati  1
      x = rau / eta2
      u = (1.0d0 - x) / x
      x2 = x * x
      ru = sqrt(u)
      rx = sqrt(x)
      tre = 1.0d0 / (u * ru * eta2)
      trb = tre * tre
      fi = (sqrt((1.0d0-x)*x) + asin(rx) - 1.570796326794897d0) * eta2
      tr1 = - 0.25d0 * log(u)
      tr2 = - ((9.0d0 * u + 6.0d0) * u + 5.0d0) / 48.0d0
      tr3 =((((-3.0d0*u-4.0d0)*u+6.0d0)*u+12.0d0)*u+5.0d0)/64.0d0
      tr4 = - ((((((u + 2.0d0)*945.0d0*u + 1395.0d0)*u + 12300.0d0)  &
     &      *u + 25191.0d0)*u + 19890.0d0)*u + 5525.0d0) / 46080.0d0
      tr5 = ((((((((-27.d0*u-72.0d0)*u-68.0d0)*u+360.0d0)*u+2190.0d0)&
     &      *u+4808.0d0)*u+5148.0d0)*u+2712.0d0)*u+565.0d0)/2048.0d0
      tr6 = -(((((((((g61*u+g62)*u+g63)*u+g64)*u+g65)*u+g66)*u+g67)  &
     &      *u+g68)*u+g69)*u+g610)*u+g611
      tr7 = ((((((((((((-81.0d0*u-324.0d0)*u-486.0d0)*u-404.0d0)     &
     &      *u+4509.0d0)*u+52344.0d0)*u+233436.0d0)*u+567864.0d0)    &
     &      *u+838521.0d0)*u+775884.0d0)*u+441450.d0)                &
     &      *u+141660.d0)*u+19675.d0)/6144.d0
      tr8 =(((((((((((((g81*u+g82)*u+g83)*u+g84)*u+g85)*u+g86)*u+g87)&
     &  *u+g88)*u+g89)*u+g810)*u+g811)*u+g812)*u+g813)*u+g814)*u+g815
      psip = psip + tra
      xxx = 138.1551055796428d0
      fi = fi + tre * (tr2 + trb * (tr4 + trb * (tr6 + trb * tr8)))
      psi = - fi
      indg = int(psi / xxx)
      idiv = 60 * indg
      tra = tr1 + trb * (tr3 + trb * (tr5 + trb * tr7))
      fi = fi + tra
      psi = psi + tra
   
      fip = ru * eta2
      tra = 1.0d0 / (x2 * u)
      tr1 = 0.25d0
      tre = tre / (x2 * x2 * u)
      trb = trb / (x2 * x2)
      tr2 = -(8.0d0 * x - 3.0d0)/32.0d0
      tr3 = ((24.0d0 * x - 12.0d0) * x + 3.0d0) / 64.0d0
      tr4 = (((-1536.0d0*x+704.0d0)*x-336.0d0)*x+63.0d0)/2048.d0
      tr5 = ((((1920.0d0*x-576.0d0)*x+504.0d0)*x-180.0d0)            &
     &    * x+27.0d0)/1024.d0
      tr6 = ((((-gp66*x+gp65)*x-gp64)*x+gp63)*x-gp62)*x+gp61
      tr7 = - ((((((-40320.d0*x-10560.d0)*x-13248.d0)*x+7560.d0)     &
     &      * x -3132.0d0)*x+756.0d0)*x-81.0d0) / 2048.0d0
      tr8 = - ((((((gp88*x+gp87)*x+gp86)*x-gp85)*x+gp84)             &
     &      * x-gp83)*x+gp82)*x-gp81
      fip = fip + tre * (tr2 + trb * (tr4 + trb * (tr6 + trb * tr8)))
      tra = tra * (tr1 + trb * (tr3 + trb * (tr5 + trb * tr7)))
      fip = fip + tra
      psip = - fip
      if (indg /= 0) then
        psi = psi - xxx * dble(indg)
        fi  = fi  + xxx * dble(indg)
      end if
      fo = 0.5d0 * exp(fi)
      go = exp(psi)
      fpo = fo *  fip / eta2
      gpo = go * psip / eta2
      if (jnd == 0) return
      rau = ro1
      go  = fo
      gpo = fpo

      x = rau0 - ro1
      x = rau0 - ro1
      x2 = x * x
      x3 = x * x2
      unr = 1.0d0 / rau0
      etr0 = 1.0d0 - 2.0d0 * eta * unr
      u0 = go
      u1 = - x * gpo
      u2 = - 0.5d0 * etr0 * x2 * u0
      s = u0 + u1 + u2
      v1 = u1 / x
      v2 = 2.0d0 * u2 / x
      t = v1 + v2
      xn = 3.0d0
!
      do n=3, 10000
!!      n=n
        xn1 = xn - 1.0d0
        xn1 = xn * xn1
        u3 = x*u2*unr*(1.0d0-2.0d0/xn)-etr0*u1*x2/xn1+x3*u0*unr/xn1
        s = s + u3
        v3 = xn * u3 / x
        t = t + v3
        if (abs(u3/s) <= 1.0d-10 .and. abs(v3/t) <= 1.0d-10) exit
        u0 = u1
        u1 = u2
        u2 = u3
        xn = xn + 1.0d0
      end do
      if (ig /= 0) then
        go = s
        gpo = - t
        fo = ho
        fpo = hpo
        return
      end if
      ho = s
      hpo = - t

   18 et0 = eta ** (0.166666666666667d0)
      etad = etac * etac
      et = eta ** (0.6666666666666667d0)
      et1 = et * et
      et2 = et1 * et1
      et3 = et2 * et
      et4 = etad * et
      et5 = et4 * et
      fo = 1.0d0- q(1)/et1-q(2)/etac- q(3)/et3- q(4)/etad- q(5)/et5
      go = 1.0d0+ q(1)/et1-q(2)/etac+ q(3)/et3- q(4)/etad+ q(5)/et5
      fpo= 1.0d0+qp(1)/et+qp(2)/etac+qp(3)/et2+qp(4)/etad+qp(5)/et4
      gpo= 1.0d0-qp(1)/et+qp(2)/etac-qp(3)/et2+qp(4)/etad-qp(5)/et4
      fo = 0.7063326373d0 * et0 * fo
      go = 1.223404016d0 * et0 * go
      fpo = 0.4086957323d0 * fpo / et0
      gpo = - 0.7078817734d0 * gpo / et0
      idiv = 0
      if (ind == 0) return
      ig = 1
      rau0 = eta2
      x = rau0 - ro1
      x = rau0 - ro1
      x2 = x * x
      x3 = x * x2
      unr = 1.0d0 / rau0
      etr0 = 1.0d0 - 2.0d0 * eta * unr
      u0 = go
      u1 = - x * gpo
      u2 = - 0.5d0 * etr0 * x2 * u0
      s = u0 + u1 + u2
      v1 = u1 / x
      v2 = 2.0d0 * u2 / x
      t = v1 + v2
      xn = 3.0d0
!
      do n=3, 10000
!!      n=n
        xn1 = xn - 1.0d0
        xn1 = xn * xn1
        u3 = x*u2*unr*(1.0d0-2.0d0/xn)-etr0*u1*x2/xn1+x3*u0*unr/xn1
        s = s + u3
        v3 = xn * u3 / x
        t = t + v3
        if (abs(u3/s) <= 1.0d-10 .and. abs(v3/t) <= 1.0d-10) exit
        u0 = u1
        u1 = u2
        u2 = u3
        xn = xn + 1.0d0
      end do
      go = s
      gpo = - t
      fo = ho
      fpo = hpo
      return

   2  x = eta2 / rau
      x2 = x * x
      u = 1.0d0 - x
      ru = sqrt(u)
      u3 = u * u * u
      trd = 1.0d0 / (u3 * eta2 * eta2)
      trc = x2 * trd
      tre = 1.0d0 / (u * ru * eta2)
      fi = - 0.25d0 * log(u)
      trb = trd / 64.d0
      tr3 = (((3.0d0*u-4.0d0)*u-6.0d0)*u+12.0d0)*u-5.0d0
      tr5 = ((((((((-27.0d0*u+72.0d0)*u-68.0d0)*u-360.0d0)*u+2190.0d0) &
     &      *u-4808.0d0)*u+5148.0d0)*u-2712.0d0)*u+565.0d0)/32.0d0
      tr7= ((((((((((((81.0d0*u-324.0d0)*u+486.0d0)*u-404.0d0)         &
     &      *u-4509.0d0)*u+52344.0d0)*u-233436.0d0)*u+567864.0d0)      &
     &      *u-838521.0d0)*u+775884.0d0)*u-441450.0d0)*u+141660.0d0)   &
     &      *u-19675.0d0)/96.0d0
      fi = fi + trb * (tr3 + trd * (tr5 + trd * tr7))
 
      fip = 0.25d0 / u
      trb = 3.0d0 * trc / (64.0d0 * u)
      tr3 = (x - 4.0d0) * x + 8.0d0
      tr5 = ((((9.0d0*x-60.0d0)*x+168.0d0)*x-192.0d0)*x+640.0d0)/16.0d0
      tr7 = ((((((-27.0d0*x+252.0d0)*x-1044.0d0)*x+2520.0d0)           &
            *x-4416.0d0)*x-3520.0d0)*x-13440.0d0)/32.0d0
      fip = fip + trb * (tr3 + trc * (tr5 + trc * tr7))
      tra = abs((ru - 1.0d0) / (ru + 1.0d0))
      psi = (0.5d0 * log(tra) + ru / x) * eta2 + 0.785398163397448d0
      tr2 = - ((9.0d0 * u - 6.0d0) * u + 5.0d0) / 48.d0
      tr4 = ((((((u-2.0d0)*945.0d0*u+1395.0d0)*u-12300.0d0)            &
            *u+25191.0d0)*u-19890.0d0)*u+5525.0d0)/46080.0d0
      tr6 = (((((((((-g61*u+g62)*u-g63)*u+g64)*u-g65)*u+g66)*u-g67)    &
     &      *u+g68)*u-g69)*u+g610)*u-g611
      tr8 = (((((((((((((g81*u-g82)*u+g83)*u-g84)*u+g85)*u-g86)*u+g87) &
     &    *u-g88)*u+g89)*u-g810)*u+g811)*u-g812)*u+g813)*u-g814)*u+g815
      psi = psi + tre * (tr2 + trd * (tr4 + trd * (tr6 + trd * tr8)))
      psip = - ru * eta2 / x2
      trb = tre * x / u
      tr2 = (3.0d0 * x - 8.0d0) / 32.d0
      tr4 = - (((63.0d0*x-336.0d0)*x+704.0d0)*x-1536.0d0)/2048.0d0
      tr6 = ((((gp61*x-gp62)*x+gp63)*x-gp64)*x+gp65)*x-gp66
      tr8 = ((((((-gp81*x+gp82)*x-gp83)*x+gp84)*x-gp85)*x+gp86)        &
     &      *x+gp87)*x+gp88
      psip = psip + trb * (tr2 + trc * (tr4 + trc * (tr6 + trc * tr8)))
      tra = exp(fi)
      fo = tra * sin(psi)
      go = tra * cos(psi)
      if (indice /= 0) then
        tra = fo
        fo = go
        go = - tra
      end if
      tra = - eta2 / rauc
      fpo = (fip * fo + psip * go) * tra
      gpo = (fip * go - psip * fo) * tra

      return
      end subroutine
!======================================================================!
      subroutine dfcz0(eta, ro, f0, fp0, g0, gp0, sigma0, iexp)
      implicit none
      integer, intent(out) :: iexp
      integer :: i, j, ii, n, m, ntruc
      real(8), intent(in) :: eta, ro
      real(8), intent(out) :: f0, fp0, g0, gp0 ,sigma0
      real(8), dimension(110) :: a1, a2, b1, b2
      real(8), parameter :: pi=3.141592653589793d0
      real(8) :: borne, tra, ro2, etap, pieta, pieta2, b, u0, u1, u2
      real(8) :: xn, u, xn1, up, h, dh, di, s, q, d, c, x1, x2, t1, t2
      real(8) :: derive, ti, reslt, z, y, x

      if (ro <= 0.0d0) then
        write (6,*) 'ro negatif ou nul **'
        f0 = 0.0d0
        g0 = 0.0d0
        fp0 = 0.0d0
        gp0 = 0.0d0
        sigma0 = 0.0d0
        iexp = 0 
        return
      end if

      if (eta > 30.0d0 .or. eta < - 8.0d0) then
        if (abs(eta) <= 5.0d2) then
          call yfrica(eta, ro, f0, fp0, g0, gp0, sigma0, iexp)
        else
          f0 = 0.0d0
          g0 = 0.0d0
          fp0 = 0.0d0
          gp0 = 0.0d0
          sigma0 = 0.0d0
          iexp = 0
          write (6,*) 'valeur absolue de eta supe-&eu-e a 500 **'
        end if
        return
      end if

      if (eta == 0.0d0) then
        f0 = sin(ro)
        g0 = cos(ro)
        fp0 = g0
        gp0 = - f0
        iexp = 0
        sigma0 = 0.0d0
        return
      end if

      borne = 1.666666666666667d0 * abs(eta) + 7.5d0

      if (ro >= borne) then
        call yfasym(eta, ro, f0, fp0, g0, gp0, sigma0, iexp)
        return
      end if

      if (eta < 10.0d0) then
        if (eta > 0.0d0 .and. ro < 2.0d0) goto 14
      else
        if (eta > (5.0d0 * ro + 6.0d1) / 7.0d0) goto 14
      end if
      call yfasym(eta, borne, f0, fp0, g0, gp0, sigma0, iexp)
      h = borne
      dh = f0 / h
!!!!
      if (eta < 0.0d0) then
        n = - 0.5d0 * eta + 5.0d0
      else
        n = 0.2d0 * eta + 5.0d0
      end if
      n = 5 * (n + 1)
      z = 4.0d0 / h
      y = 1.0d0 - (eta + eta) * z
      a1(n+2:n+4) = (/1.0d-55, 0.0d0, 1.0d-64/)
      b1(n+3:n+4) = (/1.0d-50, 1.0d-68/)
      a2(n+2:n+4) = (/0.0d0, 1.0d-74, 1.0d-53/)
      b2(n+3:n+4) = (/0.0d0, 1.0d-66/)
      m = n + 2
      di = n

      do ii=2, m
        i = m - ii + 2
        b1(i) = b1(i+2) + z * (di + 1.0d0) * a1(i+1)
        s = a1(i+2) + y * (a1(i+1) - a1(i))
        q = (di + 2.0d0) * b1(i) + (di - 1.0d0) * b1(i+1)
        a1(i-1) = s - z * q
        b2(i) = b2(i+2) + z * (di + 1.0d0) * a2(i+1)
        s = a2(i+2) + y * (a2(i+1) - a2(i))
        q = (di + 2.0d0) * b2(i) + (di - 1.0d0) * b2(i+1)
        a2(i-1) = s - z * q
        if (i >= n) then
          di = di - 1.0d0
          cycle
        end if
        d = - (b2(i+2) + b2(i)) / (b1(i+2) + b1(i))
        do j=i, m                         ! to leave as do loop 
          a2(j) = a2(j) + d * a1(j)       ! is faster ?
          b2(j) = b2(j) + d * b1(j)
        end do
        a2(i-1) = a2(i-1) + d * a1(i-1)
        di = di - 1.0d0
      end do


      c = (a1(3) - a1(1)) / (a2(3) - a2(1))
!     x1 = dh/( 0.5d0*(a1(2)-c*a2(2)) + sum(a1(3:m))-c*sum(a2(3:m)) )
      x1 = 0.5d0 * (a1(2) - c * a2(2))
      do i=3, m
        x1 = x1 + a1(i) - c * a2(i)
      end do
      x1 = dh / x1

      x2 = - c * x1
      b1(2:m) = x1 * b1(2:m) + x2 * b2(2:m)
      a1(1:m) = x1 * a1(1:m) + x2 * a2(1:m)
      b1(1) = 0.0d0
      x = ro / h       ! ?????????
      y = 2.0d0 * x - 1.0d0
      t1 = 1.0d0
      t2 = y
      reslt = 0.5d0 * a1(2) + y * a1(3)
      derive = 0.5d0 * b1(2) + y * b1(3)

      do i=2, n
        ti = 2.0d0 * y * t2 - t1
        t1 = t2
        t2 = ti
        reslt = reslt + ti * a1(i+2)
        derive = derive + ti * b1(i+2)
      end do
      f0 = reslt * ro
      fp0 = derive * ro + reslt
      go to 30                            ! ???????????

      return

!! **   serie origine reguliere
 14     call jflgam(1.0d0, eta, tra, sigma0, ntruc)
        iexp = 0
        ro2 = ro * ro
        etap = eta + eta
        pieta = pi * eta
        pieta2 = 0.5d0 * pieta
        b = exp(pieta2) * sqrt(sinh(pieta) / pieta)
        u0 = 0.0d0
        u1 = ro
        u = u0 + u1
        up = 1.0d0
        xn = 2.0d0

        do n=2, 10000
          xn1 = xn * (xn - 1.0d0)
          u2 = (etap * ro * u1 - ro2 * u0) / xn1
          u = u + u2
          up = up + xn * u2 / ro
          if (abs(u2 / u) < 1.0d-10) exit
          u0 = u1
          u1 = u2
          xn = xn + 1.0d0
        end do
        f0 = u / b
        fp0 = up / b
   30 call yfireg(eta, ro, g0, gp0)

      return
      end subroutine
!======================================================================!
      subroutine jfdelg(xd, yd, par, pai, nbchif)
      implicit none
      integer, intent(out) :: nbchif
      integer :: kr, i, k
      real(8), intent(in) :: xd, yd
      real(8), intent(out) :: par, pai
      real(8), parameter :: rac2=0.3465735902799726d0
      real(8), parameter :: pis4=0.785398163397448d0
      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: depi=6.283185307179586d0
      real(8), parameter :: supint=2147483647.0d0
      real(8) :: test(7)=(/2.9152d7, 2.2958d3, 1.4124d2, 3.9522d1,     &
     &                    19.6611d0, 12.791d0, -10.d0/)
      real(8) :: c(6)=(/8.333333333333333d-2, -8.33333333333333d-3,    &
     &                  3.968253968253968d-3, -4.166666666666667d-3,   &
     &                  7.575757575757576d-3, -2.109279609279609d-2/)
      real(8) :: x, y, u, v, tra, tra1, trai, trb, cosi, cos2i, sini
      real(8) :: sin2i, zmod, xx

      x = abs(xd)
      xx = x
      nbchif = 15

      if (yd == 0.0d0) then
!       delgamma reel
        pai = 0.0d0
        if (xd == 0.0d0) then
!       conditions d existence
          write (6,*) 'jfdelg(0) est infini'
          par = 1.0d74
          nbchif = 0
          return
        else if (xd < 0.0d0) then
          if (x - 4503599627370496.d0 >= 0.0d0) then
            write (6,*) 'argument de jfdelg trop grand'
            par = 1.0d74
            nbchif = 0
            return
          else
            y = mod(x, supint)
            if (y == 0.0d0) then
              write (6,*) 'jfdelg (-entier) est infini'
              par = 1.0d74
              nbchif = 0
              return
            end if
            if (y > 0.99d0) then
              trai = int(y + 0.1d0)
              if (abs(y - tra) <= 5.0d-115) then
                write (6,*) 'jfdelg (-entier) est infini'
                par = 1.0d74
                nbchif = 0
                return
              end if
            end if
          end if
        end if

!       translation
        i = mod(10.99d0-x, supint)
        if (i > 0) then
          tra = i
          x = x + tra
        end if
!!      developpement asymptotique
        y = log(x)
        par = y - 0.5d0 / x
        if (x <= test(1)) then
          cos2i = 1.0d0 / (x * x)
          cosi = cos2i
          k = 1
          par = par - c(k) * cosi
          k = k + 1
          do while(x <= test(k))
            cosi = cosi * cos2i
            par = par - c(k) * cosi
            k = k + 1
          end do
!!      translation inverse
          do while(i > 0)
            i = i - 1
            x = i
            x = xx + x
            par = par - 1.0d0 / x
          end do
        end if
!!      x negatif
        if (xd < 0.0d0) then
          par = par + 1.0d0 / x
          y = pi * mod(x, 2.0d0)
          par = par + pi * cos(y) / sin(y)
        end if
        return
      end if
      
      y = abs(yd)
      kr = 1
      i = mod(10.99d0-x, supint)
!      translation
      if (i > 0) then
        tra = i
        x = x + tra
      end if
!      logarithme(x+iy) (x,y positifs)
      if (x < y) then
        tra1 = x / y
        trb = 1.0d0 + tra1 * tra1
        tra = y * sqrt(trb)
        sini = 1.0d0 / (trb * y)
        cosi = sini * tra1
        tra1 = y / x
        u = log(tra)
        v = atan(tra1)
      else if (x == y) then
        u = rac2 + log(x)
        v = pis4
        sini = 0.5d0 / x
        cosi = sini
      else   
        tra1 = y / x
        trb = 1.0d0 + tra1 * tra1
        tra = x * sqrt(trb)
        cosi = 1.0d0 / (trb * x)
        sini = cosi * tra1
        u = log(tra)
        v = atan(tra1)
      end if
!     developpement asymptotique ( x superieur a 10 )
      par = u - 0.5d0 * cosi
      pai = v + 0.5d0 * sini
      zmod = x + y
      if (zmod <= test(1)) then
        sin2i = (sini * cosi) + (sini * cosi)
        cos2i = (cosi + sini) * (cosi - sini)
        sini = sin2i
        cosi = cos2i
        k = 1
        do 
          par = par - c(k) * cosi
          pai = pai + c(k) * sini
          k = k + 1
          if (zmod > test(k)) exit
          tra = cosi * cos2i - sini * sin2i
          sini = sini * cos2i + cosi * sin2i
          cosi = tra
        end do
      end if

!     translation inverse
      do
        if (i <= 0) exit
        i = i - 1
        x = i
        x = xx + x
        if (x <= y) then
          tra1 = x / y
          trb = x * tra1 + y
          sini = 1.0d0 / trb
          cosi = tra1 / trb
        else
          tra1 = y / x
          trb = x + y * tra1
          cosi = 1.0d0 / trb
          sini = tra1 / trb
        end if
        par = par - cosi
        pai = pai + sini
      end do

!!    controle du quadrant
      if (xd == 0.0d0) then
        tra = pi * y
        if (tra <= 1.0d-2) then
          trb = tra*(2.0d0+tra*(-2.0d0+tra*(1.333333333333333d0+tra*(  &
     &    - 0.6666666666666666d0+tra*(0.2666666666666666d0+tra*(       &
     &    - 0.08888888888888888d0+tra*0.02539682539682540d0))))))
          trb = (2.0d0 - trb) / trb
        else
          trb = exp(- tra - tra)
          trb = (1.0d0 + trb) / (1.0d0 - trb)
        end if
        pai = 0.5d0 * (1.0d0 / y + pi * trb)
        if (yd < 0.0d0) then
          pai = - pai
          write (6,*) 'argument de jfdelg trop grand'
          par = 1.0d74
          nbchif = 0
          return
        end if
      else if (xd > 0.0d0) then
        if (yd < 0.0d0) then
          pai = - pai
          write (6,*) 'argument de jfdelg trop grand'
          par = 1.0d74
          nbchif = 0
          return
        end if
      else
!     x+iy change en -x-iy
        tra = exp(- depi * y)
        trb = tra * tra
        cos2i = depi * mod(x, 1.0d0)
        sin2i = - 2.0d0 * tra * cos(cos2i) + 1.0d0 + trb
        par = par + cosi + depi * tra * sin(cos2i) / sin2i
        pai = pai - sini + pi * (trb - 1.0d0) / sin2i
        if (yd > 0.0d0) then
          pai = - pai
          write (6,*) 'argument de jfdelg trop grand'
          par = 1.0d74
          nbchif = 0
          return
        else
        end if
      end if

!     argument dans -pi,pi
      tra = abs(pai / depi)
      if (tra >= 1.0d15) then
        nbchif = 0
        pai = 0.0d0
      else
        if (tra > 1.0d0) then
          nbchif = log10(tra)
          nbchif = 14 - nbchif
          tra = mod(tra, supint)
          pai = mod(tra, 1.0d0) * sign(depi, pai)
        end if
        if (abs(pai) > pi) pai = pai - sign(depi, pai)
      end if

      return
      end subroutine
!**********************************************************************!
      end module
