c     ============================================
      subroutine diffuse_u(maxmx,mbc,mx,meqn,q,
     &            xlower,dx,t,dt)
c     ============================================

      implicit none
      double precision q(1-mbc:maxmx+mbc, meqn)

      integer maxmx,mbc,mx,meqn,maux
      double precision xlower,dx,t,dt

*     !local variables
      integer MXMAX,i
      parameter (MXMAX=5000)
      double precision visc,rk
      double precision r(1:MXMAX),u(1:MXMAX)

      include "digparamsdec.i"
      include "digparamscommon.i"
      return
      if (mx.gt.MXMAX) then
         write(*,*) 'diffuse_u: increase MXMAX to: ', mx
         stop
      endif

      visc = m0*mu
      rk = visc/(2.d0*dx**2)
      do i = 1,mx
         if (q(i,1).gt.dry_tol) then
            u(i) = q(i,2)/q(i,1)
         else
            u(i) = 0.d0
         endif
      enddo
      r(1)  = (1.d0-rk)*u(1) + rk*u(2)
      r(mx) = (1.d0-rk)*u(mx) + rk*u(mx-1)

      do i = 2,mx-1
         r(i) = rk*u(i-1) + (1.d0-2.d0*rk)*u(i) + rk*u(i+1)
      enddo

      call cranknicholson(u(1:mx),r(1:mx),rk,mx)

      do i = 1,mx
         q(i,2) = u(i)*q(i,1)
      enddo
      do i=1,mbc
         q(1-i,2) = q(1,2)
         q(mx+i,2) = q(mx,2)
      enddo

      return
      end

c     ------------------------------------------------------------------
      subroutine cranknicholson(u,r,rk,n)
c     ------------------------------------------------------------------

      integer n,NMAX
*     !i/o
      double precision a,b,c,rk,u(n),r(n)
*     !i/o
      PARAMETER (NMAX = 5000)
      integer j
      double precision bet,gam(NMAX)

      a = -rk
      b = 1.d0+2.d0*rk
      c = -rk

      bet = 1.d0 + rk!b
      u(1) = r(1)/bet

      do j = 2,n
         gam(j) = c/bet
         bet = b - a*gam(j)
         if (j.eq.n) bet = bet - rk
         if (bet.eq.0.) then
            write(*,*) 'tridiag failed'
            stop
         endif
         u(j) = (r(j)-a*u(j-1))/bet
      enddo

      do j = n-1,1,-1
         u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      return
      end
