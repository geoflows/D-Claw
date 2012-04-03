      subroutine cranknicholson(u,r,rk,n)

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
