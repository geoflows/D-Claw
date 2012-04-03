      subroutine tridiag(a,b,c,r,u,n)

      integer n,NMAX
*     !i/o
      double precision a(n),b(n),c(n),r(n),u(n)
*     !i/o
      PARAMETER (NMAX = 5000)
      integer j
      double precision bet,gam(NMAX)

      if (b(1).eq.0) write(*,*) 'tridiag: rewrite equations'
      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
         gam(j) = c(j-1)/bet
         bet = b(j) - a(j)*gam(j)
         if (bet.eq.0.) then
            write(*,*) 'tridiag failed'
            stop
         endif
         u(j) = (r(j)-a(j)*u(j-1))/bet
      enddo

      do j = n-1,1,-1
         u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      return
      end
