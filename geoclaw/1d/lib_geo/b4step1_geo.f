c     ============================================
      subroutine b4step1(maxmx,mbc,mx,meqn,q,
     &            xlower,dx,t,dt,maux,aux)
c     ============================================
c
c     # called from claw1 before each call to step1.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.
c
c     # b4step1_geo checks for negative depths
c     # or vanishing depths. 
c     # and sets all values in that cell to zero

 
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)

      common /geo/ grav,abs_tol 

      do i=1-mbc,maxmx+mbc
         if (q(i,1).le.abs_tol) then
            if (q(i,1).lt.-abs_tol) then
               write(*,*) 'B4STEP1: negative depth set to zero'
               write(*,*) 'B4STEP1: location i =',i
               write(*,*) 'B4STEP1: value h = ', q(i,1)
            endif

            do m=1,meqn
               q(i,m)=0.d0
            enddo
         endif
      enddo

      return
      end
