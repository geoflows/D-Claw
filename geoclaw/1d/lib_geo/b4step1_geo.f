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

      use gauges_module

      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      logical allocated_true

      common /geo/ grav,abs_tol

      do i=1-mbc,maxmx+mbc
         if (q(i,1).le.abs_tol) then
            if (q(i,1).lt.-abs_tol) then
               write(*,*) 'B4STEP1: negative depth set to zero'
               write(*,*) 'B4STEP1: location i =',i
               write(*,*) 'B4STEP1: value h = ', q(i,1)
            endif
            q(i,1) = max(q(i,1),0.d0)
            do m=2,meqn
               q(i,m)=0.d0
            enddo
         endif
      enddo

      if (allocated(igauge)) then
         mvars = meqn+maux
         allocated_true = allocated(sol)
         if (.not.allocated_true) then
            allocate(sol(1:mvars))
         endif

         do ig = 1,mgauges
            if (t.ge.t0gauge(ig).and.t.le.tFgauge(ig)) then
               call return_gauge(meqn,maux,mvars,mx,dx,xlower,
     &                   xgauge(ig),sol(1:mvars),q(1:mx,1:meqn),
     &                                             aux(1:mx,1:maux))
               write(OUTGAUGEUNIT,*) igauge(ig),1,t,(sol(j),j=1,mvars)
            endif
         enddo
      endif

      return
      end
