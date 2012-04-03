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

      implicit none
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, maux)

      integer maxmx,mbc,mx,meqn,maux
      double precision xlower,dx,t,dt

*     !local variables
      integer iflag(4)
      integer m,i,ig,mvars,j
      double precision u,em,theta,h,rho,taut

      logical allocated_true

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i=1-mbc,mx+mbc
*        !check for Nans
         do m=1,meqn
            if (q(i,m).ne.q(i,m)) then
               write(*,*) 'B4STEP: Nan in solution in component: ', m
               stop
               endif
         enddo
*        !check for admissible solution
         theta = aux(i,i_theta)
         h = q(i,1)
         call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),u,em,iflag)
      enddo

      if (pimin.eq.0.d0) p_initialized = 1
      if (init_ptype.gt.0.and.p_initialized.eq.0) then
         if (init_ptype.eq.2.and.t.le.init_ptf) then
            do i=1-mbc,mx+mbc
               h = q(i,1)
               if (h.le.dry_tol) then
                  cycle
               endif
               call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),
     &                                    u,em,iflag)
               rho = (rho_s*em + (1.d0-em)*rho_f)
               theta = aux(i,i_theta)
               if (q(i,4).lt.rho*grav*cos(theta)*h) then

                  q(i,4) = q(i,4) + (dt/init_ptf)*
     &            init_pmax_ratio*pimin*rho_f*grav*cos(theta)*h

                  call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),
     &                                    u,em,iflag)
               endif
            enddo
         elseif (init_ptype.eq.2.and.t.gt.init_ptf) then
            p_initialized = 1
         elseif (init_ptype.eq.1.and.p_initialized.eq.0) then
            do i=1-mbc,mx+mbc
               h = q(i,1)
               theta = aux(i,i_theta)
               q(i,4) = grav*cos(theta)*h*pimin*rho_f
               call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),
     &                                    u,em,iflag)
            enddo
            p_initialized = 1
         endif
      endif
      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

c      call diffuse_p(maxmx,mbc,mx,meqn,q,xlower,dx,t,dt)

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
