c
c
c ==============================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c ==============================================================

      implicit none
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx,flumegeo,hinit

      !local
      integer i
      double precision h,theta,xcell

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i=1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         if (init_htype.eq.0) then
            h = max(0.d0,flumegeo(xcell))
         elseif (init_htype.eq.1) then
            !you must supply a function hinit(xcell)
            h = max(0.d0,hinit(xcell))
         endif
         q(i,1) = h
         theta = aux(i,i_theta)
         q(i,2)= 0.d0
         q(i,3)= q(i,1)*m0
         q(i,4)= rho_f*grav*q(i,1)*cos(theta)
      enddo

      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

      p_initialized = 1
      if (init_ptype.gt.0) then
         p_initialized = 0
         pimin = 1.d0
         call calcpmin(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
         write(*,*) 'pimin', pimin
         do i = 1-mbc,mx+mbc
            theta = aux(i,i_theta)
            if (init_ptype.eq.1) then
                q(i,4) = pimin*rho_f*grav*cos(theta)*q(i,1)
            elseif (init_ptype.eq.2) then
                q(i,4) = 0.d0
            endif
         enddo
         if (init_ptype.eq.1) p_initialized = 1
      endif


      return
      end


c=======================================================================
      function flumegeo(x)
c
c     function to return the depth, given x, for the flume hopper
c     uses parameters set in digparams
c=======================================================================

      implicit none

      double precision x,h,flumegeo
      double precision pi,x0,x1,x2,x3,y0,y1,y2,y3,m

      include "digparamsdec.i"
      include "digparamscommon.i"

      pi = 3.14159
      x0 = -hopperlen
      x3 = 0.d0
      x2 = -hmax*cos(0.5d0*pi - theta1)
      x1 = x2 - hoppertop*cos(theta1 - hopperangle)

      y0 = 0.d0
      y2 = hmax*sin(0.5d0*pi - theta1)
      y1 = y2 - hoppertop*sin(theta1 - hopperangle)
      y3 = 0.d0

      if (x.ge.x0 .and. x.lt.x1) then
         m = (y1-y0)/(x1-x0)
         h = y0 + m*(x-x0)
      elseif (x.ge.x1 .and. x.lt.x2) then
         m = (y2-y1)/(x2-x1)
         h = y1 + m*(x-x1)
      elseif (x.ge.x2 .and. x.lt.x3) then
         m = (y3-y2)/(x3-x2)
         h = y2 + m*(x-x2)
      else
         h = 0.d0
      endif

      flumegeo = h
      return
      end

c ================================================================
       subroutine calcpmin(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c ================================================================

      implicit none
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx

      !local
      integer i
      double precision h,hr,hl,rhol,rhor,rho
      double precision thetal,thetar,theta,phil,phir,phi
      double precision forcemag,pcrit,drytol

      include "digparamsdec.i"
      include "digparamscommon.i"

      drytol = 2.d0*dry_tol

      do i=2-mbc,mx+mbc
         hr = q(i,1)
         hl = q(i-1,1)
         if (hl.le.dry_tol.and.hr.le.dry_tol) then
            cycle
         endif

         thetar = aux(i,i_theta)
         thetal = aux(i-1,i_theta)
         phir = aux(i,i_phi)
         phil = aux(i-1,i_phi)
         rhor = aux(i,i_rho)
         rhol = aux(i-1,i_rho)


         if (hl.gt.drytol.and.hr.gt.drytol) then
            h = 0.5d0*(hr + hl)
            phi = 0.5d0*(phir + phil)
            theta = 0.5d0*(thetar + thetal)
            rho = 0.5d0*(rhol + rhor)
         elseif (hr.gt.drytol) then
            h = hr
            phi = phir
            theta = thetar
            rho = rhor
         else
            h = hl
            phi = phil
            theta = thetal
            rho = rhol
         endif

*        !determine pressure min ratio
         forcemag = abs(grav*h*sin(theta)*dx
     &         -0.5d0*grav*(cos(thetar)*hr**2 - cos(thetal)*hl**2))

         pcrit = rho*h*grav*cos(theta) - rho*forcemag/(dx*tan(phi))
         pcrit = max(pcrit,0.0)
         pimin = min(pimin,pcrit/(rho_f*grav*cos(theta)*h))
         pimin = max(pimin,0.d0)

      enddo

      return
      end

