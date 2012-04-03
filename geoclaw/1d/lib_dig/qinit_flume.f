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
      double precision h,pimin,theta,xcell

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         if (init_htype.eq.0) then
            h = max(0.d0,flumegeo(xcell))
         elseif (init_htype = 1) then
            !you must supply a function hinit(xcell)
         endif
         q(i,1) = h
         theta = aux(i,i_theta)
         q(i,2)= 0.d0
         q(i,3)= q(i,1)*m0
         q(i,4)= rho_f*grav*q(i,1)*cos(theta)
      enddo

      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

      if (init_ptype.eq.1) then
         pimin = 1.d0
         call calcpmin(maxmx,meqn,mbc,mx,xlower,dx,pimin,q,maux,aux)
         do i = 1,mx
            theta = aux(i,i_theta)
            q(i,4) = max(pimin,0.d0)*rho_f*grav*cos(theta)*q(i,1)
         enddo
         p_initialized = 1
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
c      h = 2.d0*exp(-(x)**2)

      flumegeo = h
      return
      end

c ================================================================
       subroutine calcpmin(maxmx,meqn,mbc,mx,xlower,dx,pimin,q,maux,aux)
c ================================================================

      implicit none
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx

      !local
      integer i
      double precision h,hr,hl,rhol,rhor,rho,thetal,thetar,theta,phi
      double precision del,taumag,pimin,pcrit

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i=1,mx
         hr = q(i,1)
         hl = q(i-1,1)
         h = 0.5d0*(hr + hl)
         if (h.lt.dry_tol) then
            cycle
         endif
         phi = 0.5d0*(aux(i,i_phi) + aux(i-1,i_phi))
         thetal = aux(i-1,i_theta)
         thetar = aux(i,i_theta)
         theta = 0.5d0*(thetal + thetar)

         rhol = aux(i-1,i_rho)
         rhor = aux(i,i_rho)
         rho = 0.5d0*(rhol+rhor)
*        ! determine pressure min ratio
         del = rho*0.5d0*grav*(cos(thetar)*hr**2 - cos(thetal)*hl**2)
     &         - rho*grav*h*sin(theta)*dx
         taumag = abs(del/dx)
         pcrit = rho*h*grav - taumag/tan(phi)
         pimin = min(pimin,pcrit/(rho_f*grav*cos(theta)*h))
      enddo

      return
      end

