
c
c
c =========================================================
      subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit none

      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

c     integrate a source term for friction if ifrictiontype>0.

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx,t,dt

      integer i
      double precision h,hu,hm,p,u,D,tau,rho,theta,tanpsi,prate,peq,dp2
      double precision m,hrate,g


      include "digparamsdec.i"
      include "digparamscommon.i"

      g=grav


      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

      do i=0,mx+1

         if (q(i,1).gt.dry_tol) then
            h = q(i,1)
            hu= q(i,2)
            hm= q(i,3)
            p = q(i,4)
            u = hu/h
            D = aux(i,i_D)
            tau = aux(i,i_tau)
            rho = aux(i,i_rho)
            theta = aux(i,i_theta)
            tanpsi = aux(i,i_tanpsi)
            prate = -kappita/(h*mu)*
     &            ((2.d0/(alpha*h))+rho_f*g*cos(theta))
            peq = rho_f*g*h*cos(theta)
            dp2 = -2.d0*u*tanpsi/(alpha*h)

            h  = h + dt*D
            h = max(h,0.d0)
            if (u.gt.0.d0) then
               hu = hu + dt*(g*h*sin(theta) - tau/rho)
            elseif (u.lt.0.d0) then
               hu = hu + dt*(g*h*sin(theta) + tau/rho)
            else
c               hu = hu + dt*max(g*h*sin(theta)-tau/rho,0.d0)
                hu = hu + dt*g*h*sin(theta)
            endif

            p = peq + (p - peq)*exp(prate*dt)
            p = p + dt*dp2

            q(i,1) = h
            q(i,2) = hu
            q(i,4) = max(p,0.d0)
            q(i,4) = min(q(i,4),rho*h*g*cos(theta))
         endif
      enddo


      return
      end
