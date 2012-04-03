
c ======================================================================
      subroutine one_step(convtol,tn,dt,dtmax,theta,q0,qn,istep,iflag)
c     !find maximum possible dt and take it
c     ! dt is based on not progressing into an inadmissable state
c ======================================================================


      implicit none

*     !i/o
      integer iflag(4)
      logical istep
      double precision convtol,tn,dt,dtmax,theta
      double precision q0(4),qn(4)

*     !local
      integer md,mdmax,maxiter,miniter,iter,meq
      double precision sgnu,tau,rho,D,tanpsi,m,u
      include "digparamsdec.i"
      include "digparamscommon.i"

*     !try implicit solve

      call forwardeuler(dt,theta,q0,qn,istep,iflag)
      if (istep) then
         tn = tn + dt
         do meq = 1,4
            q0(meq) = qn(meq)
         enddo
      endif

      return
      end


c ======================================================================
      subroutine forwardeuler(dt,theta,q0,qn,istep,iflag)

c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      integer maxiter,iter
      double precision theta,convtol,dt
      double precision q0(4),qn(4)
      logical istep

*     !local
      double precision g,D,rho,h,hu,hm,u,m,P,tanpsi,tau,krate
      double precision h0,hu0,hm0,P0,dt2,krate0,kprate
      double precision Peq,lambda,kperm,compress,sig,mpos


      include "digparamsdec.i"
      include "digparamscommon.i"

      g = grav

      call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)

      if (iflag(1).gt.0) then
         istep = .false.
         return
      endif

      call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))

      h0 = q0(1)
      hu0 = q0(2)
      hm0 = q0(3)
      P0 = q0(4)
      dt2 = 0.5d0*dt
      krate0 = D*(rho - rho_f)/rho

      h = h0 + krate0*dt2
      hu = hu0*exp(dt2*krate0/h0)
      hu = hu0 + g*h0*sin(theta)*dt2
      hm = hm0*exp( -dt2*D*rho_f/(rho*h0))

      call admissibleq(theta,h,hu,hm,P0,u,m,iflag)
      call friction(theta,dt2,h,hu,hm,P0)
      call fluidfriction(theta,dt2,h,hu,hm,P0)

      call admissibleq(theta,h,hu,hm,P0,u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,P0)

      if (iflag(1).gt.0) then
c         dt = dt2
         istep = .true.
         return
      endif
      mpos = max(m,1.d-6)
      sig = max(0.d0,rho*g*h*cos(theta) - P0)
      kperm = (kappita**2*(1.d0-m)**3)/(180.d0*mpos**2)
      compress = alpha/(mpos*(sig + 1.d5 ))
      lambda = 3.d0/(compress*h) +
     &   (1.5d0*rho_f/(6.d0*rho))*(rho-rho_f)*g*cos(theta)
      Peq = rho_f*grav*h*cos(theta) -
     &            .75d0*(mu/kperm)*abs(hu)*tanpsi/(h*compress*lambda) !
      kprate = -lambda*2.d0*kperm/(compress*mu*h)
      P = Peq + (P0-Peq)*exp(kprate*dt)
      call admissibleq(theta,h,hu,hm,P,u,m,iflag)

      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,P)
      krate = D*(rho - rho_f)/rho
      h  =  h + krate*dt2
      hu = hu*exp(dt2*krate/h)
      hu = hu + g*h*sin(theta)*dt2
      hm = hm*exp( -dt2*D*rho_f/(rho*h))

      call admissibleq(theta,h,hu,hm,P,u,m,iflag)
      call friction(theta,dt2,h,hu,hm,P)
      call fluidfriction(theta,dt2,h,hu,hm,P)
      call admissibleq(theta,h,hu,hm,P,u,m,iflag)


      qn(1) = h
      qn(2) = hu
      qn(3) = hm
      qn(4) = P

      istep = .true.

      return
      end

c ======================================================================
      subroutine friction(theta,dt,h,hu,hm,p)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================


      implicit none

*     !i/o
      double precision theta,dt,h,hu,hm,p

*     !local
      integer iflag(4)
      double precision u,m,tau,rho,D,tanpsi

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,p)

      if (hu.gt.0.d0) then
         hu = max(hu - dt*tau/rho, 0.d0)
      else
         hu = min(hu + dt*tau/rho, 0.d0)
      endif
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      return
      end


c ======================================================================
      subroutine fluidfriction(theta,dt,h,hu,hm,p)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================


      implicit none

*     !i/o
      double precision theta,dt,h,hu,hm,p

*     !local
      integer iflag(4)
      double precision u,m,tau,rho,D,tanpsi,fterm

      include "digparamsdec.i"
      include "digparamscommon.i"
      return
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      if (iflag(1).gt.0) then
         return
      endif
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,p)

      fterm = 3.d0*(1.d0-m)*mu*u/(h*rho)
      if (hu.gt.0.d0) then
         hu = max(hu - dt*fterm, 0.d0)
      else
         hu = min(hu - dt*fterm, 0.d0)
      endif
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      return
      end




