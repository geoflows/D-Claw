
c
c
c =========================================================
      subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit none

*     !i/o
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx,t,dt

*     !local
      logical iflag(4)
      double precision psi(4)
      double precision p,u,m,h,hu,hm
      double precision g,t0,tf,kappa,theta,S,rho,tanpsi,D,tau,phi
      double precision sigbed,kperm,compress,pm
      double precision zeta,p_eq,p_hydro,p_litho,krate

      integer i

      include "digparamsdec.i"
      include "digparamscommon.i"

      if (p_initialized.eq.0) return
      g=grav
c     !recalculate aux due to hyperbolic step
      if (p_initialized.eq.0) return
      do i=1-mbc,mx+mbc
         call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),u,m,iflag)
         h = q(i,1)
         if (h.le.2.d0*dry_tol) then
            cycle
         endif
         hu = q(i,2)
         hm = q(i,3)
         p = q(i,4)

         call auxeval(h,u,m,p,kappa,theta,
     &              S,rho,tanpsi,D,tau,phi,sigbed,kperm,compress,pm)

c        ============pressure=====================================
         zeta = 6.d0/(compress*h*(1.d0+kappa))  +
     &        1.5d0*(rho-rho_f)*rho_f*grav*sin(theta)/(6.d0*rho)
         krate = -zeta*kperm/(h*mu)

         p_hydro = h*rho_f*grav*cos(theta)
         p_litho = (rho_s*m + (1.d0-m)*rho_f)*grav*h*cos(theta)
         p_eq = p_hydro
     &      - 3.d0*mu*abs(u)*tanpsi/(zeta*kperm*compress*(1.d0+kappa))
         p_eq = max(p_eq,p_hydro)
         p_eq = min(p_eq,p_litho)

         p = p_eq + (p-p_eq)*exp(krate*dt)
         call admissibleq(theta,h,hu,hm,p,u,m,iflag)
         q(i,4) = p
c        ============velocity======================================
c         D = (2.d0*kperm/(h*mu))*(rho_f*g*h*cos(theta)-p)
c         krate = (1.d0/(rho*h**2))*(h*D*(rho-rho_s)-(1.d0-m)*mu)
c         hu = hu*exp(krate*dt)
c         q(i,2) = hu
      enddo

      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

      return
      end


