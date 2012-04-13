c
c
c =========================================================
      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &             q,maux,aux,t,dt)
c =========================================================
      use geoclaw_module
      use digclaw_module

      implicit none

      !i/o
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      double precision xlower,ylower,dx,dy,t,dt
      integer maxmx,maxmy,meqn,mbc,mx,my,maux

      !local
      double precision gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      double precision D,tau,sigbed,kperm,compress,pm
      double precision zeta,p_hydro,p_litho,p_eq,krate
      integer i,j
c

c     # check for NANs in solution:
c     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

      if(p_initialized.eq.0) return
      gmod=grav
      pm=1.d0


      do i=1,mx
         do j=1,my
            call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),
     &            q(i,j,4),q(i,j,5),u,v,m)
            h = q(i,j,1)
            if (h.le.drytolerance) cycle
            hu = q(i,j,2)
            hv = q(i,j,3)
            hm = q(i,j,4)
            p =  q(i,j,5)
            phi = aux(i,j,i_phi)

            call auxeval(h,u,v,m,p,phi,kappa,S,rho,tanpsi,D,tau,
     &                  sigbed,kperm,compress,pm)

            !integrate pressure source terms not included in Riemann solve
            zeta = 6.d0/(compress*h*(1.d0+kperm))  +
     &        1.5d0*(rho-rho_f)*rho_f*gmod/(6.d0*rho)

            p_hydro = h*rho_f*gmod
            p_litho = (rho_s*m + (1.d0-m)*rho_f)*gmod*h
            p_eq = p_hydro
     &       - 3.d0*mu*abs(u)*tanpsi/(zeta*kperm*compress*(1.d0+kperm))
            p_eq = max(p_eq,p_hydro)
            p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dt)

         enddo
      enddo
c
      return
      end
