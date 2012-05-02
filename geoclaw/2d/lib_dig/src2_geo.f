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
      double precision g,gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      double precision D,tau,sigbed,kperm,compress,pm,coeff,tol
      double precision zeta,p_hydro,p_litho,p_eq,krate,gamma,dgamma
      double precision cx,cy,pdh,vnorm
      integer i,j
c

c     # check for NANs in solution:
c     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

      if(p_initialized.eq.0) return
      gmod=grav
      pm=1.d0

      g=grav
      coeff = coeffmanning
      tol = 1.d-30  !# to prevent divide by zero in gamma

      !Pressure relaxation ---------------------------------------------

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

            !integrate momentum source term
            call auxeval(h,u,v,m,p,phi,kappa,S,rho,tanpsi,D,tau,
     &                  sigbed,kperm,compress,pm)
            hu = hu -dsign(tau/rho,hu)*dt
            hv = hv -dsign(tau/rho,hv)*dt
            hu = hu*dexp(-(1.d0-m)*mu*dt/h**2)
            hv = hv*dexp(-(1.d0-m)*mu*dt/h**2)

            if (hu*q(i,j,2).le.0.d0) then
               hu = 0.d0
            endif
            if (hv*q(i,j,3).le.0.d0) then
               hv = 0.d0
            endif
            vnorm = dsqrt(u**2 + v**2)
            call admissibleq(h,hu,hv,hm,p,u,v,m)


            !integrate pressure source term
            call auxeval(h,u,v,m,p,phi,kappa,S,rho,tanpsi,D,tau,
     &                  sigbed,kperm,compress,pm)
            !integrate shear-induced dilatancy
            p = p - dt*3.d0*dabs(vnorm)*tanpsi/(h*compress*(1.d0+kappa))

            zeta = 6.d0/(compress*h*(1.d0+kappa))  +
     &        1.5d0*(rho-rho_f)*rho_f*gmod/(6.d0*rho)

            krate=-zeta*kperm/(h*mu)

            p_hydro = h*rho_f*gmod
            p_litho = (rho_s*m + (1.d0-m)*rho_f)*gmod*h
            p_eq = p_hydro

            !p_eq = max(p_eq,0.d0)
            !p_eq = min(p_eq,p_litho)
            p = p_eq + (p-p_eq)*dexp(krate*dt)


            !integrate solid volume source term
            !krate = -rho_f*D/(rho*h)
            !hm = hm*dexp(krate*dt)

            !call admissibleq(h,hu,hv,hm,p,u,v,m)
            !call auxeval(h,u,v,m,p,phi,kappa,S,rho,tanpsi,D,tau,
c     &                  sigbed,kperm,compress,pm)


            call admissibleq(h,hu,hv,hm,p,u,v,m)
            q(i,j,1) = h
            q(i,j,2) = hu
            q(i,j,3) = hv
            q(i,j,4) = hm
            q(i,j,5) = p!h*pdh

         enddo
      enddo

*     ! Manning friction------------------------------------------------
      if (coeffmanning.gt.0.d0.and.frictiondepth.gt.0.d0) then
         do i=1,mx
            do j=1,my

               h=q(i,j,1)
               if (h.le.frictiondepth) then
c                 # apply friction source term only in shallower water
                  hu=q(i,j,2)
                  hv=q(i,j,3)

                  if (h.lt.tol) then
                     q(i,j,2)=0.d0
                     q(i,j,3)=0.d0
                  else
                     gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7/3))
                     dgamma=1.d0 + dt*gamma
                     q(i,j,2)= q(i,j,2)/dgamma
                     q(i,j,3)= q(i,j,3)/dgamma
                  endif
               endif
            enddo
         enddo
      endif
*     ! ----------------------------------------------------------------


c
      return
      end
