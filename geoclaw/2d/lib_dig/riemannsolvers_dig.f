


c-----------------------------------------------------------------------
      subroutine riemann_dig2_conservative(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uhat,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,sw,fw,w)

      ! solve the dig Riemann problem for debris flow eqn
      ! this is for 2d version

      use geoclaw_module
      use digclaw_module

      implicit none

*     !i/o
      integer meqn,mwaves

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uhat,uL,uR,vL,vR,mL,mR
      double precision kappa,rho,kperm,compress,tanpsi,D,tau
      double precision gamma,gmod,dx

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision h,u,mbar
      double precision det1,det2,det3,detR
      double precision R(3,3),A(3,3),alph(3),delq(3),delf(3),del(3)
      double precision beta(3)



      if (hL.ge.drytolerance.and.hR.ge.drytolerance) then
         h = 0.5d0*(hL + hR)
         u = uhat
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytolerance) then
         h = hL
         u = uL
         mbar = mL
      else
         h = hR
         u = uR
         mbar = mR
      endif

*     !determine R
      R(1,1) = 1.d0
      R(2,1) = sw(1)
      R(3,1) = gamma*rho*gmod

      R(1,2) = kappa-1.d0
      R(2,2) = uhat*(kappa-1.d0)
      R(3,2) = kappa*rho*gmod

      R(1,3) = 1.d0
      R(2,3) = sw(3)
      R(3,3) = gamma*rho*gmod

      !determine delq
      delq(1) = hR-hL
      delq(2) = huR - huL
      delq(3) = pR - pL

      delf(1) = delq(2)
      delf(2) = hR*uR**2 + 0.5d0*kappa*gmod*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*gmod*hL**2)
      delf(3) = 0.d0

      del(1) = delf(1)
      del(2) = delf(2) + (1.d0-kappa)*delq(3)/rho
      del(3) = -rho*u*gmod*gamma*delq(1) + rho*gmod*gamma*delq(2)
     &         + u*delq(3)

*     !determine the source term
      do m=1,4
            psi(m)=0.d0
      enddo
      call psieval(tau,rho,D,tanpsi,kperm,compress,h,u,mbar,psi)


      del(1) = del(1) - dx*psi(1)
      del(3) = del(3) - dx*psi(4)
      if (abs(u).gt.1.d-6) then
         del(2) = del(2) - dx*psi(2) + (gmod*h*(bR-bL))
      else
         if (abs(delf(2)+gmod*h*(bR-bL)).ge.abs(dx*tau/rho)) then
            del(2)=sign(abs(abs(delf(2)+gmod*h*(bR-bL))
     &                   -abs(dx*tau/rho)),delf(2)-gmod*h*(bR-bL))
         else
            del(2)=0.d0
         endif
      endif

*     R alpha = delq
*     solve for alpha1 using Cramer's rule

      det1=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))
      det2=R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1))
      det3=R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      detR = (det1-det2+det3)


      do ks = 1,3 !solve for alpha(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = delq(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         alph(ks) = (det1-det2+det3)/detR
      enddo

      do ks = 1,3 !solve for beta(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = del(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         beta(ks) = (det1-det2+det3)/detR
      enddo

      do m=1,3
         mp=m
         if (m.eq.3) mp = 5
         do mw=1,3
            w(mp,mw)  = alph(mw)*R(m,mw)
            fw(mp,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

      !waves and fwaves for delta hum
      w(4,1) = alph(1)*mL
      w(4,3) = alph(3)*mR
      w(4,2) = hmR-hmL - mL*alph(1) - mR*alph(3)

      fw(4,1) = beta(1)*mL
      fw(4,3) = beta(3)*mR
      fw(4,2) = hmR*uR-hmL*uL - mL*beta(1)- mR*beta(3)-dx*psi(3)

      !waves and fwaves for delta huv

      w(3,1) = alph(1)*vL
      w(3,3) = alph(3)*vR
      w(3,2) = hvR-hvL - vL*alph(1) - vR*alph(3)

      fw(3,1) = beta(1)*vL
      fw(3,3) = beta(3)*vR
      fw(3,2) = hvR*uR - hvL*uL -beta(1)*vL - beta(3)*vR

      return
      end !subroutine riemann_dig2_conservative



c=============================================================================
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &             maxiter,drytol,g)

      !determine the Riemann structure (wave-type in each family)


      implicit none

      !input
      double precision hL,hR,uL,uR,drytol,g
      integer maxiter

      !output
      double precision s1m,s2m
      logical rare1,rare2

      !local
      double precision hm,u1m,u2m,um,delu
      double precision h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      integer iter



c     !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
         F_max= delu +
     &         (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            hm=(1.d0/(16.d0*g))*
     &               max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
            um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

            s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
            s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks

c           !root finding using a Newton iteration on sqrt(h)===
            h0=h_max
            do iter=1,maxiter
               gL=sqrt(.5d0*g*(1/h0 + 1/hL))
               gR=sqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+
     &                   gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*sqrt(h0)*dfdh
               h0=(sqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               s1m=u1m-sqrt(g*hm)
               s2m=u2m+sqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,maxiter
               F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max))
     &                  + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
               um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
               s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
               s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
               rare1=.true.
               rare2=.false.
            else
               s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
               s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
               um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return

      end ! subroutine riemanntype----------------------------------------------------------------


